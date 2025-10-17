use core::time::Duration;

use crate::depth_utils::{get_depth, get_depth_idx};
use crate::dive::{Stop, StopSchedule};
use crate::gas::{GasMix, TissuesLoading};
use crate::mptt::{Tissue, MVALUES, NUM_STOP_DEPTHS, TISSUES};
use crate::pressure_unit::{msw, Pa, Pressure};
use crate::setup::{initialize_model_state, initialize_profile, set_m, LAST_STOP, NUM_TISSUES};
use crate::time_utils::get_time_ms_rel;
use crate::update::first_stop_depth;
#[cfg(not(feature = "thalmann"))]
pub use crate::update_exp::{
    compute_stop_time_exp as compute_stop_time, update_model_state_exp as update_model_state,
};
#[cfg(feature = "thalmann")]
pub use crate::update_thalmann::{
    compute_stop_time_thalmann as compute_stop_time,
    update_model_state_thalmann as update_model_state,
};

#[derive(Debug, Clone, PartialEq)]
pub enum ThalmannResult {
    FinishedResult {
        iterations: usize,
        reason: &'static str,
    },
    ErrorResult {
        reason: &'static str,
    },
}

pub const MVALUES_HE9_040: MVALUES = set_m(0);

pub fn thalmann<P: Pressure, const NUM_GASES: usize>(
    loading: &mut TissuesLoading<NUM_TISSUES, Pa>,
    max_depth: P,
    gases: &[GasMix<f32>; NUM_GASES],
) -> ThalmannResult
where
    [(); NUM_GASES - 1]:,
{
    if max_depth.to_msw() < LAST_STOP {
        return ThalmannResult::ErrorResult {
            reason: "Maximum depth is shallower than the last stop.",
        };
    }

    let mut prev_ms: usize = 0;
    let mut current_ms: usize = 0;

    initialize_profile();
    initialize_model_state();

    let m_values = &MVALUES_HE9_040;

    // Current depth idx
    let mut current_maximum_allowed_depth: msw = max_depth.to_msw();
    // let mut d_delta: usize = (max_depth.to_msw().0 / DINC.to_msw().0).ceil() as usize;

    // TODO: Select and be able to change current gas
    let current_gas = 0;

    // Last times
    get_time_ms_rel(&mut current_ms);

    // Iteration count
    let mut iter_count = 0;

    loop {
        iter_count += 1;
        get_time_ms_rel(&mut current_ms);
        while prev_ms >= current_ms {
            // TODO: Better than busy waiting
            get_time_ms_rel(&mut current_ms);
        }

        let current_depth = current_maximum_allowed_depth; // TODO: Measure actual depth
        if current_maximum_allowed_depth.to_msw().to_f32() <= 0.0 {
            return ThalmannResult::FinishedResult {
                iterations: iter_count,
                reason: "Current maximum allowed depth Smaller than 0",
            };
        }

        let duration = Duration::from_millis((current_ms - prev_ms) as u64);
        prev_ms = current_ms;
        update_model_state(
            loading,
            &TISSUES,
            m_values,
            &gases[current_gas],
            current_depth,
            &duration,
        );
        let first_stop = first_stop_depth(&loading, m_values);
        if first_stop.is_none() {
            return ThalmannResult::FinishedResult {
                iterations: iter_count,
                reason: "No First Stop remaining",
            };
        }
        current_maximum_allowed_depth = first_stop.unwrap();
        let _stop_time = compute_stop_time(
            &loading,
            &TISSUES,
            &gases[current_gas],
            m_values,
            current_maximum_allowed_depth,
        );
    }
}

pub fn calc_deco_schedule<const NUM_STOPS: usize>(
    loading: &TissuesLoading<NUM_TISSUES, Pa>,
    breathing_gas: &GasMix<f32>,
) -> Result<StopSchedule<NUM_STOPS>, &'static str> {
    calc_deco_schedule_intern(loading, &TISSUES, breathing_gas, &MVALUES_HE9_040)
}

pub fn calc_deco_schedule_intern<const NUM_TS: usize, const NUM_STOPS: usize>(
    loading: &TissuesLoading<NUM_TS, Pa>,
    tissues: &[Tissue; NUM_TS],
    breathing_gas: &GasMix<f32>,
    m_values: &MVALUES,
) -> Result<StopSchedule<NUM_STOPS>, &'static str> {
    assert!(NUM_STOPS < NUM_STOP_DEPTHS);

    let mut loading = loading.clone();
    let mut stops: [Stop; NUM_STOPS] =
        [Stop::new(msw::new(0.0), Duration::from_millis(0)); NUM_STOPS];

    for i in 0..NUM_STOPS {
        stops[stop_idx_in_stops(NUM_STOPS, i)] =
            Stop::new(get_depth(i).to_msw(), Duration::from_millis(0));
    }
    while let Some(stop_depth) = first_stop_depth(&loading, m_values) {
        let depth_idx = get_depth_idx(stop_depth);
        if depth_idx > NUM_STOPS {
            return Err("Not enough space to store stops for this dive.");
        }
        let depth_idx = stop_idx_in_stops(NUM_STOPS, depth_idx);
        if !stops[depth_idx].duration().is_zero() {
            return Err("Attempting to override / repeat stop.");
        }
        assert!(stop_depth == m_values[depth_idx].depth);
        let stop_duration =
            compute_stop_time(&loading, tissues, breathing_gas, m_values, stop_depth);
        update_model_state(
            &mut loading,
            tissues,
            &m_values,
            breathing_gas,
            stop_depth,
            &stop_duration,
        );
        stops[depth_idx] = Stop::new(stop_depth, stop_duration);
    }
    Ok(StopSchedule::new(stops))
}

const fn stop_idx_in_stops(num_stops: usize, i: usize) -> usize {
    return num_stops - 1 - i;
}

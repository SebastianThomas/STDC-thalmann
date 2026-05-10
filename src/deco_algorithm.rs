use core::time::Duration;

use crate::setup::{NUM_STOP_DEPTHS, NUM_TISSUES};

use crate::depth_utils::{get_depth, get_depth_idx};
use crate::dive::{Stop, StopSchedule};
use crate::gas::{GasDensitySettings, GasMix, TissuesLoading, best_available_mix};
use crate::mptt::{self, Tissue};
#[cfg(not(feature = "lin_exp"))]
use crate::mptt_buehlmann::{self, NUM_STOP_DEPTHS_BUEHLMANN, NUM_TISSUES_BUEHLMANN};
#[cfg(feature = "lin_exp")]
use crate::mptt_thalmann::{self, NUM_STOP_DEPTHS_THALMANN, NUM_TISSUES_THALMANN};
use crate::pressure_unit::{AbsPressure, Pa, Pressure, msw};
use crate::setup::{LAST_STOP, set_m};
use crate::time_utils::get_time_ms_rel;
use crate::update::first_stop_depth;
#[cfg(not(feature = "lin_exp"))]
pub use crate::update_exp::{
    compute_stop_time_exp as compute_stop_time, update_model_state_exp as update_model_state,
};
#[cfg(feature = "lin_exp")]
pub use crate::update_exp_lin::{
    compute_stop_time_lin_exp as compute_stop_time,
    update_model_state_lin_exp as update_model_state,
};

#[derive(Debug, Clone, PartialEq)]
pub enum DecoAlgorithmResult {
    FinishedResult {
        iterations: usize,
        reason: &'static str,
    },
    ErrorResult {
        reason: &'static str,
    },
}

#[cfg(feature = "lin_exp")]
pub const MVALUES: mptt::MValues<Pa, { NUM_TISSUES_THALMANN }, { NUM_STOP_DEPTHS_THALMANN }> =
    set_m(0);
#[cfg(not(feature = "lin_exp"))]
pub const MVALUES: mptt::MValues<Pa, { NUM_TISSUES_BUEHLMANN }, { NUM_STOP_DEPTHS_BUEHLMANN }> =
    set_m(0);
#[cfg(feature = "lin_exp")]
pub const TISSUES: [Tissue; NUM_TISSUES_THALMANN] = mptt_thalmann::TISSUES;
#[cfg(not(feature = "lin_exp"))]
pub const TISSUES: [Tissue; NUM_TISSUES_BUEHLMANN] = mptt_buehlmann::TISSUES;

pub type MValues<P: const AbsPressure> = mptt::MValues<P, { NUM_TISSUES }, { NUM_STOP_DEPTHS }>;

pub fn compute_deco_algorithm<P: Pressure, const NUM_GASES: usize>(
    loading: &mut TissuesLoading<NUM_TISSUES, Pa>,
    max_depth: P,
    gases: &[GasMix<f32>; NUM_GASES],
) -> DecoAlgorithmResult
where
    [(); NUM_GASES - 1]:,
{
    if max_depth.to_msw() < LAST_STOP {
        return DecoAlgorithmResult::ErrorResult {
            reason: "Maximum depth is shallower than the last stop.",
        };
    }

    let mut prev_ms: usize = 0;
    let mut current_ms: usize = 0;

    let m_values = &MVALUES;

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
            return DecoAlgorithmResult::FinishedResult {
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
            current_depth.to_pa(),
            &duration,
        );
        let first_stop = first_stop_depth(loading, m_values);
        if first_stop.is_none() {
            return DecoAlgorithmResult::FinishedResult {
                iterations: iter_count,
                reason: "No First Stop remaining",
            };
        }
        current_maximum_allowed_depth = first_stop.unwrap();
        let _stop_time = compute_stop_time(
            loading,
            &TISSUES,
            &gases[current_gas],
            m_values,
            current_maximum_allowed_depth,
        );
    }
}

#[cfg(feature = "lin_exp")]
type TissuesLoadingNumTissues<P> = TissuesLoading<{ NUM_TISSUES_THALMANN }, P>;
#[cfg(not(feature = "lin_exp"))]
type TissuesLoadingNumTissues<P> = TissuesLoading<{ NUM_TISSUES_BUEHLMANN }, P>;

pub fn calc_deco_schedule<const NUM_STOPS: usize, const NUM_GASES: usize>(
    loading: &TissuesLoadingNumTissues<Pa>,
    gases: &[GasMix<f32>; NUM_GASES],
    deco_settings: &DecoSettings<Pa>,
) -> Result<StopSchedule<NUM_STOPS>, &'static str> {
    calc_deco_schedule_intern(loading, &TISSUES, gases, &MVALUES, deco_settings)
}

pub struct DecoSettings<P: const AbsPressure> {
    pub gas_density_settings: GasDensitySettings<P>,
    pub max_deco_po2: P,
}

fn calc_deco_schedule_intern<
    const NUM_STOPS: usize,
    const NUM_GASES: usize,
    P: const AbsPressure,
>(
    loading: &TissuesLoading<NUM_TISSUES, P>,
    tissues: &[Tissue; NUM_TISSUES],
    gases: &[GasMix<f32>; NUM_GASES],
    m_values: &MValues<P>,
    deco_settings: &DecoSettings<P>,
) -> Result<StopSchedule<NUM_STOPS>, &'static str> {
    assert!(NUM_STOPS < NUM_STOP_DEPTHS);

    let mut loading = loading.clone();
    let mut stops: [Stop; NUM_STOPS] =
        [Stop::new(msw::new(0.0), Duration::from_millis(0), None); NUM_STOPS];

    for i in 0..NUM_STOPS {
        stops[stop_idx_in_stops(NUM_STOPS, i)] =
            Stop::new(get_depth(i).to_msw(), Duration::from_millis(0), None);
    }
    while let Some(stop_depth) = first_stop_depth(&loading, m_values) {
        let mix = best_available_mix(
            deco_settings.max_deco_po2,
            stop_depth.to_pa().into(),
            gases,
            &loading,
            false,
            &deco_settings.gas_density_settings,
        );
        if mix.is_none() {
            return Err("No gas for depth.");
        }
        let (_gas_idx, breathing_gas) = mix.unwrap();

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
            m_values,
            breathing_gas,
            stop_depth.to_pa().into(),
            &stop_duration,
        );
        stops[depth_idx] = Stop::new(stop_depth, stop_duration, Some(*breathing_gas));
    }
    Ok(StopSchedule::new(stops))
}

const fn stop_idx_in_stops(num_stops: usize, i: usize) -> usize {
    num_stops - 1 - i
}

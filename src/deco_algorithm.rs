use core::time::Duration;

use crate::setup::{NUM_STOP_DEPTHS, NUM_TISSUES};

use crate::depth_utils::{get_depth, get_depth_idx};
use crate::dive::{Stop, StopSchedule};
use crate::gas::{GasDensitySettings, GasMix, TissuesLoading, best_available_mix};
use crate::mptt;
#[cfg(feature = "lin_exp")]
use crate::mptt::Tissue;
#[cfg(not(feature = "lin_exp"))]
use crate::mptt_buehlmann::{self, BuehlmannTissue as Tissue};
#[cfg(feature = "lin_exp")]
use crate::mptt_thalmann::{self, NUM_STOP_DEPTHS_THALMANN, NUM_TISSUES_THALMANN};
use crate::pressure_unit::{AbsPressure, Pa, Pressure, msw};
use crate::setup::set_m;
use crate::update::first_stop_depth_with_gf;
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
pub const MVALUES: mptt::MValues<Pa, { NUM_TISSUES }, { NUM_STOP_DEPTHS }> = set_m(0);
#[cfg(feature = "lin_exp")]
pub const TISSUES: [Tissue; NUM_TISSUES_THALMANN] = mptt_thalmann::TISSUES;
#[cfg(not(feature = "lin_exp"))]
pub const TISSUES: [Tissue; NUM_TISSUES] = mptt_buehlmann::TISSUES;

pub type MValues<P: const AbsPressure> = mptt::MValues<P, { NUM_TISSUES }, { NUM_STOP_DEPTHS }>;

const STOP_SAFETY_MARGIN: Duration = Duration::from_secs(5);

pub struct DecoSettings<P: const AbsPressure> {
    pub gas_density_settings: GasDensitySettings,
    pub max_deco_po2: P,
    pub ignore_icd: bool,
    pub gf_low: f32,
    pub gf_high: f32,
    pub last_deco_stop: msw,
}

#[derive(Debug, Clone, Copy)]
pub struct GradientFactors {
    pub low: f32,
    pub high: f32,
}

fn compute_initial_first_stop<P: const AbsPressure>(
    loading: &TissuesLoading<NUM_TISSUES, P>,
    m_values: &MValues<P>,
    gf: GradientFactors,
) -> Option<msw> {
    first_stop_depth_with_gf(loading, m_values, gf.low)
}

fn interpolate_gf_for_depth(initial_first_stop: msw, stop_depth: msw, gf: GradientFactors) -> f32 {
    if initial_first_stop.to_msw().to_f32() <= 0.0 {
        gf.high
    } else {
        let t = (initial_first_stop.to_msw().to_f32() - stop_depth.to_msw().to_f32())
            / initial_first_stop.to_msw().to_f32();
        let t = if t < 0.0 {
            0.0
        } else if t > 1.0 {
            1.0
        } else {
            t
        };
        gf.low + (gf.high - gf.low) * t
    }
}

fn add_stop_safety_margin(stop_duration: Duration) -> Duration {
    stop_duration.saturating_add(STOP_SAFETY_MARGIN)
}

fn compute_next_stop_depth<P: const AbsPressure>(
    loading: &TissuesLoading<NUM_TISSUES, P>,
    m_values: &MValues<P>,
    initial_first_stop: msw,
    current_stop_depth: msw,
    gf: GradientFactors,
    last_deco_stop: msw,
) -> Option<msw> {
    if current_stop_depth.to_msw().to_f32() <= last_deco_stop.to_msw().to_f32() {
        return None;
    }

    let current_depth_idx = get_depth_idx(current_stop_depth);
    if current_depth_idx <= 1 {
        return None;
    }

    let next_depth = get_depth(current_depth_idx - 1).to_msw();
    let next_gf = interpolate_gf_for_depth(initial_first_stop, next_depth, gf);

    match first_stop_depth_with_gf(loading, m_values, next_gf) {
        Some(depth) if depth.to_msw().to_f32() < last_deco_stop.to_msw().to_f32() => {
            Some(last_deco_stop)
        }
        Some(depth) => Some(depth),
        None => {
            // No deeper first-stop found for the interpolated GF. If the next
            // shallower rung is at or above the configured `last_deco_stop`,
            // enforce that floor as the final stop; otherwise, no further
            // stops are required.
            if next_depth.to_msw().to_f32() <= last_deco_stop.to_msw().to_f32() {
                Some(last_deco_stop)
            } else {
                None
            }
        }
    }
}

#[cfg(feature = "lin_exp")]
type TissuesLoadingNumTissues<P> = TissuesLoading<{ NUM_TISSUES_THALMANN }, P>;
#[cfg(not(feature = "lin_exp"))]
type TissuesLoadingNumTissues<P> = TissuesLoading<{ NUM_TISSUES }, P>;

pub fn calc_deco_schedule<const NUM_STOPS: usize, const NUM_GASES: usize>(
    loading: &TissuesLoadingNumTissues<Pa>,
    gases: &[GasMix<f32>; NUM_GASES],
    gases_enabled: &[bool; NUM_GASES],
    deco_settings: &DecoSettings<Pa>,
) -> Result<StopSchedule<NUM_STOPS>, &'static str> {
    let gf = GradientFactors {
        low: deco_settings.gf_low,
        high: deco_settings.gf_high,
    };
    calc_deco_schedule_intern(
        loading,
        &TISSUES,
        gases,
        gases_enabled,
        &MVALUES,
        deco_settings,
        gf,
    )
}

fn calc_deco_schedule_intern<
    const NUM_STOPS: usize,
    const NUM_GASES: usize,
    P: const AbsPressure,
>(
    loading: &TissuesLoading<NUM_TISSUES, P>,
    tissues: &[Tissue; NUM_TISSUES],
    gases: &[GasMix<f32>; NUM_GASES],
    gases_enabled: &[bool; NUM_GASES],
    m_values: &MValues<P>,
    deco_settings: &DecoSettings<P>,
    gf: GradientFactors,
) -> Result<StopSchedule<NUM_STOPS>, &'static str> {
    assert!(NUM_STOPS < NUM_STOP_DEPTHS);

    let mut loading = loading.clone();
    let mut stops: [Stop; NUM_STOPS] =
        [Stop::new(msw::new(0.0), Duration::from_millis(0), None); NUM_STOPS];

    for i in 0..NUM_STOPS {
        stops[stop_idx_in_stops(NUM_STOPS, i)] =
            Stop::new(get_depth(i).to_msw(), Duration::from_millis(0), None);
    }
    // Determine initial first stop using GFLow; if none, return empty schedule
    let initial_first_stop = match compute_initial_first_stop(&loading, m_values, gf) {
        Some(d) => d,
        None => return Ok(StopSchedule::new(stops)),
    };

    let mut iterations: usize = 0;
    const MAX_ITER: usize = 1024;
    let mut next_stop = Some(initial_first_stop);
    while let Some(stop_depth) = next_stop {
        iterations += 1;
        if iterations > MAX_ITER {
            return Err("Exceeded max iterations building schedule");
        }
        let mix = best_available_mix(
            deco_settings.max_deco_po2,
            stop_depth.to_pa().into(),
            gases,
            gases_enabled,
            &loading,
            deco_settings.ignore_icd,
            &deco_settings.gas_density_settings,
        );
        if mix.is_none() {
            return Err("No gas for depth.");
        }
        let (_gas_idx, breathing_gas) = mix.unwrap();

        let depth_idx = get_depth_idx(stop_depth);
        #[cfg(feature = "lin_exp")]
        let mvalue_idx = depth_idx
            .checked_sub(1)
            .ok_or("Thalmann stop depths start at 3m.")?;
        if depth_idx > NUM_STOPS {
            return Err("Not enough space to store stops for this dive.");
        }
        let depth_idx = stop_idx_in_stops(NUM_STOPS, depth_idx);
        #[cfg(feature = "lin_exp")]
        assert!(stop_depth == m_values[mvalue_idx].depth);
        let gf_stop = interpolate_gf_for_depth(initial_first_stop, stop_depth, gf);

        let stop_duration = compute_stop_time(
            &loading,
            tissues,
            breathing_gas,
            m_values,
            stop_depth,
            gf_stop,
            deco_settings.last_deco_stop,
        );
        if stop_duration.is_zero() {
            // nothing to add and no progress — stop scheduling further stops
            break;
        }
        let stop_duration = add_stop_safety_margin(stop_duration);
        update_model_state(
            &mut loading,
            tissues,
            m_values,
            breathing_gas,
            stop_depth.to_pa().into(),
            &stop_duration,
        );

        // Merge repeated chunks at the same depth into a single scheduled stop.
        let existing = stops[depth_idx].duration();
        if existing.is_zero() {
            stops[depth_idx] = Stop::new(stop_depth, stop_duration, Some(*breathing_gas));
        } else {
            let new_total = existing.saturating_add(stop_duration);
            stops[depth_idx] = Stop::new(stop_depth, new_total, Some(*breathing_gas));
        }

        next_stop = compute_next_stop_depth(
            &loading,
            m_values,
            initial_first_stop,
            stop_depth,
            gf,
            deco_settings.last_deco_stop,
        );
    }
    Ok(StopSchedule::new(stops))
}

const fn stop_idx_in_stops(num_stops: usize, i: usize) -> usize {
    num_stops - 1 - i
}

#[cfg(not(feature = "lin_exp"))]
use core::{f32::consts::LN_2, time::Duration};

#[cfg(not(feature = "lin_exp"))]
use num::Float;

#[cfg(not(feature = "lin_exp"))]
use crate::{
    deco_algorithm::MValues,
    depth_utils::get_depth_idx,
    gas::{Gas, GasMix, HE_IDX, N2_IDX, TissuesLoading},
    mptt_buehlmann::{BuehlmannTissue, TISSUES},
    pressure_unit::{AbsPressure, Pa, Pressure, ambient_pressure_at_depth, msw},
    time_utils::max,
    update_common::exp_pressure,
};

#[cfg(not(feature = "lin_exp"))]
type BuehlmannLoading<P> = TissuesLoading<{ TISSUES.len() }, P>;

#[cfg(not(feature = "lin_exp"))]
type BuehlmannTissues = [BuehlmannTissue; TISSUES.len()];

#[cfg(not(feature = "lin_exp"))]
/**
* Schreiner Update:
* P(t) = P_{inspired} + (P_0 - P_{inspired}) * e^{-kt}
* <=> dP = P(t) - P_0 = (P_{inspired} - P_0) + (P_0 - P_{inspired}) * e^{-kt}
* <=> dP = P(t) - P_0 = (P_{inspired} - P_0) - (P_{inspired} - P_0) * e^{-kt}
* <=> dP = P(t) - P_0 = (P_{inspired} - P_0) * (1 - e^{-kt})
*/
pub fn update_model_state_exp<P: Pressure + const AbsPressure>(
    loading: &mut BuehlmannLoading<P>,
    tissues: &BuehlmannTissues,
    _m_values: &MValues<P>,
    breathing_gas: &GasMix<f32>,
    current_depth: P,
    delta_time: &Duration,
) {
    let delta_time_minutes: f32 = delta_time.as_secs_f32() / 60.0;
    let current_depth_pa: Pa = current_depth.to_pa();

    for gas_idx in [N2_IDX, HE_IDX] {
        let mut ks: [f32; TISSUES.len()] = [0.0; TISSUES.len()];
        for tissue_idx in 0..TISSUES.len() {
            ks[tissue_idx] = match gas_idx {
                N2_IDX => LN_2 / tissues[tissue_idx].n2.half_time,
                HE_IDX => LN_2 / tissues[tissue_idx].he.half_time,
                _ => unreachable!(),
            };
        }

        let (gas_loading, p_inspired): ([P; TISSUES.len()], P) = match gas_idx {
            N2_IDX => (loading.n2, breathing_gas.pn2(current_depth)),
            HE_IDX => (loading.he, breathing_gas.phe(current_depth)),
            _ => unreachable!(),
        };
        for tissue_idx in 0..TISSUES.len() {
            let p_tissue = gas_loading[tissue_idx];
            let k = ks[tissue_idx];
            let p_new = exp_pressure(p_inspired, p_tissue, k, delta_time_minutes);
            match gas_idx {
                N2_IDX => loading.n2[tissue_idx] = p_new,
                HE_IDX => loading.he[tissue_idx] = p_new,
                _ => unreachable!(),
            }
        }
    }
}

#[cfg(not(feature = "lin_exp"))]
/**
* P_{tissue}(t) = P_{inspired} + (P_{tissue,gas} - P_{inspired}) * e^{-kt}
* t_{tissue} = -1/k * ln((M - P_{inspired}) / (P_{tissue,gas} - P_{inspired}))
*/
pub fn compute_stop_time_exp<P: const AbsPressure>(
    loading: &BuehlmannLoading<P>,
    tissues: &BuehlmannTissues,
    breathing_gas: &GasMix<f32>,
    m_values: &MValues<P>,
    stop_depth: msw,
    gf: f32,
    surface_pressure: P,
    last_deco_stop: msw,
) -> Duration {
    let stop_depth_idx = get_depth_idx(stop_depth);
    let is_last_stop = stop_depth.to_msw().to_f32() <= last_deco_stop.to_msw().to_f32();

    let mut t_stop_mins: f32 = 0.0;
    // Use total inert pressure (N2 + He) per tissue and total inspired inert.
    let stop_ambient = ambient_pressure_at_depth(surface_pressure, stop_depth);
    let inspired_inert = stop_ambient * (breathing_gas.fn2() + breathing_gas.fhe());
    for tissue_idx in 0..TISSUES.len() {
        let p_n2 = loading.n2[tissue_idx].to_pa();
        let p_he = loading.he[tissue_idx].to_pa();
        let p_tissue = p_n2 + p_he;

        // Prefer mixed Buehlmann a/b when available; otherwise use table entry
        let m_value = if is_last_stop {
            crate::update_common::mixed_buehlmann_mvalue(
                tissue_idx,
                p_n2,
                p_he,
                surface_pressure,
            )
        } else if p_tissue.to_f32() > 0.0 {
            crate::update_common::mixed_buehlmann_mvalue(tissue_idx, p_n2, p_he, stop_ambient)
        } else {
            m_values[stop_depth_idx].max_saturation[tissue_idx].to_pa()
        };

        // Apply gradient factor to derive target stopping M-value
        let p_amb = stop_ambient;
        let target_m = super::update::allowed_with_gf(p_amb, m_value, gf);

        if target_m.to_f32() - p_tissue.to_f32() >= 0.0 {
            // Tissue is safe
            continue;
        }
        // Compute desaturation rate k as weighted per-gas rate using current
        // tissue partial pressures as weights.
        let k_n2 = LN_2 / tissues[tissue_idx].n2.half_time;
        let k_he = LN_2 / tissues[tissue_idx].he.half_time;
        let k = if p_tissue.to_f32() > 0.0 {
            (k_n2 * p_n2.to_f32() + k_he * p_he.to_f32()) / p_tissue.to_f32()
        } else {
            k_n2
        };
        let t_gas_tissue = -((target_m - inspired_inert) / (p_tissue - inspired_inert)).ln() / k;
        t_stop_mins = max(t_stop_mins, t_gas_tissue);
    }

    Duration::from_secs_f32(t_stop_mins * 60.0)
}

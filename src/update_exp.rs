#[cfg(not(feature = "lin_exp"))]
use core::{f32::consts::LN_2, time::Duration};

#[cfg(not(feature = "lin_exp"))]
use crate::{
    depth_utils::get_depth_idx,
    gas::{GasMix, HE_IDX, N2_IDX, TissuesLoading},
    mptt::{MVALUES, NUM_TISSUES, Tissue},
    mptt_buehlmann::BuehlmannTissue,
    pressure_unit::{Pa, Pressure, msw},
    time_utils::max,
    update_common::exp_pressure,
};

#[cfg(not(feature = "lin_exp"))]
/**
* Schreiner Update:
* P(t) = P_{inspired} + (P_0 - P_{inspired}) * e^{-kt}
* <=> dP = P(t) - P_0 = (P_{inspired} - P_0) + (P_0 - P_{inspired}) * e^{-kt}
* <=> dP = P(t) - P_0 = (P_{inspired} - P_0) - (P_{inspired} - P_0) * e^{-kt}
* <=> dP = P(t) - P_0 = (P_{inspired} - P_0) * (1 - e^{-kt})
*/
pub fn update_model_state_exp<P: Pressure>(
    loading: &mut TissuesLoading<NUM_TISSUES, Pa>,
    tissues: &[BuehlmannTissue; NUM_TISSUES],
    breathing_gas: &GasMix<f32>,
    current_depth: P,
    delta_time: Duration,
) {
    let delta_time_minutes: f32 = delta_time.as_secs_f32() / 60.0;
    let current_depth_pa: Pa = current_depth.to_pa();

    for gas_idx in [N2_IDX, HE_IDX] {
        // Build per-tissue KSAT (and use same for desat since Buehlmann doesn't have SDR)
        let mut ks: [f32; NUM_TISSUES] = [0.0; NUM_TISSUES];
        for tissue_idx in 0..NUM_TISSUES {
            ks[tissue_idx] = match gas_idx {
                N2_IDX => LN_2 / tissues[tissue_idx].n2.half_time,
                HE_IDX => LN_2 / tissues[tissue_idx].he.half_time,
                _ => unreachable!(),
            };
        }

        let (gas_loading, p_inspired): ([Pa; NUM_TISSUES], Pa) = match gas_idx {
            N2_IDX => (loading.n2, current_depth_pa * breathing_gas.n2()),
            HE_IDX => (loading.he, current_depth_pa * breathing_gas.he()),
            _ => unreachable!(),
        };
        for tissue_idx in 0..NUM_TISSUES {
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
pub fn compute_stop_time_exp(
    loading: &TissuesLoading<NUM_TISSUES, Pa>,
    tissues: &[BuehlmannTissue; NUM_TISSUES],
    breathing_gas: &GasMix<f32>,
    m_values: &MVALUES,
    stop_depth: msw,
) -> Duration {
    let stop_depth_idx = get_depth_idx(stop_depth);

    let mut t_stop_mins: f32 = 0.0;
    // Use total inert pressure (N2 + He) per tissue and total inspired inert.
    let inspired_inert = stop_depth.to_pa() * (breathing_gas.n2() + breathing_gas.he());
    for tissue_idx in 0..NUM_TISSUES {
        let p_n2 = loading.n2[tissue_idx];
        let p_he = loading.he[tissue_idx];
        let p_tissue = p_n2 + p_he;

        // Prefer mixed Buehlmann a/b when available; otherwise use table entry
        let m_value = if p_tissue.to_f32() > 0.0 {
            crate::update_common::mixed_buehlmann_mvalue(tissue_idx, p_n2, p_he, stop_depth.to_pa())
        } else {
            m_values[stop_depth_idx].max_saturation[tissue_idx]
        };

        if m_value.to_f32() - p_tissue.to_f32() >= 0.0 {
            // Tissue is safe
            continue;
        }
        assert!(p_tissue < inspired_inert);
        // Compute desaturation rate k as weighted per-gas rate using current
        // tissue partial pressures as weights.
        let k_n2 = LN_2 / tissues[tissue_idx].n2.half_time;
        let k_he = LN_2 / tissues[tissue_idx].he.half_time;
        let k = if p_tissue.to_f32() > 0.0 {
            (k_n2 * p_n2.to_f32() + k_he * p_he.to_f32()) / p_tissue.to_f32()
        } else {
            k_n2
        };

        let t_gas_tissue = -((m_value - inspired_inert) / (p_tissue - inspired_inert)).ln() / k;
        t_stop_mins = max(t_stop_mins, t_gas_tissue);
    }

    Duration::from_secs_f32(t_stop_mins * 60.0)
}

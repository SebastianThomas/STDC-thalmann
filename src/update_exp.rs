#[cfg(not(feature = "thalmann"))]
use core::{f32::consts::LN_2, time::Duration};

#[cfg(not(feature = "thalmann"))]
use crate::{
    depth_utils::get_depth_idx,
    gas::{GasMix, TissuesLoading, HE_IDX, N2_IDX},
    mptt::{Tissue, MVALUES, NUM_TISSUES},
    pressure_unit::{msw, Pa, Pressure},
    time_utils::max,
};

#[cfg(not(feature = "thalmann"))]
/**
* Schreiner Update:
* P(t) = P_{inspired} + (P_0 - P_{inspired}) * e^{-kt}
* <=> dP = P(t) - P_0 = (P_{inspired} - P_0) + (P_0 - P_{inspired}) * e^{-kt}
* <=> dP = P(t) - P_0 = (P_{inspired} - P_0) - (P_{inspired} - P_0) * e^{-kt}
* <=> dP = P(t) - P_0 = (P_{inspired} - P_0) * (1 - e^{-kt})
*/
pub fn update_model_state_exp<P: Pressure>(
    loading: &mut TissuesLoading<NUM_TISSUES, Pa>,
    tissues: &[Tissue; NUM_TISSUES],
    breathing_gas: &GasMix<f32>,
    current_depth: P,
    delta_time: Duration,
) {
    let delta_time_minutes: f32 = delta_time.as_secs_f32() / 60.0;
    let current_depth_pa: Pa = current_depth.to_pa();

    let tissue_ks = tissues.map(|tissue| (LN_2 / tissue.half_time) * tissue.sdr); // KSAT / KDSAT
    let tissue_ks_exp = tissue_ks.map(|k| (-k * delta_time_minutes).exp());
    let one_minus_tissue_ks_exp = tissue_ks_exp.map(|exp| 1.0 - exp);
    for gas_idx in [N2_IDX, HE_IDX] {
        let (gas_loading, p_inspired): ([Pa; NUM_TISSUES], Pa) = match gas_idx {
            N2_IDX => (loading.n2, current_depth_pa * breathing_gas.n2()),
            HE_IDX => (loading.he, current_depth_pa * breathing_gas.he()),
            _ => unreachable!(),
        };
        for tissue_idx in 0..NUM_TISSUES {
            let p_tissue = gas_loading[tissue_idx];
            let delta_p = (p_inspired - p_tissue) * one_minus_tissue_ks_exp[tissue_idx]; // TEXP
            match gas_idx {
                N2_IDX => loading.n2[tissue_idx] += delta_p,
                HE_IDX => loading.he[tissue_idx] += delta_p,
                _ => unreachable!(),
            }
        }
    }
}

#[cfg(not(feature = "thalmann"))]
/**
* P_{tissue}(t) = P_{inspired} + (P_{tissue,gas} - P_{inspired}) * e^{-kt}
* t_{tissue} = -1/k * ln((M - P_{inspired}) / (P_{tissue,gas} - P_{inspired}))
*/
pub fn compute_stop_time_exp(
    loading: &TissuesLoading<NUM_TISSUES, Pa>,
    tissues: &[Tissue; NUM_TISSUES],
    breathing_gas: &GasMix<f32>,
    m_values: &MVALUES,
    stop_depth: msw,
) -> Duration {
    let stop_depth_idx = get_depth_idx(stop_depth);

    let tissue_ks = tissues.map(|tissue| (LN_2 / tissue.half_time) * tissue.sdr); // KSAT / KDSAT

    let mut t_stop_mins: f32 = 0.0;
    for gas_idx in [N2_IDX, HE_IDX] {
        let (gas_loading, p_inspired): ([Pa; NUM_TISSUES], Pa) = match gas_idx {
            N2_IDX => (loading.n2, stop_depth.to_pa() * breathing_gas.n2()),
            HE_IDX => (loading.he, stop_depth.to_pa() * breathing_gas.he()),
            _ => unreachable!(),
        };
        for (tissue_idx, &p_tissue) in gas_loading.iter().enumerate() {
            let m_value = m_values[stop_depth_idx].max_saturation[tissue_idx];
            let delta_p = m_value - p_tissue;
            if delta_p.to_f32() >= 0.0 {
                // Tissue is safe
                continue;
            }
            assert!(p_tissue < p_inspired);

            let t_gas_tissue =
                -((m_value - p_inspired) / (p_tissue - p_inspired)).ln() / tissue_ks[tissue_idx];
            t_stop_mins = max(t_stop_mins, t_gas_tissue);
        }
    }

    Duration::from_secs_f32(t_stop_mins * 60.0)
}

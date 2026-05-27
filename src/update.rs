use core::time::Duration;

use crate::{
    deco_algorithm::{MValues, update_model_state},
    dive::{DiveMeasurement, DiveProfile},
    gas::{AIR, TissuesLoading},
    pressure_unit::{AbsPressure, Pressure, msw},
    setup::NUM_TISSUES,
};

#[cfg(feature = "lin_exp")]
use crate::depth_utils::get_depth_idx;
#[cfg(feature = "lin_exp")]
use crate::mptt::Tissue;
#[cfg(not(feature = "lin_exp"))]
use crate::mptt_buehlmann::BuehlmannTissue as Tissue;

pub fn allowed_with_gf<P: const AbsPressure>(p_amb: P, target: P, gf: f32) -> P {
    p_amb + (target - p_amb) * gf
}

#[cfg(not(feature = "lin_exp"))]
pub fn tissue_mvalues_with_gf<P: const AbsPressure, const NUM_TISSUES: usize>(
    loading: &TissuesLoading<NUM_TISSUES, P>,
    current_depth: msw,
    tissue_idx: usize,
    gf: f32,
) -> (P, P) {
    #[cfg(not(feature = "lin_exp"))]
    {
        let p_n2 = loading.n2[tissue_idx].to_pa();
        let p_he = loading.he[tissue_idx].to_pa();
        let absolute = crate::update_common::mixed_buehlmann_mvalue(
            tissue_idx,
            p_n2,
            p_he,
            current_depth.to_pa(),
        );
        let gf_mvalue = allowed_with_gf(current_depth.to_pa(), absolute, gf);
        (absolute.into(), gf_mvalue.into())
    }
}

#[cfg(feature = "lin_exp")]
pub fn tissue_mvalues_with_gf<P: const AbsPressure>(
    m_values: &MValues<P>,
    current_depth: msw,
    tissue_idx: usize,
    gf: f32,
) -> (P, P) {
    {
        let depth_idx = get_depth_idx(current_depth);
        let table_idx = depth_idx.min(m_values.len().saturating_sub(1));
        let absolute = m_values[table_idx].max_saturation[tissue_idx];
        let gf_mvalue = allowed_with_gf(current_depth.to_pa().into(), absolute, gf);
        (absolute, gf_mvalue)
    }
}

pub fn first_stop_depth_with_gf<P: const AbsPressure>(
    p: &TissuesLoading<{ NUM_TISSUES }, P>,
    m_values: &MValues<P>,
    gf: f32,
) -> Option<msw> {
    for mvalues_at_depth in m_values.iter().rev() {
        #[allow(clippy::needless_range_loop)]
        for i in 0..NUM_TISSUES {
            // Prefer mixed-gas a/b when available (Buehlmann); otherwise use
            // precomputed max_saturation table entry.
            #[cfg(not(feature = "lin_exp"))]
            {
                // Buehlmann: compute allowed value via helper to avoid duplication
                let p_n2 = p.n2[i].to_pa();
                let p_he = p.he[i].to_pa();
                let total = p_n2 + p_he;
                if total.to_f32() <= 0.0 {
                    // No inert present, skip
                    continue;
                }
                let (_absolute, allowed) = tissue_mvalues_with_gf(p, mvalues_at_depth.depth, i, gf);
                let total_p: P = (p_n2 + p_he).into();
                if total_p > allowed {
                    return Some(mvalues_at_depth.depth);
                }
            }

            #[cfg(feature = "lin_exp")]
            {
                // Thalmann / precomputed M-values: compare total inert vs table
                let total_inert = p.n2[i] + p.he[i];
                let p_amb: P = mvalues_at_depth.depth.to_pa().into();
                let mval = mvalues_at_depth.max_saturation[i];
                let allowed = allowed_with_gf(p_amb, mval, gf);
                if total_inert > allowed {
                    return Some(mvalues_at_depth.depth);
                }
            }
        }
    }
    None
}

pub fn loadings_from_dive_profile<
    const NUM_GASES: usize,
    const NUM_MEASUREMENTS: usize,
    P: const AbsPressure,
>(
    tissues: &[Tissue; NUM_TISSUES],
    profile: &DiveProfile<P, f32, NUM_GASES, NUM_MEASUREMENTS>,
    m_values: &MValues<P>,
    surface: P,
) -> TissuesLoading<NUM_TISSUES, P> {
    let mut loadings = TissuesLoading::new(surface, &AIR);
    for w in profile.measurements.windows(2) {
        assert!(w.len() == 2);
        let DiveMeasurement {
            time_ms: time_ms_prev,
            depth: depth_prev,
            gas: _gas_prev,
        } = &w[0];
        let DiveMeasurement {
            time_ms,
            depth,
            gas,
        } = &w[1];
        let delta_time = Duration::from_millis((time_ms - time_ms_prev) as u64);
        let midpoint = (*depth + *depth_prev) / 2.0;
        update_model_state(
            &mut loadings,
            tissues,
            m_values,
            &profile.gases[*gas],
            midpoint,
            &delta_time,
        );
    }
    loadings
}

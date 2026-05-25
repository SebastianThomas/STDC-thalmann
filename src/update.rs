use core::time::Duration;

use crate::{
    deco_algorithm::{MValues, update_model_state},
    dive::{DiveMeasurement, DiveProfile},
    gas::{AIR, TissuesLoading},
    mptt::Tissue,
    pressure_unit::{AbsPressure, msw},
    setup::NUM_TISSUES,
};

pub fn first_stop_depth<P: const AbsPressure>(
    p: &TissuesLoading<{ NUM_TISSUES }, P>,
    m_values: &MValues<P>,
) -> Option<msw> {
    for mvalues_at_depth in m_values.iter().rev() {
        #[allow(clippy::needless_range_loop)]
        for i in 0..NUM_TISSUES {
            // Prefer mixed-gas a/b when available (Buehlmann); otherwise use
            // precomputed max_saturation table entry.
            #[cfg(not(feature = "lin_exp"))]
            {
                use crate::mptt_buehlmann::TISSUES as BUEHL_TISSUES;
                let p_n2 = p.n2[i];
                let p_he = p.he[i];
                let total = p_n2 + p_he;
                if total.to_f32() <= 0.0 {
                    // No inert present, skip
                    continue;
                }
                let a_n2 = BUEHL_TISSUES[i].n2.a.to_pa();
                let a_he = BUEHL_TISSUES[i].he.a.to_pa();
                let b_n2 = BUEHL_TISSUES[i].n2.b;
                let b_he = BUEHL_TISSUES[i].he.b;

                let a_mix = (a_n2 * p_n2 + a_he * p_he) / total;
                let b_mix = (b_n2 * p_n2 + b_he * p_he) / total;

                // Buehlmann form: P_tol = a + P_amb / b
                let p_tol = a_mix + mvalues_at_depth.depth.to_pa() / b_mix;
                let total_inert = total;
                if total_inert > p_tol {
                    return Some(mvalues_at_depth.depth);
                }
            }

            #[cfg(feature = "lin_exp")]
            {
                // Thalmann / precomputed M-values: compare total inert vs table
                let total_inert = p.n2[i] + p.he[i];
                if total_inert > mvalues_at_depth.max_saturation[i]
                {
                    return Some(mvalues_at_depth.depth);
                }
            }
        }
    }
    None
}

pub fn loadings_from_dive_profile<
    const NUM_TISSUES: usize,
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

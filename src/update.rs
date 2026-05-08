use core::time::Duration;

use crate::{
    dive::{DiveMeasurement, DiveProfile},
    gas::{AIR, HE_IDX, N2_IDX, TissuesLoading},
    mptt::{MVALUES, Tissue},
    pressure_unit::{AbsPressure, msw},
    deco_algorithm::update_model_state,
};

pub fn first_stop_depth<const NUM_TISSUES: usize, P: const AbsPressure>(
    p: &TissuesLoading<NUM_TISSUES, P>,
    m_values: &MVALUES<P>,
) -> Option<msw> {
    for mvalues_at_depth in m_values.iter().rev() {
        for &gas_idx in [N2_IDX, HE_IDX].iter() {
            let tissue_loadings = match gas_idx {
                N2_IDX => &p.n2,
                HE_IDX => &p.he,
                _ => unreachable!(),
            };
            for i in 0..NUM_TISSUES {
                if tissue_loadings[i] > mvalues_at_depth.max_saturation[i] {
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
    P: const AbsPressure
>(
    tissues: &[Tissue; NUM_TISSUES],
    profile: &DiveProfile<P, f32, NUM_GASES, NUM_MEASUREMENTS>,
    m_values: &MVALUES<P>,
    surface: P,
) -> TissuesLoading<NUM_TISSUES, P> {
    let mut loadings = TissuesLoading::new(surface, &AIR);
    for w in profile.measurements.windows(2) {
        assert!(w.len() != 2);
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
        update_model_state(
            &mut loadings,
            tissues,
            m_values,
            &profile.gases[*gas],
            *depth + ((*depth - *depth_prev) / 2.0),
            &delta_time,
        );
    }
    loadings
}

use core::{f32::consts::LN_2, time::Duration};
#[allow(unused)]
use num::Float;

use crate::{
    depth_utils::get_depth_idx,
    gas::{Gas, GasMix, TissuesLoading, HE_IDX, N2_IDX},
    mptt::{Tissue, MVALUES},
    pressure_unit::{msw, Pa, Pressure},
    time_utils::max,
};

#[cfg(feature = "thalmann")]
pub fn update_model_state_thalmann<P: Pressure, const NUM_TISSUES: usize>(
    loading: &mut TissuesLoading<NUM_TISSUES, Pa>,
    tissues: &[Tissue; NUM_TISSUES],
    &m_values: &MVALUES,
    breathing_gas: &GasMix<f32>,
    current_depth: P,
    delta_time: &Duration,
) {
    use crate::pressure_unit::Bar;

    let delta_time_minutes: f32 = delta_time.as_secs_f32() / 60.0;
    let current_depth: Pa = current_depth.to_pa();

    let p_inspired_n2 = breathing_gas.pn2(current_depth).to_pa();
    let p_inspired_he = breathing_gas.phe(current_depth).to_pa();

    let k_values: [f32; NUM_TISSUES] = tissues.map(|t| (LN_2 / t.half_time) * t.sdr);

    for (gas_idx, p_inspired) in [(N2_IDX, p_inspired_n2), (HE_IDX, p_inspired_he)] {
        let gas_loading = match gas_idx {
            N2_IDX => &mut loading.n2,
            HE_IDX => &mut loading.he,
            _ => unreachable!(),
        };

        for tissue_idx in 0..NUM_TISSUES {
            let p_old = gas_loading[tissue_idx];
            let k: f32 = k_values[tissue_idx];
            // M-value for this tissue at this depth
            let depth_idx = get_depth_idx(current_depth);
            let m_value: Bar = m_values[depth_idx].max_saturation[tissue_idx].to_bar();
            // Crossover time
            let t_x = -((m_value.to_pa() - p_inspired) / (p_old - p_inspired)).ln() / k;

            let p_crossover = get_exp_pressure(p_inspired, p_old, k, t_x);
            let r: Pa = (p_crossover - p_inspired) * k;

            let p_new: Pa = if t_x >= delta_time_minutes {
                // Exponential phase only
                get_exp_pressure(p_inspired, p_old, k, delta_time_minutes)
            } else if t_x <= 0.0 {
                // Linear only
                p_old - r.to_pa() * delta_time_minutes
            } else {
                // Split: exponential + linear
                let t_lin: f32 = delta_time_minutes - t_x;
                p_crossover.to_pa() - r.to_pa() * t_lin
            };
            gas_loading[tissue_idx] = p_new;
        }
    }
}

#[cfg(feature = "thalmann")]
pub fn compute_stop_time_thalmann<const NUM_TISSUES: usize>(
    loading: &TissuesLoading<NUM_TISSUES, Pa>,
    tissues: &[Tissue; NUM_TISSUES],
    breathing_gas: &GasMix<f32>,
    m_values: &MVALUES,
    stop_depth: msw,
) -> Duration {
    let stop_depth_pa = stop_depth.to_pa();
    let stop_idx = get_depth_idx(stop_depth);
    let mut t_stop_mins = 0.0;

    let k_values: [f32; NUM_TISSUES] = tissues.map(|t| LN_2 / t.half_time * t.sdr);

    for (gas_idx, p_inspired) in [
        (N2_IDX, stop_depth_pa * breathing_gas.fn2()),
        (HE_IDX, stop_depth_pa * breathing_gas.fhe()),
    ] {
        let gas_loading = match gas_idx {
            N2_IDX => &loading.n2,
            HE_IDX => &loading.he,
            _ => unreachable!(),
        };

        for tissue_idx in 0..NUM_TISSUES {
            let p_tissue = gas_loading[tissue_idx];
            let k = k_values[tissue_idx];
            let m_value = m_values[stop_idx].max_saturation[tissue_idx];

            if p_tissue <= m_value {
                continue; // tissue already safe
            }

            // exponential → linear
            let t_x = -((m_value - p_inspired) / (p_tissue - p_inspired)).ln() / k;

            let t_tissue = if t_x <= 0.0 {
                // linear only
                let r = (p_tissue - p_inspired) * k;
                (p_tissue - m_value) / r
            } else {
                // exponential first, then linear if needed
                let p_cross = p_inspired + (p_tissue - p_inspired) * (-k * t_x).exp();
                if p_cross <= m_value {
                    t_x // only exponential needed
                } else {
                    let r = (p_cross - p_inspired) * k;
                    t_x + (p_cross - m_value) / r
                }
            };

            t_stop_mins = max(t_stop_mins, t_tissue);
        }
    }
    Duration::from_secs_f32(t_stop_mins * 60.0)
}

fn get_exp_pressure<P: const Pressure>(p_inspired: P, p_old: P, k: f32, t_x: f32) -> P {
    let exp: f32 = (-k * t_x).exp();
    p_inspired + (p_old - p_inspired) * exp
}

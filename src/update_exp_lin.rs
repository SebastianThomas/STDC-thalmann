use core::time::Duration;
#[allow(unused)]
use num::Float;

#[cfg(feature = "lin_exp")]
use crate::deco_algorithm::MValues;
use crate::{
    depth_utils::get_depth_idx,
    gas::{Gas, GasMix, HE_IDX, N2_IDX, TissuesLoading},
    mptt::Tissue,
    pressure_unit::{AbsPressure, Pa, Pressure, msw},
    time_utils::max,
    update_common::exp_pressure,
};

#[cfg(feature = "lin_exp")]
const THALMANN_FSW_TO_PA: f32 = 3_064.305_931_38;
#[cfg(feature = "lin_exp")]
const THALMANN_PVO2: Pa = Pa::new(2.0 * THALMANN_FSW_TO_PA);
#[cfg(feature = "lin_exp")]
const THALMANN_PVCO2: Pa = Pa::new(2.3 * THALMANN_FSW_TO_PA);
#[cfg(feature = "lin_exp")]
const THALMANN_PH2O: Pa = Pa::new(0.0);
#[cfg(feature = "lin_exp")]
const THALMANN_PBOVP: Pa = Pa::new(0.0);
#[cfg(feature = "lin_exp")]
const LIN_EXP_STOP_EPSILON_PA: Pa = Pa::new(10.0);

#[cfg(feature = "lin_exp")]
fn thalmann_crossover_pressure<P: const AbsPressure>(ambient_pressure: P) -> P {
    // Thalmann defines the linear->exponential crossover at PVSAT + PBOVP.
    // PVSAT = PAMB - (PVO2 + PVCO2 + PH2O), and PBOVP is the configurable
    // gas-phase overpressure threshold. The paper's default examples use 0.
    ambient_pressure - THALMANN_PVO2.into() - THALMANN_PVCO2.into() - THALMANN_PH2O.into()
        + THALMANN_PBOVP.into()
}

#[cfg(feature = "lin_exp")]
fn thalmann_mvalue_idx(stop_depth: msw) -> usize {
    get_depth_idx(stop_depth)
        .checked_sub(1)
        .expect("Thalmann stop depths start at 3m")
}

#[cfg(feature = "lin_exp")]
pub fn update_model_state_lin_exp<P: const AbsPressure, const NUM_TISSUES: usize>(
    loading: &mut TissuesLoading<NUM_TISSUES, P>,
    tissues: &[Tissue; NUM_TISSUES],
    _m_values: &MValues<P>,
    breathing_gas: &GasMix<f32>,
    current_depth: P,
    delta_time: &Duration,
) {
    let delta_time_minutes: f32 = delta_time.as_secs_f32() / 60.0;

    let p_inspired_n2 = breathing_gas.pn2(current_depth);
    let p_inspired_he = breathing_gas.phe(current_depth);
    let crossover_pressure = thalmann_crossover_pressure(current_depth);

    // KSAT and KDSAT per tissue (KDSAT = KSAT * SDR)
    let (k_values_sat, k_values_desat) = crate::update_common::ks_arrays(tissues);

    for (gas_idx, p_inspired) in [(N2_IDX, p_inspired_n2), (HE_IDX, p_inspired_he)] {
        let gas_loading = match gas_idx {
            N2_IDX => &mut loading.n2,
            HE_IDX => &mut loading.he,
            _ => unreachable!(),
        };

        for tissue_idx in 0..NUM_TISSUES {
            let p_old = gas_loading[tissue_idx];
            let k_sat: f32 = k_values_sat[tissue_idx];
            let k_desat: f32 = k_values_desat[tissue_idx];
            if p_old <= p_inspired {
                // Ongassing pure exponential uses KSAT
                gas_loading[tissue_idx] =
                    exp_pressure(p_inspired, p_old, k_sat, delta_time_minutes);
                continue;
            }

            // For desaturation computations use KDSAT.
            // When the tissue tension is above the Thalmann crossover pressure,
            // use linear washout first. Once the tissue tension is no longer too
            // large, the remaining tail is exponential.
            let k = k_desat;

            let p_new = if crossover_pressure <= p_inspired || p_old <= crossover_pressure {
                // Exponential only: the tissue is already below the crossover pressure.
                exp_pressure(p_inspired, p_old, k, delta_time_minutes)
            } else {
                // Linear first until the crossover pressure is reached.
                let linear_rate = (p_old - p_inspired) * k;
                let t_linear = (p_old - crossover_pressure) / linear_rate;

                if delta_time_minutes <= t_linear {
                    p_old - linear_rate * delta_time_minutes
                } else {
                    let t_exp = delta_time_minutes - t_linear;
                    exp_pressure(p_inspired, crossover_pressure, k, t_exp)
                }
            };
            if p_new < Pa::new(0.0).into() {
                panic!("Illegal p_new: {:?}", p_new);
            }
            gas_loading[tissue_idx] = p_new;
        }
    }
}

#[cfg(feature = "lin_exp")]
pub fn compute_stop_time_lin_exp<const NUM_TISSUES: usize, P: const AbsPressure>(
    loading: &TissuesLoading<NUM_TISSUES, P>,
    tissues: &[Tissue; NUM_TISSUES],
    breathing_gas: &GasMix<f32>,
    m_values: &MValues<P>,
    stop_depth: msw,
    gf: f32,
    last_deco_stop: msw,
) -> Duration {
    let stop_depth_pa = stop_depth.to_pa();
    let stop_idx = thalmann_mvalue_idx(stop_depth);
    let is_last_stop = stop_depth.to_msw().to_f32() <= last_deco_stop.to_msw().to_f32();
    let mut t_stop_mins = 0.0;
    let crossover_pressure: P = thalmann_crossover_pressure(stop_depth_pa.into());

    // KDSAT values for desaturation (KSAT * SDR)
    let (_k_values_sat, k_values_desat) = crate::update_common::ks_arrays(tissues);

    // Use total inert pressure (N2 + He) per tissue and total inspired inert.
    let inspired_inert_pa = stop_depth_pa * (breathing_gas.fn2() + breathing_gas.fhe());
    let p_inspired: P = inspired_inert_pa.into();
    let surface_targets: [P; NUM_TISSUES] = core::array::from_fn(|tissue_idx| {
        let shallow = m_values[0].max_saturation[tissue_idx];
        let next = m_values[1].max_saturation[tissue_idx];
        shallow - (next - shallow)
    });

    for tissue_idx in 0..NUM_TISSUES {
        let p_tissue =
            loading.n2[tissue_idx] + loading.he[tissue_idx] + LIN_EXP_STOP_EPSILON_PA.into();
        // Use desaturation rate (KDSAT = KSAT * SDR) when computing stop times
        let k = k_values_desat[tissue_idx];
        let m_value = if is_last_stop || stop_idx == 0 {
            surface_targets[tissue_idx]
        } else {
            m_values[stop_idx - 1].max_saturation[tissue_idx]
        };
        let p_amb: P = stop_depth_pa.into();
        let target_m = super::update::allowed_with_gf(p_amb, m_value, gf);

        // TODO: Plus or minus epsilon
        if p_tissue <= target_m {
            // Tissue safe
            continue;
        }

        let t_tissue = if p_tissue <= crossover_pressure {
            // Below the crossover pressure: exponential washout to the stop limit.
            -((target_m - p_inspired) / (p_tissue - p_inspired)).ln() / k
        } else if target_m >= crossover_pressure {
            // The ceiling is still in the linear region, so stop once the tissue
            // reaches the current stop limit.
            (p_tissue - target_m) / ((p_tissue - p_inspired) * k)
        } else {
            // Above the crossover pressure and the stop limit lies below it:
            // linear first, then exponential to the stop limit.
            let linear_rate = (p_tissue - p_inspired) * k;
            let t_linear = (p_tissue - crossover_pressure) / linear_rate;
            let t_exp = -((target_m - p_inspired) / (crossover_pressure - p_inspired)).ln() / k;
            t_linear + t_exp
        };

        t_stop_mins = max(t_stop_mins, t_tissue);
    }
    Duration::from_secs_f32(t_stop_mins * 60.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        deco_algorithm::{MVALUES, TISSUES, update_model_state},
        dive::{DiveMeasurement, DiveProfile},
        gas::{AIR, NX50, NX100, TMX18_45, TissuesLoading},
        loadings_from_dive_profile,
        pressure_unit::{Pa, Pressure, msw},
        update::first_stop_depth_with_gf,
    };
    use core::f32::consts::LN_2;
    use std::println;
    #[test]
    fn compute_first_stops_from_realistic_profile() {
        let gases = [TMX18_45, NX50, NX100];

        let descent_rate_m_per_min = 20.0_f32;
        let bottom_time_min = 35.0_f32;
        // let ascent_rate_deep_m_per_min = 9.0_f32;
        // let ascent_rate_shallow_m_per_min = 3.0_f32;
        let bottom = 55.0_f32;
        // let shallow_transition = 21.0_f32;

        let descent_time_s = bottom / descent_rate_m_per_min * 60.0;
        let descent_ms = (descent_time_s * 1000.0) as usize;
        let bottom_ms = (bottom_time_min * 60.0 * 1000.0) as usize;
        // let ascent_deep_time_s = (bottom - shallow_transition) / ascent_rate_deep_m_per_min * 60.0;
        // let ascent_deep_ms = (ascent_deep_time_s * 1000.0) as usize;
        // let ascent_shallow_time_s = shallow_transition / ascent_rate_shallow_m_per_min * 60.0;
        // let ascent_shallow_ms = (ascent_shallow_time_s * 1000.0) as usize;

        let measurements = [
            DiveMeasurement {
                time_ms: 0,
                depth: msw::new(0.0).to_pa(),
                gas: 0,
            },
            DiveMeasurement {
                time_ms: descent_ms,
                depth: msw::new(bottom).to_pa(),
                gas: 0,
            },
            DiveMeasurement {
                time_ms: descent_ms + bottom_ms,
                depth: msw::new(bottom).to_pa(),
                gas: 0,
            },
            // DiveMeasurement {
            //     time_ms: descent_ms + bottom_ms + ascent_deep_ms,
            //     depth: msw::new(shallow_transition).to_pa(),
            //     gas: 1,
            // },
            // DiveMeasurement {
            //     time_ms: descent_ms + bottom_ms + ascent_deep_ms + ascent_shallow_ms,
            //     depth: msw::new(0.0).to_pa(),
            //     gas: 2,
            // },
        ];

        let profile: DiveProfile<Pa, f32, 3, 3> = DiveProfile {
            dive_id: 1,
            max_depth: msw::new(bottom).to_pa(),
            gases,
            measurements,
        };

        let mut current_loading =
            loadings_from_dive_profile(&TISSUES, &profile, &MVALUES, msw::new(0.0).to_pa());

        for _step in 0..3 {
            let Some(stop_depth) = first_stop_depth_with_gf(&current_loading, &MVALUES, 1.0) else {
                panic!("Expected at least five stops, got: {:?}", _step);
            };
            let stop_duration = compute_stop_time_lin_exp(
                &current_loading,
                &TISSUES,
                &gases[0],
                &MVALUES,
                stop_depth,
                1.0,
                msw::new(3.0),
            );
            println!("Stop {:?}: {:?}", stop_depth, stop_duration);
            assert!(!stop_duration.is_zero());
            update_model_state(
                &mut current_loading,
                &TISSUES,
                &MVALUES,
                &gases[0],
                stop_depth.to_pa().into(),
                &stop_duration,
            );
        }
    }

    #[test]
    fn update_model_state_lin_exp_uses_exponential_update_while_ongassing() {
        let mut loading: TissuesLoading<{ TISSUES.len() }, Pa> = TissuesLoading {
            n2: TISSUES.map(|_| msw::new(0.0).to_pa()),
            he: TISSUES.map(|_| msw::new(0.0).to_pa()),
        };

        update_model_state_lin_exp(
            &mut loading,
            &TISSUES,
            &MVALUES,
            &TMX18_45,
            msw::new(30.0).to_pa(),
            &Duration::from_secs(60),
        );

        assert!(loading.n2.iter().all(|value| value.to_f32().is_finite()));
        assert!(loading.he.iter().all(|value| value.to_f32().is_finite()));
        assert!(loading.n2.iter().all(|value| value.to_f32() > 0.0));
        assert!(loading.he.iter().all(|value| value.to_f32() > 0.0));
    }

    #[test]
    fn update_model_state_lin_exp_uses_kdsat_for_desaturation_branch() {
        let mut loading: TissuesLoading<{ TISSUES.len() }, Pa> = TissuesLoading {
            n2: TISSUES.map(|_| msw::new(0.0).to_pa()),
            he: TISSUES.map(|_| msw::new(0.0).to_pa()),
        };

        let depth = msw::new(10.0).to_pa();
        let tissue_idx = 2; // Thalmann compartment with half-time=20, SDR=0.67
        let p_inspired = AIR.pn2(depth);
        let crossover_pressure = thalmann_crossover_pressure(depth);

        // Force the linear-first desaturation branch:
        // p_old > crossover_pressure > p_inspired.
        let p_old = crossover_pressure + msw::new(5.0).to_pa();
        loading.n2[tissue_idx] = p_old;

        let dt = Duration::from_secs(60);
        update_model_state_lin_exp(&mut loading, &TISSUES, &MVALUES, &AIR, depth, &dt);

        let dt_min = 1.0;
        let k_desat = (LN_2 / TISSUES[tissue_idx].half_time) * TISSUES[tissue_idx].sdr;
        let linear_rate = (p_old - p_inspired) * k_desat;
        let t_linear = (p_old - crossover_pressure) / linear_rate;
        let expected = if dt_min <= t_linear {
            p_old - linear_rate * dt_min
        } else {
            let t_exp = dt_min - t_linear;
            exp_pressure(p_inspired, crossover_pressure, k_desat, t_exp)
        };
        let actual = loading.n2[tissue_idx];
        assert!((actual.to_f32() - expected.to_f32()).abs() < 1e-3);
    }

    #[test]
    fn update_model_state_lin_exp_uses_ksat_for_ongassing() {
        let mut loading: TissuesLoading<{ TISSUES.len() }, Pa> = TissuesLoading {
            n2: TISSUES.map(|_| msw::new(0.0).to_pa()),
            he: TISSUES.map(|_| msw::new(0.0).to_pa()),
        };

        let depth = msw::new(30.0).to_pa();
        let tissue_idx = 2; // half-time=20, SDR=0.67
        let p_inspired = AIR.pn2(depth);
        let p_old = msw::new(0.0).to_pa();
        loading.n2[tissue_idx] = p_old;

        let dt = Duration::from_secs(60);
        update_model_state_lin_exp(&mut loading, &TISSUES, &MVALUES, &AIR, depth, &dt);

        let dt_min = 1.0;
        let k_sat = LN_2 / TISSUES[tissue_idx].half_time;
        let expected = exp_pressure(p_inspired, p_old, k_sat, dt_min);
        let actual = loading.n2[tissue_idx];
        assert!((actual.to_f32() - expected.to_f32()).abs() < 1e-6);
    }
}

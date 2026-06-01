use std::convert::TryInto;
use std::time::Duration;

use stdc_diving_algorithms::deco_algorithm::update_model_state;
use stdc_diving_algorithms::deco_algorithm::{DecoSettings, MVALUES, TISSUES, calc_deco_schedule};
use stdc_diving_algorithms::dive::{DiveMeasurement, DiveProfile};
use stdc_diving_algorithms::gas::{AIR, TissuesLoading};
use stdc_diving_algorithms::gas::{GasDensitySettings, MAX_PO2_DECO, TMX10_80};
use stdc_diving_algorithms::pressure_unit::{Pa, Pressure, msw};
use stdc_diving_algorithms::setup::NUM_STOP_DEPTHS;

fn main() {
    const NUM_GASES: usize = 1;
    const SAMPLE_DT_S: usize = 10; // sample every 10s
    const TOTAL_S: usize = 600; // total duration: 5 min bottom + 5 min ascent
    const NUM_MEASUREMENTS: usize = TOTAL_S / SAMPLE_DT_S + 1; // inclusive

    let gases = [TMX10_80];

    // Build measurements: t=0 surface, t=300 bottom @100m, then ascent from 300..600
    let mut measurements_vec: Vec<DiveMeasurement<Pa>> = Vec::with_capacity(NUM_MEASUREMENTS);
    for i in 0..=TOTAL_S / SAMPLE_DT_S {
        let t_s = i * SAMPLE_DT_S;
        let time_ms = t_s * 1000;
        let depth_m = if t_s == 0 {
            0.0
        } else if t_s <= 300 {
            // treat as already at depth by t=300 (bottom time)
            100.0
        } else {
            // ascent: 20 m/min = 20/60 = 0.333333 m/s; elapsed since start of ascent = t_s - 300
            let ascent_elapsed_s = (t_s - 300) as f32;
            let ascent_rate_m_per_s = 20.0 / 60.0;
            let d = 100.0 - ascent_rate_m_per_s * ascent_elapsed_s;
            if d < 0.0 { 0.0 } else { d }
        };
        measurements_vec.push(DiveMeasurement {
            time_ms,
            depth: msw::new(depth_m).to_pa(),
            gas: 0,
        });
    }

    let measurements_array: [DiveMeasurement<Pa>; NUM_MEASUREMENTS] = measurements_vec
        .try_into()
        .expect("measurements length mismatch");

    let profile: DiveProfile<Pa, f32, NUM_GASES, NUM_MEASUREMENTS> = DiveProfile {
        dive_id: 1,
        max_depth: msw::new(100.0).to_pa(),
        gases,
        measurements: measurements_array,
    };

    // Build tissue loadings incrementally by applying each profile segment
    let mut loadings = TissuesLoading::new(msw::new(0.0).to_pa(), &AIR);
    for w in profile.measurements.windows(2) {
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
            &TISSUES,
            &MVALUES,
            &profile.gases[*gas],
            midpoint,
            &delta_time,
        );
    }

    let settings = DecoSettings {
        gas_density_settings: GasDensitySettings::Ignore,
        max_deco_po2: MAX_PO2_DECO.to_pa(),
        surface_pressure: msw::new(0.0).to_pa(),
        ignore_icd: false,
        gf_low: 0.50,
        gf_high: 0.85,
        last_deco_stop: msw::new(3.0),
    };

    let schedule =
        calc_deco_schedule::<{ NUM_STOP_DEPTHS - 1 }, 1>(&loadings, &gases, &[true], &settings)
            .expect("schedule");

    println!("first_stop={:?}", schedule.first_stop());
    for s in schedule.stops().iter() {
        if !s.duration().is_zero() {
            println!(
                "stop depth={:.1}m duration={:.3}s",
                s.depth().to_msw().to_f32(),
                s.duration().as_secs_f32()
            );
        }
    }
}

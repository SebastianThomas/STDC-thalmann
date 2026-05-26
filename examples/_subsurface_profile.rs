use std::convert::TryInto;
use std::time::Duration;

use stdc_diving_algorithms::deco_algorithm::update_model_state;
use stdc_diving_algorithms::deco_algorithm::{DecoSettings, MVALUES, TISSUES, calc_deco_schedule};
use stdc_diving_algorithms::dive::{DiveMeasurement, DiveProfile};
use stdc_diving_algorithms::gas::{
    AIR, GasDensitySettings, MAX_PO2_DECO, NX50, TMX10_80, TissuesLoading,
};
use stdc_diving_algorithms::pressure_unit::{Pa, Pressure, msw};
use stdc_diving_algorithms::setup::NUM_STOP_DEPTHS;
use stdc_diving_algorithms::tissue_mvalues_with_gf;

const LOG_TISSUES: [usize; 4] = [0, 1, 2, 3];
const SUBSTEP_MS: usize = 10_000;

fn push_segment(
    measurements: &mut Vec<DiveMeasurement<Pa>>,
    time_s: usize,
    depth_m: f32,
    gas: usize,
) {
    measurements.push(DiveMeasurement {
        time_ms: time_s * 1000,
        depth: msw::new(depth_m).to_pa(),
        gas,
    });
}

fn log_selected_tissues(
    label: &str,
    depth: msw,
    loadings: &TissuesLoading<{ stdc_diving_algorithms::setup::NUM_TISSUES }, Pa>,
    mvalues: &stdc_diving_algorithms::deco_algorithm::MValues<Pa>,
    gf: f32,
) {
    println!("{} at {:.1}m:", label, depth.to_msw().to_f32());
    for &idx in &LOG_TISSUES {
        let n2 = loadings.n2[idx].to_pa().to_f32();
        let he = loadings.he[idx].to_pa().to_f32();
        let (abs_mvalue, gf_mvalue) = tissue_mvalues_with_gf(loadings, mvalues, depth, idx, gf);
        println!(
            "  tissue {}: n2={:.3} Pa he={:.3} Pa total={:.3} Pa abs_mvalue={:.3} Pa gf_mvalue={:.3} Pa",
            idx,
            n2,
            he,
            n2 + he,
            abs_mvalue.to_pa().to_f32(),
            gf_mvalue.to_pa().to_f32()
        );
    }
}

fn update_segment_substepped(
    loadings: &mut TissuesLoading<{ stdc_diving_algorithms::setup::NUM_TISSUES }, Pa>,
    mvalues: &stdc_diving_algorithms::deco_algorithm::MValues<Pa>,
    gas: &stdc_diving_algorithms::gas::GasMix<f32>,
    prev_depth: Pa,
    next_depth: Pa,
    delta_time: Duration,
) {
    let total_ms = delta_time.as_millis() as usize;
    if total_ms == 0 {
        return;
    }

    let mut elapsed_ms = 0usize;
    while elapsed_ms < total_ms {
        let chunk_ms = (total_ms - elapsed_ms).min(SUBSTEP_MS);
        let start_ratio = elapsed_ms as f32 / total_ms as f32;
        let end_ratio = (elapsed_ms + chunk_ms) as f32 / total_ms as f32;
        let mid_ratio = (start_ratio + end_ratio) * 0.5;
        let interpolated_depth = prev_depth + (next_depth - prev_depth) * mid_ratio;
        update_model_state(
            loadings,
            &TISSUES,
            mvalues,
            gas,
            interpolated_depth,
            &Duration::from_millis(chunk_ms as u64),
        );
        elapsed_ms += chunk_ms;
    }
}

fn main() {
    let gases = [TMX10_80, AIR, NX50];
    let enabled = [true, true, true];

    // Approximate the Subsurface plan from the prompt:
    // 0 -> 100m in 1 min on 10/80, hold 1 min, ascend to 66m in 4 min,
    // hold 1 min on air, ascend to 21m in 5 min, hold 1 min on EAN50,
    // ascend to 6m in 2 min, hold 1 min, jump to 3m, hold 7 min, surface.
    let mut measurements = Vec::new();
    push_segment(&mut measurements, 0, 0.0, 0);
    push_segment(&mut measurements, 60, 100.0, 0);
    push_segment(&mut measurements, 120, 100.0, 0);
    push_segment(&mut measurements, 360, 66.0, 0);
    push_segment(&mut measurements, 420, 66.0, 1);
    push_segment(&mut measurements, 720, 21.0, 1);
    push_segment(&mut measurements, 780, 21.0, 2);
    push_segment(&mut measurements, 900, 6.0, 2);
    push_segment(&mut measurements, 960, 6.0, 2);
    push_segment(&mut measurements, 960, 3.0, 2);
    push_segment(&mut measurements, 1380, 3.0, 2);
    push_segment(&mut measurements, 1380 + 1, 0.0, 2);

    let measurements_array: [DiveMeasurement<Pa>; 12] = measurements.try_into().expect("len");

    let profile: DiveProfile<Pa, f32, 3, 12> = DiveProfile {
        dive_id: 42,
        max_depth: msw::new(100.0).to_pa(),
        gases,
        measurements: measurements_array,
    };

    let mut loadings = TissuesLoading::new(msw::new(0.0).to_pa(), &TMX10_80);
    log_selected_tissues("profile start", msw::new(0.0), &loadings, &MVALUES, 1.0);
    for window in profile.measurements.windows(2) {
        let prev = window[0];
        let next = window[1];
        let delta_time = Duration::from_millis((next.time_ms - prev.time_ms) as u64);
        let midpoint = (prev.depth + next.depth) / 2.0;
        let gas_idx = next.gas;
        let gas = &profile.gases[gas_idx];
        println!(
            "segment time={}s depth={:.1}m gas={} fo2={:.3} fhe={:.3}",
            next.time_ms / 1000,
            midpoint.to_msw().to_f32(),
            gas_idx,
            gas.fo2(),
            gas.fhe()
        );
        update_segment_substepped(
            &mut loadings,
            &MVALUES,
            gas,
            prev.depth,
            next.depth,
            delta_time,
        );
        log_selected_tissues("profile point", next.depth.to_msw(), &loadings, &MVALUES, 1.0);
    }

    // Debug: compute initial first stop using GFLow and print it before scheduling
    // Can't call internal helpers due to const-generic visibility; compute via schedule.

    let settings = DecoSettings {
        gas_density_settings: GasDensitySettings::Ignore,
        max_deco_po2: MAX_PO2_DECO.to_pa(),
        ignore_icd: true,
        gf_low: 0.50,
        gf_high: 0.85,
    };

    // Diagnostic: replicate first_stop_depth_with_gf logic here to see why no stop
    #[cfg(not(feature = "lin_exp"))]
    {
        use stdc_diving_algorithms::deco_algorithm::TISSUES as BUEHL_TISSUES;
        println!("Running standalone first-stop detection (Buehlmann)");
        let mut found: Option<msw> = None;
        for mvalues_at_depth in MVALUES.iter().rev() {
            let depth = mvalues_at_depth.depth;
            let mut any = false;
            // scan tissues but don't print each depth
            for i in 0..(stdc_diving_algorithms::setup::NUM_TISSUES) {
                let p_n2 = loadings.n2[i].to_pa();
                let p_he = loadings.he[i].to_pa();
                let total = p_n2 + p_he;
                if total.to_f32() <= 0.0 {
                    continue;
                }
                let a_n2 = BUEHL_TISSUES[i].n2.a.to_pa();
                let a_he = BUEHL_TISSUES[i].he.a.to_pa();
                let b_n2 = BUEHL_TISSUES[i].n2.b;
                let b_he = BUEHL_TISSUES[i].he.b;
                let n2_frac = p_n2 / total;
                let he_frac = p_he / total;
                let a_mix = a_n2 * n2_frac + a_he * he_frac;
                let b_mix = b_n2 * n2_frac + b_he * he_frac;
                let p_tol = a_mix + depth.to_pa() * b_mix;
                let p_amb = depth.to_pa();
                let allowed = p_amb + (p_tol - p_amb) * settings.gf_low;
                if total > allowed {
                    println!(
                        "first-stop candidate depth={:.1}m (GFLow={:.2})",
                        depth.to_msw().to_f32(),
                        settings.gf_low
                    );
                    log_selected_tissues("first-stop loading", depth, &loadings, &MVALUES, settings.gf_low);
                    any = true;
                    found = Some(depth);
                    break;
                }
            }
            if any {
                break;
            }
        }
        println!("diagnostic first stop: {:?}", found);
        if let Some(depth) = found {
            use stdc_diving_algorithms::deco_algorithm::compute_stop_time;
            use stdc_diving_algorithms::gas::best_available_mix;
            let mix = best_available_mix(
                settings.max_deco_po2,
                depth.to_pa().into(),
                &gases,
                &enabled,
                &loadings,
                settings.ignore_icd,
                &settings.gas_density_settings,
            );
            println!("best mix at {:?} => {:?}", depth, mix);
            if let Some((_idx, gas)) = mix {
                let dur =
                    compute_stop_time(&loadings, &TISSUES, gas, &MVALUES, depth, settings.gf_low);
                println!(
                    "computed stop duration at {:?} = {:.3}s",
                    depth,
                    dur.as_secs_f32()
                );
            }
            // done with diagnostics; continue to build schedule
        }
    }

    // Build schedule and then print per-tissue comparisons at initial first stop

    let schedule =
        calc_deco_schedule::<{ NUM_STOP_DEPTHS - 1 }, 3>(&loadings, &gases, &enabled, &settings)
            .expect("schedule");
    println!("first_stop={:?}", schedule.first_stop());
    if let Some(initial_first) = schedule.first_stop() {
        println!("initial_first_stop={:?}", initial_first);
        log_selected_tissues(
            "first-stop loading",
            initial_first.depth(),
            &loadings,
            &MVALUES,
            settings.gf_low,
        );
    } else {
        println!("no initial first stop with gf_low");
    }
    for s in schedule.stops().iter() {
        if !s.duration().is_zero() {
            println!(
                "stop depth={:.1}m duration={:.3}s gas={:?}",
                s.depth().to_msw().to_f32(),
                s.duration().as_secs_f32(),
                s.gas()
            );
        }
    }

    let mut deco_loadings = loadings.clone();
    for s in schedule.stops().iter() {
        if s.duration().is_zero() {
            continue;
        }
        if let Some(gas) = s.gas() {
            update_model_state(
                &mut deco_loadings,
                &TISSUES,
                &MVALUES,
                &gas,
                s.depth().to_pa(),
                &s.duration(),
            );
            log_selected_tissues("deco point", s.depth(), &deco_loadings, &MVALUES, settings.gf_low);
        }
    }

    // Local simulation of calc_deco_schedule_intern to compare results
    {
        use stdc_diving_algorithms::deco_algorithm::GradientFactors;
        use stdc_diving_algorithms::deco_algorithm::MVALUES as DMV;
        use stdc_diving_algorithms::depth_utils::get_depth;
        use stdc_diving_algorithms::dive::Stop;

        println!("Running local scheduler simulation");
        let gf_struct = GradientFactors {
            low: settings.gf_low,
            high: settings.gf_high,
        };
        let mut local_loading = loadings.clone();
        const NUM_STOPS_LOCAL: usize = NUM_STOP_DEPTHS - 1;
        let mut local_stops: [Stop; NUM_STOPS_LOCAL] =
            [Stop::new(msw::new(0.0), Duration::from_millis(0), None); NUM_STOPS_LOCAL];
        for i in 0..NUM_STOPS_LOCAL {
            local_stops[NUM_STOPS_LOCAL - 1 - i] =
                Stop::new(get_depth(i).to_msw(), Duration::from_millis(0), None);
        }
        // initial first
        fn first_stop_local(
            loading: &TissuesLoading<{ stdc_diving_algorithms::setup::NUM_TISSUES }, Pa>,
            gf: f32,
        ) -> Option<msw> {
            #[cfg(not(feature = "lin_exp"))]
            {
                for mv in DMV.iter().rev() {
                    for i in 0..(stdc_diving_algorithms::setup::NUM_TISSUES) {
                        let p_n2 = loading.n2[i].to_pa();
                        let p_he = loading.he[i].to_pa();
                        let total = p_n2 + p_he;
                        if total.to_f32() <= 0.0 {
                            continue;
                        }
                        let a_n2 = TISSUES[i].n2.a.to_pa();
                        let a_he = TISSUES[i].he.a.to_pa();
                        let b_n2 = TISSUES[i].n2.b;
                        let b_he = TISSUES[i].he.b;
                        let n2_frac = p_n2 / total;
                        let he_frac = p_he / total;
                        let a_mix = a_n2 * n2_frac + a_he * he_frac;
                        let b_mix = b_n2 * n2_frac + b_he * he_frac;
                        let p_tol = a_mix + mv.depth.to_pa() * b_mix;
                        let allowed = mv.depth.to_pa() + (p_tol - mv.depth.to_pa()) * gf;
                        if total > allowed {
                            return Some(mv.depth);
                        }
                    }
                }
                return None;
            }

            #[cfg(feature = "lin_exp")]
            {
                // Thalmann M-values provide a per-tissue max_saturation for each depth.
                for mv in DMV.iter().rev() {
                    for i in 0..(stdc_diving_algorithms::setup::NUM_TISSUES) {
                        let p_n2 = loading.n2[i].to_pa();
                        let p_he = loading.he[i].to_pa();
                        let total = p_n2 + p_he;
                        if total.to_f32() <= 0.0 {
                            continue;
                        }
                        let mval = mv.max_saturation[i];
                        let allowed = mv.depth.to_pa() + (mval - mv.depth.to_pa()) * gf;
                        if total > allowed {
                            return Some(mv.depth);
                        }
                    }
                }
                return None;
            }
        }
    }
}

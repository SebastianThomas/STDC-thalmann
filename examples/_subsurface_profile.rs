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
    // Plan based on user-provided schedule (times in seconds)
    push_segment(&mut measurements, 0 * 60, 0.0, 0); // 0:00 start surface
    push_segment(&mut measurements, 5 * 60, 100.0, 0); // 5:00 reach 100m on 10/80
    push_segment(&mut measurements, 12 * 60, 100.0, 0); // 12:00 hold at 100m
    push_segment(&mut measurements, 16 * 60, 66.0, 0); // 16:00 ascend to 66m
    push_segment(&mut measurements, 17 * 60, 66.0, 1); // 17:00 switch to air at 66m
    push_segment(&mut measurements, 21 * 60, 30.0, 1); // 21:00 descend/ascend to 30m on air
    push_segment(&mut measurements, 22 * 60, 30.0, 1); // 22:00 hold 30m
    push_segment(&mut measurements, 23 * 60, 24.0, 1); // 23:00 to 24m
    push_segment(&mut measurements, 25 * 60, 24.0, 1); // 25:00 hold 24m
    push_segment(&mut measurements, 25 * 60, 21.0, 1); // 25:00 move to 21m (same time)
    push_segment(&mut measurements, 26 * 60, 21.0, 2); // 26:00 switch to EAN50 at 21m
    push_segment(&mut measurements, 27 * 60, 18.0, 2); // 27:00 to 18m
    push_segment(&mut measurements, 29 * 60, 18.0, 2); // 29:00 hold 18m
    push_segment(&mut measurements, 29 * 60, 15.0, 2); // 29:00 move to 15m
    push_segment(&mut measurements, 33 * 60, 15.0, 2); // 33:00 hold 15m
    push_segment(&mut measurements, 33 * 60, 12.0, 2); // 33:00 move to 12m
    push_segment(&mut measurements, 38 * 60, 12.0, 2); // 38:00 hold 12m
    push_segment(&mut measurements, 38 * 60, 9.0, 2); // 38:00 move to 9m
    push_segment(&mut measurements, 45 * 60, 9.0, 2); // 45:00 hold 9m
    push_segment(&mut measurements, 45 * 60, 6.0, 2); // 45:00 move to 6m
    push_segment(&mut measurements, 57 * 60, 6.0, 2); // 57:00 hold 6m
    push_segment(&mut measurements, 57 * 60, 3.0, 2); // 57:00 move to 3m
    push_segment(&mut measurements, 80 * 60, 3.0, 2); // 80:00 hold 3m
    push_segment(&mut measurements, 80 * 60, 0.0, 2); // 80:00 surface

    // Ensure measurement timestamps are strictly increasing to avoid
    // zero-duration segments which produce `.5` midpoints (e.g. 22.5m).
    fn normalize_measurements(measurements: &mut Vec<DiveMeasurement<Pa>>) {
        let mut last_time: Option<usize> = None;
        for m in measurements.iter_mut() {
            if let Some(lt) = last_time {
                if m.time_ms <= lt {
                    m.time_ms = lt + 1000; // bump by 1s
                }
            }
            last_time = Some(m.time_ms);
        }
    }

    normalize_measurements(&mut measurements);

    let measurements_array: [DiveMeasurement<Pa>; 24] = measurements.try_into().expect("len");

    let profile: DiveProfile<Pa, f32, 3, 24> = DiveProfile {
        dive_id: 42,
        max_depth: msw::new(100.0).to_pa(),
        gases,
        measurements: measurements_array,
    };

    let mut loadings = TissuesLoading::new(msw::new(0.0).to_pa(), &TMX10_80);
    let settings = DecoSettings {
        gas_density_settings: GasDensitySettings::Ignore,
        max_deco_po2: MAX_PO2_DECO.to_pa(),
        surface_pressure: msw::new(0.0).to_pa(),
        ignore_icd: true,
        gf_low: 0.50,
        gf_high: 0.85,
        last_deco_stop: msw::new(6.0),
    };
    for window in profile.measurements.windows(2) {
        let prev = window[0];
        let next = window[1];
        let delta_time = Duration::from_millis((next.time_ms - prev.time_ms) as u64);
        let delta_ms = next.time_ms - prev.time_ms;
        let midpoint = if delta_ms <= 1000 {
            next.depth
        } else {
            (prev.depth + next.depth) / 2.0
        };
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

        // Before executing the next segment, check if following the profile
        // would ascend above the current deco ceiling (first-stop). We derive
        // the ceiling by building a short schedule from the current loading
        // and inspecting its first stop.
        if let Ok(schedule) =
            calc_deco_schedule::<{ NUM_STOP_DEPTHS - 1 }, 3>(&loadings, &gases, &enabled, &settings)
        {
            if let Some(first) = schedule.first_stop() {
                let first_stop = first.depth();
                if next.depth.to_msw().to_f32() < first_stop.to_msw().to_f32() {
                    println!(
                        "Enforcing full deco schedule starting at {:.1}m before continuing profile",
                        first_stop.to_msw().to_f32()
                    );
                    // Diagnostic: compute stop durations for 6m and 3m under
                    // both last_deco_stop = 6m and last_deco_stop = 3m using the
                    // current loading snapshot so we can compare behavior.
                    {
                        use stdc_diving_algorithms::deco_algorithm::compute_stop_time;
                        use stdc_diving_algorithms::deco_algorithm::TISSUES as BUEHL_TISSUES;
                        fn interp_gf(initial_first_stop: msw, stop_depth: msw, low: f32, high: f32) -> f32 {
                            if initial_first_stop.to_msw().to_f32() <= 0.0 {
                                high
                            } else {
                                let mut t = (initial_first_stop.to_msw().to_f32() - stop_depth.to_msw().to_f32()) / initial_first_stop.to_msw().to_f32();
                                if t < 0.0 { t = 0.0 } else if t > 1.0 { t = 1.0 }
                                low + (high - low) * t
                            }
                        }

                        let initial = first_stop;
                        let gf6 = interp_gf(initial, msw::new(6.0), settings.gf_low, settings.gf_high);
                        let gf3 = interp_gf(initial, msw::new(3.0), settings.gf_low, settings.gf_high);
                        // best mixes for depths
                        use stdc_diving_algorithms::gas::best_available_mix;
                        let mix6 = best_available_mix(settings.max_deco_po2, msw::new(6.0).to_pa().into(), &gases, &enabled, &loadings, settings.ignore_icd, &settings.gas_density_settings);
                        let mix3 = best_available_mix(settings.max_deco_po2, msw::new(3.0).to_pa().into(), &gases, &enabled, &loadings, settings.ignore_icd, &settings.gas_density_settings);
                        println!("Diagnostic GF: gf6={:.3} gf3={:.3}", gf6, gf3);
                            if let Some((_i, g6)) = mix6 {
                            let d6_as_final = compute_stop_time(&loadings, &BUEHL_TISSUES, g6, &MVALUES, msw::new(6.0), gf6, msw::new(0.0).to_pa(), msw::new(6.0));
                            let d6_with_3floor = compute_stop_time(&loadings, &BUEHL_TISSUES, g6, &MVALUES, msw::new(6.0), gf6, msw::new(0.0).to_pa(), msw::new(3.0));
                            println!("compute_stop_time 6m final(6m floor) = {:.1}s", d6_as_final.as_secs_f32());
                            println!("compute_stop_time 6m non-final(3m floor) = {:.1}s", d6_with_3floor.as_secs_f32());
                        }
                        if let Some((_i, g3)) = mix3 {
                            let d3_final = compute_stop_time(&loadings, &BUEHL_TISSUES, g3, &MVALUES, msw::new(3.0), gf3, msw::new(0.0).to_pa(), msw::new(3.0));
                            println!("compute_stop_time 3m final(3m floor) = {:.1}s", d3_final.as_secs_f32());
                        }
                    }
                    // Execute the computed schedule in order (this updates loadings)
                    for s in schedule.stops().iter() {
                        if s.duration().is_zero() {
                            continue;
                        }
                        let stop_depth = s.depth();
                        if let Some(gas) = s.gas() {
                            println!(
                                "  stop {:.1}m for {:.1}s gas fo2={:.3}",
                                stop_depth.to_msw().to_f32(),
                                s.duration().as_secs_f32(),
                                gas.fo2()
                            );
                            update_model_state(
                                &mut loadings,
                                &TISSUES,
                                &MVALUES,
                                &gas,
                                stop_depth.to_pa().into(),
                                &s.duration(),
                            );
                        } else {
                            // No gas provided; just simulate time at depth
                            println!(
                                "  stop {:.1}m for {:.1}s (no gas)",
                                stop_depth.to_msw().to_f32(),
                                s.duration().as_secs_f32()
                            );
                            update_model_state(
                                &mut loadings,
                                &TISSUES,
                                &MVALUES,
                                &gases[0],
                                stop_depth.to_pa().into(),
                                &s.duration(),
                            );
                        }
                    }
                }
            }
        }

        update_segment_substepped(
            &mut loadings,
            &MVALUES,
            gas,
            prev.depth,
            next.depth,
            delta_time,
        );
    }

    // Debug: compute initial first stop using GFLow and print it before scheduling
    // Can't call internal helpers due to const-generic visibility; compute via schedule.

    // let settings = DecoSettings {
    //     gas_density_settings: GasDensitySettings::Ignore,
    //     max_deco_po2: MAX_PO2_DECO.to_pa(),
    //     surface_pressure: msw::new(0.0).to_pa(),
    //     ignore_icd: true,
    //     gf_low: 0.50,
    //     gf_high: 0.85,
    //     last_deco_stop: msw::new(6.0),
    // };

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
                    compute_stop_time(
                        &loadings,
                        &TISSUES,
                        gas,
                        &MVALUES,
                        depth,
                        settings.gf_low,
                        settings.last_deco_stop,
                    );
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
        }
    }

    // Local simulation of calc_deco_schedule_intern to compare results
    {
        use stdc_diving_algorithms::depth_utils::get_depth;
        use stdc_diving_algorithms::dive::Stop;

        println!("Running local scheduler simulation");
        const NUM_STOPS_LOCAL: usize = NUM_STOP_DEPTHS - 1;
        let mut local_stops: [Stop; NUM_STOPS_LOCAL] =
            [Stop::new(msw::new(0.0), Duration::from_millis(0), None); NUM_STOPS_LOCAL];
        for i in 0..NUM_STOPS_LOCAL {
            local_stops[NUM_STOPS_LOCAL - 1 - i] =
                Stop::new(get_depth(i).to_msw(), Duration::from_millis(0), None);
        }
    }
}

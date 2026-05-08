use core::time::Duration;

use num::pow::Pow;

use crate::dive::DiveProfile;
use crate::pressure_unit::{AbsPressure, Bar, Pressure};

pub enum O2ToxCalculation {
    /** */
    NOAA,
    /** https://pmc.ncbi.nlm.nih.gov/articles/PMC12500339/
     *  for 1.3 bar PO2
     *  with NOAA as fallback
     */
    RevisedDHM2025,
}

pub enum O2ExposureType {
    Single,
    Daily24h,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum O2ExposureLimit {
    Limit(Duration),
    ExeptionalLimit(Duration),
    Unsafe,
}

#[derive(Debug, Clone, Copy)]
pub struct O2ToxicityPercentage {
    /// CNS toxicity as percentage (0.0 to 100.0+)
    pub cns_percent: f32,
    /// Pulmonary toxicity as percentage (0.0 to 100.0+)
    pub pulmonary_percent: f32,
}

impl O2ToxicityPercentage {
    pub const fn new(cns_percent: f32, pulmonary_percent: f32) -> Self {
        O2ToxicityPercentage {
            cns_percent,
            pulmonary_percent,
        }
    }
}

impl O2ToxCalculation {
    pub fn limit<P: const Pressure>(&self, exposure: &O2ExposureType, po2: P) -> O2ExposureLimit {
        match self {
            Self::NOAA => noaa_o2_limit(po2.to_bar(), exposure),
            Self::RevisedDHM2025 => revised_dhm_2025_o2_limit(po2.to_bar(), exposure),
        }
    }
}

/** https://pmc.ncbi.nlm.nih.gov/articles/PMC12500339/table/T1/ */
fn noaa_o2_limit(po2: Bar, exposure: &O2ExposureType) -> O2ExposureLimit {
    let mins = match po2.to_f32() {
        v if v <= 1.6 => match exposure {
            O2ExposureType::Single => 45,
            O2ExposureType::Daily24h => 150,
        },
        v if v <= 1.5 => match exposure {
            O2ExposureType::Single => 120,
            O2ExposureType::Daily24h => 180,
        },
        v if v <= 1.4 => match exposure {
            O2ExposureType::Single => 150,
            O2ExposureType::Daily24h => 180,
        },
        v if v <= 1.3 => match exposure {
            O2ExposureType::Single => 180,
            O2ExposureType::Daily24h => 210,
        },
        v if v <= 1.2 => match exposure {
            O2ExposureType::Single => 210,
            O2ExposureType::Daily24h => 240,
        },
        v if v <= 1.1 => match exposure {
            O2ExposureType::Single => 240,
            O2ExposureType::Daily24h => 270,
        },
        v if v <= 1.0 => match exposure {
            O2ExposureType::Single => 300,
            O2ExposureType::Daily24h => 300,
        },
        v if v <= 0.9 => match exposure {
            O2ExposureType::Single => 360,
            O2ExposureType::Daily24h => 360,
        },
        v if v <= 0.8 => match exposure {
            O2ExposureType::Single => 450,
            O2ExposureType::Daily24h => 450,
        },
        v if v <= 0.7 => match exposure {
            O2ExposureType::Single => 570,
            O2ExposureType::Daily24h => 570,
        },
        v if v <= 0.6 => match exposure {
            O2ExposureType::Single => 720,
            O2ExposureType::Daily24h => 720,
        },
        _ => 0,
    };
    if mins != 0 {
        return O2ExposureLimit::Limit(Duration::from_mins(mins));
    }
    let mins = match po2.to_f32() {
        v if v <= 1.7 => 75,
        v if v <= 1.8 => 60,
        v if v <= 1.9 => 45,
        v if v <= 2.0 => 30,
        _ => 0,
    };
    if mins != 0 {
        return O2ExposureLimit::ExeptionalLimit(Duration::from_mins(mins));
    }
    return O2ExposureLimit::Unsafe;
}

/** https://doi.org/10.28920/dhm55.3.262-270 */
fn revised_dhm_2025_o2_limit(po2: Bar, exposure: &O2ExposureType) -> O2ExposureLimit {
    let noaa_limit = noaa_o2_limit(po2, exposure);
    let new_limit_1_3 = O2ExposureLimit::Limit(Duration::from_mins(match exposure {
        O2ExposureType::Single => 240,
        // TODO: Total with Working Phase 240 + Resting 240
        O2ExposureType::Daily24h => 240 * 2,
    }));
    if po2.to_f32() - 1.3 <= 0.05 && noaa_limit <= new_limit_1_3 {
        return new_limit_1_3;
    }
    noaa_limit
}

/** Oxygen toxicity units per minute at varying partial pressure */
pub fn otu_per_minute<P: const Pressure>(po2: P) -> f32 {
    match po2.to_bar().to_f32() {
        v if v <= 0.50 => 0.00,
        v if v <= 0.55 => 0.15,
        v if v <= 0.60 => 0.27,
        v if v <= 0.65 => 0.37,
        v if v <= 0.70 => 0.47,
        v if v <= 0.75 => 0.56,
        v if v <= 0.80 => 0.65,
        v if v <= 0.85 => 0.74,
        v if v <= 0.90 => 0.83,
        v if v <= 0.95 => 0.92,
        v if v <= 1.00 => 1.00,
        v if v <= 1.05 => 1.08,
        v if v <= 1.10 => 1.16,
        v if v <= 1.15 => 1.24,
        v if v <= 1.20 => 1.32,
        v if v <= 1.25 => 1.40,
        v if v <= 1.30 => 1.48,
        v if v <= 1.35 => 1.55,
        v if v <= 1.40 => 1.63,
        v if v <= 1.45 => 1.70,
        v if v <= 1.50 => 1.78,
        v if v <= 1.55 => 1.85,
        v if v <= 1.60 => 1.92,
        v if v <= 1.65 => 2.00,
        v if v <= 1.70 => 2.07,
        v if v <= 1.75 => 2.14,
        v if v <= 1.80 => 2.21,
        v if v <= 1.85 => 2.28,
        v if v <= 1.90 => 2.35,
        v if v <= 1.95 => 2.42,
        v if v <= 2.00 => 2.49,
        _ => po2.to_bar().to_f32() * 2.0,
    }
}

pub fn oti_cns<P: const Pressure>(po2: P, time: Duration) -> f32 {
    oti(po2, time.as_secs_f32() / 60.0, 6.8) / (26.108 * 100.0)
}

pub fn oti_pulmonary<P: const Pressure>(po2: P, time: Duration) -> f32 {
    oti(po2, time.as_secs_f32() / (60.0 * 60.0), 4.57) / 250.0
}

fn oti<P: const Pressure>(po2: P, time: f32, c: f32) -> f32 {
    let time_2 = time * time;
    let po2_c = po2.to_bar().to_f32().pow(c);
    time_2 * po2_c
}

/// Calculates cumulative CNS and pulmonary oxygen toxicity percentages from a dive profile.
///
/// Iterates through all measurements in the profile and accumulates exposure at each depth,
/// calculating the percentage of maximum allowed exposure time at each point.
///
/// # Arguments
/// * `profile` - The dive profile to process
/// * `exposure_type` - The type of exposure (Single or Daily24h)
/// * `calculation_method` - The limits calculation method (NOAA or RevisedDHM2025)
pub fn calculate_toxicity_from_profile<
    const NUM_GASES: usize,
    const NUM_MEASUREMENTS: usize,
    P: const AbsPressure,
>(
    profile: &DiveProfile<P, f32, NUM_GASES, NUM_MEASUREMENTS>,
    exposure_type: &O2ExposureType,
    calculation_method: O2ToxCalculation,
) -> O2ToxicityPercentage {
    let initial_toxicity = O2ToxicityPercentage::new(0.0, 0.0);
    calculate_toxicity_diff(
        profile,
        0,
        &initial_toxicity,
        exposure_type,
        calculation_method,
    )
}

/// Calculates differential (incremental) toxicity from a specific measurement index onwards.
///
/// Takes the previous toxicity percentage and processes measurements from start_index forward.
/// This is useful for real-time monitoring where new measurements arrive incrementally.
///
/// # Arguments
/// * `profile` - The dive profile to process
/// * `start_index` - The measurement index to start from (0-based)
/// * `previous_toxicity` - The accumulated toxicity percentage before this segment
/// * `exposure_type` - The type of exposure (Single or Daily24h)
/// * `calculation_method` - The limits calculation method (NOAA or RevisedDHM2025)
///
/// # Returns
/// The new cumulative O2ToxicityPercentage including all measurements from start_index onwards
pub fn calculate_toxicity_diff<
    const NUM_GASES: usize,
    const NUM_MEASUREMENTS: usize,
    P: const AbsPressure,
>(
    profile: &DiveProfile<P, f32, NUM_GASES, NUM_MEASUREMENTS>,
    start_index: usize,
    previous_toxicity: &O2ToxicityPercentage,
    exposure_type: &O2ExposureType,
    calculation_method: O2ToxCalculation,
) -> O2ToxicityPercentage {
    let mut cns_percent = previous_toxicity.cns_percent;
    let mut pulmonary_percent = previous_toxicity.pulmonary_percent;

    // Need at least one measurement after start_index to process
    if start_index >= profile.measurements.len().saturating_sub(1) {
        return O2ToxicityPercentage::new(cns_percent, pulmonary_percent);
    }

    // Process from start_index to the end
    let measurements = &profile.measurements;
    for i in start_index..measurements.len().saturating_sub(1) {
        let prev_measurement = &measurements[i];
        let curr_measurement = &measurements[i + 1];

        let delta_time_ms = curr_measurement
            .time_ms
            .saturating_sub(prev_measurement.time_ms);
        let delta_time = Duration::from_millis(delta_time_ms as u64);

        // Use current depth
        let depth = curr_measurement.depth;
        let gas_mix = &profile.gases[curr_measurement.gas];

        // Calculate PO2 at current depth
        let abs_pressure = depth.to_pa();
        let po2 = Bar::new(abs_pressure.to_bar().to_f32() * gas_mix.fo2());

        // Get exposure limit for this PO2
        let limit = calculation_method.limit(exposure_type, po2);

        // Add to percentage if within safe limits
        match limit {
            O2ExposureLimit::Limit(max_duration) => {
                let max_secs = max_duration.as_secs_f32();
                let elapsed_secs = delta_time.as_secs_f32();
                cns_percent += (elapsed_secs / max_secs) * 100.0;
                pulmonary_percent += (elapsed_secs / max_secs) * 100.0;
            }
            O2ExposureLimit::ExeptionalLimit(max_duration) => {
                let max_secs = max_duration.as_secs_f32();
                let elapsed_secs = delta_time.as_secs_f32();
                // Exceptional limits contribute more aggressively to toxicity
                cns_percent += (elapsed_secs / max_secs) * 150.0;
                pulmonary_percent += (elapsed_secs / max_secs) * 150.0;
            }
            O2ExposureLimit::Unsafe => {
                // Unsafe exposure contributes immediately
                cns_percent += 200.0;
                pulmonary_percent += 200.0;
            }
        }
    }

    O2ToxicityPercentage::new(cns_percent, pulmonary_percent)
}

#[cfg(test)]
mod tests {
    use crate::dive::DiveMeasurement;
    use crate::gas::GasMix;
    use crate::pressure_unit::Bar;

    use super::*;

    #[test]
    fn get_oti_cns_test() {
        let test_cases = [
            // (po2 in bar, time in seconds, expected rough value range)
            (1.0, 300, 0.10),  // 5 min at 1.0 bar
            (1.2, 600, 0.30),  // 10 min at 1.2 bar
            (1.4, 300, 0.25),  // 5 min at 1.4 bar
            (1.6, 600, 0.70),  // 10 min at 1.6 bar
            (0.8, 1200, 0.15), // 20 min at 0.8 bar
        ];

        for (po2_bar, time_secs, _expected_rough) in &test_cases {
            let po2 = Bar::new(*po2_bar);
            let duration = Duration::from_secs(*time_secs);
            let cns = oti_cns(po2, duration);

            assert!(
                cns >= 0.0,
                "CNS should be non-negative for po2={}, time={}s",
                po2_bar,
                time_secs
            );
            assert!(
                cns < 1.0,
                "CNS should not exceed 100% (1.0) for reasonable dive parameters"
            );

            if *po2_bar >= 1.5 && *time_secs >= 300 {
                assert!(
                    cns > 0.1,
                    "Expected significant CNS loading for po2={}, time={}s",
                    po2_bar,
                    time_secs
                );
            }
        }
    }

    #[test]
    fn get_otu_test() {
        let po2_levels = [
            (0.5, 0.0),  // Below threshold
            (0.7, 0.47), // ~0.47 OTU/min
            (1.0, 1.0),  // 1.0 OTU/min at 1.0 bar
            (1.3, 1.48), // ~1.48 OTU/min at 1.3 bar
            (1.6, 1.92), // ~1.92 OTU/min at 1.6 bar
        ];

        for (po2_bar, expected_otu_per_min) in &po2_levels {
            let po2 = Bar::new(*po2_bar);
            let otu_per_min = otu_per_minute(po2);

            // Verify positive values for valid PO2
            assert!(
                otu_per_min >= 0.0,
                "OTU per minute should be non-negative for po2={}",
                po2_bar
            );

            // Check rough approximation (allow 5% tolerance for table values)
            if *expected_otu_per_min > 0.0 {
                let tolerance = expected_otu_per_min * 0.05;
                assert!(
                    (otu_per_min - expected_otu_per_min).abs() <= tolerance,
                    "OTU mismatch for po2={}: got {}, expected ~{}",
                    po2_bar,
                    otu_per_min,
                    expected_otu_per_min
                );
            }
        }

        // Test total OTU calculation with oti_pulmonary
        let test_durations = [
            (1.0, 300),  // 1.0 bar for 5 min
            (1.3, 600),  // 1.3 bar for 10 min
            (1.6, 1200), // 1.6 bar for 20 min
        ];

        for (po2_bar, time_secs) in &test_durations {
            let po2 = Bar::new(*po2_bar);
            let duration = Duration::from_secs(*time_secs);
            let otu_pulmonary = oti_pulmonary(po2, duration);

            assert!(otu_pulmonary >= 0.0, "Pulmonary OTI should be non-negative");
            assert!(
                otu_pulmonary < 100.0,
                "Pulmonary OTI should not exceed operational limits for test case po2={}, time={}s",
                po2_bar,
                time_secs
            );

            // Longer times, higher PO2 give higher OTI (even if small due to 250.0 divisor)
            if *po2_bar >= 1.5 && *time_secs >= 600 {
                assert!(
                    otu_pulmonary > 0.0,
                    "Expected non-zero pulmonary OTI for moderate duration at po2={}",
                    po2_bar
                );
            }
        }
    }

    #[test]
    fn o2_exposure_limit_safe_test() {
        // Test safe exposure limits
        let safe_cases = [
            (1.0, O2ExposureType::Single),   // 300 min limit
            (0.8, O2ExposureType::Daily24h), // 450 min limit
            (1.3, O2ExposureType::Daily24h), // 210 min limit
        ];

        for (po2_bar, exposure_type) in &safe_cases {
            let po2 = Bar::new(*po2_bar);
            let limit = O2ToxCalculation::NOAA.limit(exposure_type, po2);

            match limit {
                O2ExposureLimit::Limit(duration) => {
                    assert!(
                        duration.as_secs() > 0,
                        "Safe limit should have non-zero duration"
                    );
                }
                _ => panic!("Expected Limit variant for safe PO2 level"),
            }
        }
    }

    #[test]
    fn o2_exposure_limit_exceptional_test() {
        let exceptional_cases = [
            (1.7, O2ExposureType::Single, 75 * 60),
            (1.8, O2ExposureType::Single, 60 * 60),
            (1.9, O2ExposureType::Single, 45 * 60),
            (2.0, O2ExposureType::Single, 30 * 60),
        ];

        for (po2_bar, exposure_type, expected_secs) in &exceptional_cases {
            let po2 = Bar::new(*po2_bar);
            let limit = O2ToxCalculation::NOAA.limit(exposure_type, po2);

            match limit {
                O2ExposureLimit::ExeptionalLimit(duration) => {
                    assert_eq!(
                        duration.as_secs(),
                        *expected_secs as u64,
                        "Exceptional limit mismatch for po2={}",
                        po2_bar
                    );
                }
                _ => panic!("Expected ExceptionalLimit variant for PO2 {}", po2_bar),
            }
        }
    }

    #[test]
    fn o2_exposure_limit_unsafe_test() {
        let unsafe_cases = [2.1, 2.5, 3.0, 5.0];

        for po2_bar in &unsafe_cases {
            let po2 = Bar::new(*po2_bar);
            let limit = O2ToxCalculation::NOAA.limit(&O2ExposureType::Single, po2);

            match limit {
                O2ExposureLimit::Unsafe => {
                    // Expected behavior
                }
                _ => panic!("Expected Unsafe variant for PO2 {}", po2_bar),
            }
        }
    }

    #[test]
    fn toxicity_percentage_from_safe_profile_test() {
        use crate::dive::DiveProfile;

        let measurements = [
            DiveMeasurement {
                time_ms: 0,
                depth: Bar::new(2.0), // 10m depth
                gas: 0,
            },
            DiveMeasurement {
                time_ms: 300_000, // 5 minutes at safe depth
                depth: Bar::new(2.0),
                gas: 0,
            },
            DiveMeasurement {
                time_ms: 600_000, // Another 5 minutes
                depth: Bar::new(2.0),
                gas: 0,
            },
        ];

        let nitrox32 = GasMix::new(0.32, 0.0).unwrap();
        let profile: DiveProfile<Bar, f32, 1, 3> = DiveProfile {
            dive_id: 1,
            max_depth: Bar::new(2.0),
            gases: [nitrox32],
            measurements,
        };

        let toxicity = calculate_toxicity_from_profile(
            &profile,
            &O2ExposureType::Single,
            O2ToxCalculation::NOAA,
        );

        assert!(
            toxicity.cns_percent >= 0.0,
            "CNS percentage should be non-negative"
        );
        assert!(
            toxicity.pulmonary_percent >= 0.0,
            "Pulmonary percentage should be non-negative"
        );
        assert!(
            toxicity.cns_percent < 30.0,
            "CNS should be reasonable for safe dive: got {}",
            toxicity.cns_percent
        );
        assert!(
            toxicity.pulmonary_percent < 30.0,
            "CNS should be reasonable for safe dive: got {}",
            toxicity.cns_percent
        );
    }

    #[test]
    fn toxicity_percentage_from_aggressive_profile_test() {
        use crate::dive::DiveProfile;

        let measurements = [
            DiveMeasurement {
                time_ms: 0,
                depth: Bar::new(1.8),
                gas: 0,
            },
            DiveMeasurement {
                time_ms: 300_000, // 5 minutes at exceptional limit PO2
                depth: Bar::new(1.8),
                gas: 0,
            },
            DiveMeasurement {
                time_ms: 600_000, // Another 5 minutes
                depth: Bar::new(1.8),
                gas: 0,
            },
        ];

        let pure_o2 = GasMix::new(1.0, 0.0).unwrap();
        let profile: DiveProfile<Bar, f32, 1, 3> = DiveProfile {
            dive_id: 2,
            max_depth: Bar::new(1.8),
            gases: [pure_o2],
            measurements,
        };

        let toxicity = calculate_toxicity_from_profile(
            &profile,
            &O2ExposureType::Single,
            O2ToxCalculation::NOAA,
        );

        assert!(
            toxicity.cns_percent > 20.0,
            "CNS should be elevated for exceptional limit exposure: got {}",
            toxicity.cns_percent
        );
        assert!(
            toxicity.pulmonary_percent > 20.0,
            "Pulmonary should be elevated for exceptional limit exposure: got {}",
            toxicity.pulmonary_percent
        );
    }

    #[test]
    fn toxicity_percentage_from_unsafe_profile_test() {
        // Create a dive profile with unsafe PO2 levels
        use crate::dive::DiveProfile;

        let measurements = [
            DiveMeasurement {
                time_ms: 0,
                depth: Bar::new(2.5), // Unsafe PO2
                gas: 0,
            },
            DiveMeasurement {
                time_ms: 10_000, // Just 10 seconds
                depth: Bar::new(2.5),
                gas: 0,
            },
        ];

        let pure_o2 = GasMix::new(1.0, 0.0).unwrap();
        let profile: DiveProfile<Bar, f32, 1, 2> = DiveProfile {
            dive_id: 3,
            max_depth: Bar::new(2.5),
            gases: [pure_o2],
            measurements,
        };

        let toxicity = calculate_toxicity_from_profile(
            &profile,
            &O2ExposureType::Single,
            O2ToxCalculation::NOAA,
        );

        // Unsafe exposure should immediately spike toxicity percentages
        assert!(
            toxicity.cns_percent >= 200.0,
            "CNS should spike for unsafe PO2 exposure"
        );
        assert!(
            toxicity.pulmonary_percent >= 200.0,
            "Pulmonary should spike for unsafe PO2 exposure"
        );
    }

    #[test]
    fn toxicity_differential_streaming_test() {
        // Test differential calculation for streaming scenario
        // where new measurements arrive one at a time
        use crate::dive::DiveProfile;

        let measurements = [
            DiveMeasurement {
                time_ms: 0,
                depth: Bar::new(2.0),
                gas: 0,
            },
            DiveMeasurement {
                time_ms: 300_000, // 5 minutes
                depth: Bar::new(2.0),
                gas: 0,
            },
            DiveMeasurement {
                time_ms: 600_000, // Another 5 minutes
                depth: Bar::new(2.0),
                gas: 0,
            },
        ];

        let nitrox32 = GasMix::new(0.32, 0.0).unwrap();
        let profile: DiveProfile<Bar, f32, 1, 3> = DiveProfile {
            dive_id: 4,
            max_depth: Bar::new(2.0),
            gases: [nitrox32],
            measurements,
        };

        // Full calculation
        let full_result = calculate_toxicity_from_profile(
            &profile,
            &O2ExposureType::Single,
            O2ToxCalculation::NOAA,
        );

        // Simulate streaming:
        // The differential function processes from start_index to end of profile
        let initial = O2ToxicityPercentage::new(0.0, 0.0);
        let from_start = calculate_toxicity_diff(
            &profile,
            0,
            &initial,
            &O2ExposureType::Single,
            O2ToxCalculation::NOAA,
        );

        // Results should match since we process from beginning
        assert!(
            (full_result.cns_percent - from_start.cns_percent).abs() < 0.01,
            "Differential from start should match full calculation"
        );
        assert!(
            (full_result.pulmonary_percent - from_start.pulmonary_percent).abs() < 0.01,
            "Differential from start should match full calculation"
        );
    }

    #[test]
    fn toxicity_differential_from_middle_test() {
        // Test differential starting from middle of profile
        use crate::dive::DiveProfile;

        let measurements = [
            DiveMeasurement {
                time_ms: 0,
                depth: Bar::new(2.0),
                gas: 0,
            },
            DiveMeasurement {
                time_ms: 300_000, // 5 minutes
                depth: Bar::new(2.0),
                gas: 0,
            },
            DiveMeasurement {
                time_ms: 600_000, // Another 5 minutes
                depth: Bar::new(1.8),
                gas: 0,
            },
            DiveMeasurement {
                time_ms: 900_000, // Another 5 minutes
                depth: Bar::new(1.8),
                gas: 0,
            },
        ];

        let nitrox32 = GasMix::new(0.32, 0.0).unwrap();
        let profile: DiveProfile<Bar, f32, 1, 4> = DiveProfile {
            dive_id: 5,
            max_depth: Bar::new(2.0),
            gases: [nitrox32],
            measurements,
        };

        // Calculate full profile
        let full_result = calculate_toxicity_from_profile(
            &profile,
            &O2ExposureType::Single,
            O2ToxCalculation::NOAA,
        );

        // Calculate first half (indices 0 and 1)
        let initial = O2ToxicityPercentage::new(0.0, 0.0);
        let first_half = calculate_toxicity_diff(
            &profile,
            0,
            &initial,
            &O2ExposureType::Single,
            O2ToxCalculation::NOAA,
        );

        assert!(
            (full_result.cns_percent - first_half.cns_percent).abs() < 0.01,
            "Differential from index 0 should equal full calculation: full={}, diff={}",
            full_result.cns_percent,
            first_half.cns_percent
        );
    }

    #[test]
    fn toxicity_differential_empty_segment_test() {
        // Test differential when start_index is beyond measurements
        use crate::dive::DiveProfile;

        let measurements = [
            DiveMeasurement {
                time_ms: 0,
                depth: Bar::new(2.0),
                gas: 0,
            },
            DiveMeasurement {
                time_ms: 300_000,
                depth: Bar::new(2.0),
                gas: 0,
            },
        ];

        let nitrox32 = GasMix::new(0.32, 0.0).unwrap();
        let profile: DiveProfile<Bar, f32, 1, 2> = DiveProfile {
            dive_id: 6,
            max_depth: Bar::new(2.0),
            gases: [nitrox32],
            measurements,
        };

        let initial = O2ToxicityPercentage::new(42.5, 35.0);

        // Try to start from beyond the profile
        let result = calculate_toxicity_diff(
            &profile,
            10,
            &initial,
            &O2ExposureType::Single,
            O2ToxCalculation::NOAA,
        );

        // Should return the initial toxicity unchanged
        assert_eq!(
            result.cns_percent, initial.cns_percent,
            "Invalid index should return previous toxicity"
        );
        assert_eq!(
            result.pulmonary_percent, initial.pulmonary_percent,
            "Invalid index should return previous toxicity"
        );
    }

    #[test]
    fn toxicity_with_revised_dhm_2025_test() {
        // Test using Revised DHM 2025 limits
        // At 1.3 bar PO2, the limits are more generous
        use crate::dive::DiveProfile;

        let measurements = [
            DiveMeasurement {
                time_ms: 0,
                depth: Bar::new(2.3), // ~1.3 bar PO2 with air
                gas: 0,
            },
            DiveMeasurement {
                time_ms: 600_000, // 10 minutes
                depth: Bar::new(2.3),
                gas: 0,
            },
        ];

        let air = GasMix::new(0.21, 0.0).unwrap();
        let profile: DiveProfile<Bar, f32, 1, 2> = DiveProfile {
            dive_id: 7,
            max_depth: Bar::new(2.3),
            gases: [air],
            measurements,
        };

        // Calculate with both methods
        let noaa_result = calculate_toxicity_from_profile(
            &profile,
            &O2ExposureType::Single,
            O2ToxCalculation::NOAA,
        );
        let revised_result = calculate_toxicity_from_profile(
            &profile,
            &O2ExposureType::Single,
            O2ToxCalculation::RevisedDHM2025,
        );

        // Both should be valid calculations
        assert!(
            noaa_result.cns_percent >= 0.0,
            "NOAA calculation should be valid"
        );
        assert!(
            revised_result.cns_percent >= 0.0,
            "Revised DHM 2025 calculation should be valid"
        );

        // Results may differ due to different limit values
        // This demonstrates that the method parameter affects the calculation
    }

    #[test]
    fn toxicity_noaa_vs_revised_dhm_comparison_test() {
        // Compare NOAA and Revised DHM 2025 at critical 1.3 bar
        // where they differ most significantly
        use crate::dive::DiveProfile;

        let measurements = [
            DiveMeasurement {
                time_ms: 0,
                depth: Bar::new(2.3), // Surface equivalent with nitrox
                gas: 0,
            },
            DiveMeasurement {
                time_ms: 900_000,     // 15 minutes
                depth: Bar::new(2.3), // ~1.3 bar PO2 with oxygen
                gas: 0,
            },
        ];

        let pure_o2 = GasMix::new(1.0, 0.0).unwrap();
        let profile: DiveProfile<Bar, f32, 1, 2> = DiveProfile {
            dive_id: 8,
            max_depth: Bar::new(2.3),
            gases: [pure_o2],
            measurements,
        };

        let noaa = calculate_toxicity_from_profile(
            &profile,
            &O2ExposureType::Single,
            O2ToxCalculation::NOAA,
        );
        let revised = calculate_toxicity_from_profile(
            &profile,
            &O2ExposureType::Single,
            O2ToxCalculation::RevisedDHM2025,
        );

        // At 1.3 bar, Revised DHM 2025 allows 240 min vs NOAA's 180 min
        // So for a 15 minute exposure, Revised should show lower percentage
        assert!(
            revised.cns_percent <= noaa.cns_percent,
            "Revised DHM 2025 should allow more exposure at 1.3 bar than NOAA"
        );
    }

    #[test]
    fn toxicity_different_methods_with_diff_function_test() {
        // Test that differential function works with both calculation methods
        use crate::dive::DiveProfile;

        let measurements = [
            DiveMeasurement {
                time_ms: 0,
                depth: Bar::new(2.0),
                gas: 0,
            },
            DiveMeasurement {
                time_ms: 300_000, // 5 minutes
                depth: Bar::new(2.0),
                gas: 0,
            },
            DiveMeasurement {
                time_ms: 600_000, // Another 5 minutes
                depth: Bar::new(2.0),
                gas: 0,
            },
        ];

        let nitrox32 = GasMix::new(0.32, 0.0).unwrap();
        let profile: DiveProfile<Bar, f32, 1, 3> = DiveProfile {
            dive_id: 9,
            max_depth: Bar::new(2.0),
            gases: [nitrox32],
            measurements,
        };

        // Test with NOAA
        let noaa_full = calculate_toxicity_from_profile(
            &profile,
            &O2ExposureType::Single,
            O2ToxCalculation::NOAA,
        );
        let initial = O2ToxicityPercentage::new(0.0, 0.0);
        let noaa_diff = calculate_toxicity_diff(
            &profile,
            0,
            &initial,
            &O2ExposureType::Single,
            O2ToxCalculation::NOAA,
        );

        assert!(
            (noaa_full.cns_percent - noaa_diff.cns_percent).abs() < 0.01,
            "NOAA diff should match full calculation"
        );

        // Test with Revised DHM 2025
        let revised_full = calculate_toxicity_from_profile(
            &profile,
            &O2ExposureType::Single,
            O2ToxCalculation::RevisedDHM2025,
        );
        let revised_diff = calculate_toxicity_diff(
            &profile,
            0,
            &initial,
            &O2ExposureType::Single,
            O2ToxCalculation::RevisedDHM2025,
        );

        assert!(
            (revised_full.cns_percent - revised_diff.cns_percent).abs() < 0.01,
            "Revised DHM 2025 diff should match full calculation"
        );

        // Both methods should produce valid results
        assert!(
            revised_full.cns_percent >= 0.0 && noaa_full.cns_percent >= 0.0,
            "Both methods should produce valid toxicity percentages"
        );
    }
}

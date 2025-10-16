use core::time::Duration;

use crate::{
    gas,
    pressure_unit::{msw, Pa, Pressure},
};

use num;

#[derive(Debug, Clone)]
pub struct DiveMeasurement {
    pub time_ms: usize,
    pub depth: Pa,
    pub gas: usize,
}

#[derive(Debug, Clone)]
pub struct DiveProfile<F: num::Float, const G: usize, const M: usize> {
    pub dive_id: usize,
    pub max_depth: F,
    pub gases: [gas::GasMix<F>; G],
    pub measurements: [DiveMeasurement; M],
}

#[derive(Debug, Clone, Copy)]
pub struct Stop {
    depth: msw,
    duration: Duration,
}

impl Stop {
    pub fn new<P: Pressure>(depth: P, duration: Duration) -> Self {
        Stop {
            depth: depth.to_msw(),
            duration,
        }
    }

    pub fn depth(&self) -> msw {
        self.depth
    }

    pub fn duration(&self) -> Duration {
        self.duration
    }
}

#[derive(Debug, Clone)]
pub struct StopSchedule<const NUM_STOPS: usize> {
    pub stops: [Stop; NUM_STOPS],
    pub tts: Duration,
}

impl<const NUM_STOPS: usize> StopSchedule<NUM_STOPS> {
    pub fn new(stops: [Stop; NUM_STOPS]) -> Self {
        let tts: Duration = stops.iter().map(|s| s.duration()).sum();
        StopSchedule { stops, tts }
    }
}

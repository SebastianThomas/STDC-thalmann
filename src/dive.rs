use core::time::Duration;

use crate::{
    depth_utils::{get_ascent_time, get_depth},
    gas::{self, GasMix},
    pressure_unit::{msw, Pa, Pressure},
};

use num::Float;

#[derive(Debug, Clone)]
pub struct DiveMeasurement {
    pub time_ms: usize,
    pub depth: Pa,
    pub gas: usize,
}

#[derive(Debug, Clone)]
pub struct DiveProfile<P: Pressure, F: Float, const G: usize, const M: usize> {
    pub dive_id: usize,
    pub max_depth: P,
    pub gases: [gas::GasMix<F>; G],
    pub measurements: [DiveMeasurement; M],
}

#[derive(Debug, Clone, Copy)]
pub struct Stop {
    depth: msw,
    duration: Duration,
    gas: Option<GasMix<f32>>,
}

impl Stop {
    pub const fn new<P: const Pressure>(
        depth: P,
        duration: Duration,
        gas: Option<GasMix<f32>>,
    ) -> Self {
        Stop {
            depth: depth.to_msw(),
            duration,
            gas,
        }
    }

    pub fn depth(&self) -> msw {
        self.depth
    }

    pub fn duration(&self) -> Duration {
        self.duration
    }

    pub fn gas(&self) -> Option<GasMix<f32>> {
        self.gas
    }
}

#[derive(Debug, Clone)]
pub struct StopSchedule<const NUM_STOPS: usize> {
    stops: [Stop; NUM_STOPS],
}

impl<const NUM_STOPS: usize> const Default for StopSchedule<NUM_STOPS> {
    fn default() -> Self {
        let mut stops: [Stop; NUM_STOPS] =
            [Stop::new(msw::new(0.0), Duration::from_millis(0), None); NUM_STOPS];

        let mut i = 0;
        while i < NUM_STOPS {
            stops[NUM_STOPS - i - 1] =
                Stop::new(get_depth(i).to_msw(), Duration::from_millis(0), None);
            i += 1;
        }
        StopSchedule { stops }
    }
}

impl<const NUM_STOPS: usize> StopSchedule<NUM_STOPS> {
    pub fn new(stops: [Stop; NUM_STOPS]) -> Self {
        StopSchedule { stops }
    }

    pub fn first_stop(&self) -> Option<&Stop> {
        self.stops.iter().find(|stop| !stop.duration.is_zero())
    }

    pub fn get_deco_tts(&self, max_deco_ascent_rate_per_meter: &Duration) -> Duration {
        return match self.first_stop() {
            Some(first_stop) => {
                let stops_time: Duration = self.stops.iter().map(|s| s.duration()).sum();
                let stops_ascent_time: Duration =
                    get_ascent_time(first_stop.depth(), max_deco_ascent_rate_per_meter);
                stops_time + stops_ascent_time
            }
            None => Duration::ZERO,
        };
    }

    pub fn get_tt_first_stop_ascent_now<P: Pressure>(
        &self,
        current_depth: P,
        max_ascent_rate_per_meter: &Duration,
    ) -> Result<Duration, &'static str> {
        let first_stop = self.first_stop();
        if first_stop.is_none() {
            return Ok(Duration::ZERO);
        }
        let first_stop = first_stop.unwrap().depth();
        let current_depth = current_depth.to_msw();
        if current_depth < first_stop {
            return Err(
                "First stop must be still outstanding to get time to deco. Otherwise, use 0",
            );
        }
        let diff = current_depth - first_stop;
        return Ok(get_ascent_time(diff, max_ascent_rate_per_meter));
    }
}

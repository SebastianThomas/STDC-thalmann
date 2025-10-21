use core::time::Duration;

// Required for usage of float methods
#[allow(unused)]
use num::Float;

use crate::pressure_unit::{msw, Pa, Pressure};
use crate::setup::DINC;

pub const fn get_depth(d_idx: usize) -> Pa {
    DINC.to_pa() * d_idx as f32
}

/**
* If called with d in msw, no further unit conversion will be necessary.
*/
pub fn get_depth_idx<P: Pressure>(d: P) -> usize {
    assert!(d.to_f32() > 0.0);
    (d.to_msw().to_f32() / DINC.to_msw().to_f32()).ceil() as usize
}

pub fn get_ascent_time(meters: msw, max_ascent_rate_meters: &Duration) -> Duration {
    Duration::from_secs(meters.to_f32().ceil() as u64 / max_ascent_rate_meters.as_secs())
}

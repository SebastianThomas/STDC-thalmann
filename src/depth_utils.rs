use core::time::Duration;

use num::integer::mod_floor;
// Required for usage of float methods
#[allow(unused)]
use num::Float;
use num::ToPrimitive;

use crate::pressure_unit::{msw, Pa, Pressure};
use crate::setup::DINC;

pub const fn get_depth(d_idx: usize) -> Pa {
    if d_idx == 0 {
        msw::new(0.0).to_pa()
    } else {
        (DINC * d_idx as f32).to_pa()
    }
}

/**
* If called with d in msw, no further unit conversion will be necessary.
*/
pub fn get_depth_idx<P: Pressure>(d: P) -> usize {
    assert!(d.to_f32() > 0.0);
    (d.to_msw().to_f32() / DINC.to_msw().to_f32()).ceil() as usize
}

pub fn get_ascent_rate_per_meter(meters: u64) -> Duration {
    Duration::new(
        60 / meters,
        mod_floor(60_000_000_000_u64 / meters, 1_000_000_000_u64)
            .to_u32()
            .expect("Floor always smaller than 1E9"),
    )
}

pub fn get_ascent_time(meters: msw, max_ascent_rate_meters: &Duration) -> Duration {
    Duration::from_secs(meters.to_f32().ceil() as u64 / max_ascent_rate_meters.as_secs())
}

#[cfg(test)]
mod tests {
    use crate::pressure_unit::msw;

    use super::*;

    #[test]
    fn get_depth_test() {
        for i in 0..100 {
            assert_eq!(get_depth(i), msw::new((i * 3) as f32).to_pa());
        }
    }

    #[test]
    fn get_ascent_time_test() {
        assert_eq!(
            get_ascent_time(msw::new(0.0), &get_ascent_rate_per_meter(9)),
            Duration::new(0, 0)
        );
    }
}

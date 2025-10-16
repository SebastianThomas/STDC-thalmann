use num::Float;

use crate::pressure_unit::{Pa, Pressure};
use crate::setup::DINC;

pub fn get_depth(d_idx: usize) -> Pa {
    DINC.to_pa() * d_idx as f32
}

/**
* If called with d in msw, no further unit conversion will be necessary.
*/
pub fn get_depth_idx<P: Pressure>(d: P) -> usize {
    assert!(d.to_f32() > 0.0);
    (d.to_msw().to_f32() / DINC.to_msw().to_f32()).ceil() as usize
}

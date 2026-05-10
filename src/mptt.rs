use crate::pressure_unit::{AbsPressure, Pa, msw};

/** Maximum Permissible Tissue Tension */

#[derive(Copy, Clone)]
pub struct Tissue {
    pub half_time: f32,
    pub sdr: f32, /* Saturation Desaturation Ration */
}

#[derive(Copy, Clone)]
pub struct TissueRow<const TISSUES: usize, P: const AbsPressure> {
    pub depth: msw,
    pub max_saturation: [P; TISSUES],
}
impl<const TISSUES: usize, P: const AbsPressure> TissueRow<TISSUES, P> {
    pub const fn empty(val: P) -> Self {
        TissueRow {
            depth: msw(-1.0),
            max_saturation: [val; TISSUES],
        }
    }
}
impl<const TISSUES: usize> TissueRow<TISSUES, Pa> {
    pub const fn empty_pa() -> Self {
        Self::empty(Pa::new(-1.0))
    }
}
pub type MValues<P, const NUM_TISSUES: usize, const NUM_STOP_DEPTHS: usize> =
    [TissueRow<NUM_TISSUES, P>; NUM_STOP_DEPTHS];

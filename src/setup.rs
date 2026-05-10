#[cfg(not(feature = "lin_exp"))]
pub use crate::mptt_buehlmann::{
    BUEHLMANN_16C as MVALUES_TISSUES, NUM_STOP_DEPTHS_BUEHLMANN as NUM_STOP_DEPTHS,
    NUM_TISSUES_BUEHLMANN as NUM_TISSUES,
};
#[cfg(feature = "lin_exp")]
pub use crate::mptt_thalmann::{
    NUM_STOP_DEPTHS_THALMANN as NUM_STOP_DEPTHS, NUM_TISSUES_THALMANN as NUM_TISSUES,
    XVAL_HE9_040_F32 as MVALUES_TISSUES,
};
use crate::{
    mptt::{MValues, TissueRow},
    pressure_unit::{Pa, Pressure, msw},
};

pub const MSW_0_PA: Pa = msw::new(0.0).to_pa();

// Depth Increment
pub const DINC: msw = msw::new(3.0);
pub const DINC_PA: Pa = DINC.to_pa() - MSW_0_PA;
// IDX * DINC
pub const LAST_STOP: msw = msw::new(6.0);

pub const fn set_m(mode: u8) -> MValues<Pa, { NUM_TISSUES }, { NUM_STOP_DEPTHS }> {
    if mode == 0 || mode != 1 {
        return MVALUES_TISSUES;
    }
    let idx = (LAST_STOP.0 / DINC.0) as isize;
    if idx <= 1 {
        return set_m(0);
    }
    let idx = idx as usize;
    let mut result = MVALUES_TISSUES;
    // Copy surfacing MVals to IDX row // TODO: What are surfacing?
    result[idx] = result[0];
    let mut i = 0;
    while i < idx {
        // Zero MVals for stop depths shallower than IDX row
        result[i] = TissueRow {
            depth: result[i].depth,
            max_saturation: [msw::new(0.0).to_pa(); NUM_TISSUES],
        };
        i += 1;
    }
    result
}

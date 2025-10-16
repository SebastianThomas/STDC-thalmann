use crate::mptt::{TissueRow, XVAL_HE9_040_F32};
pub use crate::mptt::{MVALUES, NUM_TISSUES};
use crate::pressure_unit::{msw, Pressure};

// Depth Increment
pub const DINC: msw = msw::new(3.0);
// IDX * DINC
pub const LAST_STOP: msw = msw::new(6.0);

pub struct ModelState {
    // finished_profile: bool,            // DONE
    // first_stop_calculated: bool,       // FRSTOP
    // include_travel_time_in_stop: bool, // TTIS

    // elapsed_time: f32,                       // TIME
    // cumulative_ascent_time: f32,             // TT_SUM
    // cumulative_ascent_included_in_stop: f32, // TTSTC_SUM
    // previous_tt: f32,                        // TT0
}

pub fn set_m(mode: u8) -> MVALUES {
    if mode == 0 || mode != 1 {
        return XVAL_HE9_040_F32;
    }
    let idx = (LAST_STOP.0 / DINC.0) as isize;
    if idx <= 1 {
        return set_m(0);
    }
    let idx = idx as usize;
    let mut result = XVAL_HE9_040_F32;
    // Copy surfacing MVals to IDX row // TODO: What are surfacing?
    result[idx] = result[0];
    for i in 0..idx {
        // Zero MVals for stop depths shallower than IDX row
        result[i] = TissueRow {
            depth: result[i].depth,
            max_saturation: [msw::new(0.0).to_pa(); NUM_TISSUES],
        };
    }
    return result;
}

pub fn initialize_profile() {
    ()
}
pub fn initialize_model_state() -> ModelState {
    // --- Flags and booleans ---
    ModelState {
        // finished_profile: false,
        // first_stop_calculated: false,
        // include_travel_time_in_stop: true,

        // elapsed_time: 0.0,
        // cumulative_ascent_time: 0.0,
        // cumulative_ascent_included_in_stop: 0.0,
        // previous_tt: 0.0,
    }

    // --- Other model state ---
    // For now we assume tissue loadings (P) are already passed in
    // so we do not initialize them here

    // If needed, we could also reset repetitive group reference:
    // let reference_compartment: usize = 0;   // IREF
    // let first_profile_to_read: bool = true; // PROFL1

    // The function essentially prepares the loop to start iterating
    // without modifying the current tissue loadings
}

// Note: While the original paper on this algorithm is in imperial units, this is a trancription in
// metric units. The paper for reference can be found on this url: https://apps.dtic.mil/sti/tr/pdf/ADA549883.pdf

#![no_std]
#![feature(const_trait_impl)]
#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

mod depth_utils;
pub mod display_utils;
pub mod dive;
pub mod gas;
pub mod mptt;
pub mod pressure_unit;
mod setup;
mod thalmann;
mod time_utils;
mod update;
mod update_exp;
mod update_thalmann;

pub use setup::{DINC, LAST_STOP};
pub use thalmann::{calc_deco_schedule, thalmann, ThalmannResult, MVALUES_HE9_040};
pub use update::loadings_from_dive_profile;

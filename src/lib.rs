// Note: While the original paper on this algorithm is in imperial units, this is a trancription in
// metric units. The paper for reference can be found on this url: https://apps.dtic.mil/sti/tr/pdf/ADA549883.pdf

#![no_std]
#![feature(const_trait_impl)]
#![feature(const_default)]
#![feature(const_ops)]
#![feature(const_option_ops)]
#![feature(const_array)]
#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(const_convert)]
#![feature(const_cmp)]
#![feature(derive_const)]

mod deco_algorithm;
mod depth_utils;
pub mod display_utils;
pub mod dive;
pub mod gas;
pub mod mptt;
#[cfg(not(feature = "lin_exp"))]
mod mptt_buehlmann;
#[cfg(feature = "lin_exp")]
mod mptt_thalmann;
pub mod o2tox;
pub mod pressure_unit;
mod setup;
mod time_utils;
mod update;
mod update_exp;
mod update_exp_lin;

pub use deco_algorithm::{
    DecoAlgorithmResult, MVALUES, MValues, calc_deco_schedule, compute_deco_algorithm,
};
pub use setup::{DINC, LAST_STOP};
pub use update::loadings_from_dive_profile;

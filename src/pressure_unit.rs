use core::clone::Clone;
use core::cmp::{PartialEq, PartialOrd};
use core::convert::From;
use core::fmt::Debug;
use core::marker::Copy;
use core::ops::{Add, AddAssign, Div, Mul, Sub};

use core::prelude::rust_2024::derive;

/// Common trait for all pressure units
pub const trait Pressure:
    Copy
    + Clone
    + Debug
    + PartialEq
    + PartialOrd
    + Add
    + Sub
    + Mul<f32, Output = Self>
    + Div<f32, Output = Self>
{
    fn to_pa(self) -> Pa;
    fn to_kpa(self) -> kPa;
    fn to_hpa(self) -> hPa;
    fn to_bar(self) -> Bar;
    fn to_msw(self) -> msw;
    fn to_f32(self) -> f32;
}

/// Macro to generate pressure unit newtypes + trait impl + arithmetic
macro_rules! pressure_unit {
    ($name:ident, $to_pa_factor:expr) => {
        #[allow(non_camel_case_types)]
        #[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
        pub struct $name(f32);

        impl $name {
            pub const fn new(v: f32) -> Self {
                Self(v)
            }
        }

        impl const Pressure for $name {
            fn to_pa(self) -> Pa {
                Pa(self.0 * $to_pa_factor)
            }
            fn to_hpa(self) -> hPa {
                hPa(self.0 * $to_pa_factor / 100.0)
            }
            fn to_kpa(self) -> kPa {
                kPa(self.0 * $to_pa_factor / 1000.0)
            }
            fn to_bar(self) -> Bar {
                Bar(self.0 * $to_pa_factor / 1E5)
            }
            fn to_msw(self) -> msw {
                msw((self.0 * $to_pa_factor - 1E5) / 1.013E5)
            }
            fn to_f32(self) -> f32 {
                self.0
            }
        }

        /// Construct from f32 directly
        impl From<f32> for $name {
            fn from(v: f32) -> Self {
                $name(v)
            }
        }

        /// Arithmetic
        impl Add for $name {
            type Output = $name;
            fn add(self, rhs: $name) -> $name {
                $name(self.0 + rhs.0)
            }
        }
        impl AddAssign for $name {
            fn add_assign(&mut self, rhs: $name) {
                self.0 += rhs.0;
            }
        }
        impl Sub for $name {
            type Output = $name;
            fn sub(self, rhs: $name) -> $name {
                $name(self.0 - rhs.0)
            }
        }
        impl Mul<f32> for $name {
            type Output = $name;
            fn mul(self, rhs: f32) -> $name {
                $name(self.0 * rhs)
            }
        }
        impl Div<f32> for $name {
            type Output = $name;
            fn div(self, rhs: f32) -> $name {
                $name(self.0 / rhs)
            }
        }
        impl Div<$name> for $name {
            type Output = f32;
            fn div(self, rhs: $name) -> f32 {
                self.0 / rhs.0
            }
        }
    };
}

macro_rules! pressure_unit_relative {
    ($name:ident, $to_pa_factor:expr) => {
        #[allow(non_camel_case_types)]
        #[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
        pub struct $name(pub f32);

        impl $name {
            pub const fn new(v: f32) -> Self {
                Self(v)
            }
        }

        impl const Pressure for $name {
            fn to_pa(self) -> Pa {
                Pa((self.0 * $to_pa_factor) - 1E5)
            }
            fn to_hpa(self) -> hPa {
                hPa(self.to_pa().0 / 100.0)
            }
            fn to_kpa(self) -> kPa {
                kPa(self.to_pa().0 / 1000.0)
            }
            fn to_bar(self) -> Bar {
                Bar(self.to_pa().0 / 10E5)
            }
            fn to_msw(self) -> msw {
                msw(self.to_pa().0 / 1.013E5)
            }
            fn to_f32(self) -> f32 {
                self.0
            }
        }

        /// Construct from f32 directly
        impl From<f32> for $name {
            fn from(v: f32) -> Self {
                $name(v)
            }
        }

        /// Arithmetic
        impl Add for $name {
            type Output = $name;
            fn add(self, rhs: $name) -> $name {
                $name(self.0 + rhs.0)
            }
        }
        impl Sub for $name {
            type Output = $name;
            fn sub(self, rhs: $name) -> $name {
                $name(self.0 - rhs.0)
            }
        }
        impl Mul<f32> for $name {
            type Output = $name;
            fn mul(self, rhs: f32) -> $name {
                $name(self.0 * rhs)
            }
        }
        impl Mul<usize> for $name {
            type Output = $name;
            fn mul(self, rhs: usize) -> $name {
                $name(self.0 * rhs as f32)
            }
        }
        impl Div<f32> for $name {
            type Output = $name;
            fn div(self, rhs: f32) -> $name {
                $name(self.0 / rhs)
            }
        }
    };
}

// -------------------- Base units --------------------
pressure_unit!(Pa, 1.0);
pressure_unit!(kPa, 1000.0);
pressure_unit!(Bar, 100_000.0);
pressure_unit!(hPa, 100.0);
pressure_unit_relative!(msw, 1.013E5);
pressure_unit_relative!(fsw, 3064.30593138);
#[allow(non_camel_case_types)]
pub type mBar = hPa; // alias

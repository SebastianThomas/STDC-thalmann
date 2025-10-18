use core::ops::Mul;

use num::Float;

use crate::pressure_unit::{msw, Pressure};

pub const N2_IDX: usize = 0;
pub const HE_IDX: usize = 1;

#[derive(Debug, Clone)]
pub struct GasMix<F: Float> {
    o2: F,
    he: F,
    h2: F,
}

impl GasMix<f32> {
    pub const fn new(o2: f32, he: f32) -> Result<GasMix<f32>, &'static str> {
        if o2 + he > 1.0 {
            return Err("FO2 + FHe should be <= 1");
        }
        Ok(GasMix { o2, he, h2: 0.0 })
    }

    pub const fn o2(&self) -> f32 {
        self.o2
    }

    pub const fn he(&self) -> f32 {
        self.he
    }

    pub const fn h2(&self) -> f32 {
        self.h2
    }

    pub const fn n2(&self) -> f32 {
        1.0 - (self.o2 + self.he)
    }
}

pub const AIR: GasMix<f32> = match GasMix::new(0.79, 0.000_005_2) {
    Ok(g) => g,
    Err(_) => unreachable!(),
};

#[derive(Debug, Clone)]
pub struct TissuesLoading<const NUM_TISSUES: usize, P: Pressure> {
    pub n2: [P; NUM_TISSUES],
    pub he: [P; NUM_TISSUES],
}

impl<const NUM_TS: usize, P: const Pressure + const Mul<f32>> TissuesLoading<NUM_TS, P> {
    pub const fn new(
        ambient: P,
        breathing_gas: &GasMix<f32>,
    ) -> TissuesLoading<NUM_TS, P> {
        TissuesLoading {
            n2: [ambient * breathing_gas.n2(); NUM_TS],
            he: [ambient * breathing_gas.he(); NUM_TS],
        }
    }
}

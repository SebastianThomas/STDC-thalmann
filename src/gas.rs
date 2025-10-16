use num::Float;

use crate::pressure_unit::Pressure;

pub const N2_IDX: usize = 0;
pub const HE_IDX: usize = 1;

#[derive(Debug, Clone)]
pub struct GasMix<F: Float> {
    o2: F,
    he: F,
    h2: F,
}

impl GasMix<f32> {
    pub fn new(o2: f32, he: f32) -> Result<GasMix<f32>, &'static str> {
        if o2 + he > 1.0 {
            return Err("FO2 + FHe should be <= 1");
        }
        Ok(GasMix { o2, he, h2: 0.0 })
    }

    pub fn o2(&self) -> f32 {
        self.o2
    }

    pub fn he(&self) -> f32 {
        self.he
    }

    pub fn h2(&self) -> f32 {
        self.h2
    }

    pub const fn n2(&self) -> f32 {
        1.0 - (self.o2 + self.he)
    }
}

#[derive(Debug, Clone)]
pub struct TissuesLoading<const NUM_TISSUES: usize, P: Pressure> {
    pub n2: [P; NUM_TISSUES],
    pub he: [P; NUM_TISSUES],
}

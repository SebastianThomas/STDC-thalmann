use core::{iter::zip, ops::Mul};

#[allow(unused)]
use num::Float;

use crate::pressure_unit::{Bar, Pressure};

pub const N2_IDX: usize = 0;
pub const HE_IDX: usize = 1;

pub const MAX_PO2_CCR_WORKING: Bar = Bar::new(1.3);
pub const MAX_PO2_WORKING: Bar = Bar::new(1.4);
pub const MAX_O2_CCR_DECO: Bar = Bar::new(1.5);
pub const MAX_PO2_DECO: Bar = Bar::new(1.6);
pub const MAX_O2_DILUENT: Bar = Bar::new(1.1);

pub const DENSITY_O2: f32 = 1.43;
pub const DENSITY_N2: f32 = 1.2506;
pub const DENSITY_HE: f32 = 0.1785;
pub const DENSITY_H2: f32 = 0.0899;
pub const DENSITY_AIR: f32 = 1.205;

pub const MAX_GAS_DENSITY: Bar = Bar::new(5.2);
pub const MAX_GAS_DENSITY_LIMIT: Bar = Bar::new(6.2);

#[derive(Debug, Clone, Copy)]
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

    pub const fn o2_density(&self) -> f32 {
        self.o2() * DENSITY_O2
    }

    pub const fn n2_density(&self) -> f32 {
        self.n2() * DENSITY_N2
    }

    pub const fn he_density(&self) -> f32 {
        self.he() * DENSITY_HE
    }

    pub const fn h2_density(&self) -> f32 {
        self.h2() * DENSITY_H2
    }

    pub const fn gas_density(&self) -> f32 {
        self.o2_density() + self.n2_density() + self.he_density() + self.h2_density()
    }

    pub const fn gas_density_depth<P: Pressure>(&self, depth: P) -> P {
        depth * self.gas_density()
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

impl<const NUM_TS: usize, P: const Pressure> TissuesLoading<NUM_TS, P> {
    pub const fn new(ambient: P, breathing_gas: &GasMix<f32>) -> TissuesLoading<NUM_TS, P> {
        TissuesLoading {
            n2: [ambient * breathing_gas.n2(); NUM_TS],
            he: [ambient * breathing_gas.he(); NUM_TS],
        }
    }

    pub fn is_isobaric_counterdiffusion<D: Pressure>(
        &self,
        depth: D,
        new_gas: &GasMix<f32>,
    ) -> bool {
        let depth = depth.to_pa();
        let new_gas_n2 = depth * new_gas.n2();
        let new_gas_he = depth * new_gas.he();
        return zip(self.n2, self.he)
            .any(|(n2, he)| n2.to_pa() < new_gas_n2 && he.to_pa() > new_gas_he);
    }
}

pub fn best_mix_o2<P: Pressure>(max_po2: Bar, depth: P) -> f32 {
    max_po2 / depth.to_bar()
}

pub enum GasDensitySettings {
    Ignore,
    Limit { limit: Bar },
}

impl GasDensitySettings {
    pub fn no_violation<P: Pressure>(&self, depth: P, gas: &GasMix<f32>) -> bool {
        if let GasDensitySettings::Limit { limit } = self {
            return gas.gas_density_depth(depth).to_bar() < *limit;
        }
        return true;
    }
}

/**
* Returns None iff available_gases is empty, or no gas fits
* - o2 requirements or
* - (optional) isobaric counterdiffusion requirements or
* - (optional) gas density requirements.
*
* Best performance can be expected if called with P = D = Pa
*/
pub fn best_available_mix<
    'a,
    D: Pressure,
    P: const Pressure + const Mul<f32>,
    const G: usize,
    const NUM_TS: usize,
>(
    max_po2: Bar,
    depth: D,
    available_gases: &'a [GasMix<f32>; G],
    tissue_loading: &TissuesLoading<NUM_TS, P>,
    ignore_isobaric_counterdiffusion: bool,
    gas_density: &GasDensitySettings,
) -> Option<(usize, &'a GasMix<f32>)> {
    let best_mix_po2 = best_mix_o2(max_po2, depth);
    available_gases
        .iter()
        .enumerate()
        .filter(|(_i, g)| g.o2() <= best_mix_po2)
        .filter(|(_i, g)| {
            ignore_isobaric_counterdiffusion
                || tissue_loading.is_isobaric_counterdiffusion(depth, g)
        })
        .filter(|(_i, g)| gas_density.no_violation(depth, g))
        .reduce(|a, b| if a.1.o2() > b.1.o2() { a } else { b })
}

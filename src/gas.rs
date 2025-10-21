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

pub const trait Gas {
    fn po2<D: const Pressure>(&self, depth: D) -> Bar;
    fn pn2<D: const Pressure>(&self, depth: D) -> Bar;
    fn phe<D: const Pressure>(&self, depth: D) -> Bar;
    fn ph2<D: const Pressure>(&self, depth: D) -> Bar;
    fn pn2_phe_ph2<D: const Pressure>(&self, depth: D) -> (Bar, Bar, Bar);

    pub const fn o2_density<P: const Pressure>(&self, depth: P) -> Bar {
        self.po2(depth) * DENSITY_O2
    }

    pub const fn n2_density<P: const Pressure>(&self, depth: P) -> Bar {
        self.pn2(depth) * DENSITY_N2
    }

    pub const fn he_density<P: const Pressure>(&self, depth: P) -> Bar {
        self.phe(depth) * DENSITY_HE
    }

    pub const fn h2_density<P: const Pressure>(&self, depth: P) -> Bar {
        self.ph2(depth) * DENSITY_H2
    }

    pub const fn gas_density<P: const Pressure>(&self, depth: P) -> Bar {
        self.o2_density(depth)
            + self.n2_density(depth)
            + self.he_density(depth)
            + self.h2_density(depth)
    }
}

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

    pub const fn fo2(&self) -> f32 {
        self.o2
    }

    pub const fn fn2(&self) -> f32 {
        1.0 - (self.o2 + self.he + self.h2)
    }

    pub const fn fhe(&self) -> f32 {
        self.he
    }

    pub const fn fh2(&self) -> f32 {
        self.h2
    }
}

impl const Gas for GasMix<f32> {
    pub const fn po2<D: const Pressure>(&self, depth: D) -> Bar {
        depth.to_bar() * self.o2
    }

    pub const fn phe<D: const Pressure>(&self, depth: D) -> Bar {
        depth.to_bar() * self.he
    }

    pub const fn ph2<D: const Pressure>(&self, depth: D) -> Bar {
        depth.to_bar() * self.h2
    }

    pub const fn pn2<D: const Pressure>(&self, depth: D) -> Bar {
        depth.to_bar() * self.fn2()
    }

    pub const fn pn2_phe_ph2<D: const Pressure>(&self, depth: D) -> (Bar, Bar, Bar) {
        (self.pn2(depth), self.phe(depth), self.ph2(depth))
    }
}

pub struct CCRGas<F: Float> {
    diluent: GasMix<F>,
    set_point: Bar,
}

impl const Gas for CCRGas<f32> {
    pub const fn po2<D: const Pressure>(&self, depth: D) -> Bar {
        let depth = depth.to_bar();
        if depth.to_f32() < self.set_point.to_bar().to_f32() {
            depth
        } else {
            self.set_point
        }
    }

    pub const fn pn2_phe_ph2<D: const Pressure>(&self, depth: D) -> (Bar, Bar, Bar) {
        let depth = depth.to_bar();
        let po2 = self.po2(depth);
        let fo2_loop = po2 / depth;
        let fo2_dil = self.diluent.fo2();
        let fo2_from_dil = (fo2_loop - fo2_dil) / (1.0 - fo2_dil);
        let fn2_loop = (1.0 - fo2_from_dil) * self.diluent.fn2();
        let fhe_loop = (1.0 - fo2_from_dil) * self.diluent.fhe();
        let fh2_loop = (1.0 - fo2_from_dil) * self.diluent.fh2();
        return (depth * fn2_loop, depth * fhe_loop, depth * fh2_loop);
    }

    pub fn pn2<D: const Pressure>(&self, depth: D) -> Bar {
        self.pn2_phe_ph2(depth).0
    }

    pub fn phe<D: const Pressure>(&self, depth: D) -> Bar {
        self.pn2_phe_ph2(depth).1
    }

    pub fn ph2<D: const Pressure>(&self, depth: D) -> Bar {
        self.pn2_phe_ph2(depth).2
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
            n2: [ambient * breathing_gas.fn2(); NUM_TS],
            he: [ambient * breathing_gas.fhe(); NUM_TS],
        }
    }

    pub fn is_isobaric_counterdiffusion<D: Pressure, G: const Gas>(
        &self,
        depth: D,
        new_gas: &G,
    ) -> bool {
        let depth = depth.to_bar();
        let new_gas_n2 = new_gas.pn2(depth);
        let new_gas_he = new_gas.phe(depth);
        return zip(self.n2, self.he)
            .any(|(n2, he)| n2.to_bar() < new_gas_n2 && he.to_bar() > new_gas_he);
    }
}

pub fn best_mix_fo2<P: Pressure>(max_po2: Bar, depth: P) -> f32 {
    max_po2 / depth.to_bar()
}

pub enum GasDensitySettings {
    Ignore,
    Limit { limit: Bar },
}

impl GasDensitySettings {
    pub fn no_violation<P: const Pressure>(&self, depth: P, gas: &GasMix<f32>) -> bool {
        if let GasDensitySettings::Limit { limit } = self {
            return gas.gas_density(depth).to_bar() < *limit;
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
    D: const Pressure,
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
    let best_mix_fo2 = best_mix_fo2(max_po2, depth);
    available_gases
        .iter()
        .enumerate()
        .filter(|(_i, g)| g.fo2() <= best_mix_fo2)
        .filter(|(_i, g)| {
            ignore_isobaric_counterdiffusion
                || tissue_loading.is_isobaric_counterdiffusion(depth, g)
        })
        .filter(|(_i, g)| gas_density.no_violation(depth, g))
        .reduce(|a, b| if a.1.fo2() > b.1.fo2() { a } else { b })
}

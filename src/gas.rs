use core::iter::zip;

#[allow(unused)]
use num::Float;

use crate::pressure_unit::{AbsPressure, Bar};

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
    fn po2<D: const AbsPressure>(&self, depth: D) -> D;
    fn pn2<D: const AbsPressure>(&self, depth: D) -> D;
    fn phe<D: const AbsPressure>(&self, depth: D) -> D;
    fn ph2<D: const AbsPressure>(&self, depth: D) -> D;
    fn pn2_phe_ph2<D: const AbsPressure>(&self, depth: D) -> (D, D, D);

    fn o2_density<P: const AbsPressure>(&self, depth: P) -> P {
        self.po2(depth) * DENSITY_O2
    }

    fn n2_density<P: const AbsPressure>(&self, depth: P) -> P {
        self.pn2(depth) * DENSITY_N2
    }

    fn he_density<P: const AbsPressure>(&self, depth: P) -> P {
        self.phe(depth) * DENSITY_HE
    }

    fn h2_density<P: const AbsPressure>(&self, depth: P) -> P {
        self.ph2(depth) * DENSITY_H2
    }

    fn gas_density<P: const AbsPressure>(&self, depth: P) -> P {
        self.o2_density(depth)
            + self.n2_density(depth)
            + self.he_density(depth)
            + self.h2_density(depth)
    }

    fn fio2<D: const AbsPressure>(&self, depth: D) -> f32 {
        self.po2(depth) / depth
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
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
    fn po2<D: const AbsPressure>(&self, depth: D) -> D {
        depth * self.o2
    }

    fn phe<D: const AbsPressure>(&self, depth: D) -> D {
        depth * self.he
    }

    fn ph2<D: const AbsPressure>(&self, depth: D) -> D {
        depth * self.h2
    }

    fn pn2<D: const AbsPressure>(&self, depth: D) -> D {
        depth * self.fn2()
    }

    fn pn2_phe_ph2<D: const AbsPressure>(&self, depth: D) -> (D, D, D) {
        (self.pn2(depth), self.phe(depth), self.ph2(depth))
    }
}

pub struct CCRGas<F: Float, P: const AbsPressure> {
    pub diluent: GasMix<F>,
    pub set_point: P,
}

impl<P: const AbsPressure> CCRGas<f32, P> {
    pub fn to_fixed_gas_mix<D: const AbsPressure>(&self, depth: D) -> GasMix<f32> {
        let (current_pn2, current_phe, _) = self.pn2_phe_ph2(depth);
        let current_fn2 = current_pn2 / depth;
        let current_fhe = current_phe / depth;
        GasMix {
            o2: self.fio2(depth),
            he: current_fn2,
            h2: current_fhe,
        }
    }
}

impl<P: const AbsPressure> const Gas for CCRGas<f32, P> {
    fn po2<D: const AbsPressure>(&self, depth: D) -> D {
        let set_point = D::from(self.set_point.to_pa());
        if depth < set_point {
            depth
        } else {
            set_point
        }
    }

    fn pn2_phe_ph2<D: const AbsPressure>(&self, depth: D) -> (D, D, D) {
        let po2 = self.po2(depth);
        let fo2_loop = po2 / depth;
        let fo2_dil = self.diluent.fo2();
        let fo2_from_dil = (fo2_loop - fo2_dil) / (1.0 - fo2_dil);
        let fn2_loop = (1.0 - fo2_from_dil) * self.diluent.fn2();
        let fhe_loop = (1.0 - fo2_from_dil) * self.diluent.fhe();
        let fh2_loop = (1.0 - fo2_from_dil) * self.diluent.fh2();
        return (depth * fn2_loop, depth * fhe_loop, depth * fh2_loop);
    }

    fn pn2<D: const AbsPressure>(&self, depth: D) -> D {
        self.pn2_phe_ph2(depth).0
    }

    fn phe<D: const AbsPressure>(&self, depth: D) -> D {
        self.pn2_phe_ph2(depth).1
    }

    fn ph2<D: const AbsPressure>(&self, depth: D) -> D {
        self.pn2_phe_ph2(depth).2
    }
}

pub const AIR: GasMix<f32> = match GasMix::new(0.79, 0.000_005_2) {
    Ok(g) => g,
    Err(_) => unreachable!(),
};

#[derive(Debug, Clone)]
pub struct TissuesLoading<const NUM_TISSUES: usize, P: const AbsPressure> {
    pub n2: [P; NUM_TISSUES],
    pub he: [P; NUM_TISSUES],
}

impl<const NUM_TS: usize, P: const AbsPressure> TissuesLoading<NUM_TS, P> {
    pub const fn new(ambient: P, breathing_gas: &GasMix<f32>) -> TissuesLoading<NUM_TS, P> {
        TissuesLoading {
            n2: [ambient * breathing_gas.fn2(); NUM_TS],
            he: [ambient * breathing_gas.fhe(); NUM_TS],
        }
    }

    pub fn is_isobaric_counterdiffusion<G: Gas>(&self, depth: P, new_gas: &G) -> bool {
        let new_gas_n2 = new_gas.pn2(depth);
        let new_gas_he = new_gas.phe(depth);
        return zip(self.n2, self.he).any(|(n2, he)| n2 < new_gas_n2 && he > new_gas_he);
    }

    pub fn tick<G: Gas>(&mut self, time_delta_ms: u16, depth: P, gas: &G) {
        Self::tick_gas(time_delta_ms, gas.pn2(depth), &mut self.n2);
        Self::tick_gas(time_delta_ms, gas.phe(depth), &mut self.he);
    }

    fn tick_gas(time_delta_ms: u16, pp_insp: P, cur: &mut [P; NUM_TS]) {
        for i in 0..NUM_TS {
            let delta: f32 = f32::from(time_delta_ms) / 1000.0;
            cur[i] = cur[i] + (pp_insp - cur[i]) * delta;
        }
    }
}

pub fn best_mix_fo2<P: const AbsPressure>(max_po2: P, depth: P) -> f32 {
    max_po2 / depth
}

pub enum GasDensitySettings<P: const AbsPressure> {
    Ignore,
    Limit { limit: P },
}

impl<P: const AbsPressure> GasDensitySettings<P> {
    pub fn no_violation(&self, depth: P, gas: &GasMix<f32>) -> bool {
        if let GasDensitySettings::Limit { limit } = self {
            return gas.gas_density(depth) < *limit;
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
pub fn best_available_mix<'a, P: const AbsPressure, const G: usize, const NUM_TS: usize>(
    max_po2: P,
    depth: P,
    available_gases: &'a [GasMix<f32>; G],
    tissue_loading: &TissuesLoading<NUM_TS, P>,
    ignore_isobaric_counterdiffusion: bool,
    gas_density: &GasDensitySettings<P>,
) -> Option<(usize, &'a GasMix<f32>)> {
    let best_mix_fo2 = best_mix_fo2(max_po2, depth);
    available_gases
        .iter()
        .enumerate()
        .filter(|(_i, g)| g.fo2() <= best_mix_fo2)
        .filter(|(_i, g)| {
            ignore_isobaric_counterdiffusion
                || tissue_loading.is_isobaric_counterdiffusion(depth, *g)
        })
        .filter(|(_i, g)| gas_density.no_violation(depth, g))
        .reduce(|a, b| {
            let better_fo2 = a.1.fo2() > b.1.fo2();
            let same_fo2_better_he = a.1.fo2() == b.1.fo2() && a.1.fhe() > b.1.fhe();
            if better_fo2 || same_fo2_better_he {
                a
            } else {
                b
            }
        })
}

#[cfg(test)]
mod tests {
    use crate::pressure_unit::{msw, Pressure};

    use super::*;

    #[test]
    fn best_mix_fo2_test() {
        assert_eq!(
            best_mix_fo2(Bar::new(1.6).to_pa(), msw::new(0.0).to_pa()),
            1.6
        );
        assert_eq!(
            best_mix_fo2(Bar::new(1.6).to_pa(), msw::new(6.0).to_pa()) - 1.0 < 0.01,
            true
        );
        assert_eq!(
            best_mix_fo2(Bar::new(1.6).to_pa(), msw::new(21.0).to_pa()) - 0.5 < 0.1,
            true
        );
        assert_eq!(
            best_mix_fo2(Bar::new(1.4).to_pa(), msw::new(0.0).to_pa()),
            1.4
        );
        assert_eq!(
            best_mix_fo2(Bar::new(1.4).to_pa(), msw::new(4.0).to_pa()) - 1.0 < 0.01,
            true
        );
        assert_eq!(
            best_mix_fo2(Bar::new(1.4).to_pa(), msw::new(18.0).to_pa()) - 0.5 < 0.1,
            true
        );
    }

    #[test]
    fn best_available_mix_test() {
        let gases = [
            AIR,
            GasMix::new(0.21, 0.35).expect("21 + 35 < 100"),
            GasMix::new(0.5, 0.0).expect("50 < 100"),
        ];
        let empty_tissues = TissuesLoading {
            n2: [msw::new(0.0).to_pa()],
            he: [msw::new(0.0).to_pa()],
        };
        assert_eq!(
            best_available_mix(
                Bar::new(1.6).to_pa(),
                msw::new(21.0).to_pa(),
                &gases,
                &empty_tissues,
                true,
                &GasDensitySettings::Ignore
            )
            .expect("There are gases, so reduce should return a result"),
            (2_usize, &gases[2])
        );
        assert_eq!(
            best_available_mix(
                Bar::new(1.4).to_pa(),
                msw::new(21.0).to_pa(),
                &gases,
                &empty_tissues,
                true,
                &GasDensitySettings::Ignore
            )
            .expect("There are gases, so reduce should return a result"),
            (1_usize, &gases[1])
        );
        assert_eq!(
            best_available_mix(
                Bar::new(1.6).to_pa(),
                msw::new(22.0).to_pa(),
                &gases,
                &empty_tissues,
                true,
                &GasDensitySettings::Ignore
            )
            .expect("There are gases, so reduce should return a result"),
            (1_usize, &gases[1])
        );
    }
}

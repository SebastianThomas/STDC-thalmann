use core::iter::zip;

#[allow(unused)]
use num::Float;

use crate::pressure_unit::{AbsPressure, Bar, Pressure};

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

#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct gL(f32);

impl gL {
    pub const fn new(v: f32) -> Self {
        Self(v)
    }

    pub const fn to_f32(self) -> f32 {
        self.0
    }
}

const impl core::ops::Add for gL {
    type Output = gL;

    fn add(self, rhs: gL) -> gL {
        gL(self.0 + rhs.0)
    }
}

const impl core::ops::Mul<f32> for gL {
    type Output = gL;

    fn mul(self, rhs: f32) -> gL {
        gL(self.0 * rhs)
    }
}

pub const MAX_GAS_DENSITY: gL = gL::new(5.2);
pub const MAX_GAS_DENSITY_LIMIT: gL = gL::new(6.2);

pub const fn gas_density_limit_from_air_multiplier(air_multiplier: f32) -> gL {
    gL::new(DENSITY_AIR * air_multiplier)
}

pub fn air_density_at_depth<P: const AbsPressure>(depth: P) -> gL {
    gL::new(depth.to_bar().to_f32() * DENSITY_AIR)
}

pub fn gas_density_limit_at_depth_from_air_multiplier<P: const AbsPressure>(
    depth: P,
    air_multiplier: f32,
) -> gL {
    air_density_at_depth(depth) * air_multiplier
}

pub const trait Gas {
    fn po2<D: const AbsPressure>(&self, depth: D) -> D;
    fn pn2<D: const AbsPressure>(&self, depth: D) -> D;
    fn phe<D: const AbsPressure>(&self, depth: D) -> D;
    fn ph2<D: const AbsPressure>(&self, depth: D) -> D;
    fn pn2_phe_ph2<D: const AbsPressure>(&self, depth: D) -> (D, D, D);

    fn o2_density<P: const AbsPressure>(&self, depth: P) -> gL {
        gL::new(self.po2(depth).to_bar().to_f32() * DENSITY_O2)
    }

    fn n2_density<P: const AbsPressure>(&self, depth: P) -> gL {
        gL::new(self.pn2(depth).to_bar().to_f32() * DENSITY_N2)
    }

    fn he_density<P: const AbsPressure>(&self, depth: P) -> gL {
        gL::new(self.phe(depth).to_bar().to_f32() * DENSITY_HE)
    }

    fn h2_density<P: const AbsPressure>(&self, depth: P) -> gL {
        gL::new(self.ph2(depth).to_bar().to_f32() * DENSITY_H2)
    }

    fn gas_density<P: const AbsPressure>(&self, depth: P) -> gL {
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

const impl Gas for GasMix<f32> {
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

const impl<P: const AbsPressure> Gas for CCRGas<f32, P> {
    fn po2<D: const AbsPressure>(&self, depth: D) -> D {
        let set_point = D::from(self.set_point.to_pa());
        if depth < set_point { depth } else { set_point }
    }

    fn pn2_phe_ph2<D: const AbsPressure>(&self, depth: D) -> (D, D, D) {
        let po2 = self.po2(depth);
        let fo2_loop = po2 / depth;
        let fo2_dil = self.diluent.fo2();
        let fo2_from_dil = (fo2_loop - fo2_dil) / (1.0 - fo2_dil);
        let fn2_loop = (1.0 - fo2_from_dil) * self.diluent.fn2();
        let fhe_loop = (1.0 - fo2_from_dil) * self.diluent.fhe();
        let fh2_loop = (1.0 - fo2_from_dil) * self.diluent.fh2();
        (depth * fn2_loop, depth * fhe_loop, depth * fh2_loop)
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

pub const AIR: GasMix<f32> = match GasMix::new(0.209, 0.000_005_2) {
    Ok(g) => g,
    Err(_) => unreachable!(),
};
pub const NX100: GasMix<f32> = match GasMix::new(0.99, 0.0) {
    Ok(g) => g,
    Err(_) => unreachable!(),
};
pub const NX50: GasMix<f32> = match GasMix::new(0.50, 0.000_005_2) {
    Ok(g) => g,
    Err(_) => unreachable!(),
};
pub const TMX21_35: GasMix<f32> = match GasMix::new(0.21, 0.35) {
    Ok(g) => g,
    Err(_) => unreachable!(),
};
pub const TMX18_45: GasMix<f32> = match GasMix::new(0.18, 0.45) {
    Ok(g) => g,
    Err(_) => unreachable!(),
};
pub const TMX15_55: GasMix<f32> = match GasMix::new(0.15, 0.55) {
    Ok(g) => g,
    Err(_) => unreachable!(),
};
pub const TMX12_65: GasMix<f32> = match GasMix::new(0.12, 0.65) {
    Ok(g) => g,
    Err(_) => unreachable!(),
};
pub const TMX10_80: GasMix<f32> = match GasMix::new(0.10, 0.80) {
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
        zip(self.n2, self.he).any(|(n2, he)| n2 < new_gas_n2 && he > new_gas_he)
    }

    pub fn tick<G: Gas>(&mut self, time_delta_ms: u16, depth: P, gas: &G) {
        Self::tick_gas(time_delta_ms, gas.pn2(depth), &mut self.n2);
        Self::tick_gas(time_delta_ms, gas.phe(depth), &mut self.he);
    }

    fn tick_gas(time_delta_ms: u16, pp_insp: P, cur: &mut [P; NUM_TS]) {
        for item in cur.iter_mut().take(NUM_TS) {
            let delta: f32 = f32::from(time_delta_ms) / 1000.0;
            *item = *item + (pp_insp - *item) * delta;
        }
    }
}

pub fn best_mix_fo2<P: const AbsPressure>(max_po2: P, depth: P) -> f32 {
    max_po2 / depth
}

pub enum GasDensitySettings {
    Ignore,
    Limit { limit_g_l: gL },
}

impl GasDensitySettings {
    pub const fn limit_g_l(limit_g_l: gL) -> Self {
        Self::Limit { limit_g_l }
    }

    pub const fn limit_from_air_multiplier(air_multiplier: f32) -> Self {
        Self::Limit {
            limit_g_l: gas_density_limit_from_air_multiplier(air_multiplier),
        }
    }

    pub fn no_violation<P: const AbsPressure>(&self, depth: P, gas: &GasMix<f32>) -> bool {
        if let GasDensitySettings::Limit { limit_g_l } = self {
            return gas.gas_density(depth) < *limit_g_l;
        }
        true
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
    gases_enabled: &[bool; G],
    tissue_loading: &TissuesLoading<NUM_TS, P>,
    ignore_isobaric_counterdiffusion: bool,
    gas_density: &GasDensitySettings,
) -> Option<(usize, &'a GasMix<f32>)> {
    let best_mix_fo2 = best_mix_fo2(max_po2, depth);
    available_gases
        .iter()
        .enumerate()
        .filter(|(i, _g)| gases_enabled[*i])
        .filter(|(_i, g)| g.fo2() <= best_mix_fo2)
        .filter(|(_i, g)| {
            ignore_isobaric_counterdiffusion
                || !tissue_loading.is_isobaric_counterdiffusion(depth, *g)
        })
        .filter(|(_i, g)| gas_density.no_violation(depth, g))
        .reduce(|(ai, ag), (bi, bg)| {
            let better_fo2 = ag.fo2() > bg.fo2();
            let same_fo2_better_he = ag.fo2() == bg.fo2() && ag.fhe() > bg.fhe();
            if better_fo2 || same_fo2_better_he {
                (ai, ag)
            } else {
                (bi, bg)
            }
        })
}

#[cfg(test)]
mod tests {
    use crate::pressure_unit::{Pa, Pressure, msw};

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

    fn best_available_mix_fixture() -> ([GasMix<f32>; 4], [bool; 4], TissuesLoading<1, Pa>) {
        let gases = [
            AIR,
            GasMix::new(0.21, 0.35).expect("21 + 35 < 100"),
            GasMix::new(0.5, 0.0).expect("50 < 100"),
            GasMix::new(0.10, 0.80).expect("10 + 80 < 100"),
        ];
        let gases_enabled = [true; 4];
        let empty_tissues = TissuesLoading {
            n2: [msw::new(0.0).to_pa()],
            he: [msw::new(0.0).to_pa()],
        };
        (gases, gases_enabled, empty_tissues)
    }

    #[test]
    fn best_available_mix_at_21m_ppo2_1_6_selects_gas_2() {
        let (gases, gases_enabled, empty_tissues) = best_available_mix_fixture();
        assert_eq!(
            best_available_mix(
                Bar::new(1.6).to_pa(),
                msw::new(21.0).to_pa(),
                &gases,
                &gases_enabled,
                &empty_tissues,
                true,
                &GasDensitySettings::Ignore
            )
            .expect("There are gases, so reduce should return a result"),
            (2_usize, &gases[2])
        );
    }

    #[test]
    fn best_available_mix_at_21m_ppo2_1_4_selects_gas_1() {
        let (gases, gases_enabled, empty_tissues) = best_available_mix_fixture();
        assert_eq!(
            best_available_mix(
                Bar::new(1.4).to_pa(),
                msw::new(21.0).to_pa(),
                &gases,
                &gases_enabled,
                &empty_tissues,
                true,
                &GasDensitySettings::Ignore
            )
            .expect("There are gases, so reduce should return a result"),
            (1_usize, &gases[1])
        );
    }

    #[test]
    fn best_available_mix_at_22m_ppo2_1_6_selects_gas_1() {
        let (gases, gases_enabled, empty_tissues) = best_available_mix_fixture();
        assert_eq!(
            best_available_mix(
                Bar::new(1.6).to_pa(),
                msw::new(22.0).to_pa(),
                &gases,
                &gases_enabled,
                &empty_tissues,
                true,
                &GasDensitySettings::Ignore
            )
            .expect("There are gases, so reduce should return a result"),
            (1_usize, &gases[1])
        );
    }

    #[test]
    fn best_available_mix_at_90m_ppo2_1_4_selects_gas_3() {
        let (gases, gases_enabled, empty_tissues) = best_available_mix_fixture();
        assert_eq!(
            best_available_mix(
                Bar::new(1.4).to_pa(),
                msw::new(90.0).to_pa(),
                &gases,
                &gases_enabled,
                &empty_tissues,
                true,
                &GasDensitySettings::Ignore
            ),
            Some((3_usize, &gases[3]))
        );
    }

    #[test]
    fn best_available_mix_with_air_multiplier_limit_returns_none() {
        let (gases, gases_enabled, empty_tissues) = best_available_mix_fixture();
        assert_eq!(
            best_available_mix(
                Bar::new(1.4).to_pa(),
                msw::new(130.0).to_pa(),
                &gases,
                &gases_enabled,
                &empty_tissues,
                true,
                &GasDensitySettings::limit_from_air_multiplier(3.0)
            ),
            None,
            // Some((3_usize, &gases[3]))
        );
    }

    #[test]
    fn best_available_mix_with_6gl_limit_selects_gas_3() {
        let (gases, gases_enabled, empty_tissues) = best_available_mix_fixture();
        assert_eq!(
            best_available_mix(
                Bar::new(1.4).to_pa(),
                msw::new(119.0).to_pa(),
                &gases,
                &gases_enabled,
                &empty_tissues,
                true,
                &GasDensitySettings::limit_g_l(gL::new(6.3))
            ),
            Some((3_usize, &gases[3]))
        );
    }

    #[test]
    fn best_available_mix_with_2gl_limit_selects_gas_3() {
        let (gases, gases_enabled, empty_tissues) = best_available_mix_fixture();
        assert_eq!(
            best_available_mix(
                Bar::new(1.4).to_pa(),
                msw::new(119.0).to_pa(),
                &gases,
                &gases_enabled,
                &empty_tissues,
                true,
                &GasDensitySettings::limit_g_l(gL::new(2.0))
            ),
            None,
        );
    }

    #[test]
    fn best_available_mix_no_valid_gas_returns_none() {
        let (gases, gases_enabled, empty_tissues) = best_available_mix_fixture();
        assert_eq!(
            best_available_mix(
                Bar::new(1.0).to_pa(),
                msw::new(200.0).to_pa(),
                &gases,
                &gases_enabled,
                &empty_tissues,
                true,
                &GasDensitySettings::Ignore
            ),
            None
        );
    }

    #[test]
    fn is_isobaric_counterdiffusion_true() {
        let depth = msw::new(30.0).to_pa();
        let tissues: TissuesLoading<1, Pa> = TissuesLoading {
            n2: [depth * 0.1],
            he: [depth * 0.7],
        };
        // AIR has relatively high N2 and negligible He compared to the tissue above
        let new_gas = AIR;
        assert_eq!(tissues.is_isobaric_counterdiffusion(depth, &new_gas), true);
    }

    #[test]
    fn is_isobaric_counterdiffusion_false_when_no_match() {
        let depth = msw::new(30.0).to_pa();
        let tissues: TissuesLoading<1, Pa> = TissuesLoading {
            n2: [depth * 0.5],
            he: [depth * 0.1],
        };
        // TMX10_80 is helium rich; this should not trigger the check (he > new_he false)
        let new_gas = TMX10_80;
        assert_eq!(tissues.is_isobaric_counterdiffusion(depth, &new_gas), false);
    }

    #[test]
    fn is_isobaric_counterdiffusion_any_tissue_true() {
        let depth = msw::new(30.0).to_pa();
        let tissues: TissuesLoading<2, Pa> = TissuesLoading {
            n2: [depth * 0.5, depth * 0.1],
            he: [depth * 0.1, depth * 0.7],
        };
        let new_gas = AIR;
        // second tissue should trigger the condition
        assert_eq!(tissues.is_isobaric_counterdiffusion(depth, &new_gas), true);
    }

    #[test]
    fn gas_density_no_violation_and_violation() {
        let depth = msw::new(50.0).to_pa();
        // TMX10_80 is helium rich and generally light
        let light_gas = TMX10_80;
        let settings_ok = GasDensitySettings::limit_g_l(gL::new(10.0));
        assert_eq!(settings_ok.no_violation(depth, &light_gas), true);

        let deep = msw::new(100.0).to_pa();
        let heavy_gas = AIR;
        let settings_strict = GasDensitySettings::limit_g_l(gL::new(1.0));
        assert_eq!(settings_strict.no_violation(deep, &heavy_gas), false);
    }
}

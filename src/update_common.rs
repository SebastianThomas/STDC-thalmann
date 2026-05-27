#[cfg(feature = "lin_exp")]
use core::f32::consts::LN_2;

#[cfg(not(test))]
use num::Float;

#[cfg(feature = "lin_exp")]
use crate::mptt::Tissue;
use crate::pressure_unit::AbsPressure;
#[cfg(not(feature = "lin_exp"))]
use crate::pressure_unit::{Pa, Pressure};

pub(crate) fn exp_pressure<P: const AbsPressure>(p_inspired: P, p_old: P, k: f32, t: f32) -> P {
    let exp = (-k * t).exp();
    p_inspired + (p_old - p_inspired) * exp
}

/// Return (KSAT_array, KDSAT_array) for the given tissues.
#[cfg(feature = "lin_exp")]
pub(crate) fn ks_arrays<const N: usize>(tissues: &[Tissue; N]) -> ([f32; N], [f32; N]) {
    let ks_sat: [f32; N] = tissues.map(|t| LN_2 / t.half_time);
    let ks_desat: [f32; N] = tissues.map(|t| (LN_2 / t.half_time) * t.sdr);
    (ks_sat, ks_desat)
}

/// Compute mixed Buehlmann a/b based M-value for a tissue using current
/// tissue partial pressures as weights. Falls back to zero if total inert is zero.
#[cfg(not(feature = "lin_exp"))]
pub(crate) fn mixed_buehlmann_mvalue(
    tissue_idx: usize,
    p_n2: Pa,
    p_he: Pa,
    stop_depth_pa: Pa,
) -> Pa {
    use crate::mptt_buehlmann::TISSUES as BUEHL_TISSUES;
    let total = p_n2 + p_he;
    if total.to_f32() <= 0.0 {
        return BUEHL_TISSUES[tissue_idx].n2.a.to_pa()
            + stop_depth_pa * BUEHL_TISSUES[tissue_idx].n2.b;
    }
    let a_n2 = BUEHL_TISSUES[tissue_idx].n2.a.to_pa();
    let a_he = BUEHL_TISSUES[tissue_idx].he.a.to_pa();
    let b_n2 = BUEHL_TISSUES[tissue_idx].n2.b;
    let b_he = BUEHL_TISSUES[tissue_idx].he.b;

    let n2_frac = p_n2 / total;
    let he_frac = p_he / total;
    let a_mix = a_n2 * n2_frac + a_he * he_frac;
    let b_mix = b_n2 * n2_frac + b_he * he_frac;
    a_mix + stop_depth_pa * b_mix
}

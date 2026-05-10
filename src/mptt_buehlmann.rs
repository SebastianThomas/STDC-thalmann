use crate::{
    mptt::{MValues, TissueRow},
    pressure_unit::{AbsPressure, Pa, Pressure, msw},
};

pub const NUM_TISSUES_BUEHLMANN: usize = 16;
pub const NUM_STOP_DEPTHS_BUEHLMANN: usize = 32;
pub const DEFAULT_BUEHLMANN_HE_RATIO: f32 = 0.0;

pub struct BuehlmannTissue {
    n2: BuehlmannTissueGas,
    he: BuehlmannTissueGas,
}
pub struct BuehlmannTissueGas {
    half_time: f32,
    a: msw,
    b: f32,
}

pub const TISSUES: [BuehlmannTissue; NUM_TISSUES_BUEHLMANN] = [
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 4.0,
            a: msw::new(1.1696),
            b: 0.5578,
        },
        he: BuehlmannTissueGas {
            half_time: 4.0,
            a: msw::new(1.6189),
            b: 0.4770,
        },
    },
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 8.0,
            a: msw::new(1.0000),
            b: 0.6514,
        },
        he: BuehlmannTissueGas {
            half_time: 8.0,
            a: msw::new(1.3830),
            b: 0.5747,
        },
    },
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 12.5,
            a: msw::new(0.8618),
            b: 0.7222,
        },
        he: BuehlmannTissueGas {
            half_time: 12.5,
            a: msw::new(1.1919),
            b: 0.6527,
        },
    },
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 18.5,
            a: msw::new(0.7562),
            b: 0.7825,
        },
        he: BuehlmannTissueGas {
            half_time: 18.5,
            a: msw::new(1.0458),
            b: 0.7223,
        },
    },
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 27.0,
            a: msw::new(0.6200),
            b: 0.8126,
        },
        he: BuehlmannTissueGas {
            half_time: 27.0,
            a: msw::new(0.9220),
            b: 0.7582,
        },
    },
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 38.3,
            a: msw::new(0.5034),
            b: 0.8434,
        },
        he: BuehlmannTissueGas {
            half_time: 38.3,
            a: msw::new(0.8205),
            b: 0.7957,
        },
    },
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 54.3,
            a: msw::new(0.4410),
            b: 0.8693,
        },
        he: BuehlmannTissueGas {
            half_time: 54.3,
            a: msw::new(0.7305),
            b: 0.8279,
        },
    },
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 77.0,
            a: msw::new(0.4000),
            b: 0.8910,
        },
        he: BuehlmannTissueGas {
            half_time: 77.0,
            a: msw::new(0.6502),
            b: 0.8553,
        },
    },
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 109.0,
            a: msw::new(0.3750),
            b: 0.9092,
        },
        he: BuehlmannTissueGas {
            half_time: 109.0,
            a: msw::new(0.5950),
            b: 0.8757,
        },
    },
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 146.0,
            a: msw::new(0.3500),
            b: 0.9222,
        },
        he: BuehlmannTissueGas {
            half_time: 146.0,
            a: msw::new(0.5545),
            b: 0.8903,
        },
    },
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 187.0,
            a: msw::new(0.3295),
            b: 0.9319,
        },
        he: BuehlmannTissueGas {
            half_time: 187.0,
            a: msw::new(0.5333),
            b: 0.8997,
        },
    },
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 239.0,
            a: msw::new(0.3065),
            b: 0.9403,
        },
        he: BuehlmannTissueGas {
            half_time: 239.0,
            a: msw::new(0.5189),
            b: 0.9073,
        },
    },
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 305.0,
            a: msw::new(0.2835),
            b: 0.9477,
        },
        he: BuehlmannTissueGas {
            half_time: 305.0,
            a: msw::new(0.5181),
            b: 0.9122,
        },
    },
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 390.0,
            a: msw::new(0.2610),
            b: 0.9544,
        },
        he: BuehlmannTissueGas {
            half_time: 390.0,
            a: msw::new(0.5176),
            b: 0.9171,
        },
    },
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 498.0,
            a: msw::new(0.2480),
            b: 0.9602,
        },
        he: BuehlmannTissueGas {
            half_time: 498.0,
            a: msw::new(0.5172),
            b: 0.9217,
        },
    },
    BuehlmannTissue {
        n2: BuehlmannTissueGas {
            half_time: 635.0,
            a: msw::new(0.2327),
            b: 0.9653,
        },
        he: BuehlmannTissueGas {
            half_time: 635.0,
            a: msw::new(0.5119),
            b: 0.9267,
        },
    },
];

#[derive(Clone, Copy)]
pub struct TissueRowBuehlmann<const NUM_TISSUES: usize, P: const AbsPressure> {
    depth: msw,
    total: TissueRow<NUM_TISSUES, P>,
    he: TissueRow<NUM_TISSUES, P>,
    n2: TissueRow<NUM_TISSUES, P>,
}

const fn buehlmann_16c_mvalues_table(
    he_ratio: Option<f32>,
) -> MValues<Pa, { NUM_TISSUES_BUEHLMANN }, { NUM_STOP_DEPTHS_BUEHLMANN }> {
    let mut result: [TissueRowBuehlmann<NUM_TISSUES_BUEHLMANN, Pa>; NUM_STOP_DEPTHS_BUEHLMANN] =
        [TissueRowBuehlmann {
            depth: msw(-1.0),
            total: TissueRow::empty_pa(),
            n2: TissueRow::empty_pa(),
            he: TissueRow::empty_pa(),
        }; NUM_STOP_DEPTHS_BUEHLMANN];
    let mut i: usize = 0;
    let mut depth: f32 = 0.0;
    while i < NUM_STOP_DEPTHS_BUEHLMANN {
        let d = msw(depth);
        result[i] = buehlmann_16c_depth(d, he_ratio);

        i += 1;
        depth += 3.0;
    }
    const fn get_total(
        b: TissueRowBuehlmann<NUM_TISSUES_BUEHLMANN, Pa>,
    ) -> TissueRow<{ NUM_TISSUES_BUEHLMANN }, Pa> {
        b.total
    }
    result.map(get_total)
}

const fn buehlmann_16c_depth(
    depth: msw,
    he_ratio: Option<f32>,
) -> TissueRowBuehlmann<{ NUM_TISSUES_BUEHLMANN }, Pa> {
    let r = he_ratio.unwrap_or(0.0);
    let mut max_saturation = [Pa::new(0.0); NUM_TISSUES_BUEHLMANN];
    let mut max_saturation_n2 = [Pa::new(0.0); NUM_TISSUES_BUEHLMANN];
    let mut max_saturation_he = [Pa::new(0.0); NUM_TISSUES_BUEHLMANN];
    let mut i = 0;
    let p_amb = depth.to_pa();
    while i < NUM_TISSUES_BUEHLMANN {
        let tissue = &TISSUES[i];
        (
            max_saturation[i],
            max_saturation_n2[i],
            max_saturation_he[i],
        ) = buehlmann_saturation_tolerance(p_amb, r, tissue);
        i += 1;
    }
    TissueRowBuehlmann {
        depth,
        total: TissueRow {
            depth,
            max_saturation,
        },
        n2: TissueRow {
            depth,
            max_saturation: max_saturation_n2,
        },
        he: TissueRow {
            depth,
            max_saturation: max_saturation_he,
        },
    }
}

/** P_{igtol} = Inert Gas Tolerance */
const fn buehlmann_saturation_tolerance(
    p_amb: Pa,
    r: f32,
    tissue: &BuehlmannTissue,
) -> (Pa, Pa, Pa) {
    let a_n2_pa = tissue.n2.a.to_pa();
    let a_he_pa = tissue.he.a.to_pa();
    let b_n2 = tissue.n2.b;
    let b_he = tissue.he.b;
    let a_combined = a_n2_pa * (1.0 - r) + a_he_pa * r;
    let b_combined = b_n2 * (1.0 - r) + b_he * r;
    (
        a_combined + p_amb / b_combined,
        tissue.n2.a.to_pa() + p_amb / tissue.n2.b,
        tissue.he.a.to_pa() + p_amb / tissue.he.b,
    )
}

pub const BUEHLMANN_16C: MValues<Pa, { NUM_TISSUES_BUEHLMANN }, { NUM_STOP_DEPTHS_BUEHLMANN }> =
    buehlmann_16c_mvalues_table(Some(DEFAULT_BUEHLMANN_HE_RATIO));

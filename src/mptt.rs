use crate::pressure_unit::{fsw, msw, Pa, Pressure};

#[derive(Copy, Clone)]
pub struct Tissue {
    pub half_time: f32,
    pub sdr: f32, /* Saturation Desaturation Ration */
}

#[derive(Copy, Clone)]
pub struct TissueRow<const TISSUES: usize> {
    pub depth: msw,
    pub max_saturation: [Pa; TISSUES],
}

pub const NUM_TISSUES: usize = 5;
pub const NUM_STOP_DEPTHS: usize = 32;
pub type MVALUES = [TissueRow<NUM_TISSUES>; NUM_STOP_DEPTHS];

pub const TISSUES: [Tissue; NUM_TISSUES] = [
    Tissue {
        half_time: 10.0,
        sdr: 1.0,
    },
    Tissue {
        half_time: 20.0,
        sdr: 2.0,
    },
    Tissue {
        half_time: 20.0,
        sdr: 0.67,
    },
    Tissue {
        half_time: 120.0,
        sdr: 1.0,
    },
    Tissue {
        half_time: 200.0,
        sdr: 1.0,
    },
];

// TODO: Maybe allow f64?
pub const XVAL_HE9_040_F32: MVALUES = [
    // XVAL-HE-9_040 (MSW)
    // Half-times (mins)
    //      10 20 20 120 200
    // Stop Depth SDR
    // (msw) 1 2 0.67 1 1
    TissueRow {
        depth: msw(3.0),
        max_saturation: [
            fsw(85.000).to_pa(),
            fsw(64.000).to_pa(),
            fsw(83.000).to_pa(),
            fsw(41.731).to_pa(),
            fsw(38.274).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(6.0),
        max_saturation: [
            fsw(94.843).to_pa(),
            fsw(73.843).to_pa(),
            fsw(92.843).to_pa(),
            fsw(61.416).to_pa(),
            fsw(49.969).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(9.0),
        max_saturation: [
            fsw(104.685).to_pa(),
            fsw(83.685).to_pa(),
            fsw(102.685).to_pa(),
            fsw(81.101).to_pa(),
            fsw(61.664).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(12.0),
        max_saturation: [
            fsw(114.528).to_pa(),
            fsw(93.528).to_pa(),
            fsw(112.528).to_pa(),
            fsw(100.786).to_pa(),
            fsw(73.359).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(15.0),
        max_saturation: [
            fsw(124.370).to_pa(),
            fsw(103.370).to_pa(),
            fsw(122.370).to_pa(),
            fsw(120.471).to_pa(),
            fsw(85.054).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(18.0),
        max_saturation: [
            fsw(134.213).to_pa(),
            fsw(113.213).to_pa(),
            fsw(132.213).to_pa(),
            fsw(140.156).to_pa(),
            fsw(96.749).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(21.0),
        max_saturation: [
            fsw(144.055).to_pa(),
            fsw(123.055).to_pa(),
            fsw(142.055).to_pa(),
            fsw(159.841).to_pa(),
            fsw(108.444).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(24.0),
        max_saturation: [
            fsw(153.898).to_pa(),
            fsw(132.898).to_pa(),
            fsw(151.898).to_pa(),
            fsw(179.526).to_pa(),
            fsw(120.139).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(27.0),
        max_saturation: [
            fsw(163.740).to_pa(),
            fsw(142.740).to_pa(),
            fsw(161.740).to_pa(),
            fsw(199.211).to_pa(),
            fsw(131.834).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(30.0),
        max_saturation: [
            fsw(173.583).to_pa(),
            fsw(152.583).to_pa(),
            fsw(171.583).to_pa(),
            fsw(218.896).to_pa(),
            fsw(143.529).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(33.0),
        max_saturation: [
            fsw(183.425).to_pa(),
            fsw(162.425).to_pa(),
            fsw(181.425).to_pa(),
            fsw(238.581).to_pa(),
            fsw(155.224).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(36.0),
        max_saturation: [
            fsw(193.268).to_pa(),
            fsw(172.268).to_pa(),
            fsw(191.268).to_pa(),
            fsw(258.266).to_pa(),
            fsw(166.919).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(39.0),
        max_saturation: [
            fsw(203.110).to_pa(),
            fsw(182.110).to_pa(),
            fsw(201.110).to_pa(),
            fsw(277.951).to_pa(),
            fsw(178.614).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(42.0),
        max_saturation: [
            fsw(212.953).to_pa(),
            fsw(191.953).to_pa(),
            fsw(210.953).to_pa(),
            fsw(297.637).to_pa(),
            fsw(190.309).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(45.0),
        max_saturation: [
            fsw(222.795).to_pa(),
            fsw(201.795).to_pa(),
            fsw(220.795).to_pa(),
            fsw(317.322).to_pa(),
            fsw(202.004).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(48.0),
        max_saturation: [
            fsw(232.638).to_pa(),
            fsw(211.638).to_pa(),
            fsw(230.638).to_pa(),
            fsw(337.007).to_pa(),
            fsw(213.699).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(51.0),
        max_saturation: [
            fsw(242.480).to_pa(),
            fsw(221.480).to_pa(),
            fsw(240.480).to_pa(),
            fsw(356.692).to_pa(),
            fsw(225.394).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(54.0),
        max_saturation: [
            fsw(252.323).to_pa(),
            fsw(231.323).to_pa(),
            fsw(250.323).to_pa(),
            fsw(376.377).to_pa(),
            fsw(237.089).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(57.0),
        max_saturation: [
            fsw(262.165).to_pa(),
            fsw(241.165).to_pa(),
            fsw(260.165).to_pa(),
            fsw(396.062).to_pa(),
            fsw(248.784).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(60.0),
        max_saturation: [
            fsw(272.008).to_pa(),
            fsw(251.008).to_pa(),
            fsw(270.008).to_pa(),
            fsw(415.747).to_pa(),
            fsw(260.479).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(63.0),
        max_saturation: [
            fsw(281.850).to_pa(),
            fsw(260.850).to_pa(),
            fsw(279.850).to_pa(),
            fsw(435.432).to_pa(),
            fsw(272.173).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(66.0),
        max_saturation: [
            fsw(291.693).to_pa(),
            fsw(270.693).to_pa(),
            fsw(289.693).to_pa(),
            fsw(455.117).to_pa(),
            fsw(283.868).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(69.0),
        max_saturation: [
            fsw(301.535).to_pa(),
            fsw(280.535).to_pa(),
            fsw(299.535).to_pa(),
            fsw(474.802).to_pa(),
            fsw(295.563).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(72.0),
        max_saturation: [
            fsw(311.378).to_pa(),
            fsw(290.378).to_pa(),
            fsw(309.378).to_pa(),
            fsw(494.487).to_pa(),
            fsw(307.258).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(75.0),
        max_saturation: [
            fsw(321.220).to_pa(),
            fsw(300.220).to_pa(),
            fsw(319.220).to_pa(),
            fsw(514.172).to_pa(),
            fsw(318.953).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(78.0),
        max_saturation: [
            fsw(331.063).to_pa(),
            fsw(310.063).to_pa(),
            fsw(329.063).to_pa(),
            fsw(533.857).to_pa(),
            fsw(330.648).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(81.0),
        max_saturation: [
            fsw(340.906).to_pa(),
            fsw(319.906).to_pa(),
            fsw(338.906).to_pa(),
            fsw(553.542).to_pa(),
            fsw(342.343).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(84.0),
        max_saturation: [
            fsw(350.748).to_pa(),
            fsw(329.748).to_pa(),
            fsw(348.748).to_pa(),
            fsw(573.227).to_pa(),
            fsw(354.038).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(87.0),
        max_saturation: [
            fsw(360.591).to_pa(),
            fsw(339.591).to_pa(),
            fsw(358.591).to_pa(),
            fsw(592.912).to_pa(),
            fsw(365.733).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(90.0),
        max_saturation: [
            fsw(370.433).to_pa(),
            fsw(349.433).to_pa(),
            fsw(368.433).to_pa(),
            fsw(612.597).to_pa(),
            fsw(377.428).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(93.0),
        max_saturation: [
            fsw(380.276).to_pa(),
            fsw(359.276).to_pa(),
            fsw(378.276).to_pa(),
            fsw(632.282).to_pa(),
            fsw(389.123).to_pa(),
        ],
    },
    TissueRow {
        depth: msw(96.0),
        max_saturation: [
            fsw(390.118).to_pa(),
            fsw(369.118).to_pa(),
            fsw(388.118).to_pa(),
            fsw(651.967).to_pa(),
            fsw(400.818).to_pa(),
        ],
    },
];

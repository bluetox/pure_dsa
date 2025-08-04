use crate::polyvec::{
    Polyveck, 
    PolyveckStruct, 
    Polyvecl, 
    PolyveclStruct
};
pub const SHAKE256_RATE: usize = 136;
pub const STREAM256_BLOCKBYTES: usize = SHAKE256_RATE;

pub enum PolyUniformGamma1Buffer {
    Mode2(PolyUniformGamma1BufferStruct<{(576 + 136 - 1) / 136 * 136}>),
    Mode3(PolyUniformGamma1BufferStruct<{(640 + 136 - 1) / 136 * 136}>),
    Mode5(PolyUniformGamma1BufferStruct<{(640 + 136 - 1) / 136 * 136}>),
}

impl PolyUniformGamma1Buffer {
    pub fn buf(&mut self) -> &mut [u8] {
        match self {
            PolyUniformGamma1Buffer::Mode2(v) => &mut v.buf[..],
            PolyUniformGamma1Buffer::Mode3(v) => &mut v.buf[..],
            PolyUniformGamma1Buffer::Mode5(v) => &mut v.buf[..],
        }
    }
}
pub(crate) struct PolyUniformGamma1BufferStruct<const N: usize> {
    pub buf: [u8; N],
}

#[derive(Clone)]
pub enum Mat {
    Mode2([Polyvecl; 4]),
    Mode3([Polyvecl; 6]),
    Mode5([Polyvecl; 8]),
}




pub trait DilithiumParams {
    const SEEDBYTES: usize = 32;
    const CRHBYTES: usize = 64;
    const TRBYTES: usize = 64;
    const RNDBYTES: usize = 32;
    const N: usize = 256;
    const Q: usize = 8380417;
    const D: usize = 13;
    const ROOT_OF_UNITY: usize = 1753;

    const K: usize;
    const L: usize;
    const ETA: usize;
    const TAU: usize;
    const BETA: usize;
    const GAMMA1: usize;
    const GAMMA2: usize;
    const OMEGA: usize;
    const CTILDEBYTES: usize;

    const POLYT1_PACKEDBYTES: usize = 320;
    const POLYT0_PACKEDBYTES: usize = 416;
    const POLYETA_PACKEDBYTES: usize;
    const POLYZ_PACKEDBYTES: usize;
    const POLYVECH_PACKEDBYTES: usize = Self::OMEGA + Self::K;
    const SIGNBYTES: usize =
        Self::SEEDBYTES + Self::L * Self::POLYZ_PACKEDBYTES + Self::POLYVECH_PACKEDBYTES;
    const PUBLIC_KEY_BYTES: usize = Self::SEEDBYTES + Self::K * Self::POLYT1_PACKEDBYTES;
    const SECRET_KEY_BYTES: usize = 2 * Self::SEEDBYTES
            + Self::TRBYTES
            + Self::L * Self::POLYETA_PACKEDBYTES
            + Self::K * Self::POLYETA_PACKEDBYTES
            + Self::K * Self::POLYT0_PACKEDBYTES;
    const POLY_UNIFORM_GAMMA1_NBLOCKS: usize = (Self::POLYZ_PACKEDBYTES + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES;
    fn polyvecknew() -> Polyveck;
    fn polyveclnew() -> Polyvecl;
     fn poly_uniform_gamma1_buffer() -> PolyUniformGamma1Buffer {
        match Self::K {
            4 => PolyUniformGamma1Buffer::Mode2(PolyUniformGamma1BufferStruct { buf: [0u8; { (576 + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES * STREAM256_BLOCKBYTES }] }),
            6 => PolyUniformGamma1Buffer::Mode3(PolyUniformGamma1BufferStruct { buf: [0u8; { (640 + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES * STREAM256_BLOCKBYTES}] }),
            8 => PolyUniformGamma1Buffer::Mode5(PolyUniformGamma1BufferStruct { buf: [0u8; { (640 + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES * STREAM256_BLOCKBYTES}] }),
            _ => panic!("Invalid GAMMA1 or POLYZ_PACKEDBYTES for gamma1 buffer"),
        }
    }
    fn mat() -> Mat {
    match Self::K {
        4 => Mat::Mode2(core::array::from_fn(|_| Self::polyveclnew())),
        6 => Mat::Mode3(core::array::from_fn(|_| Self::polyveclnew())),
        8 => Mat::Mode5(core::array::from_fn(|_| Self::polyveclnew())),
        _ => panic!("Invalid K value"),
    }
}

    fn polyz_packedbytes() -> usize {
        if Self::GAMMA1 == (1 << 17) {
            576
        } else if Self::GAMMA1 == (1 << 19) {
            640
        } else {
            panic!("Invalid GAMMA1");
        }
    }

    fn polyw1_packedbytes() -> usize {
        if Self::GAMMA2 == (Self::Q - 1) / 88 {
            192
        } else if Self::GAMMA2 == (Self::Q - 1) / 32 {
            128
        } else {
            panic!("Invalid GAMMA2");
        }
    }

    fn crypto_publickeybytes() -> usize {
        Self::SEEDBYTES + Self::K * Self::POLYT1_PACKEDBYTES
    }
}

pub struct Mode2;
impl DilithiumParams for Mode2 {
    const K: usize = 4;
    const L: usize = 4;
    const ETA: usize = 2;
    const POLYETA_PACKEDBYTES: usize = 96;
    const POLYZ_PACKEDBYTES: usize = 576;
    const TAU: usize = 39;
    const BETA: usize = 78;
    const GAMMA1: usize = 1 << 17;
    const GAMMA2: usize = (Self::Q - 1) / 88;
    const OMEGA: usize = 80;
    const CTILDEBYTES: usize = 32;
    fn polyvecknew() -> Polyveck {
        Polyveck::Mode2(PolyveckStruct::default())
    }
    fn polyveclnew() -> Polyvecl {
        Polyvecl::Mode2(PolyveclStruct::default())
    }
}

pub struct Mode3;
impl DilithiumParams for Mode3 {
    const K: usize = 6;
    const L: usize = 5;
    const ETA: usize = 4;
    const POLYETA_PACKEDBYTES: usize = 128;
    const POLYZ_PACKEDBYTES: usize = 640;
    const TAU: usize = 49;
    const BETA: usize = 196;
    const GAMMA1: usize = 1 << 19;
    const GAMMA2: usize = (Self::Q - 1) / 32;
    const OMEGA: usize = 55;
    const CTILDEBYTES: usize = 48;
        fn polyvecknew() -> Polyveck {
        Polyveck::Mode3(PolyveckStruct::default())
    }
    fn polyveclnew() -> Polyvecl {
        Polyvecl::Mode3(PolyveclStruct::default())
    }
}

pub struct Mode5;
impl DilithiumParams for Mode5 {
    const K: usize = 8;
    const L: usize = 7;
    const ETA: usize = 2;
    const POLYETA_PACKEDBYTES: usize = 96;
    const POLYZ_PACKEDBYTES: usize = 640;
    const TAU: usize = 60;
    const BETA: usize = 120;
    const GAMMA1: usize = 1 << 19;
    const GAMMA2: usize = (Self::Q - 1) / 32;
    const OMEGA: usize = 75;
    const CTILDEBYTES: usize = 64;
    fn polyvecknew() -> Polyveck {
        Polyveck::Mode5(PolyveckStruct::default())
    }
    fn polyveclnew() -> Polyvecl {
        Polyvecl::Mode5(PolyveclStruct::default())
    }
}
use crate::polyvec::{
    Polyveck, 
    PolyveckStruct, 
    Polyvecl, 
    PolyveclStruct
};

#[derive(Copy, Clone)]
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
    fn polyvecknew() -> Polyveck;
    fn polyveclnew() -> Polyvecl;

    fn mat() -> Mat {
        match Self::K {
            4 => Mat::Mode2([Self::polyveclnew(); 4]),
            6 => Mat::Mode3([Self::polyveclnew(); 6]),
            8 => Mat::Mode5([Self::polyveclnew(); 8]),
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

use crate::params::DilithiumParams;

pub fn power2round<P: DilithiumParams>(a0: &mut i32, a: i32) -> i32 {
    let a1: i32 = (a + (1 << (P::D - 1)) - 1) >> P::D;
    *a0 = a - (a1 << P::D);
    a1
}

pub fn decompose<P: DilithiumParams>(a0: &mut i32, a: i32) -> i32 {
    let mut a1 = (a + 127) >> 7;

    let q = P::Q as i32;
    let gamma2 = P::GAMMA2 as i32;

    if P::GAMMA2 == (P::Q - 1) / 32 {
        a1 = (a1 * 1025 + (1 << 21)) >> 22;
        a1 &= 15;
    } else if P::GAMMA2 == (P::Q - 1) / 88 {
        a1 = (a1 * 11275 + (1 << 23)) >> 24;
        a1 ^= ((43 - a1) >> 31) & a1;
    }

    *a0 = a - a1 * 2 * gamma2;
    *a0 -= (((q - 1) / 2 - *a0) >> 31) & q;

    a1
}

pub fn make_hint<P: DilithiumParams>(a0: i32, a1: i32) -> u32 {
    let gamma2 = P::GAMMA2 as i32;

    if a0 > gamma2 || a0 < -gamma2 || (a0 == -gamma2 && a1 != 0) {
        1
    } else {
        0
    }
}

pub fn use_hint<P: DilithiumParams>(a: i32, hint: u32) -> i32 {
    let q = P::Q as i32;
    let gamma2 = P::GAMMA2 as i32;
    let mut a0: i32 = 0;
    let a1: i32;

    a1 = decompose::<P>(&mut a0, a);

    if hint == 0 {
        return a1;
    }

    if gamma2 == (q - 1) / 32 {
        if a0 > 0 {
            (a1 + 1) & 15
        } else {
            (a1 - 1) & 15
        }
    } else if gamma2 == (q - 1) / 88 {
        if a0 > 0 {
            if a1 == 43 { 0 } else { a1 + 1 }
        } else {
            if a1 == 0 { 43 } else { a1 - 1 }
        }
    } else {
        panic!("Unsupported GAMMA2 value");
    }
}

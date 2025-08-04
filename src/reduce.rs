#![allow(dead_code)]

use crate::params::DilithiumParams;

const MONT: i32 = -4186625;
const QINV: i32 =  58728449;

pub fn montgomery_reduce<P: DilithiumParams>(a: i64) -> i32 {
    // Cast a to i32 truncates the lower 32 bits (like (int32_t)a in C)
    let t: i32 = (a as i32).wrapping_mul(QINV);
    // Now compute: (a - t * Q) >> 32
    let result = ((a - (t as i64) * (P::Q as i64)) >> 32) as i32;
    result
}

pub fn reduce32<P: DilithiumParams>(a: i32) -> i32 {
    let t = (a + (1 << 22)) >> 23;
    a - t * P::Q as i32
}

pub fn caddq<P: DilithiumParams>(a: i32) -> i32 {
    a + ((a >> 31) & P::Q as i32)
}


pub fn freeze<P: DilithiumParams>(a: i32) -> i32 {
    let a = reduce32::<P>(a);
    caddq::<P>(a)
}
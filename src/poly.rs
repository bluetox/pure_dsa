#![allow(dead_code)]

use crate::{
  reduce::*, 
  params::DilithiumParams, 
  ntt::*, 
  rounding::*,
  fips202::*
};

pub const SEEDBYTES: usize = 32;
pub const CRHBYTES: usize = 64;
pub const TRBYTES: usize = 64;
pub const RNDBYTES: usize = 32;
pub const N: usize = 256;
pub const Q: usize = 8380417;
pub const D: usize = 13;
pub const ROOT_OF_UNITY: usize = 1753;
pub const POLYT1_PACKEDBYTES: usize =  320;
pub const POLYT0_PACKEDBYTES: usize =  416;
pub const STREAM128_BLOCKBYTES: usize = SHAKE128_RATE;
pub const STREAM256_BLOCKBYTES: usize = SHAKE256_RATE;
pub const D_SHL: i32 = 1i32 << (D - 1);

#[derive(Copy, Clone, Debug)]
pub struct Poly {
  pub coeffs: [i32; N],
}

impl Default for Poly {
  fn default() -> Self {
    Poly { coeffs: [0i32; N] }
  }
}

pub fn poly_reduce<P: DilithiumParams>(a: &mut Poly) {
  for i in 0..N {
    a.coeffs[i] = reduce32::<P>(a.coeffs[i]);
  }
}

pub fn poly_caddq<P: DilithiumParams>(a: &mut Poly) {
  for i in 0..N {
    a.coeffs[i] = caddq::<P>(a.coeffs[i]);
  }
}

pub fn poly_add(c: &mut Poly, a: &Poly, b: &Poly) {
  for i in 0..N {
    c.coeffs[i] = a.coeffs[i] + b.coeffs[i];
  }
}

pub fn poly_sub(c: &mut Poly, a: &Poly, b: &Poly) {
  for i in 0..N {
    c.coeffs[i] = a.coeffs[i] - b.coeffs[i];
  }
}

pub fn poly_shiftl(a: &mut Poly) {
  for i in 0..N {
    a.coeffs[i] <<= D;
  }
}

pub fn poly_ntt<P: DilithiumParams>(a: &mut Poly) {
  ntt::<P>(&mut a.coeffs);
}

pub fn poly_invntt_tomont<P: DilithiumParams>(a: &mut Poly) {
  invntt_tomont::<P>(&mut a.coeffs);
}

pub fn poly_pointwise_montgomery<P: DilithiumParams>(c: &mut Poly, a: &Poly, b: &Poly) {
  for i in 0..N {
    c.coeffs[i] = montgomery_reduce::<P>((a.coeffs[i] as i64) * b.coeffs[i] as i64);
  }
}

pub fn poly_power2round<P: DilithiumParams>(a1: &mut Poly, a0: &mut Poly, _a: &Poly) {
  for i in 0..N {
    a1.coeffs[i] = power2round::<P>(&mut a0.coeffs[i], a1.coeffs[i]);
  }
}


pub fn poly_decompose<P: DilithiumParams>(a1: &mut Poly, a0: &mut Poly, &a: &Poly) {
  for i in 0..N {
    a1.coeffs[i] = decompose::<P>(&mut a0.coeffs[i], a.coeffs[i]);
  }
}

pub fn poly_make_hint<P: DilithiumParams>(h: &mut Poly, a0: &Poly, a1: &Poly) -> i32 {
  let mut s = 0i32;
  for i in 0..N {
    h.coeffs[i] = make_hint::<P>(a0.coeffs[i], a1.coeffs[i]) as i32;
    s += h.coeffs[i];
  }
  s
}

pub fn poly_use_hint<P: DilithiumParams>(b: &mut Poly, a: &Poly, h: &Poly) {
  for i in 0..N {
    b.coeffs[i] = use_hint::<P>(a.coeffs[i], h.coeffs[i] as u32);
  }
}

pub fn poly_chknorm<P: DilithiumParams>(a: &Poly, b: i32) -> u8 {
  let mut t;

  if b > ((P::Q - 1) / 8 ) as i32 {
    return 1;
  }
  for i in 0..N {
    t = a.coeffs[i] >> 31;
    t = a.coeffs[i] - (t & 2 * a.coeffs[i]);

    if t >= b {
      return 1;
    }
  }
  return 0;
}

#[inline(always)]
pub fn rej_uniform(a: &mut [i32], len: u32, buf: &[u8], buflen: usize) -> u32 {
  let (mut ctr, mut pos) = (0usize, 0usize);
  let mut t;
  while ctr < len as usize && pos + 3 <= buflen {
    t = buf[pos] as u32;
    pos += 1;
    t |= (buf[pos] as u32) << 8;
    pos += 1;
    t |= (buf[pos] as u32) << 16;
    pos += 1;
    t &= 0x7FFFFF;

    if t < Q as u32 {
      a[ctr] = t as i32;
      ctr += 1;
    }
  }
  ctr as u32
}

const POLY_UNIFORM_NBLOCKS: usize =
  (768 + STREAM128_BLOCKBYTES - 1) / STREAM128_BLOCKBYTES;

#[inline(always)]
pub fn poly_uniform<P: DilithiumParams>(a: &mut Poly, seed: &[u8], nonce: u16) {
  
  let mut buflen = POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES;
  let mut buf = [0u8; POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES + 2];
  let mut state = KeccakState::default();

  state.dilithium_shake128_stream_init(seed, nonce);

  state.shake128_squeezeblocks(&mut buf, POLY_UNIFORM_NBLOCKS);

  let mut ctr = rej_uniform(&mut a.coeffs, P::N as u32, &mut buf, buflen);
  let mut off;
  while ctr < P::N as u32 {
    off = buflen % 3;
    for i in 0..off {
      buf[i] = buf[buflen - off + i];
    }
    buflen = STREAM128_BLOCKBYTES + off;
    state.shake128_squeezeblocks(&mut buf[off..], 1);
    ctr += rej_uniform(
      &mut a.coeffs[(ctr as usize)..],
      P::N as u32 - ctr,
      &mut buf,
      buflen,
    );
  }
}

pub fn rej_eta<P: DilithiumParams>(a: &mut [i32], len: usize, buf: &[u8], buflen: usize) -> u32 {
  let (mut ctr, mut pos) = (0usize, 0usize);
  let (mut t0, mut t1);
  while ctr < len && pos < buflen {
    t0 = (buf[pos] & 0x0F) as u32;
    t1 = (buf[pos] >> 4) as u32;
    pos += 1;

    if P::ETA == 2 {
      if t0 < 15 {
        t0 = t0 - (205 * t0 >> 10) * 5;
        a[ctr] = 2 - t0 as i32;
        ctr += 1;
      }
      if t1 < 15 && ctr < len {
        t1 = t1 - (205 * t1 >> 10) * 5;
        a[ctr] = 2 - t1 as i32;
        ctr += 1;
      }
    } else if P::ETA == 4 {
      if t0 < 9 {
        a[ctr] = 4 - t0 as i32;
        ctr += 1;
      }
      if t1 < 9 && ctr < len {
        a[ctr] = 4 - t1 as i32;
        ctr += 1;
      }
    }
  }
  ctr as u32
}
pub const fn poly_uniform_eta_nblocks(eta: usize) -> usize {
    if eta == 2 {
        (136 + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES
    } else {
        (227 + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES
    }
}

#[allow(non_snake_case)]
pub fn poly_uniform_eta<P: DilithiumParams>(a: &mut Poly, seed: &[u8], nonce: u16) {
  
  let mut max = [0u8; 362];
  
  let mut buf = &mut max[..poly_uniform_eta_nblocks(P::ETA) * STREAM256_BLOCKBYTES];
  let buf_len = buf.len();
  let mut state = KeccakState::default();
  state.dilithium_shake256_stream_init(seed, nonce);
  state.shake256_squeezeblocks(
    buf,
    poly_uniform_eta_nblocks(P::ETA)
  );

  let mut ctr = rej_eta::<P>(&mut a.coeffs, N, &buf, buf_len);

  while ctr < P::N as u32 {
    state.shake256_squeezeblocks(&mut buf, 1);
    ctr += rej_eta::<P>(
      &mut a.coeffs[ctr as usize..],
      N - ctr as usize,
      &buf,
      STREAM256_BLOCKBYTES,
    );
  }
}


pub fn poly_uniform_gamma1<P: DilithiumParams>(a: &mut Poly, seed: &[u8], nonce: u16) {
  let mut poly_buf = P::poly_uniform_gamma1_buffer();
  
  let mut buf = poly_buf.buf();

  let mut state = KeccakState::default();

  state.dilithium_shake256_stream_init(seed, nonce);
  state.shake256_squeezeblocks(
    &mut buf,
    P::POLY_UNIFORM_GAMMA1_NBLOCKS
  );
  polyz_unpack::<P>(a, &mut buf);
}

pub fn poly_challenge<P: DilithiumParams>(c: &mut Poly, seed: &[u8]) {
  let mut _signs = 0u64;
  let mut buf = [0u8; SHAKE256_RATE];
  let mut state = KeccakState::default();

  state.shake256_absorb(seed, SEEDBYTES);
  state.shake256_finalize();
  state.shake256_squeezeblocks(&mut buf, 1);

  for i in 0..8 {
    _signs |= (buf[i] as u64) << 8 * i;
  }
  let mut pos: usize = 8;

  let mut b;
  c.coeffs.fill(0);
  for i in N - P::TAU..N {
    loop {
      if pos >= SHAKE256_RATE {
        state.shake256_squeezeblocks(&mut buf, 1);
        pos = 0;
      }
      b = buf[pos] as usize;
      pos += 1;
      if b <= i {
        break;
      }
    }
    c.coeffs[i] = c.coeffs[b as usize];
    c.coeffs[b as usize] = 1i32 - 2 * (_signs & 1) as i32;
    _signs >>= 1;
  }
}

pub fn polyeta_pack<P: DilithiumParams>(r: &mut [u8], a: &Poly) {
  let mut t = [0u8; 8];
  if P::ETA == 2 {
    for i in 0..N / 8 {
      t[0] = (P::ETA as i32 - a.coeffs[8 * i + 0]) as u8;
      t[1] = (P::ETA as i32 - a.coeffs[8 * i + 1]) as u8;
      t[2] = (P::ETA as i32 - a.coeffs[8 * i + 2]) as u8;
      t[3] = (P::ETA as i32 - a.coeffs[8 * i + 3]) as u8;
      t[4] = (P::ETA as i32 - a.coeffs[8 * i + 4]) as u8;
      t[5] = (P::ETA as i32 - a.coeffs[8 * i + 5]) as u8;
      t[6] = (P::ETA as i32 - a.coeffs[8 * i + 6]) as u8;
      t[7] = (P::ETA as i32 - a.coeffs[8 * i + 7]) as u8;

      r[3 * i + 0] = (t[0] >> 0) | (t[1] << 3) | (t[2] << 6);
      r[3 * i + 1] = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
      r[3 * i + 2] = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
    }
  } else {
    for i in 0..N / 2 {
      t[0] = (P::ETA as i32 - a.coeffs[2 * i + 0]) as u8;
      t[1] = (P::ETA as i32 - a.coeffs[2 * i + 1]) as u8;
      r[i] = t[0] | (t[1] << 4);
    }
  }
}

pub fn polyeta_unpack<P: DilithiumParams>(r: &mut Poly, a: &[u8]) {
  if P::ETA == 2 {
    for i in 0..N / 8 {
      r.coeffs[8 * i + 0] = (a[3 * i + 0] & 0x07) as i32;
      r.coeffs[8 * i + 1] = ((a[3 * i + 0] >> 3) & 0x07) as i32;
      r.coeffs[8 * i + 2] =
        (((a[3 * i + 0] >> 6) | (a[3 * i + 1] << 2)) & 0x07) as i32;
      r.coeffs[8 * i + 3] = ((a[3 * i + 1] >> 1) & 0x07) as i32;
      r.coeffs[8 * i + 4] = ((a[3 * i + 1] >> 4) & 0x07) as i32;
      r.coeffs[8 * i + 5] =
        (((a[3 * i + 1] >> 7) | (a[3 * i + 2] << 1)) & 0x07) as i32;
      r.coeffs[8 * i + 6] = ((a[3 * i + 2] >> 2) & 0x07) as i32;
      r.coeffs[8 * i + 7] = ((a[3 * i + 2] >> 5) & 0x07) as i32;

      r.coeffs[8 * i + 0] = (P::ETA as i32 - r.coeffs[8 * i + 0]) as i32;
      r.coeffs[8 * i + 1] = (P::ETA as i32 - r.coeffs[8 * i + 1]) as i32;
      r.coeffs[8 * i + 2] = (P::ETA as i32 - r.coeffs[8 * i + 2]) as i32;
      r.coeffs[8 * i + 3] = (P::ETA as i32 - r.coeffs[8 * i + 3]) as i32;
      r.coeffs[8 * i + 4] = (P::ETA as i32 - r.coeffs[8 * i + 4]) as i32;
      r.coeffs[8 * i + 5] = (P::ETA as i32 - r.coeffs[8 * i + 5]) as i32;
      r.coeffs[8 * i + 6] = (P::ETA as i32 - r.coeffs[8 * i + 6]) as i32;
      r.coeffs[8 * i + 7] = (P::ETA as i32 - r.coeffs[8 * i + 7]) as i32;
    }
  } else {
    for i in 0..N / 2 {
      r.coeffs[2 * i + 0] = (a[i] & 0x0F) as i32;
      
      r.coeffs[2 * i + 1] = (a[i] >> 4) as i32;
      
      r.coeffs[2 * i + 0] = (P::ETA as i32 - r.coeffs[2 * i + 0]) as i32;
      r.coeffs[2 * i + 1] = (P::ETA as i32 - r.coeffs[2 * i + 1]) as i32;
    }
  }
}

pub fn polyt1_pack(r: &mut [u8], a: &Poly) {
  for i in 0..N / 4 {
    r[5 * i + 0] = (a.coeffs[4 * i + 0] >> 0) as u8;
    r[5 * i + 1] =
      ((a.coeffs[4 * i + 0] >> 8) | (a.coeffs[4 * i + 1] << 2)) as u8;
    r[5 * i + 2] =
      ((a.coeffs[4 * i + 1] >> 6) | (a.coeffs[4 * i + 2] << 4)) as u8;
    r[5 * i + 3] =
      ((a.coeffs[4 * i + 2] >> 4) | (a.coeffs[4 * i + 3] << 6)) as u8;
    r[5 * i + 4] = (a.coeffs[4 * i + 3] >> 2) as u8;
  }
}

pub fn polyt1_unpack(r: &mut Poly, a: &[u8]) {
  for i in 0..N / 4 {
    r.coeffs[4 * i + 0] = (((a[5 * i + 0] >> 0) as u32
      | (a[5 * i + 1] as u32) << 8)
      & 0x3FF) as i32;
    r.coeffs[4 * i + 1] = (((a[5 * i + 1] >> 2) as u32
      | (a[5 * i + 2] as u32) << 6)
      & 0x3FF) as i32;
    r.coeffs[4 * i + 2] = (((a[5 * i + 2] >> 4) as u32
      | (a[5 * i + 3] as u32) << 4)
      & 0x3FF) as i32;
    r.coeffs[4 * i + 3] = (((a[5 * i + 3] >> 6) as u32
      | (a[5 * i + 4] as u32) << 2)
      & 0x3FF) as i32;
  }
}

pub fn polyt0_pack(r: &mut [u8], a: &Poly) {
  let mut t = [0i32; 8];

  for i in 0..N / 8 {
    t[0] = D_SHL - a.coeffs[8 * i + 0];
    t[1] = D_SHL - a.coeffs[8 * i + 1];
    t[2] = D_SHL - a.coeffs[8 * i + 2];
    t[3] = D_SHL - a.coeffs[8 * i + 3];
    t[4] = D_SHL - a.coeffs[8 * i + 4];
    t[5] = D_SHL - a.coeffs[8 * i + 5];
    t[6] = D_SHL - a.coeffs[8 * i + 6];
    t[7] = D_SHL - a.coeffs[8 * i + 7];

    r[13 * i + 0] = (t[0]) as u8;
    r[13 * i + 1] = (t[0] >> 8) as u8;
    r[13 * i + 1] |= (t[1] << 5) as u8;
    r[13 * i + 2] = (t[1] >> 3) as u8;
    r[13 * i + 3] = (t[1] >> 11) as u8;
    r[13 * i + 3] |= (t[2] << 2) as u8;
    r[13 * i + 4] = (t[2] >> 6) as u8;
    r[13 * i + 4] |= (t[3] << 7) as u8;
    r[13 * i + 5] = (t[3] >> 1) as u8;
    r[13 * i + 6] = (t[3] >> 9) as u8;
    r[13 * i + 6] |= (t[4] << 4) as u8;
    r[13 * i + 7] = (t[4] >> 4) as u8;
    r[13 * i + 8] = (t[4] >> 12) as u8;
    r[13 * i + 8] |= (t[5] << 1) as u8;
    r[13 * i + 9] = (t[5] >> 7) as u8;
    r[13 * i + 9] |= (t[6] << 6) as u8;
    r[13 * i + 10] = (t[6] >> 2) as u8;
    r[13 * i + 11] = (t[6] >> 10) as u8;
    r[13 * i + 11] |= (t[7] << 3) as u8;
    r[13 * i + 12] = (t[7] >> 5) as u8;
  }
}

pub fn polyt0_unpack(r: &mut Poly, a: &[u8]) {
  for i in 0..N / 8 {
    r.coeffs[8 * i + 0] = a[13 * i + 0] as i32;
    r.coeffs[8 * i + 0] |= (a[13 * i + 1] as i32) << 8;
    r.coeffs[8 * i + 0] &= 0x1FFF;

    r.coeffs[8 * i + 1] = (a[13 * i + 1] as i32) >> 5;
    r.coeffs[8 * i + 1] |= (a[13 * i + 2] as i32) << 3;
    r.coeffs[8 * i + 1] |= (a[13 * i + 3] as i32) << 11;
    r.coeffs[8 * i + 1] &= 0x1FFF;

    r.coeffs[8 * i + 2] = (a[13 * i + 3] as i32) >> 2;
    r.coeffs[8 * i + 2] |= (a[13 * i + 4] as i32) << 6;
    r.coeffs[8 * i + 2] &= 0x1FFF;

    r.coeffs[8 * i + 3] = (a[13 * i + 4] as i32) >> 7;
    r.coeffs[8 * i + 3] |= (a[13 * i + 5] as i32) << 1;
    r.coeffs[8 * i + 3] |= (a[13 * i + 6] as i32) << 9;
    r.coeffs[8 * i + 3] &= 0x1FFF;

    r.coeffs[8 * i + 4] = (a[13 * i + 6] as i32) >> 4;
    r.coeffs[8 * i + 4] |= (a[13 * i + 7] as i32) << 4;
    r.coeffs[8 * i + 4] |= (a[13 * i + 8] as i32) << 12;
    r.coeffs[8 * i + 4] &= 0x1FFF;

    r.coeffs[8 * i + 5] = (a[13 * i + 8] as i32) >> 1;
    r.coeffs[8 * i + 5] |= (a[13 * i + 9] as i32) << 7;
    r.coeffs[8 * i + 5] &= 0x1FFF;

    r.coeffs[8 * i + 6] = (a[13 * i + 9] as i32) >> 6;
    r.coeffs[8 * i + 6] |= (a[13 * i + 10] as i32) << 2;
    r.coeffs[8 * i + 6] |= (a[13 * i + 11] as i32) << 10;
    r.coeffs[8 * i + 6] &= 0x1FFF;

    r.coeffs[8 * i + 7] = (a[13 * i + 11] as i32) >> 3;
    r.coeffs[8 * i + 7] |= (a[13 * i + 12] as i32) << 5;
    r.coeffs[8 * i + 7] &= 0x1FFF; // TODO: Unnecessary mask?

    r.coeffs[8 * i + 0] = D_SHL - r.coeffs[8 * i + 0];
    r.coeffs[8 * i + 1] = D_SHL - r.coeffs[8 * i + 1];
    r.coeffs[8 * i + 2] = D_SHL - r.coeffs[8 * i + 2];
    r.coeffs[8 * i + 3] = D_SHL - r.coeffs[8 * i + 3];
    r.coeffs[8 * i + 4] = D_SHL - r.coeffs[8 * i + 4];
    r.coeffs[8 * i + 5] = D_SHL - r.coeffs[8 * i + 5];
    r.coeffs[8 * i + 6] = D_SHL - r.coeffs[8 * i + 6];
    r.coeffs[8 * i + 7] = D_SHL - r.coeffs[8 * i + 7];
  }
}

pub fn polyz_pack<P: DilithiumParams>(r: &mut [u8], a: &Poly) {
  let mut t = [0i32; 4];
  if P::GAMMA1 == (1 << 17) {
    for i in 0..N / 4 {
      t[0] = P::GAMMA1 as i32 - a.coeffs[4 * i + 0];
      t[1] = P::GAMMA1 as i32 - a.coeffs[4 * i + 1];
      t[2] = P::GAMMA1 as i32 - a.coeffs[4 * i + 2];
      t[3] = P::GAMMA1 as i32 - a.coeffs[4 * i + 3];

      r[9 * i + 0] = (t[0]) as u8;
      r[9 * i + 1] = (t[0] >> 8) as u8;
      r[9 * i + 2] = (t[0] >> 16) as u8;
      r[9 * i + 2] |= (t[1] << 2) as u8;
      r[9 * i + 3] = (t[1] >> 6) as u8;
      r[9 * i + 4] = (t[1] >> 14) as u8;
      r[9 * i + 4] |= (t[2] << 4) as u8;
      r[9 * i + 5] = (t[2] >> 4) as u8;
      r[9 * i + 6] = (t[2] >> 12) as u8;
      r[9 * i + 6] |= (t[3] << 6) as u8;
      r[9 * i + 7] = (t[3] >> 2) as u8;
      r[9 * i + 8] = (t[3] >> 10) as u8;
    }
  } else if P::GAMMA1 == 1 << 19 {
    for i in 0..N / 2 {
      t[0] = P::GAMMA1 as i32 - a.coeffs[2 * i + 0];
      t[1] = P::GAMMA1 as i32 - a.coeffs[2 * i + 1];

      r[5 * i + 0] = (t[0]) as u8;
      r[5 * i + 1] = (t[0] >> 8) as u8;
      r[5 * i + 2] = (t[0] >> 16) as u8;
      r[5 * i + 2] |= (t[1] << 4) as u8;
      r[5 * i + 3] = (t[1] >> 4) as u8;
      r[5 * i + 4] = (t[1] >> 12) as u8;
    }
  }
}

pub fn polyz_unpack<P: DilithiumParams>(r: &mut Poly, a: &[u8]) {
  if P::GAMMA1 == (1 << 17) {
    for i in 0..N / 4 {
      r.coeffs[4 * i + 0] = a[9 * i + 0] as i32;
      r.coeffs[4 * i + 0] |= (a[9 * i + 1] as i32) << 8;
      r.coeffs[4 * i + 0] |= (a[9 * i + 2] as i32) << 16;
      r.coeffs[4 * i + 0] &= 0x3FFFF;

      r.coeffs[4 * i + 1] = (a[9 * i + 2] as i32) >> 2;
      r.coeffs[4 * i + 1] |= (a[9 * i + 3] as i32) << 6;
      r.coeffs[4 * i + 1] |= (a[9 * i + 4] as i32) << 14;
      r.coeffs[4 * i + 1] &= 0x3FFFF;

      r.coeffs[4 * i + 2] = (a[9 * i + 4] as i32) >> 4;
      r.coeffs[4 * i + 2] |= (a[9 * i + 5] as i32) << 4;
      r.coeffs[4 * i + 2] |= (a[9 * i + 6] as i32) << 12;
      r.coeffs[4 * i + 2] &= 0x3FFFF;

      r.coeffs[4 * i + 3] = (a[9 * i + 6] as i32) >> 6;
      r.coeffs[4 * i + 3] |= (a[9 * i + 7] as i32) << 2;
      r.coeffs[4 * i + 3] |= (a[9 * i + 8] as i32) << 10;
      r.coeffs[4 * i + 3] &= 0x3FFFF; // TODO: Unnecessary mask?

      r.coeffs[4 * i + 0] = P::GAMMA1 as i32 - r.coeffs[4 * i + 0];
      r.coeffs[4 * i + 1] = P::GAMMA1 as i32 - r.coeffs[4 * i + 1];
      r.coeffs[4 * i + 2] = P::GAMMA1 as i32 - r.coeffs[4 * i + 2];
      r.coeffs[4 * i + 3] = P::GAMMA1 as i32 - r.coeffs[4 * i + 3];
    }
  } else if P::GAMMA1 == 1 << 19 {
    for i in 0..N / 2 {
      r.coeffs[2 * i + 0] = a[5 * i + 0] as i32;
      r.coeffs[2 * i + 0] |= (a[5 * i + 1] as i32) << 8;
      r.coeffs[2 * i + 0] |= (a[5 * i + 2] as i32) << 16;
      r.coeffs[2 * i + 0] &= 0xFFFFF;

      r.coeffs[2 * i + 1] = (a[5 * i + 2] as i32) >> 4;
      r.coeffs[2 * i + 1] |= (a[5 * i + 3] as i32) << 4;
      r.coeffs[2 * i + 1] |= (a[5 * i + 4] as i32) << 12;
      r.coeffs[2 * i + 0] &= 0xFFFFF; // TODO: Unnecessary mask?

      r.coeffs[2 * i + 0] = P::GAMMA1 as i32 - r.coeffs[2 * i + 0];
      r.coeffs[2 * i + 1] = P::GAMMA1 as i32 - r.coeffs[2 * i + 1];
    }
  }
}

pub fn polyw1_pack<P: DilithiumParams>(r: &mut [u8], a: &Poly) {
  if P::GAMMA2 == (Q - 1) / 88 {
    for i in 0..N / 4 {
      r[3 * i + 0] = a.coeffs[4 * i + 0] as u8;
      r[3 * i + 0] |= (a.coeffs[4 * i + 1] << 6) as u8;
      r[3 * i + 1] = (a.coeffs[4 * i + 1] >> 2) as u8;
      r[3 * i + 1] |= (a.coeffs[4 * i + 2] << 4) as u8;
      r[3 * i + 2] = (a.coeffs[4 * i + 2] >> 4) as u8;
      r[3 * i + 2] |= (a.coeffs[4 * i + 3] << 2) as u8;
    }
  } else {
    for i in 0..N / 2 {
      r[i] = (a.coeffs[2 * i + 0] | (a.coeffs[2 * i + 1] << 4)) as u8;
    }
  }
}

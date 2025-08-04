use crate::{
  poly::*, 
  params::DilithiumParams
};

#[derive(Copy, Clone, Debug)]
pub enum Polyveck {
    Mode2(PolyveckStruct<4>),
    Mode3(PolyveckStruct<6>),
    Mode5(PolyveckStruct<8>)
}

#[derive(Copy, Clone, Debug)]
pub struct PolyveckStruct<const K: usize> {
    pub vec: [Poly; K],
}

impl<const K: usize> Default for PolyveckStruct<K> {
  fn default() -> Self {
    PolyveckStruct {
      vec: [Poly::default(); K],
    }
  }
}
#[derive(Copy, Clone, Debug)]
pub enum Polyvecl {
    Mode2(PolyveclStruct<4>),
    Mode3(PolyveclStruct<5>),
    Mode5(PolyveclStruct<7>)
}


#[derive(Copy, Clone, Debug)]
pub struct PolyveclStruct<const K: usize> {
    pub vec: [Poly; K],
}

impl<const L: usize> Default for PolyveclStruct<L> {
  fn default() -> Self {
    PolyveclStruct {
      vec: [Poly::default(); L],
    }
  }
}

#[inline(always)]
pub fn vec_from_polyveck_mut(v: &mut Polyveck) -> &mut [Poly] {
    match v {
        Polyveck::Mode2(v) => &mut v.vec,
        Polyveck::Mode3(v) => &mut v.vec,
        Polyveck::Mode5(v) => &mut v.vec,
    }
}
#[inline(always)]
pub fn vec_from_polyveck(v: &Polyveck) -> &[Poly] {
    match v {
        Polyveck::Mode2(inner) => &inner.vec,
        Polyveck::Mode3(inner) => &inner.vec,
        Polyveck::Mode5(inner) => &inner.vec,
    }
}
#[inline(always)]
pub fn vec_from_polyvecl_mut(v: &mut Polyvecl) -> &mut [Poly] {
    match v {
        Polyvecl::Mode2(v) => &mut v.vec,
        Polyvecl::Mode3(v) => &mut v.vec,
        Polyvecl::Mode5(v) => &mut v.vec,
    }
}
#[inline(always)]
pub fn vec_from_polyvecl(v: &Polyvecl) -> &[Poly] {
    match v {
        Polyvecl::Mode2(inner) => &inner.vec,
        Polyvecl::Mode3(inner) => &inner.vec,
        Polyvecl::Mode5(inner) => &inner.vec,
    }
}

pub fn polyvec_matrix_expand<P: DilithiumParams>(mat: &mut [Polyvecl], rho: &[u8]) {
    for (i, row) in mat.iter_mut().enumerate() {
      let row_vec = vec_from_polyvecl_mut(row);
      for (j, poly) in row_vec.iter_mut().enumerate() {
          poly_uniform::<P>(poly, rho, ((i << 8) + j) as u16);
      }
  }
}

pub fn polyvec_matrix_pointwise_montgomery<P: DilithiumParams>(
    t: &mut Polyveck,
    mat: &[Polyvecl],
    v: &Polyvecl,
) {
    let t_vec = vec_from_polyveck_mut(t);
    for i in 0..t_vec.len() {
        polyvecl_pointwise_acc_montgomery::<P>(&mut t_vec[i], &mat[i], v);
    }
}

//*********** Vectors of polynomials of length L ****************************

pub fn polyvecl_uniform_eta<P: DilithiumParams>(
  v: &mut Polyvecl, 
  seed: &[u8], 
  mut nonce: u16
) {
    let v_vec = vec_from_polyvecl_mut(v);
    for i in 0..P::L {
      poly_uniform_eta::<P>(&mut v_vec[i], seed, nonce);
    nonce += 1;
  }
}

pub fn polyvecl_uniform_gamma1<P: DilithiumParams>(v: &mut Polyvecl, seed: &[u8], nonce: u16) {
  let v_vec = vec_from_polyvecl_mut(v);
  for i in 0..v_vec.len() {
    poly_uniform_gamma1::<P>(&mut v_vec[i], seed, P::L as u16 * nonce + i as u16);
  }
}
pub fn polyvecl_reduce<P: DilithiumParams>(v: &mut Polyvecl) {
  let v_vec = vec_from_polyvecl_mut(v);
  for i in 0..v_vec.len() {
    poly_reduce::<P>(&mut v_vec[i]);
  }
}

/// Add vectors of polynomials of length L.
/// No modular reduction is performed.
pub fn polyvecl_add(w: &mut Polyvecl, u: &Polyvecl, v: &Polyvecl) {
  let w_vec = vec_from_polyvecl_mut(w);
  let u_vec = vec_from_polyvecl(u);
  let v_vec = vec_from_polyvecl(v);
  for i in 0..u_vec.len() {
    poly_add(&mut w_vec[i], &u_vec[i],&v_vec[i]);
  }
}

/// Forward NTT of all polynomials in vector of length L. Output
/// coefficients can be up to 16*Q larger than input coefficients.*
pub fn polyvecl_ntt<P: DilithiumParams>(v: &mut Polyvecl) {
    let v_vec = vec_from_polyvecl_mut(v);
    for i in 0..v_vec.len() {
    poly_ntt::<P>(&mut v_vec[i]);
  }
}

pub fn polyvecl_invntt_tomont<P: DilithiumParams>(v: &mut Polyvecl) {
  let v_vec = vec_from_polyvecl_mut(v);
  for i in 0..v_vec.len() {
    poly_invntt_tomont::<P>(&mut v_vec[i]);
  }
}

pub fn polyvecl_pointwise_poly_montgomery<P: DilithiumParams>(
  r: &mut Polyvecl,
  a: &Poly,
  v: &Polyvecl,
) {
  let r_vec = vec_from_polyvecl_mut(r);
  let v_vec = vec_from_polyvecl(v);
  for i in 0..r_vec.len() {
    poly_pointwise_montgomery::<P>(&mut r_vec[i], a, &v_vec[i]);
  }
}

/// Pointwise multiply vectors of polynomials of length L, multiply
/// resulting vector by 2^{-32} and add (accumulate) polynomials
/// in it. Input/output vectors are in NTT domain representation.
/// Input coefficients are assumed to be less than 22*Q. Output
/// coeffcient are less than 2*L*Q.
pub fn polyvecl_pointwise_acc_montgomery<P: DilithiumParams>(
  w: &mut Poly,
  u: &Polyvecl,
  v: &Polyvecl,
) {
  let u_vec = vec_from_polyvecl(u);
  let v_vec = vec_from_polyvecl(v);
  let mut t = Poly::default();
  poly_pointwise_montgomery::<P>(w, &u_vec[0], &v_vec[0]);
  for i in 1..u_vec.len() {
    poly_pointwise_montgomery::<P>(&mut t, &u_vec[i], &v_vec[i]);
    let tmp = w.clone();
    poly_add(w, &tmp,&t);
  }
}

/// Check infinity norm of polynomials in vector of length L.
/// Assumes input coefficients to be standard representatives.
/// Returns 0 if norm of all polynomials is strictly smaller than B and 1
/// otherwise.
pub fn polyvecl_chknorm<P: DilithiumParams>(v: &Polyvecl, bound: i32) -> u8 {
  let v_vec = vec_from_polyvecl(v);
  for i in 0..v_vec.len() {
    if poly_chknorm::<P>(&v_vec[i], bound) > 0 {
      return 1;
    }
  }
  return 0;
}

//*********** Vectors of polynomials of length K ****************************

pub fn polyveck_uniform_eta<P: DilithiumParams>(v: &mut Polyveck, seed: &[u8], mut nonce: u16) {
  let v_vec = vec_from_polyveck_mut(v);
  for i in 0..v_vec.len() {
    poly_uniform_eta::<P>(&mut v_vec[i], seed, nonce);
    nonce += 1
  }
}

/// Reduce coefficients of polynomials in vector of length K
/// to representatives in [0,2*Q].
pub fn polyveck_reduce<P: DilithiumParams>(v: &mut Polyveck) {
  let v_vec = vec_from_polyveck_mut(v);
  for i in 0..v_vec.len() {
    poly_reduce::<P>(&mut v_vec[i]);
  }
}

/// For all coefficients of polynomials in vector of length K
/// add Q if coefficient is negative.
pub fn polyveck_caddq<P: DilithiumParams>(v: &mut Polyveck) {
  let v_vec = vec_from_polyveck_mut(v);
  
  for i in 0..P::K {
    poly_caddq::<P>(&mut v_vec[i]);
  }
}

/// Add vectors of polynomials of length K.
/// No modular reduction is performed.
pub fn polyveck_add(w: &mut Polyveck, u: &Polyveck,v: &Polyveck) {
  let w_vec = vec_from_polyveck_mut(w);
  let v_vec = vec_from_polyveck(v);
  let u_vec = vec_from_polyveck(u);
  for i in 0..v_vec.len() {
    poly_add(&mut w_vec[i], &u_vec[i], &v_vec[i]);
  }
}

/// Subtract vectors of polynomials of length K.
/// Assumes coefficients of polynomials in second input vector
/// to be less than 2*Q. No modular reduction is performed.
pub fn polyveck_sub(w: &mut Polyveck, u: &Polyveck, v: &Polyveck) {
  let w_vec = vec_from_polyveck_mut(w);
  let v_vec = vec_from_polyveck(v);
  let u_vec = vec_from_polyveck(u);

  for i in 0..w_vec.len() {
    poly_sub(&mut w_vec[i], &u_vec[i], &v_vec[i]);
  }
}

/// Multiply vector of polynomials of Length K by 2^D without modular
/// reduction. Assumes input coefficients to be less than 2^{32-D}.
pub fn polyveck_shiftl(v: &mut Polyveck) {
  let v_vec = vec_from_polyveck_mut(v);
  for i in 0..v_vec.len() {
    poly_shiftl(&mut v_vec[i]);
  }
}

/// Forward NTT of all polynomials in vector of length K. Output
/// coefficients can be up to 16*Q larger than input coefficients.
pub fn polyveck_ntt<P: DilithiumParams>(v: &mut Polyveck) {
  let v_vec = vec_from_polyveck_mut(v);
  for i in 0..v_vec.len() {
    poly_ntt::<P>(&mut v_vec[i]);
  }
}

/// Inverse NTT and multiplication by 2^{32} of polynomials
/// in vector of length K. Input coefficients need to be less
/// than 2*Q.
pub fn polyveck_invntt_tomont<P: DilithiumParams>(v: &mut Polyveck) {
  let v_vec = vec_from_polyveck_mut(v);
  for i in 0..v_vec.len() {
    poly_invntt_tomont::<P>(&mut v_vec[i]);
  }
}

pub fn polyveck_pointwise_poly_montgomery<P: DilithiumParams>(
  r: &mut Polyveck,
  a: &Poly,
  v: &Polyveck,
) {
  let r_vec = vec_from_polyveck_mut(r);
  let v_vec = vec_from_polyveck(v);
  for i in 0..r_vec.len() {
    poly_pointwise_montgomery::<P>(&mut r_vec[i], a, &v_vec[i]);
  }
}

/// Check infinity norm of polynomials in vector of length K.
/// Assumes input coefficients to be standard representatives.
//
/// Returns 0 if norm of all polynomials are strictly smaller than B and 1
/// otherwise.
pub fn polyveck_chknorm<P: DilithiumParams>(v: &Polyveck, bound: i32) -> u8 {
  let v_vec = vec_from_polyveck(v);
  for i in 0..v_vec.len() {
    if poly_chknorm::<P>(&v_vec[i], bound) > 0 {
      return 1;
    }
  }
  return 0;
}

/// For all coefficients a of polynomials in vector of length K,
/// compute a0, a1 such that a mod Q = a1*2^D + a0
/// with -2^{D-1} < a0 <= 2^{D-1}. Assumes coefficients to be
/// standard representatives.
pub fn polyveck_power2round<P: DilithiumParams>(v1: &mut Polyveck, v0: &mut Polyveck, v: &Polyveck) {
  let v1_vec = vec_from_polyveck_mut(v1);
  let v0_vec = vec_from_polyveck_mut(v0);
  let v_vec = vec_from_polyveck(v);
  for i in 0.. v1_vec.len() {
    poly_power2round::<P>(&mut v1_vec[i], &mut v0_vec[i], &v_vec[i]);
  }
}

/// For all coefficients a of polynomials in vector of length K,
/// compute high and low bits a0, a1 such a mod Q = a1*ALPHA + a0
/// with -ALPHA/2 < a0 <= ALPHA/2 except a1 = (Q-1)/ALPHA where we
/// set a1 = 0 and -ALPHA/2 <= a0 = a mod Q - Q < 0.
/// Assumes coefficients to be standard representatives.
pub fn polyveck_decompose<P: DilithiumParams>(v1: &mut Polyveck, v0: &mut Polyveck, v: &Polyveck) {
  let v1_vec = vec_from_polyveck_mut(v1);
  let v0_vec = vec_from_polyveck_mut(v0);
  let v_vec = vec_from_polyveck(v);
  for i in 0..v1_vec.len() {
    poly_decompose::<P>(&mut v1_vec[i], &mut v0_vec[i], &v_vec[i]);
  }
}

/// Compute hint vector.
///
/// Returns number of 1 bits.
pub fn polyveck_make_hint<P: DilithiumParams>(
  h: &mut Polyveck,
  v0: &Polyveck,
  v1: &Polyveck,
) -> i32 {
  let v0_vec = vec_from_polyveck(v0);
  let v1_vec = vec_from_polyveck(v1);
  let h_vec = vec_from_polyveck_mut(h);
  let mut s = 0i32;
  for i in 0.. h_vec.len() {
    s += poly_make_hint::<P>(&mut h_vec[i], &v0_vec[i], &v1_vec[i]);
  }
  s
}

/// Use hint vector to correct the high bits of input vector.
pub fn polyveck_use_hint<P: DilithiumParams>(w: &mut Polyveck, u: &Polyveck, h: &Polyveck) {
  let w_vec = vec_from_polyveck_mut(w);
  let u_vec = vec_from_polyveck(u);
  let h_vec = vec_from_polyveck(h);
  for i in 0..w_vec.len() {
    poly_use_hint::<P>(&mut w_vec[i], &u_vec[i], &h_vec[i]);
  }
}

pub fn polyveck_pack_w1<P: DilithiumParams>(r: &mut [u8], w1: &Polyveck) {
  let w1_vec = vec_from_polyveck(w1);
  for i in 0..P::K {
    polyw1_pack::<P>(&mut r[i * P::polyw1_packedbytes()..], &w1_vec[i]);
  }
}

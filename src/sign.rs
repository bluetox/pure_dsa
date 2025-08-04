#![allow(dead_code)]

use rand::RngCore;

use crate::{
    randombytes::randombytes, 
    params::{
      DilithiumParams,
      Mat
    },
    polyvec::*,
    fips202::*,
    packing::*,
    poly::*
};

fn mat_from(v: &mut Mat) -> &mut [Polyvecl] {
    match v {
        Mat::Mode2(v) => &mut v[..],
        Mat::Mode3(v) => &mut v[..],
        Mat::Mode5(v) => &mut v[..],
    }
}
const SEEDBYTES: usize = 32;
const CRHBYTES: usize = 64;
const RNDBYTES: usize = 32;
const N: usize = 256;
const Q: usize = 8380417;
const D: usize = 13;
const ROOT_OF_UNITY: usize = 1753;

pub fn crypto_sign_keypair<P: DilithiumParams, R: RngCore>(
  pk: &mut [u8],
  sk: &mut [u8],
  rng: &mut R,
) -> u8 {
    let mut seedbuf = [0u8; 2 * SEEDBYTES + CRHBYTES];
    let mut tr = [0u8; SEEDBYTES];
    let mut mat_data = P::mat();
   
    let mut seed: [u8; SEEDBYTES] = [0; SEEDBYTES];

    randombytes(&mut seed, rng);

    shake256(&mut seedbuf, &seed, 2 * SEEDBYTES + CRHBYTES, SEEDBYTES);
  
    let rho_slice = &seedbuf[..SEEDBYTES];
    let rhoprime_slice = &seedbuf[SEEDBYTES..SEEDBYTES + CRHBYTES];
    let key_slice = &seedbuf[SEEDBYTES + CRHBYTES..];

    let mut mat = mat_from(&mut mat_data);

    polyvec_matrix_expand::<P>(&mut mat, rho_slice);

  let mut s1 = P::polyveclnew();
  let mut s2 = P::polyvecknew();
  
  polyvecl_uniform_eta::<P>(&mut s1, rhoprime_slice, 0);
  
  polyveck_uniform_eta::<P>(&mut s2, rhoprime_slice, P::L as u16);
  let mut s1hat = s1;
  polyvecl_ntt::<P>(&mut s1hat);
  
  

  let mut t1 = P::polyvecknew();
  polyvec_matrix_pointwise_montgomery::<P>(&mut t1, &mat, &s1hat);
  
  polyveck_reduce::<P>(&mut t1);
  polyveck_invntt_tomont::<P>(&mut t1);

  let t1_clone = t1.clone();
  polyveck_add(&mut t1, &t1_clone, &s2);

  polyveck_caddq::<P>(&mut t1);
  
  let mut t0 = P::polyvecknew();
  let t1_clone = t1.clone();
  polyveck_power2round::<P>(&mut t1, &mut t0, &t1_clone);

  pack_pk::<P>(pk, rho_slice, &t1);
  

  shake256(&mut tr, &pk[..P::crypto_publickeybytes()], SEEDBYTES, P::PUBLIC_KEY_BYTES);

  pack_sk::<P>(sk, rho_slice, &tr, key_slice, &t0, &s1, &s2);
  return 0;
}

pub fn crypto_sign_signature<P: DilithiumParams, R: RngCore>(sig: &mut [u8], m: &[u8], sk: &[u8], rng: &mut R) {

  let mut keymu = [0u8; SEEDBYTES + CRHBYTES];

  let mut nonce = 0u16;
  let mut mat_data = P::mat();
   let mut mat = mat_from(&mut mat_data);
  let (mut s1, mut y) = (P::polyveclnew(), P::polyveclnew());
  let (mut s2, mut t0) = (P::polyvecknew(), P::polyvecknew());
  let (mut w1, mut w0) = (P::polyvecknew(), P::polyvecknew());
  let mut h = P::polyvecknew();
  let mut cp = Poly::default();
  let mut state = KeccakState::default();
  let mut rho = [0u8; SEEDBYTES];
  let mut tr = [0u8; SEEDBYTES];

  let mut rhoprime = [0u8; CRHBYTES];

  unpack_sk::<P>(
    &mut rho,
    &mut tr,
    &mut keymu[..SEEDBYTES],
    &mut t0,
    &mut s1,
    &mut s2,
    &sk,
  );
  
  state.shake256_absorb(&tr, SEEDBYTES);
  state.shake256_absorb(m, m.len());
  state.shake256_finalize();
  state.shake256_squeeze(&mut keymu[SEEDBYTES..], CRHBYTES);
  
  #[cfg(feature = "random")]
  randombytes(&mut rhoprime, rng);

  #[cfg(not(feature = "random"))]
  shake256(&mut rhoprime, &keymu, CRHBYTES, SEEDBYTES + CRHBYTES);

  // Expand matrix and transform vectors
  polyvec_matrix_expand::<P>(&mut mat, &rho);
  polyvecl_ntt::<P>(&mut s1);
  polyveck_ntt::<P>(&mut s2);
  polyveck_ntt::<P>(&mut t0);

  loop {
    // Sample intermediate vector y
    polyvecl_uniform_gamma1::<P>(&mut y, &rhoprime, nonce);
    nonce += 1;

    // Matrix-vector multiplication
    let mut z = y;
    polyvecl_ntt::<P>(&mut z);
    polyvec_matrix_pointwise_montgomery::<P>(&mut w1, &mat, &z);
    polyveck_reduce::<P>(&mut w1);
    polyveck_invntt_tomont::<P>(&mut w1);

    // Decompose w and call the random oracle
    polyveck_caddq::<P>(&mut w1);
    let w1_clone = w1.clone();
    polyveck_decompose::<P>(&mut w1, &mut w0, &w1_clone);
    polyveck_pack_w1::<P>(sig, &w1);

    state.init();
    state.shake256_absorb(&keymu[SEEDBYTES..], CRHBYTES);
    state.shake256_absorb(&sig, P::K * P::polyw1_packedbytes());
    state.shake256_finalize();
    state.shake256_squeeze(sig, SEEDBYTES);

    poly_challenge::<P>(&mut cp, sig);

    poly_ntt::<P>(&mut cp);

    // Compute z, reject if it reveals secret
    polyvecl_pointwise_poly_montgomery::<P>(&mut z, &cp, &s1);
    polyvecl_invntt_tomont::<P>(&mut z);
    let z_clone = z.clone();
    polyvecl_add(&mut z, &z_clone, &y);
    polyvecl_reduce::<P>(&mut z);
    if polyvecl_chknorm::<P>(&z, (P::GAMMA1 - P::BETA) as i32) > 0 {
      continue;
    }

    /* Check that subtracting cs2 does not change high bits of w and low bits
     * do not reveal secret information */
    polyveck_pointwise_poly_montgomery::<P>(&mut h, &cp, &s2);
    polyveck_invntt_tomont::<P>(&mut h);
    let w0_clone = w0.clone();
    polyveck_sub(&mut w0, &w0_clone, &h);
    polyveck_reduce::<P>(&mut w0);
    if polyveck_chknorm::<P>(&w0, (P::GAMMA2 - P::BETA) as i32) > 0 {
      continue;
    }

    // Compute hints for w1
    polyveck_pointwise_poly_montgomery::<P>(&mut h, &cp, &t0);
    polyveck_invntt_tomont::<P>(&mut h);
    polyveck_reduce::<P>(&mut h);
    if polyveck_chknorm::<P>(&h, P::GAMMA2 as i32) > 0 {
      continue;
    }
    let w0_clone = w0.clone();
    polyveck_add(&mut w0, &w0_clone, &h);
    let n = polyveck_make_hint::<P>(&mut h, &w0, &w1);
    if n > P::OMEGA as i32 {
      continue;
    }

    // Write signature
    pack_sig::<P>(sig, None, &z, &h);
    return;
  }
}

pub fn crypto_sign_verify<P: DilithiumParams>(
  sig: &[u8],
  m: &[u8],
  pk: &[u8],
) -> Result<(), String> {
  let mut buf = vec![0u8; P::K * P::polyw1_packedbytes()];
  let mut rho = [0u8; SEEDBYTES];
  let mut mu = [0u8; CRHBYTES];
  let mut c = [0u8; SEEDBYTES];
  let mut c2 = [0u8; SEEDBYTES];
  let mut cp = Poly::default();
  let mut mat_data = P::mat();
  let mut mat = mat_from(&mut mat_data);
  let mut z = P::polyveclnew();
  let (mut t1, mut w1, mut h) = (
    P::polyvecknew(),
    P::polyvecknew(),
    P::polyvecknew()
  );
  let mut state = KeccakState::default(); // shake256_init()

  if sig.len() != P::SIGNBYTES {
    return Err("Signature length mismatch".to_string());
  }

  unpack_pk::<P>(&mut rho, &mut t1, pk);
  if let Err(e) = unpack_sig::<P>(&mut c, &mut z, &mut h, sig) {
    return Err(e);
  }
  if polyvecl_chknorm::<P>(&z, (P::GAMMA1 - P::BETA) as i32) > 0 {
    return Err("Invalid z".to_string());
  }

  // Compute CRH(CRH(rho, t1), msg)
  shake256(&mut mu, pk, SEEDBYTES, P::PUBLIC_KEY_BYTES);
  state.shake256_absorb(&mu, SEEDBYTES);
  state.shake256_absorb(m, m.len());
  state.shake256_finalize();
  state.shake256_squeeze(&mut mu, CRHBYTES);

  // Matrix-vector multiplication; compute Az - c2^dt1
  poly_challenge::<P>(&mut cp, &c);
  polyvec_matrix_expand::<P>(&mut mat, &rho);

  polyvecl_ntt::<P>(&mut z);
  polyvec_matrix_pointwise_montgomery::<P>(&mut w1, &mat, &z);

  poly_ntt::<P>(&mut cp);
  polyveck_shiftl(&mut t1);
  polyveck_ntt::<P>(&mut t1);
  let t1_2 = t1.clone();
  polyveck_pointwise_poly_montgomery::<P>(&mut t1, &cp, &t1_2);
  let w1_clone = w1.clone();
  polyveck_sub(&mut w1, &w1_clone, &t1);
  polyveck_reduce::<P>(&mut w1);
  polyveck_invntt_tomont::<P>(&mut w1);

  // Reconstruct w1
  polyveck_caddq::<P>(&mut w1);
  let w1_clone = w1.clone();
  polyveck_use_hint::<P>(&mut w1, &w1_clone, &h);
  polyveck_pack_w1::<P>(&mut buf, &w1);

  // Call random oracle and verify challenge
  state.init();
  state.shake256_absorb(&mu, CRHBYTES);
  state.shake256_absorb(&buf, P::K * P::polyw1_packedbytes());
  state.shake256_finalize();
  state.shake256_squeeze(&mut c2, SEEDBYTES);
  // Doesn't require constant time equality check
  if c != c2 {
    Err("Invalid signature".to_string())
  } else {
    Ok(())
  }
}

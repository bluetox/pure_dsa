use crate::params::DilithiumParams;
use crate::polyvec::*;
use crate::poly::*;

pub fn pack_pk<P: DilithiumParams>(pk: &mut [u8], rho: &[u8], t1: &Polyveck) {

    // Copy rho into beginning of pk
    pk[..P::SEEDBYTES].copy_from_slice(rho);

    // Offset pointer
    let mut offset = P::SEEDBYTES;

    // Pack each t1[i]
    let t1_vec = vec_from_polyveck(t1);
    for i in 0..P::K {
        polyt1_pack(
            &mut pk[offset..offset + POLYT1_PACKEDBYTES],
            &t1_vec[i],
        );
        offset += POLYT1_PACKEDBYTES;
    }
}
pub fn pack_sk<P: DilithiumParams>(
  sk: &mut [u8],
  rho: &[u8],
  tr: &[u8],
  key: &[u8],
  t0: &Polyveck,
  s1: &Polyvecl,
  s2: &Polyveck,
) {
  let mut idx = 0usize;

    let s1_vec = vec_from_polyvecl(s1);
    let s2_vec = vec_from_polyveck(s2);
    let t0_vec = vec_from_polyveck(t0);

  sk[idx..SEEDBYTES].copy_from_slice(&rho[0..SEEDBYTES]);
  idx += SEEDBYTES;

  sk[idx..idx + SEEDBYTES].copy_from_slice(&key[0..SEEDBYTES]);
  idx += SEEDBYTES;

  sk[idx..idx + SEEDBYTES].copy_from_slice(&tr[0..SEEDBYTES]);

  idx += SEEDBYTES;

  for i in 0..P::L {
    polyeta_pack::<P>(&mut sk[idx + i * P::POLYETA_PACKEDBYTES..], &s1_vec[i]);
  }
  idx += P::L * P::POLYETA_PACKEDBYTES;

  for i in 0..P::K {
    polyeta_pack::<P>(&mut sk[idx + i * P::POLYETA_PACKEDBYTES..], &s2_vec[i]);
  }
  idx += P::K * P::POLYETA_PACKEDBYTES;

  for i in 0..P::K {
    polyt0_pack(&mut sk[idx + i * POLYT0_PACKEDBYTES..], &t0_vec[i]);
  }
}

pub fn unpack_sk<P: DilithiumParams>(
  rho: &mut [u8],
  tr: &mut [u8],
  key: &mut [u8],
  t0: &mut Polyveck,
  s1: &mut Polyvecl,
  s2: &mut Polyveck,
  sk: &[u8],
) {
  let mut idx = 0usize;
let s1_vec = vec_from_polyvecl_mut(s1);
let s2_vec = vec_from_polyveck_mut(s2);
let t0_vec = vec_from_polyveck_mut(t0);
  rho[..SEEDBYTES].copy_from_slice(&sk[..SEEDBYTES]);
  idx += SEEDBYTES;

  key[..SEEDBYTES].copy_from_slice(&sk[idx..idx + SEEDBYTES]);
  idx += SEEDBYTES;

  tr[..SEEDBYTES].copy_from_slice(&sk[idx..idx + SEEDBYTES]);
  idx += SEEDBYTES;

for i in 0..P::L {
    let start = idx + i * P::POLYETA_PACKEDBYTES;
    let end = start + P::POLYETA_PACKEDBYTES;
    let input_slice = &sk[start..end];
    polyeta_unpack::<P>(&mut s1_vec[i], input_slice);
}
  idx += P::L * P::POLYETA_PACKEDBYTES;

  for i in 0..P::K {
    polyeta_unpack::<P>(&mut s2_vec[i], &sk[idx + i * P::POLYETA_PACKEDBYTES..]);
  }
  idx += P::K * P::POLYETA_PACKEDBYTES;

  for i in 0..P::K {
    polyt0_unpack(&mut t0_vec[i], &sk[idx + i * POLYT0_PACKEDBYTES..]);
  }
}


/// Bit-pack signature sig = (c, z, h).
pub fn pack_sig<P: DilithiumParams>(sig: &mut [u8], c: Option<&[u8]>, z: &Polyvecl, h: &Polyveck) {
    let z_vec = vec_from_polyvecl(z);
    let h_vec = vec_from_polyveck(h);
  let mut idx = 0usize;

  if let Some(challenge) = c {
    sig[..SEEDBYTES].copy_from_slice(&challenge[..SEEDBYTES]);
  }

  idx += SEEDBYTES;

  for i in 0..P::L {
    polyz_pack::<P>(&mut sig[idx + i * P::polyz_packedbytes()..], &z_vec[i]);
  }
  idx += P::L * P::polyz_packedbytes();
  // Encode H
  sig[idx..idx + P::OMEGA + P::K].copy_from_slice(&vec![0u8; P::OMEGA + P::K]);

  let mut k = 0;
  for i in 0..P::K {
    for j in 0..N {
      if h_vec[i].coeffs[j] != 0 {
        sig[idx + k] = j as u8;
        k += 1;
      }
    }
    sig[idx + P::OMEGA + i] = k as u8;
  }
}

/// Unpack signature sig = (z, h, c).
pub fn unpack_sig<P: DilithiumParams>(
  c: &mut [u8],
  z: &mut Polyvecl,
  h: &mut Polyveck,
  sig: &[u8],
) -> Result<(), String> {
  let mut idx = 0usize;
    let z_vec = vec_from_polyvecl_mut(z);
    let h_vec = vec_from_polyveck_mut(h);
  c[..SEEDBYTES].copy_from_slice(&sig[..SEEDBYTES]);
  idx += SEEDBYTES;

  for i in 0..P::L {
    polyz_unpack::<P>(&mut z_vec[i], &sig[idx + i * P::polyz_packedbytes()..]);
  }
  idx += P::L * P::polyz_packedbytes();

  // Decode h
  let mut k = 0usize;
  for i in 0..P::K {
    if sig[idx + P::OMEGA + i] < k as u8 || sig[idx + P::OMEGA + i] > P::OMEGA as u8 {
      return Err("INVALID OMEGA".to_string());
    }
    for j in k..sig[idx + P::OMEGA + i] as usize {
      // Coefficients are ordered for strong unforgeability
      if j > k && sig[idx + j as usize] <= sig[idx + j as usize - 1] {
        return Err("INVALID H".to_string());
      }
      h_vec[i].coeffs[sig[idx + j] as usize] = 1;
    }
    k = sig[idx + P::OMEGA + i] as usize;
  }

  // Extra indices are zero for strong unforgeability
  for j in k..P::OMEGA {
    if sig[idx + j as usize] > 0 {
      return Err("INVALID H".to_string());
    }
  }

  Ok(())
}

pub fn unpack_pk<P: DilithiumParams>(rho: &mut [u8], t1: &mut Polyveck, pk: &[u8]) {
  rho[..SEEDBYTES].copy_from_slice(&pk[..SEEDBYTES]);
  let t1_vec = vec_from_polyveck_mut(t1);
  for i in 0..P::K {
    polyt1_unpack(&mut t1_vec[i], &pk[SEEDBYTES + i * POLYT1_PACKEDBYTES..])
  }
}
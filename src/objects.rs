use crate::sign::*;
use crate::params::*;
use rand::{
    rngs::OsRng,
    RngCore
};

#[cfg(feature = "zeroize")]
use zeroize::{Zeroize, ZeroizeOnDrop};


pub enum Algorithm {
    Mode2,
    Mode3,
    Mode5,
}

impl Algorithm {
    pub fn generate(&self) -> Keypair {
        match self {
            Algorithm::Mode2 => {
                let mut pk = [0u8; Mode2::PUBLIC_KEY_BYTES];
                let mut sk = [0u8; Mode2::SECRET_KEY_BYTES];
                let mut rng = OsRng;
                crypto_sign_keypair::<Mode2, OsRng>(&mut pk, &mut sk, &mut rng);
                Keypair::Mode2(pk, sk)
            }
            Algorithm::Mode3 => {
                let mut pk = [0u8; Mode3::PUBLIC_KEY_BYTES];
                let mut sk = [0u8; Mode3::SECRET_KEY_BYTES];
                let mut rng = OsRng;
                crypto_sign_keypair::<Mode3, _>(&mut pk, &mut sk, &mut rng);
                Keypair::Mode3(pk, sk)
            }
            Algorithm::Mode5 => {
                let mut pk = [0u8; Mode5::PUBLIC_KEY_BYTES];
                let mut sk = [0u8; Mode5::SECRET_KEY_BYTES];
                let mut rng = OsRng;
                crypto_sign_keypair::<Mode5, _>(&mut pk, &mut sk, &mut rng);
                Keypair::Mode5(pk, sk)
            }
        }
    }

    pub fn load_from_bytes(&self, sk: &[u8], pk: &[u8]) -> Result<Keypair, String> {
        match self {
            Algorithm::Mode2 => {
                if pk.len() != Mode2::PUBLIC_KEY_BYTES || sk.len() != Mode2::SECRET_KEY_BYTES {
                    return Err("Invalid key lengths for Mode2".into());
                }
                let mut pk_buf = [0u8; Mode2::PUBLIC_KEY_BYTES];
                let mut sk_buf = [0u8; Mode2::SECRET_KEY_BYTES];
                pk_buf.copy_from_slice(pk);
                sk_buf.copy_from_slice(sk);
                Ok(Keypair::Mode2(pk_buf, sk_buf))
            }
            Algorithm::Mode3 => {
                if pk.len() != Mode3::PUBLIC_KEY_BYTES || sk.len() != Mode3::SECRET_KEY_BYTES {
                    return Err("Invalid key lengths for Mode3".into());
                }
                let mut pk_buf = [0u8; Mode3::PUBLIC_KEY_BYTES];
                let mut sk_buf = [0u8; Mode3::SECRET_KEY_BYTES];
                pk_buf.copy_from_slice(pk);
                sk_buf.copy_from_slice(sk);
                Ok(Keypair::Mode3(pk_buf, sk_buf))
            }
            Algorithm::Mode5 => {
                if pk.len() != Mode5::PUBLIC_KEY_BYTES || sk.len() != Mode5::SECRET_KEY_BYTES {
                    return Err("Invalid key lengths for Mode5".into());
                }
                let mut pk_buf = [0u8; Mode5::PUBLIC_KEY_BYTES];
                let mut sk_buf = [0u8; Mode5::SECRET_KEY_BYTES];
                pk_buf.copy_from_slice(pk);
                sk_buf.copy_from_slice(sk);
                Ok(Keypair::Mode5(pk_buf, sk_buf))
            }
        }
    }

    pub fn generate_with_rng<R: RngCore>(&self, rng: &mut R) -> Keypair {
        match self {
            Algorithm::Mode2 => {
                let mut pk = [0u8; Mode2::PUBLIC_KEY_BYTES];
                let mut sk = [0u8; Mode2::SECRET_KEY_BYTES];
                crypto_sign_keypair::<Mode2, R>(&mut pk, &mut sk, rng);
                Keypair::Mode2(pk, sk)
            }
            Algorithm::Mode3 => {
                let mut pk = [0u8; Mode3::PUBLIC_KEY_BYTES];
                let mut sk = [0u8; Mode3::SECRET_KEY_BYTES];
                crypto_sign_keypair::<Mode3, R>(&mut pk, &mut sk, rng);
                Keypair::Mode3(pk, sk)
            }
            Algorithm::Mode5 => {
                let mut pk = [0u8; Mode5::PUBLIC_KEY_BYTES];
                let mut sk = [0u8; Mode5::SECRET_KEY_BYTES];
                crypto_sign_keypair::<Mode5, R>(&mut pk, &mut sk, rng);
                Keypair::Mode5(pk, sk)
            }
        }
    }

    pub fn verify(&self, signature: &Signature, msg: &[u8], public_key: &[u8]) -> Result<(), &'static str> {
        match self {
            Algorithm::Mode2 => {
                crypto_sign_verify::<Mode2>(signature.bytes(), msg, public_key)
            },
            Algorithm::Mode3 => {
                crypto_sign_verify::<Mode3>(signature.bytes(), msg, public_key)
            },
            Algorithm::Mode5 => {
                crypto_sign_verify::<Mode5>(signature.bytes(), msg, public_key)
            },
        }
    }

    pub fn verify_raw(&self, signature: &[u8], msg: &[u8], public_key: &[u8]) -> Result<(), &'static str> {
        match self {
            Algorithm::Mode2 => crypto_sign_verify::<Mode2>(signature, msg, public_key),
            Algorithm::Mode3 => crypto_sign_verify::<Mode3>(signature, msg, public_key),
            Algorithm::Mode5 => crypto_sign_verify::<Mode5>(signature, msg, public_key),
        }
    }
}

#[cfg_attr(feature = "zeroize", derive(Zeroize, ZeroizeOnDrop))]
#[derive(Clone, Debug)]
pub enum Keypair {
    Mode2([u8; Mode2::PUBLIC_KEY_BYTES], [u8; Mode2::SECRET_KEY_BYTES]),
    Mode3([u8; Mode3::PUBLIC_KEY_BYTES], [u8; Mode3::SECRET_KEY_BYTES]),
    Mode5([u8; Mode5::PUBLIC_KEY_BYTES], [u8; Mode5::SECRET_KEY_BYTES]),
}

impl Keypair {
    pub fn sign(&self, msg: &[u8]) -> Signature {
        match self {
            Keypair::Mode2(_, sk) => {
                let mut sig = [0u8; Mode2::SIGNBYTES];
                crypto_sign_signature::<Mode2, _>(&mut sig, msg, sk, &mut OsRng);
                Signature {
                    bytes: SignType::SignMode2(sig)
                }
            }
            Keypair::Mode3(_, sk) => {
                let mut sig = [0u8; Mode3::SIGNBYTES];
                crypto_sign_signature::<Mode3, _>(&mut sig, msg, sk, &mut OsRng);
                Signature {
                    bytes: SignType::SignMode3(sig)
                }
            }
            Keypair::Mode5(_, sk) => {
                let mut sig = [0u8; Mode5::SIGNBYTES];
                crypto_sign_signature::<Mode5, _>(&mut sig, msg, sk, &mut OsRng);
                Signature {
                    bytes: SignType::SignMode5(sig),
                }
            }
        }
    }
    #[cfg(not(feature = "no_std"))]
    pub fn sign_to_slice(&self, msg: &[u8], sk: &[u8]) -> Vec<u8> {
        match self {
            Keypair::Mode2(_, _) => {
                let mut sig = [0u8; Mode2::SIGNBYTES];
                crypto_sign_signature::<Mode2, _>(&mut sig, msg, sk, &mut OsRng);
                sig.to_vec()
            }
            Keypair::Mode3(_, _) => {
                let mut sig = [0u8; Mode3::SIGNBYTES];
                crypto_sign_signature::<Mode3, _>(&mut sig, msg, sk, &mut OsRng);
                sig.to_vec()
            }
            Keypair::Mode5(_, _) => {
                let mut sig = [0u8; Mode5::SIGNBYTES];
                crypto_sign_signature::<Mode5, _>(&mut sig, msg, sk, &mut OsRng);
                sig.to_vec()
            }
        }
    }
        
    pub fn public(&self) -> &[u8] {
        match self {
            Keypair::Mode2(pk, _) => pk,
            Keypair::Mode3(pk, _) => pk,
            Keypair::Mode5(pk, _) => pk,
        }
    }

    pub fn secret(&self) -> &[u8] {
        match self {
            Keypair::Mode2(_, sk) => sk,
            Keypair::Mode3(_, sk) => sk,
            Keypair::Mode5(_, sk) => sk,
        }
    }
}
#[cfg_attr(feature = "zeroize", derive(Zeroize, ZeroizeOnDrop))]
pub struct Signature {
    bytes: SignType
}

impl Signature {
    #[inline(always)]
    pub fn bytes(&self) -> &[u8] {
        match &self.bytes {
            SignType::SignMode2(arr) => &arr[..],
            SignType::SignMode3(arr) => &arr[..],
            SignType::SignMode5(arr) => &arr[..],
        }
    }
}

#[cfg_attr(feature = "zeroize", derive(Zeroize, ZeroizeOnDrop))]
enum SignType {
    SignMode2([u8; 2420]),
    SignMode3([u8; 3293]),
    SignMode5([u8; 4595])
}

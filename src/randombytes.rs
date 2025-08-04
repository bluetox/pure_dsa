use rand::{
    RngCore
};


pub fn randombytes<R: RngCore>(out: &mut [u8], rng: &mut R) {
    R::fill_bytes(rng, out)
}

use crate::fips202::KeccakState;


impl KeccakState {
  pub fn dilithium_shake256_stream_init(
    &mut self,
    seed: &[u8],
    nonce: u16,
  ) {
    let t = [nonce as u8, (nonce >> 8) as u8];
    self.init();
    self.shake256_absorb(seed, 64);
    self.shake256_absorb( &t, 2);
    self.shake256_finalize();
  }
  pub fn dilithium_shake128_stream_init(
    &mut self,
    seed: &[u8],
    nonce: u16,
  ) {
    let t = [nonce as u8, (nonce >> 8) as u8];
    self.init();
    self.shake128_absorb(seed);
    self.shake128_absorb( &t);
    self.shake128_finalize();
  }
  
}
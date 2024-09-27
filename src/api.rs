#[derive(Copy, Clone, PartialEq, Eq, Hash)]
pub struct Keypair<const PUBLICKEYBYTES: usize, const SECRETKEYBYTES: usize> {
    pub public: [u8; PUBLICKEYBYTES],
    secret: [u8; SECRETKEYBYTES],
}

/// Secret key elided
impl<const PUBLICKEYBYTES: usize, const SECRETKEYBYTES: usize> std::fmt::Debug for Keypair<PUBLICKEYBYTES, SECRETKEYBYTES> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "public: {:?}\nsecret: <elided>", self.public)
    }
}

pub enum SignError {
    Input,
    Verify,
}

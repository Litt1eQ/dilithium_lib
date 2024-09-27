pub const SEEDBYTES: usize = 32;
pub const TRBYTES: usize = 64;
pub const RNDBYTES: usize = 32;
pub const CRHBYTES: usize = 64;
pub const N: usize = 256;
pub const Q: usize = 8380417;
pub const D: usize = 13;
// pub const ROOT_OF_UNITY: usize = 1753;

pub const POLYT1_PACKEDBYTES: usize = 320;
pub const POLYT0_PACKEDBYTES: usize = 416;

// Concise types to avoid cast cluttering
pub const Q_I32: i32 = Q as i32;
pub const N_U32: u32 = N as u32;
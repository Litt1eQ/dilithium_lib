use crate::poly::*;

#[derive(Copy, Clone)]
pub struct Polyveck<const K: usize> {
    pub vec: [Poly; K],
}

impl<const K: usize> Default for Polyveck<K> {
    fn default() -> Self {
        Polyveck {
            vec: [Poly::default(); K],
        }
    }
}

#[derive(Copy, Clone)]
pub struct Polyvecl<const L: usize> {
    pub vec: [Poly; L],
}

impl<const L: usize> Default for Polyvecl<L> {
    fn default() -> Self {
        Polyvecl {
            vec: [Poly::default(); L],
        }
    }
}

/// Implementation of ExpandA. Generates matrix A with uniformly
/// random coefficients a_{i,j} by performing rejection
/// sampling on the output stream of SHAKE128(rho|j|i)
/// or AES256CTR(rho,j|i).
pub fn polyvec_matrix_expand<const K: usize, const L: usize, const POLY_UNIFORM_NBLOCKS: usize, const STREAM128_BLOCKBYTES: usize>(mat: &mut [Polyvecl<L>], rho: &[u8]) {
    for i in 0..K {
        for j in 0..L {
            poly_uniform::<POLY_UNIFORM_NBLOCKS, STREAM128_BLOCKBYTES>(&mut mat[i].vec[j], rho, ((i << 8) + j) as u16);
        }
    }
}

pub fn polyvec_matrix_pointwise_montgomery<const K: usize, const L: usize>(
    t: &mut Polyveck<K>,
    mat: &[Polyvecl<L>],
    v: &Polyvecl<L>,
) {
    for i in 0..K {
        polyvecl_pointwise_acc_montgomery(&mut t.vec[i], &mat[i], v);
    }
}

//*********** Vectors of polynomials of length L ****************************

pub fn polyvecl_uniform_eta<const L: usize, const POLY_UNIFORM_ETA_NBLOCKS: usize, const ETA: usize>(v: &mut Polyvecl<L>, seed: &[u8], mut nonce: u16) {
    for i in 0..L {
        poly_uniform_eta::<POLY_UNIFORM_ETA_NBLOCKS, ETA>(&mut v.vec[i], seed, nonce);
        nonce += 1;
    }
}

pub fn polyvecl_uniform_gamma1<const L: usize, const POLY_UNIFORM_GAMMA1_NBLOCKS: usize, const GAMMA1: usize>(v: &mut Polyvecl<L>, seed: &[u8], nonce: u16) {
    for i in 0..L {
        poly_uniform_gamma1::<POLY_UNIFORM_GAMMA1_NBLOCKS, GAMMA1>(&mut v.vec[i], seed, (L as u16) * nonce + i as u16);
    }
}
pub fn polyvecl_reduce<const L: usize>(v: &mut Polyvecl<L>) {
    for i in 0..L {
        poly_reduce(&mut v.vec[i]);
    }
}

/// Add vectors of polynomials of length L.
/// No modular reduction is performed.
pub fn polyvecl_add<const L: usize>(w: &mut Polyvecl<L>, v: &Polyvecl<L>) {
    for i in 0..L {
        poly_add(&mut w.vec[i], &v.vec[i]);
    }
}

/// Forward NTT of all polynomials in vector of length L. Output
/// coefficients can be up to 16*Q larger than input coefficients.*
pub fn polyvecl_ntt<const L: usize>(v: &mut Polyvecl<L>) {
    for i in 0..L {
        poly_ntt(&mut v.vec[i]);
    }
}

pub fn polyvecl_invntt_tomont<const L: usize>(v: &mut Polyvecl<L>) {
    for i in 0..L {
        poly_invntt_tomont(&mut v.vec[i]);
    }
}

pub fn polyvecl_pointwise_poly_montgomery<const L: usize>(
    r: &mut Polyvecl<L>,
    a: &Poly,
    v: &Polyvecl<L>,
) {
    for i in 0..L {
        poly_pointwise_montgomery(&mut r.vec[i], a, &v.vec[i]);
    }
}

/// Pointwise multiply vectors of polynomials of length L, multiply
/// resulting vector by 2^{-32} and add (accumulate) polynomials
/// in it. Input/output vectors are in NTT domain representation.
/// Input coefficients are assumed to be less than 22*Q. Output
/// coeffcient are less than 2*L*Q.
pub fn polyvecl_pointwise_acc_montgomery<const L: usize>(
    w: &mut Poly,
    u: &Polyvecl<L>,
    v: &Polyvecl<L>,
) {
    let mut t = Poly::default();
    poly_pointwise_montgomery(w, &u.vec[0], &v.vec[0]);
    for i in 1..L {
        poly_pointwise_montgomery(&mut t, &u.vec[i], &v.vec[i]);
        poly_add(w, &t);
    }
}

/// Check infinity norm of polynomials in vector of length L.
/// Assumes input coefficients to be standard representatives.
/// Returns 0 if norm of all polynomials is strictly smaller than B and 1
/// otherwise.
pub fn polyvecl_chknorm<const L: usize>(v: &Polyvecl<L>, bound: i32) -> u8 {
    for i in 0..L {
        if poly_chknorm(&v.vec[i], bound) > 0 {
            return 1;
        }
    }
    return 0;
}

//*********** Vectors of polynomials of length K ****************************

pub fn polyveck_uniform_eta<const K: usize, const POLY_UNIFORM_ETA_NBLOCKS: usize, const ETA: usize>(v: &mut Polyveck<K>, seed: &[u8], mut nonce: u16) {
    for i in 0..K {
        poly_uniform_eta::<POLY_UNIFORM_ETA_NBLOCKS, ETA>(&mut v.vec[i], seed, nonce);
        nonce += 1
    }
}

/// Reduce coefficients of polynomials in vector of length K
/// to representatives in [0,2*Q].
pub fn polyveck_reduce<const K: usize>(v: &mut Polyveck<K>) {
    for i in 0..K {
        poly_reduce(&mut v.vec[i]);
    }
}

/// For all coefficients of polynomials in vector of length K
/// add Q if coefficient is negative.
pub fn polyveck_caddq<const K: usize>(v: &mut Polyveck<K>) {
    for i in 0..K {
        poly_caddq(&mut v.vec[i]);
    }
}

/// Add vectors of polynomials of length K.
/// No modular reduction is performed.
pub fn polyveck_add<const K: usize>(w: &mut Polyveck<K>, v: &Polyveck<K>) {
    for i in 0..K {
        poly_add(&mut w.vec[i], &v.vec[i]);
    }
}

/// Subtract vectors of polynomials of length K.
/// Assumes coefficients of polynomials in second input vector
/// to be less than 2*Q. No modular reduction is performed.
pub fn polyveck_sub<const K: usize>(w: &mut Polyveck<K>, v: &Polyveck<K>) {
    for i in 0..K {
        poly_sub(&mut w.vec[i], &v.vec[i]);
    }
}

/// Multiply vector of polynomials of Length K by 2^D without modular
/// reduction. Assumes input coefficients to be less than 2^{32-D}.
pub fn polyveck_shiftl<const K: usize>(v: &mut Polyveck<K>) {
    for i in 0..K {
        poly_shiftl(&mut v.vec[i]);
    }
}

/// Forward NTT of all polynomials in vector of length K. Output
/// coefficients can be up to 16*Q larger than input coefficients.
pub fn polyveck_ntt<const K: usize>(v: &mut Polyveck<K>) {
    for i in 0..K {
        poly_ntt(&mut v.vec[i]);
    }
}

/// Inverse NTT and multiplication by 2^{32} of polynomials
/// in vector of length K. Input coefficients need to be less
/// than 2*Q.
pub fn polyveck_invntt_tomont<const K: usize>(v: &mut Polyveck<K>) {
    for i in 0..K {
        poly_invntt_tomont(&mut v.vec[i]);
    }
}

pub fn polyveck_pointwise_poly_montgomery<const K: usize>(
    r: &mut Polyveck<K>,
    a: &Poly,
    v: &Polyveck<K>,
) {
    for i in 0..K {
        poly_pointwise_montgomery(&mut r.vec[i], a, &v.vec[i]);
    }
}

/// Check infinity norm of polynomials in vector of length K.
/// Assumes input coefficients to be standard representatives.
//
/// Returns 0 if norm of all polynomials are strictly smaller than B and 1
/// otherwise.
pub fn polyveck_chknorm<const K: usize>(v: &Polyveck<K>, bound: i32) -> u8 {
    for i in 0..K {
        if poly_chknorm(&v.vec[i], bound) > 0 {
            return 1;
        }
    }
    return 0;
}

/// For all coefficients a of polynomials in vector of length K,
/// compute a0, a1 such that a mod Q = a1*2^D + a0
/// with -2^{D-1} < a0 <= 2^{D-1}. Assumes coefficients to be
/// standard representatives.
pub fn polyveck_power2round<const K: usize>(v1: &mut Polyveck<K>, v0: &mut Polyveck<K>) {
    for i in 0..K {
        poly_power2round(&mut v1.vec[i], &mut v0.vec[i]);
    }
}

/// For all coefficients a of polynomials in vector of length K,
/// compute high and low bits a0, a1 such a mod Q = a1*ALPHA + a0
/// with -ALPHA/2 < a0 <= ALPHA/2 except a1 = (Q-1)/ALPHA where we
/// set a1 = 0 and -ALPHA/2 <= a0 = a mod Q - Q < 0.
/// Assumes coefficients to be standard representatives.
pub fn polyveck_decompose<const K: usize, const GAMMA2: usize>(v1: &mut Polyveck<K>, v0: &mut Polyveck<K>) {
    for i in 0..K {
        poly_decompose::<GAMMA2>(&mut v1.vec[i], &mut v0.vec[i]);
    }
}

/// Compute hint vector.
///
/// Returns number of 1 bits.
pub fn polyveck_make_hint<const K: usize, const GAMMA2: usize>(
    h: &mut Polyveck<K>,
    v0: &Polyveck<K>,
    v1: &Polyveck<K>,
) -> i32 {
    let mut s = 0i32;
    for i in 0..K {
        s += poly_make_hint::<GAMMA2>(&mut h.vec[i], &v0.vec[i], &v1.vec[i]);
    }
    s
}

/// Use hint vector to correct the high bits of input vector.
pub fn polyveck_use_hint<const K: usize, const GAMMA2: usize>(w: &mut Polyveck<K>, h: &Polyveck<K>) {
    for i in 0..K {
        poly_use_hint::<GAMMA2>(&mut w.vec[i], &h.vec[i]);
    }
}

pub fn polyveck_pack_w1<const K: usize, const POLYW1_PACKEDBYTES: usize, const GAMMA2: usize>(r: &mut [u8], w1: &Polyveck<K>) {
    for i in 0..K {
        polyw1_pack::<GAMMA2>(&mut r[i * POLYW1_PACKEDBYTES..], &w1.vec[i]);
    }
}

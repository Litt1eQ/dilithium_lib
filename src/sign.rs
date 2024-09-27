use crate::{
    fips202::*, packing::*, params::*, poly::*, polyvec::*, randombytes::*,
};
use crate::api::SignError;

pub fn crypto_sign_keypair<const K: usize, const L: usize, const PUBLICKEYBYTES: usize, const POLY_UNIFORM_NBLOCKS: usize, const STREAM128_BLOCKBYTES: usize, const POLY_UNIFORM_ETA_NBLOCKS: usize, const ETA: usize, const POLYETA_PACKEDBYTES: usize>(
    pk: &mut [u8],
    sk: &mut [u8],
    seed: Option<&[u8]>,
) -> u8 {
    let mut init_seed = [0u8; SEEDBYTES];
    match seed {
        Some(x) => init_seed.copy_from_slice(x),
        None => randombytes(&mut init_seed, SEEDBYTES),
    };
    let mut seedbuf = [0u8; 2 * SEEDBYTES + CRHBYTES];

    let mut new_init_seed = [0u8; SEEDBYTES + 2];
    new_init_seed[..SEEDBYTES].copy_from_slice(&init_seed);
    new_init_seed[SEEDBYTES + 0] = K as u8;
    new_init_seed[SEEDBYTES + 1] = L as u8;

    let mut tr = [0u8; TRBYTES];
    let (mut rho, mut rhoprime, mut key) =
        ([0u8; SEEDBYTES], [0u8; CRHBYTES], [0u8; SEEDBYTES]);
    let mut mat: [Polyvecl<L>; K] = [Polyvecl::<L>::default(); K];
    let mut s1 = Polyvecl::default();
    let (mut s2, mut t1, mut t0) = (
        Polyveck::<K>::default(),
        Polyveck::<K>::default(),
        Polyveck::<K>::default(),
    );

    // Get randomness for rho, rhoprime and key
    shake256(
        &mut seedbuf,
        2 * SEEDBYTES + CRHBYTES,
        &new_init_seed,
        SEEDBYTES + 2,
    );
    rho.copy_from_slice(&seedbuf[..SEEDBYTES]);
    rhoprime.copy_from_slice(&seedbuf[SEEDBYTES..SEEDBYTES + CRHBYTES]);
    key.copy_from_slice(&seedbuf[SEEDBYTES + CRHBYTES..]);

    // Expand matrix
    polyvec_matrix_expand::<K, L, POLY_UNIFORM_NBLOCKS, STREAM128_BLOCKBYTES>(&mut mat, &rho);
    // Sample short vectors s1 and s2
    polyvecl_uniform_eta::<L, POLY_UNIFORM_ETA_NBLOCKS, ETA>(&mut s1, &rhoprime, 0);
    polyveck_uniform_eta::<K, POLY_UNIFORM_ETA_NBLOCKS, ETA>(&mut s2, &rhoprime, L as u16);

    // Matrix-vector multiplication
    let mut s1hat = s1;
    polyvecl_ntt(&mut s1hat);

    polyvec_matrix_pointwise_montgomery(&mut t1, &mat, &s1hat);
    polyveck_reduce(&mut t1);
    polyveck_invntt_tomont(&mut t1);

    // Add error vector s2
    polyveck_add(&mut t1, &s2);
    // Extract t1 and write public key
    polyveck_caddq(&mut t1);
    polyveck_power2round(&mut t1, &mut t0);
    pack_pk(pk, &rho, &t1);

    // Compute H(rho, t1) and write secret key
    shake256(&mut tr, TRBYTES, pk, PUBLICKEYBYTES);
    pack_sk::<K, L, POLYETA_PACKEDBYTES, ETA>(sk, &rho, &tr, &key, &t0, &s1, &s2);

    0
}

pub fn crypto_sign_signature<const K: usize, const L: usize, const POLYW1_PACKEDBYTES: usize, const CTILDEBYTES: usize, const GAMMA1: usize, const GAMMA2: usize, const BETA: usize, const OMEGA: usize, const POLYETA_PACKEDBYTES: usize, const POLY_UNIFORM_NBLOCKS: usize, const STREAM128_BLOCKBYTES: usize, const POLY_UNIFORM_GAMMA1_NBLOCKS: usize, const TAU: usize, const POLYZ_PACKEDBYTES: usize, const ETA: usize>
(sig: &mut [u8], m: &[u8], sk: &[u8], ctx: Option<&[u8]>, using_randomized_signing: bool) {
    // `key` and `mu` are concatenated
    let mut keymu = [0u8; SEEDBYTES + RNDBYTES + CRHBYTES];

    let mut nonce = 0u16;
    let mut mat: [Polyvecl<L>; K] = [Polyvecl::default(); K];
    let (mut s1, mut y) = (Polyvecl::default(), Polyvecl::default());
    let (mut s2, mut t0) = (Polyveck::<K>::default(), Polyveck::default());
    let (mut w1, mut w0) = (Polyveck::<K>::default(), Polyveck::default());
    let mut h = Polyveck::default();
    let mut cp = Poly::default();
    let mut state = KeccakState::default(); //shake256_init()
    let mut rho = [0u8; SEEDBYTES];
    let mut tr = [0u8; TRBYTES];
    let mut rhoprime = [0u8; CRHBYTES];

    let ctx = match ctx {
        Some(x) => x,
        None => &[0u8; 0],
    };

    if ctx.len() > 255 {
        panic!("ctx is too long");
    }

    unpack_sk::<K, L, POLYETA_PACKEDBYTES, ETA>(
        &mut rho,
        &mut tr,
        &mut keymu[..SEEDBYTES],
        &mut t0,
        &mut s1,
        &mut s2,
        &sk,
    );

    // Compute CRH(tr, 0, ctxlen, ctx, msg)
    // keymu[0] = 0;
    // keymu[1] = ctxLen;
    let mu = [0u8, ctx.len() as u8];
    shake256_absorb(&mut state, &tr, TRBYTES);
    shake256_absorb(&mut state, &mu, 2);
    shake256_absorb(&mut state, &ctx, ctx.len());
    shake256_absorb(&mut state, m, m.len());
    shake256_finalize(&mut state);
    shake256_squeeze(&mut keymu[SEEDBYTES + RNDBYTES..], CRHBYTES, &mut state);

    if using_randomized_signing {
        randombytes(&mut rhoprime, CRHBYTES);
    } else {
        shake256(&mut rhoprime, CRHBYTES, &keymu, SEEDBYTES + RNDBYTES + CRHBYTES);
    }

    // Expand matrix and transform vectors
    polyvec_matrix_expand::<K, L, POLY_UNIFORM_NBLOCKS, STREAM128_BLOCKBYTES>(&mut mat, &rho);
    polyvecl_ntt(&mut s1);
    polyveck_ntt(&mut s2);
    polyveck_ntt(&mut t0);

    loop {
        // Sample intermediate vector y
        polyvecl_uniform_gamma1::<L, POLY_UNIFORM_GAMMA1_NBLOCKS, GAMMA1>(&mut y, &rhoprime, nonce);
        nonce += 1;

        // Matrix-vector multiplication
        let mut z = y;
        polyvecl_ntt(&mut z);
        polyvec_matrix_pointwise_montgomery(&mut w1, &mat, &z);
        polyveck_reduce::<K>(&mut w1);
        polyveck_invntt_tomont::<K>(&mut w1);

        // Decompose w and call the random oracle
        polyveck_caddq::<K>(&mut w1);
        polyveck_decompose::<K, GAMMA2>(&mut w1, &mut w0);
        polyveck_pack_w1::<K, POLYW1_PACKEDBYTES, GAMMA2>(sig, &w1);

        state.init();
        shake256_absorb(&mut state, &keymu[SEEDBYTES + RNDBYTES..], CRHBYTES);
        shake256_absorb(&mut state, &sig, K * POLYW1_PACKEDBYTES);
        shake256_finalize(&mut state);
        shake256_squeeze(sig, CTILDEBYTES, &mut state);
        poly_challenge::<CTILDEBYTES, TAU>(&mut cp, sig);
        poly_ntt(&mut cp);

        // Compute z, reject if it reveals secret
        polyvecl_pointwise_poly_montgomery(&mut z, &cp, &s1);
        polyvecl_invntt_tomont(&mut z);
        polyvecl_add(&mut z, &y);
        polyvecl_reduce(&mut z);
        if polyvecl_chknorm::<L>(&z, (GAMMA1 - BETA) as i32) > 0 {
            continue;
        }

        /* Check that subtracting cs2 does not change high bits of w and low bits
         * do not reveal secret information */
        polyveck_pointwise_poly_montgomery(&mut h, &cp, &s2);
        polyveck_invntt_tomont(&mut h);
        polyveck_sub(&mut w0, &h);
        polyveck_reduce(&mut w0);
        if polyveck_chknorm::<K>(&w0, (GAMMA2 - BETA) as i32) > 0 {
            continue;
        }

        // Compute hints for w1
        polyveck_pointwise_poly_montgomery(&mut h, &cp, &t0);
        polyveck_invntt_tomont(&mut h);
        polyveck_reduce(&mut h);
        if polyveck_chknorm::<K>(&h, GAMMA2 as i32) > 0 {
            continue;
        }

        polyveck_add(&mut w0, &h);
        let n = polyveck_make_hint::<K, GAMMA2>(&mut h, &w0, &w1);
        if n > OMEGA as i32 {
            continue;
        }

        // Write signature
        pack_sig::<K, L, CTILDEBYTES, POLYZ_PACKEDBYTES, OMEGA, GAMMA1>(sig, None, &z, &h);
        return;
    }
}

pub fn crypto_sign_verify<const K: usize, const L: usize, const POLYW1_PACKEDBYTES: usize, const CTILDEBYTES: usize, const GAMMA1: usize, const GAMMA2: usize, const BETA: usize, const OMEGA: usize, const SIGNBYTES: usize, const PUBLICKEYBYTES: usize, const POLYZ_PACKEDBYTES: usize, const TAU: usize, const POLY_UNIFORM_NBLOCKS: usize, const STREAM128_BLOCKBYTES: usize>(
    sig: &[u8],
    m: &[u8],
    pk: &[u8],
    ctx: Option<&[u8]>,
) -> Result<(), SignError> {
    let mut buf = vec![0u8; K * POLYW1_PACKEDBYTES];
    let mut rho = [0u8; SEEDBYTES];
    let mut mu = [0u8; CRHBYTES];
    let mut c = [0u8; CTILDEBYTES];
    let mut c2 = [0u8; CTILDEBYTES];
    let mut cp = Poly::default();
    let (mut mat, mut z) = ([Polyvecl::<L>::default(); K], Polyvecl::<L>::default());
    let (mut t1, mut w1, mut h) = (
        Polyveck::<K>::default(),
        Polyveck::<K>::default(),
        Polyveck::<K>::default(),
    );
    let mut state = KeccakState::default(); // shake256_init()

    let ctx = match ctx {
        Some(x) => x,
        None => &[0u8; 0],
    };

    if ctx.len() > 255 {
        return Err(SignError::Input);
    }

    if sig.len() != SIGNBYTES {
        return Err(SignError::Input);
    }

    unpack_pk::<K>(&mut rho, &mut t1, pk);
    if let Err(e) = unpack_sig::<K, L, CTILDEBYTES, POLYZ_PACKEDBYTES, OMEGA, GAMMA1>(&mut c, &mut z, &mut h, sig) {
        return Err(e);
    }
    if polyvecl_chknorm(&z, (GAMMA1 - BETA) as i32) > 0 {
        return Err(SignError::Input);
    }

    // Compute CRH(CRH(rho, t1), msg)
    shake256(&mut mu, TRBYTES, pk, PUBLICKEYBYTES);
    shake256_absorb(&mut state, &mu, TRBYTES);
    let mm = [0u8, ctx.len() as u8];
    shake256_absorb(&mut state, &mm, 2);
    shake256_absorb(&mut state, ctx, ctx.len());
    shake256_absorb(&mut state, m, m.len());
    shake256_finalize(&mut state);
    shake256_squeeze(&mut mu, CRHBYTES, &mut state);

    // Matrix-vector multiplication; compute Az - c2^dt1
    poly_challenge::<CTILDEBYTES, TAU>(&mut cp, &c);
    polyvec_matrix_expand::<K, L, POLY_UNIFORM_NBLOCKS, STREAM128_BLOCKBYTES>(&mut mat, &rho);

    polyvecl_ntt(&mut z);
    polyvec_matrix_pointwise_montgomery(&mut w1, &mat, &z);

    poly_ntt(&mut cp);
    polyveck_shiftl(&mut t1);
    polyveck_ntt(&mut t1);
    let t1_2 = t1.clone();
    polyveck_pointwise_poly_montgomery(&mut t1, &cp, &t1_2);

    polyveck_sub(&mut w1, &t1);
    polyveck_reduce(&mut w1);
    polyveck_invntt_tomont(&mut w1);

    // Reconstruct w1
    polyveck_caddq(&mut w1);
    polyveck_use_hint::<K, GAMMA2>(&mut w1, &h);
    polyveck_pack_w1::<K, POLYW1_PACKEDBYTES, GAMMA2>(&mut buf, &w1);

    // Call random oracle and verify challenge
    state.init();
    shake256_absorb(&mut state, &mu, CRHBYTES);
    shake256_absorb(&mut state, &buf, K * POLYW1_PACKEDBYTES);
    shake256_finalize(&mut state);
    shake256_squeeze(&mut c2, CTILDEBYTES, &mut state);
    // Doesn't require constant time equality check
    if c != c2 {
        Err(SignError::Verify)
    } else {
        Ok(())
    }
}

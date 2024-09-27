use crate::{params::*, poly::*, polyvec::*};
use crate::api::SignError;

/// Bit-pack public key pk = (rho, t1).
pub fn pack_pk<const K: usize>(pk: &mut [u8], rho: &[u8], t1: &Polyveck<K>) {
    pk[..SEEDBYTES].copy_from_slice(&rho[..SEEDBYTES]);
    for i in 0..K {
        polyt1_pack(&mut pk[SEEDBYTES + i * POLYT1_PACKEDBYTES..], &t1.vec[i]);
    }
}

/// Unpack public key pk = (rho, t1).
pub fn unpack_pk<const K: usize>(rho: &mut [u8], t1: &mut Polyveck<K>, pk: &[u8]) {
    rho[..SEEDBYTES].copy_from_slice(&pk[..SEEDBYTES]);
    for i in 0..K {
        polyt1_unpack(&mut t1.vec[i], &pk[SEEDBYTES + i * POLYT1_PACKEDBYTES..])
    }
}

/// Bit-pack secret key sk = (rho, key, tr, s1, s2, t0).
pub fn pack_sk<const K: usize, const L: usize, const POLYETA_PACKEDBYTES: usize, const ETA: usize>(
    sk: &mut [u8],
    rho: &[u8],
    tr: &[u8],
    key: &[u8],
    t0: &Polyveck<K>,
    s1: &Polyvecl<L>,
    s2: &Polyveck<K>,
) {
    let mut idx = 0usize;

    sk[idx..SEEDBYTES].copy_from_slice(&rho[0..SEEDBYTES]);
    idx += SEEDBYTES;

    sk[idx..idx + SEEDBYTES].copy_from_slice(&key[0..SEEDBYTES]);
    idx += SEEDBYTES;

    sk[idx..idx + TRBYTES].copy_from_slice(&tr[0..TRBYTES]);
    idx += TRBYTES;

    for i in 0..L {
        polyeta_pack::<ETA>(&mut sk[idx + i * POLYETA_PACKEDBYTES..], &s1.vec[i]);
    }
    idx += L * POLYETA_PACKEDBYTES;

    for i in 0..K {
        polyeta_pack::<ETA>(&mut sk[idx + i * POLYETA_PACKEDBYTES..], &s2.vec[i]);
    }
    idx += K * POLYETA_PACKEDBYTES;

    for i in 0..K {
        polyt0_pack(&mut sk[idx + i * POLYT0_PACKEDBYTES..], &t0.vec[i]);
    }
}

/// Unpack secret key sk = (rho, key, tr, s1, s2, t0).
pub fn unpack_sk<const K: usize, const L: usize, const POLYETA_PACKEDBYTES: usize, const ETA: usize>(
    rho: &mut [u8],
    tr: &mut [u8],
    key: &mut [u8],
    t0: &mut Polyveck<K>,
    s1: &mut Polyvecl<L>,
    s2: &mut Polyveck<K>,
    sk: &[u8],
) {
    let mut idx = 0usize;

    rho[..SEEDBYTES].copy_from_slice(&sk[..SEEDBYTES]);
    idx += SEEDBYTES;

    key[..SEEDBYTES].copy_from_slice(&sk[idx..idx + SEEDBYTES]);
    idx += SEEDBYTES;

    tr[..TRBYTES].copy_from_slice(&sk[idx..idx + TRBYTES]);
    idx += TRBYTES;

    for i in 0..L {
        polyeta_unpack::<ETA>(&mut s1.vec[i], &sk[idx + i * POLYETA_PACKEDBYTES..]);
    }
    idx += L * POLYETA_PACKEDBYTES;

    for i in 0..K {
        polyeta_unpack::<ETA>(&mut s2.vec[i], &sk[idx + i * POLYETA_PACKEDBYTES..]);
    }
    idx += K * POLYETA_PACKEDBYTES;

    for i in 0..K {
        polyt0_unpack(&mut t0.vec[i], &sk[idx + i * POLYT0_PACKEDBYTES..]);
    }
}

/// Bit-pack signature sig = (c, z, h).
pub fn pack_sig<const K: usize, const L: usize, const CTILDEBYTES: usize, const POLYZ_PACKEDBYTES: usize, const OMEGA: usize, const GAMMA1: usize>
(sig: &mut [u8], c: Option<&[u8]>, z: &Polyvecl<L>, h: &Polyveck<K>) {
    let mut idx = 0usize;

    if let Some(challenge) = c {
        sig[..CTILDEBYTES].copy_from_slice(&challenge[..CTILDEBYTES]);
    }

    idx += CTILDEBYTES;

    for i in 0..L {
        polyz_pack::<GAMMA1>(&mut sig[idx + i * POLYZ_PACKEDBYTES..], &z.vec[i]);
    }
    idx += L * POLYZ_PACKEDBYTES;
    // Encode H
    sig[idx..idx + OMEGA + K].copy_from_slice(&vec![0u8; OMEGA + K]);

    let mut k = 0;
    for i in 0..K {
        for j in 0..N {
            if h.vec[i].coeffs[j] != 0 {
                sig[idx + k] = j as u8;
                k += 1;
            }
        }
        sig[idx + OMEGA + i] = k as u8;
    }
}

/// Unpack signature sig = (z, h, c).
pub fn unpack_sig<const K: usize, const L: usize, const CTILDEBYTES: usize, const POLYZ_PACKEDBYTES: usize, const OMEGA: usize, const GAMMA1: usize>(
    c: &mut [u8],
    z: &mut Polyvecl<L>,
    h: &mut Polyveck<K>,
    sig: &[u8],
) -> Result<(), SignError> {
    let mut idx = 0usize;

    c[..CTILDEBYTES].copy_from_slice(&sig[..CTILDEBYTES]);
    idx += CTILDEBYTES;

    for i in 0..L {
        polyz_unpack::<GAMMA1>(&mut z.vec[i], &sig[idx + i * POLYZ_PACKEDBYTES..]);
    }
    idx += L * POLYZ_PACKEDBYTES;

    // Decode h
    let mut k = 0usize;
    for i in 0..K {
        if sig[idx + OMEGA + i] < k as u8 || sig[idx + OMEGA + i] > (OMEGA as u8) {
            return Err(SignError::Input);
        }
        for j in k..sig[idx + OMEGA + i] as usize {
            // Coefficients are ordered for strong unforgeability
            if j > k && sig[idx + j as usize] <= sig[idx + j as usize - 1] {
                return Err(SignError::Input);
            }
            h.vec[i].coeffs[sig[idx + j] as usize] = 1;
        }
        k = sig[idx + OMEGA + i] as usize;
    }

    // Extra indices are zero for strong unforgeability
    for j in k..OMEGA {
        if sig[idx + j as usize] > 0 {
            return Err(SignError::Input);
        }
    }

    Ok(())
}

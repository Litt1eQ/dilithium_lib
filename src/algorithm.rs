use crate::api::SignError;
use crate::sign::{crypto_sign_keypair, crypto_sign_signature, crypto_sign_verify};

pub struct Dilithium<
    const K: usize,
    const L: usize,
    const ETA: usize,
    const TAU: usize,
    const BETA: usize,
    const GAMMA1: usize,
    const GAMMA2: usize,
    const OMEGA: usize,
    const PUBLICKEYBYTES: usize,
    const SECRETKEYBYTES: usize,
    const POLY_UNIFORM_NBLOCKS: usize,
    const STREAM128_BLOCKBYTES: usize,
    const POLY_UNIFORM_ETA_NBLOCKS: usize,
    const POLYETA_PACKEDBYTES: usize,
    const SIGNBYTES: usize,
    const POLYW1_PACKEDBYTES: usize,
    const CTILDEBYTES: usize,
    const POLYZ_PACKEDBYTES: usize,
    const POLYVECH_PACKEDBYTES: usize,
    const POLY_UNIFORM_GAMMA1_NBLOCKS: usize,
>;

impl<
    const K: usize,
    const L: usize,
    const ETA: usize,
    const TAU: usize,
    const BETA: usize,
    const GAMMA1: usize,
    const GAMMA2: usize,
    const OMEGA: usize,
    const PUBLICKEYBYTES: usize,
    const SECRETKEYBYTES: usize,
    const POLY_UNIFORM_NBLOCKS: usize,
    const STREAM128_BLOCKBYTES: usize,
    const POLY_UNIFORM_ETA_NBLOCKS: usize,
    const POLYETA_PACKEDBYTES: usize,
    const SIGNBYTES: usize,
    const POLYW1_PACKEDBYTES: usize,
    const CTILDEBYTES: usize,
    const POLYZ_PACKEDBYTES: usize,
    const POLYVECH_PACKEDBYTES: usize,
    const POLY_UNIFORM_GAMMA1_NBLOCKS: usize,
>
Dilithium<
    K,
    L,
    ETA,
    TAU,
    BETA,
    GAMMA1,
    GAMMA2,
    OMEGA,
    PUBLICKEYBYTES,
    SECRETKEYBYTES,
    POLY_UNIFORM_NBLOCKS,
    STREAM128_BLOCKBYTES,
    POLY_UNIFORM_ETA_NBLOCKS,
    POLYETA_PACKEDBYTES,
    SIGNBYTES,
    POLYW1_PACKEDBYTES,
    CTILDEBYTES,
    POLYZ_PACKEDBYTES,
    POLYVECH_PACKEDBYTES,
    POLY_UNIFORM_GAMMA1_NBLOCKS,
>
{
    pub fn key_gen(zeta: Option<Vec<u8>>) -> anyhow::Result<(Vec<u8>, Vec<u8>)> {
        let mut pk = [0u8; PUBLICKEYBYTES];
        let mut sk = [0u8; SECRETKEYBYTES];
        match zeta {
            None => {
                crypto_sign_keypair::<K, L, PUBLICKEYBYTES, POLY_UNIFORM_NBLOCKS, STREAM128_BLOCKBYTES, POLY_UNIFORM_ETA_NBLOCKS, ETA, POLYETA_PACKEDBYTES>(&mut pk, &mut sk, None);
            }
            Some(zeta) => {
                crypto_sign_keypair::<K, L, PUBLICKEYBYTES, POLY_UNIFORM_NBLOCKS, STREAM128_BLOCKBYTES, POLY_UNIFORM_ETA_NBLOCKS, ETA, POLYETA_PACKEDBYTES>(&mut pk, &mut sk, Some(&zeta));
            }
        }
        Ok((sk.to_vec(), pk.to_vec()))
    }

    pub fn sign(data: Vec<u8>, sk: Vec<u8>, ctx: Option<Vec<u8>>, using_randomized_signing: bool) -> anyhow::Result<Vec<u8>> {
        let mut sig = [0u8; SIGNBYTES];
        match ctx {
            None => {
                crypto_sign_signature::<K, L, POLYW1_PACKEDBYTES, CTILDEBYTES, GAMMA1, GAMMA2, BETA, OMEGA, POLYETA_PACKEDBYTES, POLY_UNIFORM_NBLOCKS, STREAM128_BLOCKBYTES, POLY_UNIFORM_GAMMA1_NBLOCKS, TAU, POLYZ_PACKEDBYTES, ETA>(&mut sig, &data, &sk, None, using_randomized_signing);
            }
            Some(ctx) => {
                crypto_sign_signature::<K, L, POLYW1_PACKEDBYTES, CTILDEBYTES, GAMMA1, GAMMA2, BETA, OMEGA, POLYETA_PACKEDBYTES, POLY_UNIFORM_NBLOCKS, STREAM128_BLOCKBYTES, POLY_UNIFORM_GAMMA1_NBLOCKS, TAU, POLYZ_PACKEDBYTES, ETA>(&mut sig, &data, &sk, Some(&ctx), using_randomized_signing);
            }
        }
        Ok(sig.to_vec())
    }

    pub fn verify(data: Vec<u8>, sig: Vec<u8>, pk: Vec<u8>, zeta: Option<Vec<u8>>) -> Result<(), SignError> {
        match zeta {
            None => {
                crypto_sign_verify::<K, L, POLYW1_PACKEDBYTES, CTILDEBYTES, GAMMA1, GAMMA2, BETA, OMEGA, SIGNBYTES, PUBLICKEYBYTES, POLYZ_PACKEDBYTES, TAU, POLY_UNIFORM_NBLOCKS, STREAM128_BLOCKBYTES>(&sig, &data, &pk, None)
            }
            Some(zeta) => {
                crypto_sign_verify::<K, L, POLYW1_PACKEDBYTES, CTILDEBYTES, GAMMA1, GAMMA2, BETA, OMEGA, SIGNBYTES, PUBLICKEYBYTES, POLYZ_PACKEDBYTES, TAU, POLY_UNIFORM_NBLOCKS, STREAM128_BLOCKBYTES>(&sig, &data, &pk, Some(&zeta))
            }
        }
    }
}

pub type Dilithium2 = Dilithium::<
    4, 4, 2, 39, 78, 0x20000, 95232, 80, 1312, 2560, 5, 168, 1, 96, 2420, 192, 32, 576, 84, 5
>;

pub type Dilithium3 = Dilithium::<
    6, 5, 4, 49, 196, 0x80000, 261888, 55, 1952, 4032, 5, 168, 2, 128, 3309, 128, 48, 640, 71, 5
>;

pub type Dilithium5 = Dilithium::<
    8, 7, 2, 60, 120, 0x80000, 261888, 75, 2592, 4896, 5, 168, 1, 96, 4627, 128, 64, 640, 83, 5
>;

#[test]
fn test_dilithium_mode2() {
    let seed = vec![0u8; 32];
    let (sk, pk) = Dilithium2::key_gen(Some(seed)).unwrap();
    let data = vec![49u8; 10];
    let sign = Dilithium2::sign(data, sk, None, false).unwrap();
    let data = vec![49u8; 10];
    let x = Dilithium2::verify(data, sign, pk, None);
    assert!(x.is_ok());
}


#[test]
fn test_dilithium_mode3() {
    type D3 = Dilithium::<
        6, 5, 4, 49, 196, 0x80000, 261888, 55, 1952, 4032, 5, 168, 2, 128, 3309, 128, 48, 640, 71, 5
    >;
    let (sk, pk) = D3::key_gen(None).unwrap();
    let data = vec![31u8; 10];
    let sign = D3::sign(data, sk, None, false).unwrap();
    let data = vec![31u8; 10];
    let x = D3::verify(data, sign, pk, None);
    assert!(x.is_ok());
}

#[test]
fn test_dilithium_mode5() {
    type D5 = Dilithium::<
        8, 7, 2, 60, 120, 0x80000, 261888, 75, 2592, 4896, 5, 168, 1, 96, 4627, 128, 64, 640, 83, 5
    >;
    let (sk, pk) = D5::key_gen(None).unwrap();
    let data = vec![31u8; 10];
    let sign = D5::sign(data, sk, None, false).unwrap();
    let data = vec![31u8; 10];
    let x = D5::verify(data, sign, pk, None);
    assert!(x.is_ok());
}
mod distributions;
mod fft;
mod gbfv;
mod lwe_old;
mod montgomery;
mod ntt;
mod polynomial;
mod rlwe;

use byteorder::{LittleEndian, WriteBytesExt};
use fft::m_ladder;
use lwe_old::delta_polyn;
//use lwe_old::delta_polyn;
//use lwe_new::mult_ct;
use montgomery::U1024;
use ntt::Ntt;
use num_bigint::{BigInt, BigUint};
use num_traits::one;
use polynomial::{single_term, Parameters, Polynomial, Sampling};
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::fmt::Display;
use std::fmt::Write;
use std::fs;
use std::rc::Rc;
use std::time::Instant;
use std::usize;
use uint::construct_uint;

use crate::rlwe::RLWE;

fn decrypt_vs_3decrypt(t: i64) {
    let parameters: Parameters<i64> = Parameters {
        degree: 4,
        q: 768100,
        t: vec![t],
        t_relin: 32,
        p: 4,
        log_range: 768100_i64.ilog(32) as usize + 1, //println!("{}", (self.q().checked_ilog(self.t_relin()).unwrap()) + 1);,
    };
    let _p = &Rc::new(RefCell::new(parameters));
    let m1 = Polynomial::uniform_sample(_p).mod_t();
    let m2 = Polynomial::uniform_sample(_p).mod_t();

    //let m1 = Polynomial::new(vec![1, 2, 3, 4], _p);
    //let m2 = Polynomial::new(vec![5, 6, 7, 8], _p);
    let mut secret = RLWE::new(_p);
    let relin = secret.evaluate_key_gen();

    let c1 = RLWE::encrypt(&secret.public, &m1);
    let c2 = RLWE::encrypt(&secret.public, &m2);

    let c = RLWE::mult_ct(&c1, &c2);

    let p = secret.decrypt_degree_3(&c);
    let c_relin = RLWE::relin_ct(c, relin);
    let p_relin = secret.decrypt(&c_relin);
    println!("{}", t);
    println!("combined message: {}", (&m1 * &m2).get_mod().mod_t());
    //println!("??? {}", (&m1 * &m2).get_mod().mod_t());
    //println!("??? {}", secret.decrypt(&c_relin));
    println!("Test: [{}]", p);
    println!("Relin: [{}]", p_relin);
}
fn example_decrypt() {
    let mut parameters = polynomial::Parameters {
        degree: 4,
        q: 100,
        t: vec![20],
        t_relin: 32,
        p: 4,
        log_range: 0,
    };
    let _p = &Rc::new(RefCell::new(parameters));

    let mut secret = RLWE::_new(vec![1, 0, 1, 0], _p);

    let s = Polynomial::new(vec![1, 0, 1, 0], _p);
    let c0 = Polynomial::new(vec![8, 84, 86, 98], _p);
    let c1 = Polynomial::new(vec![14, 74, 16, 49], _p);

    let res_1 = (&c1 * &s).get_mod().mod_q();

    println!("res_1, {}", res_1);

    let res_2 = (&c0 + &res_1).mod_q();

    println!("res_2, {}", res_2);
    let vals = RLWE::delta_polyn(res_2);
    println!("vals, {}", vals.mod_t());
    println!("MAYBE, {}", secret.decrypt(&(c0, c1)));
}
fn example_encrypt() {
    let mut parameters = polynomial::Parameters {
        degree: 4,
        q: 100,
        t: vec![20],
        t_relin: 32,
        p: 4,
        log_range: 0,
    };
    let _p = &Rc::new(RefCell::new(parameters));

    let pk_0 = Polynomial::new(vec![74, 71, 77, 98], _p);
    let u = Polynomial::new(vec![0, 1, 0, 0], _p);

    let res_1 = &pk_0 * &u;
    println!("result {}", res_1);
    println!("result with mod {}", res_1.get_mod());
    let res_1 = res_1.get_mod().mod_q();
    println!("result with mod {}", res_1);
    let e_1 = Polynomial::new(vec![1, 0, 0, 1], _p);
    let res_2 = &res_1 + &e_1;
    println!("second result {}", res_2);
    let res_2 = res_2.get_mod().mod_q();
    println!("second result with mod {}", res_2);
    let scaled_message = Polynomial::new(vec![5, 10, 15, 20], _p);
    let res_3 = &res_2 + &scaled_message;
    println!("third result {}", res_3);
    let res_3 = res_3.get_mod().mod_q();
    println!("third result with mod {}", res_3);

    println!("PK_1:");
    let pk_1 = Polynomial::new(vec![74, 15, 49, 86], _p);
    let res_1 = &pk_1 * &u;
    println!("result {}", res_1);
    println!("result with mod {}", res_1.get_mod().mod_q());
}
fn example_pk() {
    let mut parameters = polynomial::Parameters {
        degree: 4,
        q: 100,
        t: vec![20],
        t_relin: 32,
        p: 4,
        log_range: 0,
    };
    let _p = &Rc::new(RefCell::new(parameters));
    let a = Polynomial::new(vec![74, 15, 49, 86], _p);
    let e = Polynomial::new(vec![1, 0, 0, 1], _p);
    let secret = Polynomial::new(vec![1, 0, 1, 0], _p);

    let pk = (&Polynomial::new(vec![], _p) - &(&(&a * &secret) + &e), a);
    println!(
        "PK: {} : {}",
        pk.0.get_mod().mod_q(),
        pk.1.get_mod().mod_q()
    );

    let pk0 = pk.0; //.get_mod().mod_q();
    let pk1 = pk.1; //.get_mod().mod_q();
    let ct = RLWE::compute_encryption(
        &(pk0, pk1),
        &Polynomial::new(vec![1, 2, 3, 4], _p),
        Polynomial::new(vec![0, 1, 0, 0], _p),
        Polynomial::new(vec![1, 0, 0, 1], _p),
        Polynomial::new(vec![0, 0, 1, 0], _p),
    );
    println!(
        "CT: {} : {}",
        ct.0.get_mod().mod_q(),
        ct.1.get_mod().mod_q()
    );

    let mut secret = RLWE::_new(vec![1, 0, 1, 0], _p);

    let s = Polynomial::new(vec![1, 0, 1, 0], _p);
    let c0 = Polynomial::new(vec![8, 84, 86, 98], _p);
    let c1 = Polynomial::new(vec![14, 74, 16, 49], _p);

    let res_1 = (&c1 * &s);

    println!("res_1, {}", res_1);

    let res_2 = &c0 + &res_1;

    println!("res_2, {}", res_2);
    let vals = RLWE::delta_polyn(res_2);
    println!("vals, {}", vals.mod_t());
    println!("MAYBE, {}", secret.decrypt(&(c0, c1)));

    return;
    let x = Polynomial::new(vec![1, 0, 1, 0], _p);
    let y = Polynomial::new(vec![74, 15, 49, 86], _p);
    let res = &x * &y;
    println!("{}", res);
    println!("{}", res.get_mod());
    let testx = Polynomial::new(vec![74, 15, 49, 86], _p);
    let testy = Polynomial::new(vec![0, 0, 74, 15, 49, 86], _p);
    println!("{}", &testx + &testy);

    let a_s = Polynomial::new(vec![25, 29, 23, 1], _p);
    let e = Polynomial::new(vec![1, 0, 0, 1], _p);
    let _part1 = &a_s + &e;
    println!("{}", _part1);
    let part1 = &Polynomial::new(vec![], _p) - &_part1;
    println!("{}", part1.mod_q());
    println!("{}", Polynomial::new(vec![-26, -29, -23, -2], _p).mod_q());
}
fn p_q_gradient() {
    let mut parameters = polynomial::Parameters {
        degree: 16,
        q: 7681,
        t: vec![6],
        t_relin: 32,
        p: 4,
        log_range: 0,
    };
    let mut buf = String::new();
    for i in 1..(parameters.q) {
        parameters.t = vec![i];
        let _p = &Rc::new(RefCell::new(parameters.clone()));
        let m = Polynomial::new(vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], _p);

        let mut secret = RLWE::new(_p);
        let rlk = secret.evaluate_key_gen();

        let c = RLWE::encrypt(&secret.public, &m);
        let p = secret.decrypt(&c);
        write!(&mut buf, "{}\n", p);
    }
    fs::write(
        "/Users/rensweerman/Documents/PythonScripts/Q_vs_T/vectors.txt",
        buf,
    )
    .expect("Unable to write file");
}
fn timing_test() {
    let parameters: Parameters<BigInt> = Parameters {
        degree: 32,
        q: BigInt::try_from(1329227995784915872903807060280344576_i128).unwrap(),
        t: vec![BigInt::try_from(536903681).unwrap()],
        t_relin: BigInt::try_from(16).unwrap(),
        p: BigInt::try_from(4).unwrap(),
        log_range: 0,
    };
    let _p = &Rc::new(RefCell::new(parameters));

    let p1 = Polynomial::uniform_sample(_p).mod_t(); //Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], _p);
    let p2 = Polynomial::uniform_sample(_p).mod_t(); //Polynomial::new(vec![3, 2, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], _p);
    let p3 = Polynomial::uniform_sample(_p).mod_t(); //Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], _p);

    let mut secret = RLWE::new(_p);
    let rlk = secret.evaluate_key_gen();

    let c1 = RLWE::encrypt(&secret.public, &p1);
    let c2 = RLWE::encrypt(&secret.public, &p2);
    let c3 = RLWE::encrypt(&secret.public, &p3);

    let start = Instant::now();

    let add_ct = RLWE::add_ct(&RLWE::add_ct(&c1, &c2), &c3);
    //let add_ct = c1.clone();

    let mult_ct_1_2 = RLWE::relin_ct(RLWE::mult_ct(&c1, &c2), rlk.clone());
    let mult_ct = RLWE::relin_ct(RLWE::mult_ct(&mult_ct_1_2, &c3), rlk);

    let add_res = secret.decrypt(&add_ct);
    let mult_res = secret.decrypt(&mult_ct);
    let duration = start.elapsed().as_micros();
    println!("Time taken for Timing Test: {:?}", duration);
    print_vec(&(&p1 * &(&p2 * &p3)).get_mod().val);
    //print_vec(&p1.val);
    print_vec(&mult_res.val);
}

fn test_iter_cooley_tukey_2() {
    let parameters = polynomial::Parameters {
        degree: 4,
        q: 7681,
        t: vec![6],
        t_relin: 32,
        p: 4,
        log_range: 0,
    };
    let _p = &Rc::new(RefCell::new(parameters));

    let psi = &Ntt::second_primitive_root(&7681, 8);
    let q = &7681;
    let ls = Ntt::iter_cooley_tukey(
        Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7, 8], _p),
        psi,
        &7681,
    );
    println!("Cooley-Tukey result: ");
    print_vec(&ls.val);
    let test_val = Ntt::second_ntt(Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7, 8], _p), psi, q);
    println!("naive result: ");
    print_vec(&test_val.val);
    let ls = Ntt::second_inv_ntt(ls, psi, q);
    assert_eq!(ls, Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7, 8], _p));
}
fn test_multiplying_speeds() {
    let parameters = polynomial::Parameters {
        degree: 256,
        q: 7681,
        t: vec![6],
        t_relin: 32,
        p: 4,
        log_range: 0,
    };
    let _p = &Rc::new(RefCell::new(parameters));
    let q = &_p.borrow().q;

    let psi = &Ntt::second_primitive_root(q, _p.borrow().degree as i64);

    let a = Polynomial::normal_sample(_p);
    let b = Polynomial::normal_sample(_p);
    print_vec(&a.val);
    print_vec(&b.val);
    let a1 = a.clone();
    let b1 = b.clone();
    let start = Instant::now();
    let ls = Ntt::ntt_mult(a1, b1, q, psi);
    let duration = start.elapsed();
    println!("Time taken ntt: {:?}", duration);

    let start = Instant::now();
    let rs = (&a * &b).get_mod();
    let duration = start.elapsed();
    println!("Time taken simple: {:?}", duration);
    //print_vec(&ls.val);
    //print_vec(&rs.val);
    assert_eq!(ls, rs);
    println!("Worked");
}
fn testing_big_ints() {
    let parameters: Parameters<BigInt> = Parameters {
        degree: 16,
        q: BigInt::try_from(7681).unwrap(),
        t: vec![BigInt::try_from(6).unwrap()],
        t_relin: BigInt::try_from(32).unwrap(),
        p: BigInt::try_from(4).unwrap(),
        log_range: 0,
    };
    let _p = &Rc::new(RefCell::new(parameters));
    let q = &_p.borrow().q.clone();

    let psi = &Ntt::second_primitive_root(q, _p.borrow().degree as i64);

    let a = Polynomial::uniform_sample(_p);
    let b = Polynomial::uniform_sample(_p);
    print_vec(&a.val);
    print_vec(&b.val);
    let a1 = a.clone();
    let b1 = b.clone();
    let start = Instant::now();
    let ls = Ntt::ntt_mult(a1, b1, q, psi).get_mod().mod_q();
    let duration = start.elapsed();
    println!("Time taken ntt: {:?}", duration);

    let start = Instant::now();
    let rs = (&a * &b).get_mod().mod_q();
    //let rs = Polynomial::new(vec![one()], _p);
    let duration = start.elapsed();
    println!("Time taken simple: {:?}", duration);
    //print_vec(&ls.val);
    //print_vec(&rs.val);
    assert_eq!(ls, rs);
    println!("Worked");
}
fn main() {
    for i in 2..300 {
        if (i % 40 == 0) {
            decrypt_vs_3decrypt(i);
        }
    }
    println!("ENCRYPT BELOW------");
    example_encrypt();
    println!("DECRYPT BELOW------");
    example_decrypt();
    println!("PK BELOW------");
    p_q_gradient();
    example_pk();
    return;
    //p_q_gradient();
    //timing_test();
    //test_multiplying_speeds();
    //testing_big_ints();
    let parameters = polynomial::Parameters {
        degree: 4,
        q: 7681,
        t: vec![6],
        t_relin: 32,
        p: 4,
        log_range: 0,
    };
    let _p = &Rc::new(RefCell::new(parameters));
    let secret = Polynomial::new(vec![0, 1, 1, 0], _p);
    let a: Polynomial<i64> = Polynomial::new(vec![1280, -1, -1945, -2597, -2315, -1851, -1468], _p);
    let b = Polynomial::new(vec![0, 1946, 1932, 383, 1468], _p);
    let c = Polynomial::new(vec![-1, -1945, -2597, -2314, -1851, -1468], _p);
    let d = Polynomial::new(vec![-1, 1946, 1932, 383, 1468], _p);

    let ct0 = &a * &c;
    print_vec(&a.val);
    print_vec(&c.val);
    println!("Old Mult, then new");
    print_vec(&a.old_mul(&c).get_mod().val);
    print_vec(&ct0.get_mod().val);

    let ct1 = &(&a * &d) + &(&b * &c);
    print_vec(&(ct1).val);
    let ct2 = &b * &d;
    print_vec(&(ct2).val);

    println!("Testing delta_polyn");
    let ct0 = RLWE::delta_polyn(ct0.clone()).mod_q();
    let ct1 = RLWE::delta_polyn(ct1.clone()).mod_q();
    let ct2 = RLWE::delta_polyn(ct2.clone()).mod_q();
    print_vec(&ct0.val);
    print_vec(&ct1.val);
    print_vec(&ct2.val);

    // let mut first_part = &ct0 * &single_term(_p, 6, 0);
    // print!("EEE: ");
    // print_vec(&ct0.val);
    // print_vec(&first_part.val);
    // first_part.val = first_part.val.into_iter().map(|i| (i / 7681)).collect();
    // print_vec(&(&(&a * &single_term(_p, 6, 0)) / &single_term(_p, 7681, 0)).val);

    // let c0 = &ct1.0 * &ct2.0;
    // let c1 = &(&ct1.0 * &ct2.1) + &(&ct1.1 * &ct2.0);
    // let c2 = &ct1.1 * &ct2.1;
    //KEY_GEN:

    let a1 = Polynomial::new(vec![0, 1, 0, 0], _p);
    let e1 = Polynomial::new(vec![0, 0, 1, 0], _p);

    let a2 = Polynomial::new(vec![1, 0, 0, 0], _p);
    let e2 = Polynomial::new(vec![0, 0, 1, 0], _p);

    let a3 = Polynomial::new(vec![1, 0, 0, 1], _p);
    let e3 = Polynomial::new(vec![0, 1, 0, 0], _p);

    let squared_key = &secret * &secret;
    let empty = Polynomial::new(vec![], _p);
    let squared_1 = &Polynomial::new(vec![32_i64.pow(1)], _p) * &squared_key;
    let squared_2 = &Polynomial::new(vec![32_i64.pow(2)], _p) * &squared_key;
    let squared_3 = &Polynomial::new(vec![32_i64.pow(3)], _p) * &squared_key;

    let rlk1 = ((&(&empty - &(&(&a1 * &secret) + &e1))) + &squared_1, a1);
    let rlk2 = (
        &(&empty - &(&(&a2 * &secret) + &e2)) + &squared_2,
        a2.clone(),
    );
    let rlk3 = (&(&empty - &(&(&a3 * &secret) + &e3)) + &squared_3, a3);

    println!("rlk1:");
    print_vec(&rlk1.0.val);
    print_vec(&rlk1.1.val);
    println!("rlk2:");
    print_vec(&rlk2.0.val);
    print_vec(&rlk2.1.val);
    println!("EFFE");
    print_vec(&(&empty - &(&a2 * &secret)).val);
    print_vec(&(&empty - &(&(&a2 * &secret) + &e2)).val);
    print_vec(&(squared_key).val);
    print_vec(&(squared_2).val);
    println!("....");
    println!("rlk3:");
    print_vec(&rlk3.0.val);
    print_vec(&rlk3.1.val);

    //RELIN:
    println!("RELIN TEST:");
    print_vec(&ct2.val);
    let c_2r = &ct2 % &Polynomial::new(vec![rlk1.0.t_relin()], rlk1.0.parameters());
    print_vec(&c_2r.val);

    let ct0 = &ct0 + &(&c_2r * &rlk1.0);
    let ct1 = &ct1 + &(&c_2r * &rlk1.1);

    println!("Completed relin:");
    print_vec(&ct0.val);
    print_vec(&ct1.val);

    //DECRYPT:
    println!("Both before:");
    print_vec(&ct0.val);
    print_vec(&(&ct1 * &secret).val);
    let pt0 = (&ct0 + &(&ct1 * &secret));
    print_vec(&pt0.val);
    let pt0 = pt0.get_mod();
    print_vec(&pt0.val);
    let pt1 = RLWE::delta_polyn(pt0);
    print_vec(&pt1.val);
    println!("Result:");
    print_vec(&pt1.mod_t().val);
    let real_m = &Polynomial::new(vec![1, 0, 0, 1], _p) * &Polynomial::new(vec![0, 0, 1, 1], _p);
    print_vec(&real_m.get_mod().mod_q().val);

    /*
    let a = /*8Polynomial::new(vec![], self.parameters()); */Polynomial::normal_sample(self.parameters());
    let e = /*Polynomial::new(vec![], self.parameters()); */Polynomial::normal_sample(self.parameters());
    let mut k_0 = Polynomial::new(vec![], self.parameters());
    let k_1 = &(&a * &self.secret) + &e;
    k_0 = &k_0 - &k_1;
    k_0 = &k_0
        + &(&secret_2
            * &Polynomial::new(
                vec![square_multiply(self.t_relin(), i, self.q())],
                self.parameters(),
            ));
    k_0.normal_mod(self.q());
    k_0 = &k_0 % &make_fx(self.parameters());
    rlk.push((k_0, a)); */
    return;
    /*loop {
    match range.next() {
        Some(x) => {
            println!("{}", x);
        }
        None => break,
    }
    }*/
    /*let mut count = 0;
    distributions::Poly::new(2, 2).for_each(|x| {
        count += 1;
        print_vec(&x)
    });
    println!("{}", count);*/

    // let example1 = Polynomial::normal_sample();
    // let example2 = Polynomial::uniform_sample();
    // print_vec(&example1.val);
    // print_vec(&example2.val);
    // print_vec(&(&example1 + &example2).val);
    // let mut med = &example1 + &example2;
    // let example1 = &example1 + &example2;
    // print!("Polynomial: ");
    // print_vec(&example1.val);
    // print!("Normal Mult: ");
    // let start = Instant::now();
    // for _ in 0..2 - 1 {
    //     med = &med * &example1;
    //     println!("EAFA");
    //     med = med.get_mod();
    //     println!("eadaed");
    // }
    // let duration = start.elapsed();
    // print_vec(&med.val);
    // println!("Time taken: {:?}", duration);
    // let start = Instant::now();
    // //let res = polynomial::m_ladder(example1, 2);
    // let duration = start.elapsed();
    // //let _fx = polynomial::make_fx(res.fx.clone());
    // //res = res / _fx;
    // print!("Result:      ");
    // //print_vec(&res.val);
    // println!("Time taken: {:?}", duration);

    // let message = Polynomial::new(vec![1, 1, 0, 5, 1, 18, 41, 99, 1, 1]);
    // //let mut keys = lwe_new::LWE::new();
    // let mut keys = lwe_new::LWE::new_given_secret(Polynomial::new(vec![1, 0, 1]));
    // print!("BE Message: ");
    // print_vec(&message.get_mod().val);
    // let cipher = lwe_new::encrypt(&keys.public, &message);
    // //print!("Ciphertext: ");
    // //print_vec(&cipher.val);
    // let decrypted_message = keys.decrypt(&cipher);
    // print!("AE Message: ");
    // print_vec(&decrypted_message.val);
    // //return;
    // let message1 = Polynomial::new(vec![1, 0, 1]);
    // let message2 = Polynomial::new(vec![1, 0, 1]);
    // print_vec(&message1.val);
    // print_vec(&message2.val);
    // print!("Multiplied Message: ");
    // print_vec(&(&message1 * &message2).get_mod().val);
    // let ct1 = lwe_new::encrypt(&keys.public, &message1);
    // let ct2 = lwe_new::encrypt(&keys.public, &message2);
    // //let ct3 = lwe_new::add_ct(&ct1, &ct2);
    // //let plaintext = keys.decrypt(&ct3);
    // //print!("Decrypted Plaintext: ");
    // //print_vec(&plaintext.val); //combined_ct, keys.evaluate_key_gen(10), 10);
    // let combined_ct = mult_ct(&ct1, &ct2);
    // let relin_ct = lwe_new::relinearization_1(combined_ct, keys.relin_key_gen_1());
    // let combined_ct = mult_ct(&ct1, &ct2);
    // let relin_ct_2 = lwe_new::relin_ct_2(combined_ct, keys.evaluate_key_gen_2());
    // print!("RELIN CT.0: ");
    // print_vec(&relin_ct.0.val);
    // print!("RELIN CT.1: ");
    // print_vec(&relin_ct.1.val);
    // print!("2RELIN CT.0: ");
    // print_vec(&relin_ct_2.0.val);
    // print!("2RELIN CT.1: ");
    // print_vec(&relin_ct_2.1.val);
    // //let to_add = lwe_new::encrypt(&keys.public, &message1, size);
    // //relin_ct = (&relin_ct.0 + &to_add.0, &relin_ct.1 + &to_add.1);
    // let plaintext = keys.decrypt(&relin_ct);
    // let plaintext_2 = keys.decrypt(&relin_ct_2);
    // print!("Decrypted Plaintext: ");
    // print_vec(&plaintext.val); //combined_ct, keys.evaluate_key_gen(10), 10);
    // print_vec(&plaintext_2.val);
    // let plaintext = keys.decrypt(&ct1);
    // let plaintext_2 = keys.decrypt(&ct2);
    // print_vec(&plaintext.val);
    // return;
    // print_vec(&plaintext_2.val);
    // println!("MULT TEST: ");
    // let mult_test = mult_ct(
    //     &(Polynomial::new(vec![400, 1]), Polynomial::new(vec![0, 1])),
    //     &(Polynomial::new(vec![0, 1]), Polynomial::new(vec![0, 200])),
    // );
    // println!("{}", 3 / 991);
    // //ct1.0 * ct2.0
    // print_vec(&mult_test.0.val);
    // //(ct1.0 * ct2.1) + (ct1.1 * ct2.0)
    // //x^2 + x^2
    // //2x^2
    // print_vec(&mult_test.1.val);
    // //ct1.1 * ct2.1
    // print_vec(&mult_test.2.val);

    // println!("DECOMPOSING TEST");
    // let test_value = lwe_new::encrypt(&keys.public, &message1);
    // let _test_value = test_value.clone();
    // print_vec(&test_value.1.val);
    // let value_1 = test_value.1.decompose();
    // let mut testing = Polynomial::new(vec![]);
    // for i in 0..value_1.len() {
    //     let part = value_1.get(i).unwrap();
    //     print_vec(&part.val);
    //     testing = &testing + &(part * &Polynomial::new(vec![T_RELIN.pow(i as u32)]));
    // }
    // print_vec(&testing.val);
    // println!("-------");
    // _test_value.1.other_decompose();
    // return;

    // println!("{}", m_ladder(2, 11));
    // println!("{}", m_ladder(6, 9));

    // distributions::parallel_test(
    //     2,
    //     2,
    //     vec![1, 0, 1],
    //     false,
    //     "/Users/rensweerman/PycharmProjects/pythonProject/valuesp.txt",
    //     "/Users/rensweerman/PycharmProjects/pythonProject/labelsp.txt",
    // );
    // distributions::test_no_e();
    // distributions::parallel_test(
    //     2,
    //     2,
    //     vec![1, 0, 1],
    //     false,
    //     "/Users/rensweerman/PycharmProjects/pythonProject/valuesnoe.txt",
    //     "/Users/rensweerman/PycharmProjects/pythonProject/labelsnoe.txt",
    // );
    //distributions::test3();
    return;
    // for x in ffff {
    //     print_vec(&x);
    //     //print!("!");
    // }
    // return;
    // let q = 1000;
    // let degree = 8;

    // let mut fx = Vec::new();
    // fx.push(1);
    // for _ in 0..(degree - 1) {
    //     fx.push(0);
    // }
    // fx.push(1);
    // let _p = &Rc::new(RefCell::new(PARAMETERS));

    // let mut a: RLWE<i64> = RLWE::new(_p);
    // let fx = vec![1, 0, 0, 0, 0, 0, 0, 0, 1];

    // let _ = String::from("Hello World!") //Dit is een faketekst. Alles wat hier staat is slechts om een indruk te geven van het grafische effect van tekst op deze plek. Wat u hier leest is een voorbeeldtekst. Deze wordt later vervangen door de uiteindelijke tekst, die nu nog niet bekend is. De faketekst is dus een tekst die eigenlijk nergens over gaat. Het grappige is, dat mensen deze toch vaak lezen. Zelfs als men weet dat het om een faketekst gaat, lezen ze toch door.")
    //     .into_bytes()
    //     .into_iter()
    //     .map(|x| into_bits(x))
    //     .collect::<Vec<_>>();
    // //print!("TEXT: ");
    // //print_vec(&m);
    // let pk = &a.public.clone();
    // let secret = vec![-1, 0, 0, 0, 0, 0, 0, 1, 0];
    // let public = (
    //     vec![
    //         863, 966, 593, 826, 374, 791, 450, 373, 964, 407, 175, 625, 209, 549, 765, 70,
    //     ],
    //     vec![863, 966, 593, 825, 375, 791, 451, 235, 930],
    // );
    // let size = 20;
    // let encoded: Vec<u8> = bincode::serialize(&secret).unwrap();
    // print!("Secret key: ");
    // //print_vec(&encoded);
    // let s = base64::encode(encoded);
    // println!("{}", s);

    // let encoded: Vec<u8> = bincode::serialize(&pk).unwrap();
    // print!("Public key: ");
    // //print_vec(&encoded);
    // let s = base64::encode(encoded);
    // println!("{}", s);
    // let mut pair = RLWE::_new(secret, public, degree, q, fx.clone());

    // let pt = pair.decrypt_encoding("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", size);
    // print!("Decrypted Message: {}", pt);

    // let m1 = vec![1, 0, 1, 1];
    // let m2 = vec![1, 1, 1, 0, 1];
    // let c1 = encrypt(&pair.public, &m1, degree, q, &fx, size);
    // let c2 = encrypt(&pair.public, &m2, degree, q, &fx, size);
    // let mut c3 = add_ct(&c1, &c2, q);
    // for i in 0..(20 - 1) {
    //     c3 = add_ct(&c3, &c2, q);
    //     let pt = pair.decrypt(&c3, size);
    //     print!("Plaintext after {}: ", i);
    //     print_vec(&pt);
    // }
    //let ct = lwe::mult_ct(&c1, &c2, q);
    //let mpt = pair.decrypt(&ct, size);
    //print!("MULT PT: ");
    //print_vec(&mpt);

    /*let pt = pair.decrypt(&c3, size);
    print!("Plaintext: ");
    print_vec(&pt);
    return;*/
}

fn into_bits(byte: u8) -> Vec<i32> {
    let mut bits = Vec::new();
    for i in 0..8 {
        let mask = 1 << i;
        let bit_is_set = (mask & byte) > 0;
        bits.push(bit_is_set as i32);
    }
    bits
}

#[derive(Clone, Copy, Serialize, Deserialize)]
struct Tup2<A, B>(A, B);
impl<A, B> std::fmt::Display for Tup2<A, B>
where
    A: std::fmt::Display,
    B: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}:{})", self.0, self.1)
    }
}
impl Tup2<i32, i32> {
    fn encode(&self) -> Vec<u8> {
        let mut wtr = vec![];
        wtr.write_i32::<LittleEndian>(self.0).unwrap();
        wtr.write_i32::<LittleEndian>(self.1).unwrap();
        wtr
    }
}

struct CipherText<'r> {
    val: &'r Vec<ArrC<Tup2<i32, i32>, 8>>,
    //val: ArrC<Tup2<i32, i32>, 8>,
}
impl<'r> CipherText<'r> {
    fn new(val: &'r Vec<ArrC<Tup2<i32, i32>, 8>>) -> Self {
        CipherText { val }
    }
    fn encode(&self) -> Vec<u8> {
        //let testing_encoding = my_message1.into_iter().map(|x| x.0).collect::<Vec<_>>();
        let enc = self.val.into_iter().flat_map(|x| x.encode()).collect();
        enc
    }
}

struct ArrC<A, const C: usize>([A; C]);
impl ArrC<Tup2<i32, i32>, 8> {
    fn encode(&self) -> Vec<u8> {
        let mut wtr = vec![];
        for i in 0..8 {
            wtr.append(&mut self.0[i].encode());
        }
        wtr
    }
    /*fn encode(&self) -> Vec<u8> {
    let mut wtr = vec![];
    wtr.write_i32::<LittleEndian>(self.0).unwrap();
    wtr.write_i32::<LittleEndian>(self.1).unwrap();
    wtr
    }*/
}
impl<A, const C: usize> std::fmt::Display for ArrC<A, C>
where
    A: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[");
        for i in 0..C {
            write!(f, "{}", self.0[i]);
            if i < C - 1 {
                write!(f, ", ");
            }
        }
        write!(f, "]")
    }
}

struct Pair {
    secret: i32,
    q: i32,
    rand_v: Vec<i32>,
    public_v: Vec<i32>,
}
impl Pair {
    pub fn new(secret: i32, q: i32, rand_v: Vec<i32>, public_v: Vec<i32>) -> Self {
        Self {
            secret,
            q,
            rand_v,
            public_v,
        }
    }
    fn combine_ct(&mut self, ct1: Tup2<i32, i32>, ct2: Tup2<i32, i32>) -> Tup2<i32, i32> {
        Tup2(
            (ct1.0 + ct2.0).rem_euclid(self.q),
            (ct1.1 + ct2.1).rem_euclid(self.q),
        )
    }
    fn combine_ct_bytes(
        &mut self,
        byte1: &ArrC<Tup2<i32, i32>, 8>,
        byte2: &ArrC<Tup2<i32, i32>, 8>,
    ) -> ArrC<Tup2<i32, i32>, 8> {
        let mut combined = [Tup2(0, 0); 8];
        for i in 0..8 {
            combined[i] = self.combine_ct(byte1.0[i], byte2.0[i]);
        }
        ArrC(combined)
    }
    fn encrypt(&mut self, m: i32) -> Tup2<i32, i32> {
        let new_vec = self.sample(10);
        let init = Tup2(0, 0);
        let f = |Tup2(x1, y1), &Tup2(x2, y2)| Tup2((x1 + x2) % self.q, (y1 + y2) % self.q);
        let val: Tup2<i32, i32> = new_vec.iter().fold(init, f);
        // [(1, 3), (4, 5), (3, 7)]
        // (1 + 4 + 3, 3 + 5 + 7 + q / 2 * m)
        // (8, 15 + q / 2 * m)
        Tup2(val.0, (val.1 + self.q / 2 * m) % self.q)
    }
    fn encrypt_byte(&mut self, m: u8) -> ArrC<Tup2<i32, i32>, 8> {
        let mut byte = [Tup2(0, 0); 8];
        for i in 0..8 {
            let mask = 1 << i;
            let bit_is_set = (mask & m) > 0;
            byte[i] = self.encrypt(bit_is_set as i32)
        }
        ArrC(byte)
    }
    fn decrypt(&mut self, ct: Tup2<i32, i32>) -> i32 {
        if (ct.1 - self.secret * ct.0).rem_euclid(self.q) < self.q / 2 {
            return 0;
        }
        1
    }
    fn decrypt_byte(&mut self, byte: ArrC<Tup2<i32, i32>, 8>) -> u8 {
        let mut val: u8 = 0;
        for i in 0..8 {
            val += (self.decrypt(byte.0[i]) * 2_i32.pow(i as u32)) as u8;
        }
        val
    }
    fn sample(&mut self, val: i32) -> Vec<Tup2<i32, i32>> {
        let mut new_vec: Vec<Tup2<i32, i32>> = Vec::new();
        let mut rng = rand::thread_rng();
        let lim = self.rand_v.len();
        // random numbers: [1,5,10,3]
        // public key:     [4,1,5,8]
        // sample(2)
        // [(5,1), (3,8)]
        for _ in 0..val {
            let index = rng.gen_range(0..lim);
            new_vec.push(Tup2(
                self.rand_v.get(index).unwrap().clone(),
                self.public_v.get(index).unwrap().clone(),
            ));
        }
        new_vec
    }
}
pub fn print_vec(some_vec: &Vec<impl std::fmt::Display>) {
    print!("[");
    for i in 0..some_vec.len() {
        print!("{}", some_vec.get(i).unwrap());
        if i < some_vec.len() - 1 {
            print!(", ");
        }
    }
    println!("]");
}
/*

let secret: i32 = String::from("Re")
    .into_bytes()
    .iter()
    .fold(0, |acc, digit| (acc << 8) + *digit as i32); //99999;
let q: i32 = 7607;
let e: i32 = 200;
let seed: u64 = 200;

println!("Secret Key: {}", secret);
let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

let mut rand_v: Vec<i32> = Vec::new();
let mut public_v: Vec<i32> = Vec::new();
let between = Uniform::try_from(0..q).unwrap();
for i in 0..100 {
    rand_v.push(between.sample(&mut rng));
    public_v.push(rand_v.get(i).unwrap() * secret + rng.gen_range(0..e))
}

print!("Random Values: ");
print_vec(&rand_v);
print!("Public Key: ");
print_vec(&public_v);

let mut my_pair = Pair::new(secret, q, rand_v, public_v);
let some_vec = my_pair.sample(5);
print!("Sample Values: ");
print_vec(&some_vec);
let m1 = 1;
let ct1 = my_pair.encrypt(m1);
println!("Ciphertext for Message '{}': {}", m1, ct1);

let m2 = 1;
let ct2 = my_pair.encrypt(m2);
println!("Ciphertext for Message '{}': {}", m2, ct2);

let combined = my_pair.combine_ct(ct1, ct2);
println!("Ciphertext for Message '{}': {}", (m2 + m1) % 2, combined);

let pt = my_pair.decrypt(combined);
println!("Plaintext (Test): {}", pt);
let my_message1 = String::from("Hello") //Dit is een faketekst. Alles wat hier staat is slechts om een indruk te geven van het grafische effect van tekst op deze plek. Wat u hier leest is een voorbeeldtekst. Deze wordt later vervangen door de uiteindelijke tekst, die nu nog niet bekend is. De faketekst is dus een tekst die eigenlijk nergens over gaat. Het grappige is, dat mensen deze toch vaak lezen. Zelfs als men weet dat het om een faketekst gaat, lezen ze toch door.")
    .into_bytes()
    .into_iter()
    .map(|x| my_pair.encrypt_byte(x as u8))
    .collect::<Vec<_>>();
print!("Encrypted message: ");
let some_test = CipherText { val: &my_message1 };
let encode_test = some_test.encode();
println!("Encoding Length: {}", encode_test.len());
let s = base64::encode(encode_test);
println!("------ TEST BELOW ------");
println!("{}", s);

/*let testing_encoding = my_message1.into_iter().map(|x| x.0).collect::<Vec<_>>();
let encoded: Vec<u8> = bincode::serialize(&testing_encoding).unwrap();
print!("WUT");
print_vec(&encoded);
let s = base64::encode(encoded);
println!("{}", s);*/
let my_m = base64::decode("BQAAAAAAAAD2CgAAAAAAAG4UAAAAAAAA/RwAAAAAAAD7CAAAAAAAAMwOAAAAAAAAcwgAAAAAAAD8DQAAAAAAAKkKAAAAAAAAzg8AAAAAAACBHQAAAAAAAD8JAAAAAAAAfgAAAAAAAAC6AgAAAAAAANQEAAAAAAAA7BQAAAAAAAA5GgAAAAAAAMEbAAAAAAAA3g0AAAAAAABUAwAAAAAAAK4FAAAAAAAAzBcAAAAAAABpBwAAAAAAAPEXAAAAAAAADw0AAAAAAADSBgAAAAAAAEcGAAAAAAAAQx0AAAAAAAC6CgAAAAAAAA4IAAAAAAAAUAAAAAAAAACkDQAAAAAAAAYIAAAAAAAAIQ4AAAAAAAALDwAAAAAAALEaAAAAAAAA3xsAAAAAAAAoCAAAAAAAAMIJAAAAAAAAahAAAAAAAAANCgAAAAAAACUZAAAAAAAAZgwAAAAAAAAGDgAAAAAAAGIAAAAAAAAAoBoAAAAAAACyAwAAAAAAABAGAAAAAAAA/wIAAAAAAABpAwAAAAAAAGkFAAAAAAAAAAoAAAAAAADgHAAAAAAAAOkPAAAAAAAAdA4AAAAAAACyFwAAAAAAADUdAAAAAAAA/AkAAAAAAACJCAAAAAAAAEQFAAAAAAAAnxcAAAAAAAC7FQAAAAAAALYFAAAAAAAA5AAAAAAAAABMGQAAAAAAANEPAAAAAAAAHQ8AAAAAAADVEQAAAAAAAFccAAAAAAAAjw0AAAAAAACgGAAAAAAAADYDAAAAAAAAORYAAAAAAACeCAAAAAAAAD4KAAAAAAAAOxsAAAAAAABQCgAAAAAAAHACAAAAAAAAuBsAAAAAAAC4GQAAAAAAAOsGAAAAAAAA").unwrap();
let my_mm: Vec<[Tup2<i32, i32>; 8]> = bincode::deserialize(&my_m).unwrap();
let my_message1 = my_mm.into_iter().map(|x| ArrC(x)).collect::<Vec<_>>();
print_vec(&my_message1);

let my_message2 = String::from("")
    .into_bytes()
    .into_iter()
    .map(|x| my_pair.encrypt_byte(x as u8))
    .collect::<Vec<_>>();
print!("Encrypted message: ");

//print_vec(&my_message2);
//let mut combined_message/*Vec<Arr8<Tup2<i32,i32>>*/ = Vec::new();
/*for i in 0..my_message1.len() {
combined_message.push(
    my_pair.combine_ct_bytes(my_message1.get(i).unwrap(), my_message2.get(i).unwrap()),
)
}*/
let plain_text = my_message1
    .into_iter()
    .map(|x| my_pair.decrypt_byte(x))
    .collect::<Vec<_>>();
print!("decrypted message: ");
//print_vec(&plain_text);

//plain_text.iter().map(|x| )
let s = unsafe { String::from_utf8_unchecked(plain_text) }; //{
                                                            /*Ok(v) => v,
                                                            Err(_) => rng
                                                                .sample_iter(&Alphanumeric)
                                                                .take(1)
                                                                .map(char::from)
                                                                .collect(), //panic!("Invalid UTF-8 sequence: {}", e),
                                                                };*/

println!("result: {}", s);
*/

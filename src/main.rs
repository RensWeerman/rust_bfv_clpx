mod distributions;
mod fft;
mod fractions;
mod gbfv;
mod lwe_old;
mod lwe_testing;
mod montgomery;
mod ntt;
mod polynomial;
mod rlwe;

use lwe_testing::print_vec;
//use lwe_old::delta_polyn;
//use lwe_new::mult_ct;
use ntt::Ntt;
use num_bigint::BigInt;
use polynomial::{Parameters, Polynomial};
use std::cell::RefCell;
use std::fmt::Write;
use std::fs;
use std::rc::Rc;
use std::str::FromStr;
use std::time::Instant;
use std::usize;

use crate::rlwe::RLWE;

fn decrypt_vs_3decrypt(t: i64) {
    let parameters: Parameters<i64> = Parameters {
        degree: 4,
        q: 768100,
        t: vec![t],
        t_relin: 32,
        p: 4,
        log_range: 768100_i64.ilog(32) as usize + 1, //println!("{}", (self.q().checked_ilog(self.t_relin()).unwrap()) + 1);,
        root: 0,
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
    let parameters = polynomial::Parameters {
        degree: 4,
        q: 100,
        t: vec![20],
        t_relin: 32,
        p: 4,
        log_range: 0,
        root: 0,
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
    let parameters = polynomial::Parameters {
        degree: 4,
        q: 100,
        t: vec![20],
        t_relin: 32,
        p: 4,
        log_range: 0,
        root: 0,
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
    let parameters = polynomial::Parameters {
        degree: 4,
        q: 100,
        t: vec![20],
        t_relin: 32,
        p: 4,
        log_range: 0,
        root: 0,
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

    let res_1 = &c1 * &s;

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
        degree: 4,
        q: 100,
        t: vec![6],
        t_relin: 32,
        p: 4,
        log_range: 0,
        root: 0,
    };
    let mut buf = String::new();
    for i in 1..(parameters.q) {
        for _ in 1..40 {
            parameters.t = vec![i];
            let _p = &Rc::new(RefCell::new(parameters.clone()));
            let m = Polynomial::new(vec![5, 5, 5, 5], _p);

            let mut secret = RLWE::new(_p);

            let c = RLWE::encrypt(&secret.public, &m);
            let p = secret.decrypt(&c);
            write!(&mut buf, "{}\n", p).unwrap();
        }
    }
    fs::write(
        "/Users/rensweerman/Documents/PythonScripts/Q_vs_T/vectors.txt",
        buf,
    )
    .expect("Unable to write file");
}
fn timing_addition() {
    let parameters: Parameters<BigInt> = Parameters {
        degree: 16384,
        q: BigInt::from_str("2707685248005711112098097238309689728616229107667371065645293317992434822133490690704680432116610806966309052648878173576364033").unwrap(),
        t: vec![BigInt::try_from(65537).unwrap()],
        t_relin: BigInt::try_from(256).unwrap(),
        p: BigInt::try_from(4).unwrap(),
        log_range: 0,
        root: BigInt::from_str("160122687026568703260160820745333646557040369694442191388605703127148838809696464960419374948428326604331490224104058569398").unwrap(),
    };
    let _p = &Rc::new(RefCell::new(parameters));

    let p1 = Polynomial::uniform_sample(_p).mod_t(); //Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], _p);
    let p2 = Polynomial::uniform_sample(_p).mod_t(); //Polynomial::new(vec![3, 2, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], _p);
    let p3 = Polynomial::uniform_sample(_p).mod_t(); //Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], _p);

    let mut secret = RLWE::new(_p);
    //let rlk = secret.evaluate_key_gen();

    let c1 = RLWE::encrypt(&secret.public, &p1);
    let c2 = RLWE::encrypt(&secret.public, &p2);

    let start = Instant::now();

    for _ in 0..10 {
        let _ = RLWE::add_ct(&c1, &c2);
    }
    let duration = start.elapsed().as_micros();
    println!("Time taken for Addition: {:?}", duration);
}
fn timing_multiplication() {
    let parameters: Parameters<BigInt> = Parameters {
        degree: 32,
        q: BigInt::from_str("1532495540865518635130821056977027158796330141975560193").unwrap(),
        t: vec![BigInt::try_from(65537).unwrap()],
        t_relin: BigInt::try_from(256).unwrap(),
        p: BigInt::try_from(4).unwrap(),
        log_range: 0,
        root: BigInt::from_str("354363855417910436849378356982599500396178680480945").unwrap(),
    };
    let _p = &Rc::new(RefCell::new(parameters));

    let p1 = Polynomial::uniform_sample(_p).mod_t(); //Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], _p);
    let p2 = Polynomial::uniform_sample(_p).mod_t(); //Polynomial::new(vec![3, 2, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], _p);
    let p3 = Polynomial::uniform_sample(_p).mod_t(); //Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], _p);

    let mut secret = RLWE::new(_p);
    let rlk = secret.evaluate_key_gen();

    let c1 = RLWE::encrypt(&secret.public, &p1);
    let c2 = RLWE::encrypt(&secret.public, &p2);

    let start = Instant::now();

    for _ in 0..1 {
        let temp = RLWE::mult_ct(&c1, &c2);
        let _ = RLWE::relin_ct(temp, rlk.clone());
    }
    let duration = start.elapsed().as_micros();
    println!("Time taken for Multiplication: {:?}", duration);
}

fn timing_encryption() {
    let parameters: Parameters<BigInt> = Parameters {
        degree: 16384,
        q: BigInt::from_str("2707685248005711112098097238309689728616229107667371065645293317992434822133490690704680432116610806966309052648878173576364033").unwrap(),
        t: vec![BigInt::try_from(65537).unwrap()],
        t_relin: BigInt::try_from(256).unwrap(),
        p: BigInt::try_from(4).unwrap(),
        log_range: 0,
        root: BigInt::from_str("160122687026568703260160820745333646557040369694442191388605703127148838809696464960419374948428326604331490224104058569398").unwrap(),
    };
    let _p = &Rc::new(RefCell::new(parameters));

    let p1 = Polynomial::uniform_sample(_p).mod_t();

    let secret = RLWE::new(_p);

    let start = Instant::now();
    for _ in 0..10 {
        let _ = RLWE::encrypt(&secret.public, &p1);
    }
    let duration = start.elapsed().as_micros();
    println!("Time taken for encryption: {:?}", duration);
}
fn timing_decryption() {
    let parameters: Parameters<BigInt> = Parameters {
        degree: 16384,
        q: BigInt::from_str("2707685248005711112098097238309689728616229107667371065645293317992434822133490690704680432116610806966309052648878173576364033").unwrap(),
        t: vec![BigInt::try_from(65537).unwrap()],
        t_relin: BigInt::try_from(256).unwrap(),
        p: BigInt::try_from(4).unwrap(),
        log_range: 0,
        root: BigInt::from_str("160122687026568703260160820745333646557040369694442191388605703127148838809696464960419374948428326604331490224104058569398").unwrap(),
    };
    let _p = &Rc::new(RefCell::new(parameters));

    let p1 = Polynomial::uniform_sample(_p).mod_t();

    let mut secret = RLWE::new(_p);
    //let rlk = secret.evaluate_key_gen();
    let c1 = RLWE::encrypt(&secret.public, &p1);
    let start = Instant::now();
    for _ in 0..10 {
        let _ = secret.decrypt(&c1);
    }
    let duration = start.elapsed().as_micros();
    println!("Time taken for decryption: {:?}", duration);
}

fn test_iter_cooley_tukey_2() {
    let parameters = polynomial::Parameters {
        degree: 4,
        q: 7681,
        t: vec![6],
        t_relin: 32,
        p: 4,
        log_range: 0,
        root: 0,
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
        root: 0,
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
        root: BigInt::try_from(0).unwrap(),
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
fn broken_example() {
    let mut parameters = polynomial::Parameters {
        degree: 4,
        q: 100,
        t: vec![5],
        t_relin: 32,
        p: 4,
        log_range: 0,
        root: 0,
    };
    let _p = &Rc::new(RefCell::new(parameters.clone()));
    let m = Polynomial::new(vec![0, 0, 0, 0], _p);

    let mut secret = RLWE::new(_p);

    for i in 0..100 {
        let c = RLWE::encrypt(&secret.public, &m);
        let p = secret.decrypt(&c);
        if m != p {
            println!("FAILED on {i}, {m}, {p}");
            return;
        }
    }
}
fn main() {
    // for i in 2..300 {
    //     if (i % 40 == 0) {
    //         decrypt_vs_3decrypt(i);
    //     }
    // }
    // println!("ENCRYPT BELOW------");
    // example_encrypt();
    // println!("DECRYPT BELOW------");
    // example_decrypt();
    // println!("PK BELOW------");
    // p_q_gradient();
    // example_pk();
    // broken_example();
    // timing_addition();
    // timing_encryption();
    // timing_decryption();
    timing_multiplication();
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
        root: 0,
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
    let pt0 = &ct0 + &(&ct1 * &secret);
    print_vec(&pt0.val);
    let pt0 = pt0.get_mod();
    print_vec(&pt0.val);
    let pt1 = RLWE::delta_polyn(pt0);
    print_vec(&pt1.val);
    println!("Result:");
    print_vec(&pt1.mod_t().val);
    let real_m = &Polynomial::new(vec![1, 0, 0, 1], _p) * &Polynomial::new(vec![0, 0, 1, 1], _p);
    print_vec(&real_m.get_mod().mod_q().val);
}

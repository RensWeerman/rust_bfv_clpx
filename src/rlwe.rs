use std::{cell::RefCell, rc::Rc};

use crate::{
    fractions::Fraction,
    montgomery::{square_multiply, square_multiply_no_mod},
    ntt::Ntt,
    polynomial::{make_fx, Parameters, Polynomial, Value},
    print_vec,
};

pub type Params<Values> = Rc<Parameters<Values>>;
pub struct RLWE<Val> {
    secret: Polynomial<Val>,
    pub public: (Polynomial<Val>, Polynomial<Val>),
    pub parameters: Params<Val>,
}
impl<Val> RLWE<Val>
where
    Val: Value,
{
    pub fn _new(sec: Vec<Val>, parameters: &Params<Val>) -> Self {
        RLWE {
            secret: Polynomial::new(sec, parameters),
            public: (
                Polynomial::new(vec![], parameters),
                Polynomial::new(vec![], parameters),
            ),
            parameters: Rc::clone(parameters),
        }
    }
    pub fn t_relin(&self) -> Val {
        return self.parameters.t_relin.clone();
    }
    pub fn q(&self) -> Val {
        return self.parameters.q.clone();
    }
    pub fn p(&self) -> Val {
        return self.parameters.p.clone();
    }
    pub fn t(&self) -> Val {
        return self.parameters.t.clone().get(0).unwrap().clone();
    }
    pub fn t_full(&self) -> Vec<Val> {
        return self.parameters.t.clone();
    }
    pub fn log_range(&self) -> usize {
        return self.parameters.log_range.clone();
    }
    pub fn parameters(&self) -> &Params<Val> {
        return &self.parameters;
    }
    pub fn new(parameters: &Params<Val>) -> Self {
        let secret = Polynomial::normal_sample(parameters);
        let e = Polynomial::normal_sample(parameters);
        let a = Polynomial::uniform_sample(parameters);
        let test = &e + &a;
        println!("????");
        let test = &e * &a;
        println!("????");
        let pk = (-&(&(&a * &secret) + &e), a);
        RLWE {
            secret,
            public: pk,
            parameters: Rc::clone(parameters),
        }
    }
    pub fn new_given_secret(secret: Polynomial<Val>, parameters: &Params<Val>) -> Self {
        //let rng = rand::thread_rng();
        let e = Polynomial::normal_sample(parameters);
        let a = Polynomial::uniform_sample(parameters); //mod_polyn(uniform_sample(degree + 1, q), &fx, q);
        let pk = (
            &Polynomial::new(vec![], parameters) - &(&(&a * &secret) + &e),
            (&a + &Polynomial::new(vec![], parameters)),
        );
        /*print!("Secret RWLE: ");
        print_vec(&secret);
        print!("Noise: ");
        print_vec(&e);
        print!("Public RWLE: ");
        print_vec(&pk.0);
        print_vec(&pk.1);*/
        RLWE {
            secret,
            public: pk,
            parameters: Rc::clone(parameters),
        }
    }
    /*pub fn decrypt_encoding(&mut self, ct: &str, size: i32) -> String {
    let my_m = base64::decode(ct).unwrap();
    let cct: Vec<(Vec<i32>, Vec<i32>)> = bincode::deserialize(&my_m).unwrap();
    //let ct = rlwe::encrypt(&a.public, &m, degree, q, &fx, size);
    let pt = cct
        .into_iter()
        .map(|x| self.decrypt(&x, size))
        .collect::<Vec<_>>();
    print!("Decrypted Message: ");
    let mut buf = String::new();
    for chunk in pt {
        let mut val: u8 = 0;
        for i in 0..8 {
            //println!("{}", chunk.get(i).unwrap_or(&0));
            val += (chunk.get(i).unwrap_or(&0) * 2_i32.pow(i as u32)) as u8;
        }
        write!(&mut buf, "{}", val as char);
    }
    buf
    }*/
    pub fn decrypt(&mut self, ct: &(Polynomial<Val>, Polynomial<Val>)) -> Polynomial<Val> {
        let mut pt = &ct.0 + &(&ct.1 * &self.secret);
        pt = &pt + &Polynomial::new(vec![], &self.parameters);
        //println!("PT before mod {}", pt);
        pt = pt.mod_q();
        //println!("PT going into delta_poly {}", pt);
        pt = RLWE::delta_polyn(pt);
        pt.mod_t()
    }
    pub fn decrypt_degree_3(
        &mut self,
        ct: &(Polynomial<Val>, Polynomial<Val>, Polynomial<Val>),
    ) -> Polynomial<Val> {
        let mut pt = &(&ct.0 + &(&ct.1 * &self.secret)) + &(&ct.2 * &(&self.secret * &self.secret));
        pt = &pt + &Polynomial::new(vec![], &self.parameters);
        pt = pt.mod_q();
        pt = RLWE::delta_polyn(pt);
        pt = pt.mod_t();
        return pt;
    }
    pub fn evaluate_key_gen(&self) -> Vec<(Polynomial<Val>, Polynomial<Val>)> {
        let mut rlk = Vec::new();
        let secret_2 = &self.secret * &self.secret;
        let x = self.log_range();
        //println!("{}", (self.q().checked_ilog(self.t_relin()).unwrap()) + 1);
        for i in 0..x {
            let a = /*8Polynomial::new(vec![], self.parameters()); */Polynomial::normal_sample(self.parameters());
            let e = /*Polynomial::new(vec![], self.parameters()); */Polynomial::normal_sample(self.parameters());
            let mut k_0 = Polynomial::new(vec![], self.parameters());
            let k_1 = &(&a * &self.secret) + &e;
            k_0 = &k_0 - &k_1;
            k_0 = &k_0
                + &(&secret_2
                    * &Polynomial::new(
                        vec![square_multiply_no_mod(self.t_relin(), i, self.q())],
                        self.parameters(),
                    ));
            k_0.mod_q();
            k_0 = &k_0 % &make_fx(self.parameters());
            rlk.push((k_0, a));
        }
        rlk
    }
    // pub fn relin_key_gen_1(&self) -> Vec<(Polynomial<Val>, Polynomial<Val>)> {
    //     let degree = self.parameters().borrow().degree;
    //     let s = self.secret.clone();
    //     // l is the number of levels to decompose s^2 and c_2 into.
    //     // l is a function of base (T in the paper): l = floor(log_T(q)).
    //     let l = (self.q() as f64).log(self.t_relin() as f64).floor() as usize + 1;
    //     println!("{}", l);

    //     let val = (0..l)
    //         .map(|i| {
    //             let a = Polynomial::new(vec![], self.parameters()); //Polynomial<Val>::normal_sample();
    //             let e = Polynomial::new(vec![], self.parameters()); //Polynomial<Val>::normal_sample();
    //             let base_i = self.t().pow(i as u32);
    //             let rlk_i_raw = &(&Polynomial::new(vec![], self.parameters()) - &(&(&a * &s) + &e))
    //                 + &(&(&s * &s) * &Polynomial::new(vec![base_i], self.parameters()));
    //             let rlk_i = &rlk_i_raw % &make_fx(self.parameters());
    //             (rlk_i, a)
    //         })
    //         .collect();
    //     val
    // }
    // pub fn evaluate_key_gen_2(&self) -> (Polynomial<Val>, Polynomial<Val>) {
    //     let secret_2 = &self.secret * &self.secret;
    //     let a = Polynomial::uniform_sample(self.parameters());
    //     let e = Polynomial::normal_sample(self.parameters());
    //     let mut k_0 = Polynomial::new(vec![], self.parameters());
    //     let k_1 = &(&a * &self.secret);
    //     k_0 = &k_0 - &k_1;
    //     k_0 = &k_0 - &e;
    //     k_0 = &k_0 + &(&secret_2 * &Polynomial::new(vec![self.p()], self.parameters()));
    //     (k_0, a)
    // }
    pub fn compute_encryption(
        pk: &(Polynomial<Val>, Polynomial<Val>),
        m: &Polynomial<Val>,
        u: Polynomial<Val>,
        e1: Polynomial<Val>,
        e2: Polynomial<Val>,
    ) -> (Polynomial<Val>, Polynomial<Val>) {
        let _ct0 = &(&pk.0 * &u) + &e1;
        let __ct0 = &_ct0 + &RLWE::delta_mult(m.clone());
        //println!("WHAT?1: {}", &RLWE::mult_delta(m.clone()));
        //println!("WHAT?2: {}", &RLWE::delta_mult(m.clone()));
        let __ct1 = &(&pk.1 * &u) + &e2;
        (__ct0, __ct1)
    }
    pub fn encrypt(
        pk: &(Polynomial<Val>, Polynomial<Val>),
        m: &Polynomial<Val>,
    ) -> (Polynomial<Val>, Polynomial<Val>) {
        let u = Polynomial::normal_sample(m.parameters());
        //print!("u polynomial: ");
        //print_vec(&u.val);
        let e1 = Polynomial::normal_sample(m.parameters());
        let e2 = Polynomial::normal_sample(m.parameters());
        RLWE::compute_encryption(pk, m, u, e1, e2)
    }
    pub fn flatten_m(m: Polynomial<Val>) -> Polynomial<Val> {
        let params = m.parameters();
        let t = Polynomial::new(m.t_full(), params);
        let f_params = Fraction::frac_params(params);
        let _fp = &Rc::new(f_params);
        let t_frac = Fraction::frac_poly(t.clone(), _fp);
        let fx = make_fx(params);
        let fx_frac = Fraction::frac_poly(fx.clone(), _fp);
        let old_s = Polynomial::extended_gcd(&t_frac, &fx_frac);
        println!("CURRENT T: {} and FX: {}", t, fx);
        println!("???? OLD S: {}", old_s);
        let mut temp = vec![];
        for x in old_s.val {
            temp.push(x.reduce());
        }
        let delta = Polynomial::new(temp, _fp);
        println!("THIS IS 1/T: {}", delta);
        let m_frac = Fraction::frac_poly(m.clone(), _fp);
        //let res = Fraction::round_multiply(&m_frac, &delta, m.parameters());
        let res = m_frac.old_mul(&delta);
        println!("m/t IS: {}", res);
        let mut temp = vec![];
        for mut x in res.val {
            x.flatten();
            temp.push(x.reduce());
        }
        let before_mult_t = Polynomial::new(temp, _fp);
        println!("before final t: {}, {}", before_mult_t, t_frac);
        let full_res = Fraction::round_multiply(&before_mult_t, &t_frac, m.parameters());
        println!("done, {}", full_res);
        return full_res;
    }
    pub fn delta_mult(c: Polynomial<Val>) -> Polynomial<Val> {
        //println!("DELTA MULT");
        let params = c.parameters();
        let t = Polynomial::new(c.t_full(), params);
        let f_params = Fraction::frac_params(params);
        let _fp = &Rc::new(f_params);
        let t_frac = Fraction::frac_poly(t, _fp);
        let fx = make_fx(params);
        let fx_frac = Fraction::frac_poly(fx, _fp);
        let old_s = Polynomial::extended_gcd(&t_frac, &fx_frac);
        let mut temp = vec![];
        for x in old_s.val {
            temp.push(_fp.q.clone() * x.reduce());
        }
        let delta = Polynomial::new(temp, _fp);
        let c_frac = Fraction::frac_poly(c.clone(), _fp);
        let res = Fraction::round_multiply(&c_frac, &delta, c.parameters());
        //println!("C_FRAC WAS: {}, RES IS: {}", c_frac, res);
        return res;
    }
    pub fn mult_delta(c: Polynomial<Val>) -> Polynomial<Val> {
        return RLWE::delta_mult(c);
        let params = c.parameters();
        let t = Polynomial::new(c.t_full(), params);
        if (t.degree() == 1) {
            let q = c.q();
            //println!("Incoming poly: {}", c);
            let temp = Polynomial::new(vec![q], c.parameters());
            //println!("q: {}", temp);
            let next_temp = &temp / &Polynomial::new(temp.t_full(), c.parameters());
            //println!("q/t: {}", next_temp);
            let test = next_temp.old_mul(&c);
            //println!("end result: {}", test);
            return test;
        }

        let t = Polynomial::new(c.t_full(), params);
        let range = params.degree as usize / t.degree();
        println!("LOOP AMOUNT: {}", range);
        let mut x = vec![Val::zero(); params.degree as usize];
        let val = (t.val.get(0).unwrap()).clone();

        for i in 1..(range + 1) {
            let res = square_multiply(Val::zero() - val.clone(), i - 1, params.q.clone());
            println!("Val : {} ", res);
            x[params.degree as usize - (t.degree() * i)] = res;
        }
        let delta = Polynomial::new(x.to_vec(), params);
        println!("Delta before scaling: {}", delta);
        let div =
            square_multiply((t.val.get(0).unwrap()).clone(), range, c.q().clone()) + Val::one();
        let scaler = (Val::zero() - c.q().clone()) / div.clone();
        println!("Num : {}, Den: {}", c.q(), div);
        println!("Scaler: {}", scaler);
        let delta = &delta * &Polynomial::new(vec![scaler], params);
        //let delta = &delta % &t;
        println!("Testing other delta: {}", delta);
        let temp = &c * &delta;
        println!("Before: {}", c);
        println!("Result: {}", temp);
        return temp;

        // let temp = &c * &Polynomial::new(vec![c.q().clone()], c.parameters()); // * &Polynomial::new(vec![q], c.parameters());
        // let t = temp.t_full();
        // println!("Delta stuff now: {}", temp);
        // let res = &temp / &Polynomial::new(t, temp.parameters());
        // println!("Delta stuff now: {}", res);
        // return res;
        //

        // let t = Polynomial::new(c.t_full(), c.parameters());
        // let q = Polynomial::new(vec![c.q()], c.parameters());
        // let temp = &t / &q;
        // let temp = &c / &temp;
        // return temp;
    }
    pub fn delta_polyn(c: Polynomial<Val>) -> Polynomial<Val> {
        let t = Polynomial::new(c.t_full(), c.parameters());
        let mut temp = c.clone().old_mul(&t);
        println!("HERE");
        //println!("Hmm: {}", c.clone());
        //println!("T is: {}", t);
        //println!("Temp is: {}", temp);
        temp.val = temp
            .val
            .into_iter()
            .map(|i| (i + c.q() / (Val::one() + Val::one())) / c.q())
            .collect();
        println!("WOW DAS SLOOM");
        //println!("Temp is: {}", temp);
        temp
        //&temp % &make_fx(c.parameters())
    }

    pub fn add_ct(
        a: &(Polynomial<Val>, Polynomial<Val>),
        b: &(Polynomial<Val>, Polynomial<Val>),
    ) -> (Polynomial<Val>, Polynomial<Val>) {
        let c0 = (&a.0 + &b.0).mod_q();
        let c1 = (&a.1 + &b.1).mod_q();
        (c0, c1)
    }

    pub fn mult_ct(
        ct1: &(Polynomial<Val>, Polynomial<Val>),
        ct2: &(Polynomial<Val>, Polynomial<Val>),
    ) -> (Polynomial<Val>, Polynomial<Val>, Polynomial<Val>) {
        // let c0 = &ct1.0 * &ct2.0;
        // let c1 = &(&ct1.0 * &ct2.1) + &(&ct1.1 * &ct2.0);
        // let c2 = &ct1.1 * &ct2.1;
        let c0 = ct1.0.old_mul(&ct2.0);
        println!("1");
        let c1 = &(ct1.0.old_mul(&ct2.1)) + &(ct1.1.old_mul(&ct2.0));
        println!("2");
        let c2 = ct1.1.old_mul(&ct2.1);
        println!("333");
        //(c0, c1, c2)
        (
            RLWE::delta_polyn(c0).get_mod().mod_q(),
            RLWE::delta_polyn(c1).get_mod().mod_q(),
            RLWE::delta_polyn(c2).get_mod().mod_q(),
        )
    }
    pub fn new_mult_ct(
        ct1: &(Polynomial<Val>, Polynomial<Val>),
        ct2: &(Polynomial<Val>, Polynomial<Val>),
    ) -> (Polynomial<Val>, Polynomial<Val>, Polynomial<Val>) {
        let degree = ct1.0.parameters().degree;
        let q = &ct1.0.parameters().q.clone();
        let psi = &Ntt::second_primitive_root(q, degree);
        let ct0 = Ntt::iter_cooley_tukey(ct1.0.clone(), psi, q);
        let ct1 = Ntt::iter_cooley_tukey(ct1.1.clone(), psi, q);
        let dt0 = Ntt::iter_cooley_tukey(ct2.0.clone(), psi, q);
        let dt1 = Ntt::iter_cooley_tukey(ct2.1.clone(), psi, q);

        let mut c0 = Ntt::element_wise_mult(ct0.clone(), dt0.clone(), q);
        c0 = RLWE::delta_polyn(c0);

        let mut c2 = Ntt::element_wise_mult(ct1.clone(), dt1.clone(), q);
        c2 = RLWE::delta_polyn(c2);

        let c10 = Ntt::element_wise_mult(ct0, dt1, q);
        let c11 = Ntt::element_wise_mult(ct1, dt0, q);
        let mut c1 = &c10 + &c11;
        c1 = RLWE::delta_polyn(c1);

        c0 = Ntt::iter_gentleman_sande(c0, psi, q);
        c1 = Ntt::iter_gentleman_sande(c1, psi, q);
        c2 = Ntt::iter_gentleman_sande(c2, psi, q);

        (
            c0.get_mod().mod_q(),
            c1.get_mod().mod_q(),
            c2.get_mod().mod_q(),
        )
    }
    pub fn relin_ct(
        ct: (Polynomial<Val>, Polynomial<Val>, Polynomial<Val>),
        rlk: Vec<(Polynomial<Val>, Polynomial<Val>)>,
    ) -> (Polynomial<Val>, Polynomial<Val>) {
        let mut c_0 = Polynomial::new(vec![], ct.0.parameters());
        let mut c_1 = Polynomial::new(vec![], ct.0.parameters());
        let mut c_2 = &Polynomial::new(vec![], ct.0.parameters()) + &ct.2;
        let mut testing = Polynomial::new(vec![], ct.0.parameters());
        // print!("c_2: ");
        // print_vec(&ct.2.val);
        for i in 0..ct.0.log_range() {
            //println!("{}", i);
            let c_2r = &c_2 % &Polynomial::new(vec![ct.0.t_relin()], ct.0.parameters());
            let c_2r = c_2r.mod_t_relin();
            //let c_2r = &c_2 % &Polynomial<Val>::new(vec![t.pow(i)]);
            testing = &testing
                + &(&c_2r
                    * &Polynomial::new(
                        vec![square_multiply_no_mod(ct.0.t_relin(), i, ct.0.q())],
                        ct.0.parameters(),
                    ));
            print_vec(&c_2r.val);
            c_2 = &c_2 - &c_2r; //* &Polynomial<Val>::new(vec![t.pow(i)]));
            c_2 = &c_2 / &Polynomial::new(vec![ct.0.t_relin()], ct.0.parameters());
            let val = rlk.get(i as usize).unwrap();
            //let val = rlk.pop().unwrap();
            c_0 = &c_0 + &(&c_2r * &val.0);
            c_1 = &c_1 + &(&c_2r * &val.1);
        }
        // print!("TEST: ");
        // print_vec(&testing.get_mod().val);
        // print!("ct.0: ");
        // print_vec(&ct.0.get_mod().val);
        // print!("c_0: ");
        // print_vec(&c_0.get_mod().val);

        // print!("ct.1: ");
        // print_vec(&ct.1.get_mod().val);
        // print!("c_1: ");
        // print_vec(&c_1.get_mod().val);

        ((&ct.0 + &c_0), (&ct.1 + &c_1))

        /*
        let mut rlk = Vec::new();
        for i in 0..(Q.checked_ilog(t).unwrap()) {
            let a = Polynomial<Val>::normal_sample();
            let e = Polynomial<Val>::normal_sample();
            let mut k_0 = &(&a * &self.secret) + &e;
            k_0 = &k_0 + &(&Polynomial<Val>::new(vec![t.pow(i)]) * &(&self.secret * &self.secret));
            rlk.push((k_0, a));
        }
        rlk
        */
    }

    // pub fn other_relin(
    //     ct: (Polynomial<Val>, Polynomial<Val>, Polynomial<Val>),
    //     rlk: Vec<(Polynomial<Val>, Polynomial<Val>)>,
    // ) -> (Polynomial<Val>, Polynomial<Val>) {
    //     let degree = ct.0.degree();

    //     // Decompose c_2 in base T (rlk_base), such that:
    //     // $ c_2 = \sum_{i=0}^l c_2^(i) T^i $ with $ c_2^(i) \in R_T $
    //     let c_2_dec: Vec<Polynomial<Val>> = ct.2.decompose();

    //     // Calculate the contributions of the decomposed c_2 for c_0 and c_1.
    //     let mut c_2_0 = Polynomial::new(vec![zero(); degree], ct.0.parameters());
    //     let mut c_2_1 = Polynomial::new(vec![zero(); degree], ct.0.parameters());
    //     for i in 0..(rlk.len() as usize) {
    //         // Calculate the sum of the first entry of the relinearization key and decomposed c_2:
    //         // $ \sum_{i=0}^l rlk[i][0] * c_2^(i) $
    //         c_2_0 = &c_2_0 + &(&rlk[i].0 * &c_2_dec[i]);

    //         // Calculate the sum of the second entry of the relinearization key and decomposed c_2:
    //         // $ \sum_{i=0}^l rlk[i][1] * c_2^(i) $
    //         c_2_1 = &c_2_1 + &(&rlk[i].1 * &c_2_dec[i]);
    //     }
    //     ((&ct.0 + &c_2_0), (&ct.1 + &c_2_1))
    // }

    pub fn div_p(c: Polynomial<Val>) -> Polynomial<Val> {
        let mut temp = c.clone();
        temp.val = temp
            .val
            .into_iter()
            .map(|i| ((i + (c.p() / (Val::one() + Val::one()))) / c.p()))
            .collect();
        temp.rm_trailing_zeroes();
        &temp % &make_fx(c.parameters())
    }
    pub fn relin_ct_2(
        ct: (Polynomial<Val>, Polynomial<Val>, Polynomial<Val>),
        rlk: (Polynomial<Val>, Polynomial<Val>),
    ) -> (Polynomial<Val>, Polynomial<Val>) {
        let mut c_0 = &ct.2 * &rlk.0;
        let mut c_1 = &ct.2 * &rlk.1;
        c_0 = RLWE::div_p(c_0);
        c_1 = RLWE::div_p(c_1);
        (&ct.0 + &c_0, &ct.1 + &c_1)
        /*
        let mut rlk = Vec::new();
        for i in 0..(Q.checked_ilog(t).unwrap()) {
            let a = Polynomial<Val>::normal_sample();
            let e = Polynomial<Val>::normal_sample();
            let mut k_0 = &(&a * &self.secret) + &e;
            k_0 = &k_0 + &(&Polynomial<Val>::new(vec![t.pow(i)]) * &(&self.secret * &self.secret));
            rlk.push((k_0, a));
        }
        rlk
        */
    }
    // pub fn relinearization_1(
    //     ct: (Polynomial<Val>, Polynomial<Val>, Polynomial<Val>),
    //     rlk: Vec<(Polynomial<Val>, Polynomial<Val>)>,
    // ) -> (Polynomial<Val>, Polynomial<Val>) {
    //     let degree = ct.1.degree();

    //     // Decompose c_2 in base T (rlk_base), such that:
    //     // $ c_2 = \sum_{i=0}^l c_2^(i) T^i $ with $ c_2^(i) \in R_T $
    //     let c_2_dec: Vec<Polynomial<Val>> = ct.2.decompose();

    //     // Calculate the contributions of the decomposed c_2 for c_0 and c_1.
    //     let mut c_2_0 = Polynomial::new(vec![zero(); degree], ct.0.parameters());
    //     let mut c_2_1 = Polynomial::new(vec![zero(); degree], ct.0.parameters());
    //     for i in 0..ct.0.log_range() {
    //         // Calculate the sum of the first entry of the relinearization key and decomposed c_2:
    //         // $ \sum_{i=0}^l rlk[i][0] * c_2^(i) $
    //         c_2_0 = &c_2_0 + &(&rlk[i as usize].0.clone() * &c_2_dec[i as usize].clone());

    //         // Calculate the sum of the second entry of the relinearization key and decomposed c_2:
    //         // $ \sum_{i=0}^l rlk[i][1] * c_2^(i) $
    //         c_2_1 = &c_2_1 + &(&rlk[i as usize].1.clone() * &c_2_dec[i as usize].clone());
    //     }

    //     (&ct.0 + &c_2_0, &ct.1 + &c_2_1)
    // }
}
pub fn get_def_params() -> Parameters<i64> {
    return Parameters {
        degree: 4,
        q: 7681,
        t: vec![6],
        t_relin: 32,
        p: 4,
        log_range: 7681_i64.ilog(32) as usize + 1, //println!("{}", (self.q().checked_ilog(self.t_relin()).unwrap()) + 1);,
        root: 0,
    };
}

#[cfg(test)]
mod tests {
    use crate::rlwe;

    use super::*;

    #[test]
    fn encrypt_decrypt() {
        let _p = &Rc::new(get_def_params());
        let mut keys = RLWE::new(_p);
        let m = Polynomial::new(vec![1, 5, 0, 0], _p);
        let ct = RLWE::encrypt(&keys.public, &m);
        let pt = keys.decrypt(&ct);
        assert_eq!(m, pt);
    }
    #[test]
    fn encrypt_add_decrypt() {
        let _p = &Rc::new(get_def_params());
        let mut keys = RLWE::new(_p);
        let m1 = Polynomial::new(vec![1, 0, 0, 1], _p);
        let m2 = Polynomial::new(vec![0, 0, 1, 1], _p);
        let m = &m1 + &m2;
        let ct1 = RLWE::encrypt(&keys.public, &m1);
        let ct2 = RLWE::encrypt(&keys.public, &m2);
        let ct = RLWE::add_ct(&ct1, &ct2);
        let pt = keys.decrypt(&ct);
        assert_eq!(m, pt);
    }
    #[test]
    fn encrypt_mult_decrypt() {
        let parameters: Parameters<i64> = Parameters {
            degree: 4,
            q: 7681,
            t: vec![10],
            t_relin: 32,
            p: 4,
            log_range: 7681_i64.ilog(32) as usize + 1, //println!("{}", (self.q().checked_ilog(self.t_relin()).unwrap()) + 1);,
            root: 0,
        };
        let _p = &Rc::new(parameters);
        println!("log_range: {}", _p.log_range);
        let mut keys = RLWE::new(_p);
        let m1 = Polynomial::new(vec![4, 19, 0, 7], _p);
        let m2 = Polynomial::new(vec![1, 1, 8, 2], _p);
        let m = &m1.old_mul(&m2);
        let ct1 = RLWE::encrypt(&keys.public, &m1);
        let ct2 = RLWE::encrypt(&keys.public, &m2);
        let ct = RLWE::mult_ct(&ct1, &ct2);
        println!("Next three are the ciphertext before relin");
        print_vec(&ct.0.val);
        print_vec(&ct.1.val);
        print_vec(&ct.2.val);
        let rlk = keys.evaluate_key_gen();
        // for i in &rlk {
        //     print_vec(&i.0.val);
        //     print_vec(&i.1.val);
        // }
        let ct = RLWE::relin_ct(ct, rlk);
        println!("DECRYPTING");
        let pt = keys.decrypt(&ct);
        print!("THE FINAL VALUE: ");
        print_vec(&pt.val);
        println!("message: {}", m.get_mod());
        print_vec(&m.get_mod().mod_t().val);
        //print_vec
        assert_eq!(m.get_mod().mod_t(), pt);
    }
    #[test]
    fn is_mult_working() {
        let parameters: Parameters<i64> = Parameters {
            degree: 8,
            q: 7681,
            t: vec![3],
            t_relin: 32,
            p: 4,
            log_range: 7681_i64.ilog(32) as usize + 1, //println!("{}", (self.q().checked_ilog(self.t_relin()).unwrap()) + 1);,
            root: 0,
        };
        let _p = &Rc::new(parameters);
        println!("log_range: {}", _p.log_range);
        let mut keys = RLWE::new(_p);
        let m1 = Polynomial::new(vec![1, 0, 0, 1, 0, 1, 1, 0], _p);
        let m2 = Polynomial::new(vec![1, 1, 0, 0, 1, 1, 0, 1], _p);
        let m = m1.old_mul(&m2);
        let ct1 = RLWE::encrypt(&keys.public, &m1);
        let ct2 = RLWE::encrypt(&keys.public, &m2);
        let ct = RLWE::mult_ct(&ct1, &ct2);
        println!("Next three are the ciphertext before relin");
        print_vec(&ct.0.val);
        print_vec(&ct.1.val);
        print_vec(&ct.2.val);
        let rlk = keys.evaluate_key_gen();
        // for i in &rlk {
        //     print_vec(&i.0.val);
        //     print_vec(&i.1.val);
        // }
        let ct = RLWE::relin_ct(ct, rlk);
        println!("DECRYPTING");
        let pt = keys.decrypt(&ct);
        print!("THE FINAL VALUE: ");
        print_vec(&pt.val);
        print_vec(&m.get_mod().mod_t().val);
        //print_vec
        assert_eq!(m.get_mod().mod_t(), pt);
    }
    #[test]
    fn encrypt_mult_decrypt_random() {
        let parameters: Parameters<i64> = Parameters {
            degree: 4,
            q: 7681,
            t: vec![10],
            t_relin: 32,
            p: 4,
            log_range: 7681_i64.ilog(32) as usize + 1, //println!("{}", (self.q().checked_ilog(self.t_relin()).unwrap()) + 1);,
            root: 0,
        };
        let _p = &Rc::new(parameters);
        println!("log_range: {}", _p.log_range);
        let mut keys = RLWE::new(_p);
        let m1 = Polynomial::uniform_sample(_p).mod_t();
        let m2 = Polynomial::uniform_sample(_p).mod_t();
        let m = &m1.old_mul(&m2);
        let ct1 = RLWE::encrypt(&keys.public, &m1);
        let ct2 = RLWE::encrypt(&keys.public, &m2);
        let ct = RLWE::mult_ct(&ct1, &ct2);
        println!("Next three are the ciphertext before relin");
        print_vec(&ct.0.val);
        print_vec(&ct.1.val);
        print_vec(&ct.2.val);
        let rlk = keys.evaluate_key_gen();
        // for i in &rlk {
        //     print_vec(&i.0.val);
        //     print_vec(&i.1.val);
        // }
        let ct = RLWE::relin_ct(ct, rlk);
        println!("DECRYPTING");
        let pt = keys.decrypt(&ct);
        print!("THE FINAL VALUE: ");
        print_vec(&pt.val);
        print_vec(&m.get_mod().mod_t().val);
        //print_vec
        assert_eq!(m.get_mod().mod_t(), pt);
    }
    #[test]
    fn testing_mod_q() {
        let _p = &Rc::new(get_def_params());
        let q = _p.q;
        assert_eq!(
            Polynomial::new(vec![1, (-2_i64).rem_euclid(q), 3, 6], _p),
            Polynomial::new(vec![1 + q, -2 + q, 3 + q, 6 + q], _p).mod_q()
        );
    }
    // fn relin_gen() {
    //     let _p = &Rc::new(RefCell::new(PARAMETERS));
    //     let secret = RLWE::new(_p);
    //     let relin_mine = secret.evaluate_key_gen();
    //     let relin_other = secret.relin_key_gen_1();
    //     assert_eq!(relin_mine, relin_other);
    // }
    //
    //#[test]
    fn comparing_old_new_mul() {
        let _p = &Rc::new(get_def_params());
        let secret = RLWE::new(_p);
        let ciphertext_1 = RLWE::encrypt(
            &secret.public,
            &Polynomial::new(vec![1, 999, 4, 5, 1, 0, 3, 4], _p),
        );
        let ciphertext_2 = RLWE::encrypt(
            &secret.public,
            &Polynomial::new(vec![6, 4, 2, 1, 0, 1, 2, 304, 463, 342, 325], _p),
        );
        let res_1 = RLWE::mult_ct(&ciphertext_1, &ciphertext_2);
        let res_2 = RLWE::new_mult_ct(&ciphertext_1, &ciphertext_2);

        print_vec(&res_1.0.get_mod().val);
        print_vec(&res_2.0.get_mod().val);
        assert_eq!(res_1, res_2);
    }
    //#[test]
    fn testing_flatten() {
        let parameters: Parameters<i64> = Parameters {
            degree: 4,
            q: 7681,
            t: vec![20, 1],
            t_relin: 32,
            p: 1,
            log_range: 7681_i64.ilog(32) as usize + 1, //println!("{}", (self.q().checked_ilog(self.t_relin()).unwrap()) + 1);,
            root: 0,
        };
        let _p = &Rc::new(parameters);

        println!("make_fx : {}", make_fx(_p));
        let mut keys = RLWE::new(_p);
        let test_t = Polynomial::new(_p.t.clone(), _p);
        let m = Polynomial::new(vec![18], _p);
        let flat_m = RLWE::flatten_m(m.clone());
        let ct = RLWE::encrypt(&keys.public, &flat_m);
        let pt = keys.decrypt(&ct);
        let mut pt = &pt % &test_t;
        //let mut pt = pt.mod_t();
        pt.rm_trailing_zeroes();
        println!("Plaintext: {}", pt);
        println!("FLATTEN RES: {}", flat_m);
        assert_eq!(m, pt);
    }
    #[test]
    fn clpx_encrypt_decrypt() {
        let parameters: Parameters<i64> = Parameters {
            degree: 4,
            q: 7681,
            t: vec![4, 1],
            t_relin: 32,
            p: 1,
            log_range: 7681_i64.ilog(32) as usize + 1, //println!("{}", (self.q().checked_ilog(self.t_relin()).unwrap()) + 1);,
            root: 0,
        };
        let _p = &Rc::new(parameters);

        println!("make_fx : {}", make_fx(_p));
        let mut keys = RLWE::new(_p);
        let test_t = Polynomial::new(_p.t.clone(), _p);
        //let m = Polynomial::new(vec![2, 2, 3, 1, 2], _p);
        let m = Polynomial::new(vec![1], _p);
        //let m = &m % &test_t;
        let ct = RLWE::encrypt(&keys.public, &m);
        let pt = keys.decrypt(&ct);
        let mut pt = &pt % &test_t;
        //let mut pt = pt.mod_t();
        pt.rm_trailing_zeroes();
        assert_eq!(m, pt);
    }
    #[test]
    fn gbfv_encrypt_decrypt() {
        let parameters: Parameters<i64> = Parameters {
            degree: 4,
            q: 7681,
            t: vec![4, 3, 1],
            t_relin: 32,
            p: 1,
            log_range: 7681_i64.ilog(32) as usize + 1, //println!("{}", (self.q().checked_ilog(self.t_relin()).unwrap()) + 1);,
            root: 0,
        };
        let _p = &Rc::new(parameters);

        println!("make_fx : {}", make_fx(_p));
        let mut keys = RLWE::new(_p);
        let test_t = Polynomial::new(_p.t.clone(), _p);
        //let m = Polynomial::new(vec![2, 2, 3, 1, 2], _p);
        let m = Polynomial::new(vec![10], _p);
        let flat_m = RLWE::flatten_m(m.clone());
        //let m = &m % &test_t;
        let ct = RLWE::encrypt(&keys.public, &m);
        let flat_ct = RLWE::encrypt(&keys.public, &flat_m);
        let mut pt = &keys.decrypt(&ct) % &test_t;
        let mut flat_pt = &keys.decrypt(&flat_ct) % &test_t;
        //let mut pt = pt.mod_t();
        flat_pt.rm_trailing_zeroes();
        pt.rm_trailing_zeroes();
        println!("FLATTEN: {}", flat_m);
        assert_eq!(m, pt);
        //assert_eq!(1, 0);
        assert_eq!(pt, flat_pt);
    }
}

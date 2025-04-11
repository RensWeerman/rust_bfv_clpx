use num_bigint::{BigInt, RandBigInt};

use crate::{fractions::Fraction, rlwe::Params};
use num_traits::{one, zero, Euclid, FromPrimitive, One, Zero};
use rand_distr::{Distribution, Normal, Uniform};
use std::{
    cell::RefCell,
    cmp::max,
    fmt::Display,
    ops::{Add, Div, Mul, Neg, Rem, Sub},
    rc::Rc,
};
#[derive(PartialEq, Debug, Clone)]
pub struct Polynomial<Value> {
    pub val: Vec<Value>,
    pub parameters: Rc<RefCell<Parameters<Value>>>,
}

use crate::ntt::Ntt;

#[derive(PartialEq, Debug, Clone)]
pub struct Parameters<Val> {
    pub degree: i64,
    pub q: Val,
    pub t: Vec<Val>,
    pub t_relin: Val,
    pub p: Val,
    pub log_range: usize,
}

pub trait Value:
    Clone
    + Add<Output = Self>
    + Mul<Output = Self>
    + Sub<Output = Self>
    + PartialOrd
    + Euclid
    + Sampling<Value = Self>
    + Display
    + FromPrimitive
    + Zero
    + One
    + TryInto<usize>
    + PartialEq
    + Ord
{
}

impl Value for i64 {}

impl Value for BigInt {}

pub trait Sampling {
    type Value;
    fn normal_sample(parameters: Self::Value, size: usize) -> Vec<Self::Value>;
    fn uniform_sample(parameters: Self::Value, size: usize) -> Vec<Self::Value>;
}

impl Sampling for i64 {
    type Value = Self;
    fn normal_sample(q: i64, degree: usize) -> Vec<i64> {
        //let mut rng = rand::rngs::StdRng::seed_from_u64(100);
        let normal = Normal::new(0.0, 1.0).unwrap();
        let mut sample = Vec::new();
        for _ in 0..degree {
            let x: f32 = normal.sample(&mut rand::thread_rng());
            let x = x.trunc();
            sample.push(x as i64);
        }
        //println!("SAMPLE MEAN: {}", val / 100.0);
        sample
    }
    fn uniform_sample(q: i64, degree: usize) -> Vec<i64> {
        //let mut rng = rand::rngs::StdRng::seed_from_u64(100);
        let uniform = Uniform::try_from(-q..q + 1).unwrap();
        let mut sample = Vec::new();
        for _ in 0..degree {
            let x: i64 = uniform.sample(&mut rand::thread_rng());
            //print!("{}, ", x);
            sample.push(x as i64);
        }
        //println!("SAMPLE MEAN: {}", val / 100.0);
        sample
    }
}
impl Sampling for f64 {
    type Value = Self;
    fn normal_sample(q: f64, degree: usize) -> Vec<f64> {
        //let mut rng = rand::rngs::StdRng::seed_from_u64(100);
        let normal = Normal::new(0.0, 1.0).unwrap();
        let mut sample = Vec::new();
        for _ in 0..degree {
            let x: f32 = normal.sample(&mut rand::thread_rng());
            let x = x.trunc();
            sample.push(x as f64);
        }
        //println!("SAMPLE MEAN: {}", val / 100.0);
        sample
    }
    fn uniform_sample(q: f64, degree: usize) -> Vec<f64> {
        //let mut rng = rand::rngs::StdRng::seed_from_u64(100);
        let uniform = Uniform::try_from(-q..q + 1.0).unwrap();
        let mut sample = Vec::new();
        for _ in 0..degree {
            let x: f64 = uniform.sample(&mut rand::thread_rng());
            //print!("{}, ", x);
            sample.push(x as f64);
        }
        //println!("SAMPLE MEAN: {}", val / 100.0);
        sample
    }
}

impl Sampling for BigInt {
    type Value = Self;
    fn normal_sample(q: BigInt, size: usize) -> Vec<BigInt> {
        //let mut rng = rand::rngs::StdRng::seed_from_u64(100);
        let normal = Normal::new(0.0, 1.0).unwrap();
        let mut sample = Vec::new();
        for _ in 0..size {
            let x: f32 = normal.sample(&mut rand::thread_rng());
            let x = x.trunc();
            sample.push(one());
        }
        //println!("SAMPLE MEAN: {}", val / 100.0);
        sample
    }
    fn uniform_sample(q: BigInt, size: usize) -> Vec<BigInt> {
        //let mut rng = rand::rngs::StdRng::seed_from_u64(100);
        //let uniform = Uniform::try_from(-q..q + 1).unwrap();
        let mut sample = Vec::new();
        let mut rng = rand::thread_rng();
        for _ in 0..size {
            //let x: i64 = uniform.sample(&mut rand::thread_rng());
            let x = rng.gen_bigint_range(&zero(), &q);
            //print!("{}, ", x);
            sample.push(x);
        }
        //println!("SAMPLE MEAN: {}", val / 100.0);
        sample
    }
}
impl<Val> Display for Polynomial<Val>
where
    Val: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        //write!(f, "[");
        for i in 0..self.val.len() {
            write!(f, "{}", self.val.get(i).unwrap());
            if i < self.val.len() - 1 {
                write!(f, " ");
            }
        }
        write!(f, "")
        //write!(f, "]")
    }
}

impl<Val> Polynomial<Val>
where
    Val: Value,
{
    pub fn new(val: Vec<Val>, parameters: &Params<Val>) -> Self {
        Self {
            val,
            parameters: Rc::clone(parameters),
        }
    }
    pub fn uniform_sample(parameters: &Params<Val>) -> Self {
        Polynomial {
            val: Val::uniform_sample(
                parameters.borrow().q.clone(),
                parameters.borrow().degree as usize,
            ),
            parameters: Rc::clone(parameters),
        }
    }
    pub fn normal_sample(parameters: &Params<Val>) -> Self {
        Polynomial {
            val: Val::normal_sample(
                parameters.borrow().q.clone(),
                parameters.borrow().degree as usize,
            ),
            parameters: Rc::clone(parameters),
        }
    }
    pub fn t_relin(&self) -> Val {
        return self.parameters.borrow().t_relin.clone();
    }
    pub fn q(&self) -> Val {
        return self.parameters.borrow().q.clone();
    }
    pub fn p(&self) -> Val {
        return self.parameters.borrow().p.clone();
    }
    pub fn t(&self) -> Val {
        return self.parameters.borrow().t.clone().get(0).unwrap().clone();
    }
    pub fn t_full(&self) -> Vec<Val> {
        return self.parameters.borrow().t.clone();
    }
    pub fn log_range(&self) -> usize {
        return self.parameters.borrow().log_range.clone();
    }
    pub fn parameters(&self) -> &Params<Val> {
        return &self.parameters;
    }
    pub fn delta(&self) -> Polynomial<Val> {
        let c = Polynomial::new(self.t_full(), self.parameters());
        let q = c.q();
        let temp = Polynomial::new(
            c.val.clone().into_iter().map(|i| i * q.clone()).collect(),
            c.parameters(),
        );
        let next_temp = &temp / &Polynomial::new(temp.t_full(), c.parameters());
        next_temp
    }
    // pub fn other_decompose(self) {
    //     let c_2 = self.clone();
    //     let mut testing = Polynomial::new(vec![], self.parameters());
    //     let mut c_2 = &Polynomial::new(vec![], self.parameters()) + &c_2;
    //     //print_vec(&c_2.val);
    //     //print_vec(&(&c_2 % &*FX).val);
    //     for i in 0..(self.q().checked_ilog(self.t_relin()).unwrap()) + 1 {
    //         println!("{}", i);
    //         let c_2r = &c_2 % &Polynomial::new(vec![self.t_relin()], self.parameters());
    //         let c_2r = c_2r.normal_mod(self.t_relin());
    //         testing = &testing
    //             + &(&c_2r
    //                 * &Polynomial::new(
    //                     vec![montgomery::more_secure_mont_ladder(
    //                         self.t_relin(),
    //                         i as i64,
    //                         self.q(),
    //                     )],
    //                     self.parameters(),
    //                 )); //self.t_relin.pow(i)]));
    //                     //print_vec(&c_2r.val);
    //         c_2 = &c_2 - &c_2r; //* &Polynomial::new(vec![t.pow(i)]));
    //         c_2 = &c_2 / &Polynomial::new(vec![self.t_relin()], self.parameters());
    //     }
    //     //print_vec(&testing.val);
    // }

    //NOT WORKING YET
    /*pub fn decompose(self) -> Vec<Polynomial<Value>> {
        let base = self.t_relin;
        let mut mut_poly = self.clone();
        //mut_poly = &mut_poly % &*FX;

        // Iterate i: from highest to lowest level, starting with l
        let out_polys: Vec<Polynomial<Value>> = (0..(Q.checked_ilog(T_RELIN).unwrap()) + 1)
            .rev()
            .map(|i| {
                // T^i, which is the multiplier for that level i
                //let base_i = base.pow(i as u32);
                let base_i = m_ladder(base, i as i64, self.q);

                // Iterate j: through the coefficients in poly, to decompose for level i
                let dec_val_i: Vec<i32> = mut_poly
                    .val
                    .iter_mut()
                    .map(|val_j| {
                        // Calculate how many times T^i divides the coefficient, to get decomposition
                        let fl_div = *val_j / base_i;
                        let int_div = fl_div; // = if fl_div > zero() { fl_div } else { fl_div } as i32;
                                              // Update the coefficient by subtracting T^i * the decomposed value
                        *val_j = *val_j - base_i * int_div;
                        // Return the decomposed value for that coefficient for level i
                        int_div
                    })
                    .collect();
                Polynomial::new(dec_val_i)
            })
            .collect();
        // We can't reverse within the original expression because the two "rev" calls cancel each other out
        // and we get the wrong decomposition answer (decomposing starting from the smallest levels).
        out_polys.into_iter().rev().collect()
    }*/
    pub fn degree(&self) -> usize {
        let no_zeroes = rm_trailing_zeroes(self.val.clone());
        return degree(&no_zeroes);
    }
    pub fn get_mod(&self) -> Polynomial<Val> {
        let mut x = self % &make_fx(self.parameters());
        x.val
            .resize(self.parameters().borrow().degree as usize, Val::zero());
        x
    }
    pub fn set_back(&self) -> Polynomial<Val> {
        let mut c = Vec::new();
        for i in 0..self.parameters().borrow().degree as usize {
            let mut val = mod_coeff((self.val.get(i).unwrap_or(&Val::zero())).clone(), self.q());
            if val > self.q() / (Val::one() + Val::one()) {
                val = val - self.q();
            }
            c.push(val);
        }
        //c = rm_trailing_zeroes(c);
        Polynomial::new(c, self.parameters())
    }
    pub fn mod_q(&self) -> Polynomial<Val> {
        let mut c = Vec::new();
        for i in 0..self.parameters().borrow().degree as usize {
            c.push(mod_coeff(
                (self.val.get(i).unwrap_or(&Val::zero())).clone(),
                self.q(),
            ));
        }
        //c = rm_trailing_zeroes(c);
        Polynomial::new(c, self.parameters())
    }
    pub fn mod_t(&self) -> Polynomial<Val> {
        let mut c = Vec::new();
        let t = &Polynomial::new(self.t_full(), self.parameters());
        if t.degree() == 1 {
            for i in 0..self.parameters().borrow().degree as usize {
                c.push(mod_coeff(
                    (self.val.get(i).unwrap_or(&Val::zero())).clone(),
                    self.t(),
                ));
            }
            return Polynomial::new(c, self.parameters());
        }
        let mut temp = &self.clone() % &Polynomial::new(self.t_full(), self.parameters());
        temp.val
            .truncate(self.parameters().borrow().degree as usize);
        println!("WHAT? {}", self);
        //c = rm_trailing_zeroes(c);
        Polynomial::new(temp.val, self.parameters())
    }
    pub fn mod_t_relin(&self) -> Polynomial<Val> {
        let mut c = Vec::new();
        for i in 0..self.parameters().borrow().degree as usize {
            c.push(mod_coeff(
                (self.val.get(i).unwrap_or(&Val::zero())).clone(),
                self.t_relin(),
            ));
        }
        //c = rm_trailing_zeroes(c);
        Polynomial::new(c, self.parameters())
    }
    pub fn new_mul(&self, rhs: &Self) -> Self {
        let degree = self.parameters().borrow().degree;
        let q = self.parameters().borrow().q.clone();
        let psi = &Ntt::second_primitive_root(&q, degree);
        let res = Ntt::ntt_mult(self.clone(), rhs.clone(), &q, psi);
        res
        //res.set_back()
    }
    pub fn old_mul(&self, rhs: &Self) -> Self {
        let mut c = Vec::new();
        let zero = Val::zero();
        for i in 0..(self.degree() + rhs.degree() + 1) {
            let mut s = Val::zero();
            for f in 0..(i + 1) {
                let _a = self.val.get(f).unwrap_or(&zero);
                let _b = rhs.val.get(i - f).unwrap_or(&zero);
                s = s + (_a.clone() * _b.clone());
            }
            c.push(s);
        }
        //c = rm_trailing_zeroes(c);
        Polynomial {
            val: c,
            parameters: Rc::clone(self.parameters()),
        }
    }
    pub fn rm_trailing_zeroes(&mut self) {
        //return c;
        let first = self.val.iter().rev().position(|x| *x != Val::zero());
        match first {
            Some(v) => {
                let mut _c = self
                    .val
                    .clone()
                    .into_iter()
                    .rev()
                    .collect::<Vec<_>>()
                    .split_off(v);
                self.val = _c.into_iter().rev().collect::<Vec<_>>();
            }
            None => self.val = Vec::new(),
            //None -> return,
            //Some(v) -> v,
        }
    }
    pub fn extended_gcd(a: &Self, b: &Self) -> Self {
        let (mut old_r, mut r) = (a.clone(), b.clone());
        let (mut old_s, mut s) = (
            Polynomial::new(vec![one()], a.parameters()),
            Polynomial::new(vec![zero()], a.parameters()),
        );
        let (mut old_t, mut t) = (
            Polynomial::new(vec![zero()], a.parameters()),
            Polynomial::new(vec![one()], a.parameters()),
        );

        while r.degree() > 0 {
            let quotient = &old_r / &r;
            (old_r, r) = (r.clone(), &old_r % &r);
            (old_s, s) = (s.clone(), &old_s - &(&quotient * &s));
            (old_t, t) = (t.clone(), &old_t - &(&quotient * &t));
            println!("????: {}", old_s);
            println!("WHAT IS HIER GAANDE: {}", r);
        }
        old_s = &old_s / &old_r;
        old_t = &old_t / &old_r;
        println!("Bezout coefficients: {}, {}", old_s, old_t);
        println!("GCD: {}", old_r);
        println!("quotients by the gcd: {}, {}", t, s);
        return old_s;
    }
}
impl<Val> Add for &Polynomial<Val>
where
    Val: Value,
{
    type Output = Polynomial<Val>;

    fn add(self, rhs: Self) -> Self::Output {
        let mut c = Vec::new();
        let zero = Val::zero();
        for i in 0..max(self.degree(), rhs.degree()) {
            let a = self.val.get(i).unwrap_or(&zero);
            let b = rhs.val.get(i).unwrap_or(&zero);
            c.push(a.clone() + b.clone());
        }
        //c = rm_trailing_zeroes(c);
        Polynomial::new(c, self.parameters())
    }
}
impl<Val> Mul for &Polynomial<Val>
where
    Val: Value,
{
    type Output = Polynomial<Val>;

    fn mul(self, rhs: Self) -> Self::Output {
        self.old_mul(rhs)
        // let degree = self.parameters().borrow().degree;
        // let q = self.parameters().borrow().q.clone();
        // let psi = &Ntt::second_primitive_root(&q, degree);
        // let res = Ntt::ntt_mult(self.clone(), rhs.clone(), &q, psi);
        // res.set_back()
    }
}
impl<Val> Neg for &Polynomial<Val>
where
    Val: Value,
{
    type Output = Polynomial<Val>;
    fn neg(self) -> Self::Output {
        return &Polynomial::new(vec![], self.parameters()) - self;
    }
}

impl<Val> Sub for &Polynomial<Val>
where
    Val: Value,
{
    type Output = Polynomial<Val>;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut c = Vec::new();
        let zero = Val::zero();
        for i in 0..max(self.degree(), rhs.degree()) {
            let a = self.val.get(i).unwrap_or(&zero);
            let b = rhs.val.get(i).unwrap_or(&zero);
            c.push(a.clone() - b.clone());
        }
        //c = rm_trailing_zeroes(c);
        Polynomial {
            val: c,
            parameters: Rc::clone(self.parameters()),
        }
    }
}
impl<Val> Rem for &Polynomial<Val>
where
    Val: Value,
{
    type Output = Polynomial<Val>;

    fn rem(self, rhs: Self) -> Self::Output {
        let mut quotient = Polynomial::new(vec![], self.parameters());
        let mut remainder = Polynomial::new(self.val.clone(), self.parameters());
        //println!("WHAT {}, {}", remainder.degree(), rhs.degree());
        /*println!(
            "Wattisdis {}",
            ((&remainder.val).into_iter().any(|x| *x != zero()))
        );*/
        while (remainder.degree() >= rhs.degree())
            && ((&remainder.val).into_iter().any(|x| *x != Val::zero()))
        {
            /*println!("EEEE");
            print_vec(&remainder.val);
            print_vec(&rhs.val);*/
            /*
            We need to get the leading terms of both the remainder, and the divisor. Then we
            divide them, to get the number for our quotient. We do this by getting both the
            coefficients and the degree.
            */
            //print_vec(&remainder);
            let rem_c = remainder.val.get(remainder.degree() - 1).unwrap();
            /*println!(
                "Rem_c should not be zero, but it is: {}. Degree is: {}. below is entire poly:",
                rem_c,
                remainder.degree() - 1,
            );
            print_vec(&remainder.val);*/
            let div_c = rhs.val.get(rhs.degree() - 1).unwrap();
            let t = single_term(
                self.parameters(),
                rem_c.clone() / div_c.clone(),
                remainder.degree() - rhs.degree(),
            );
            //println!("{}, {}", rem_c, div_c);
            if rem_c.clone() / div_c.clone() == Val::zero() {
                //println!("??");
                return remainder;
            }
            quotient = &quotient + &t;
            //print!("%Quotient: ");
            //print_vec(&quotient.val);
            quotient = &quotient % &make_fx(self.parameters());
            /*print!("%remainder and other stuff: ");
            print_vec(&remainder.val);
            print_vec(&(&t * &rhs).val);*/
            remainder = &remainder - &(&t.old_mul(&rhs));
            /*print!("%remainder: ");
            print_vec(&remainder.val);*/
            remainder = &remainder % &make_fx(self.parameters());
        }
        remainder
    }
}
impl<Val> Div for &Polynomial<Val>
where
    Val: Value,
{
    type Output = Polynomial<Val>;

    fn div(self, rhs: Self) -> Self::Output {
        let mut quotient: Polynomial<Val> = Polynomial::new(vec![], self.parameters());
        let mut remainder = Polynomial::new(self.val.clone(), self.parameters());
        while (remainder.degree() >= rhs.degree())
            && ((&remainder.val).into_iter().any(|x| !x.is_zero()))
        {
            println!("WHAT {}, {}", remainder.degree(), rhs.degree());

            /*let t = Vec<i32> = Vec::new();
            for i in 0..degree(remainder) {
                t.push(0);
            }
            let t =*/
            /*
            We need to get the leading terms of both the remainder, and the divisor. Then we
            divide them, to get the number for our quotient. We do this by getting both the
            coefficients and the degree.
            */
            //print_vec(&remainder);
            let rem_c = remainder.val.get(remainder.degree() - 1).unwrap();
            let div_c = rhs.val.get(rhs.degree() - 1).unwrap();
            let t = single_term(
                self.parameters(),
                rem_c.clone() / div_c.clone(),
                remainder.degree() - rhs.degree(),
            );
            if rem_c.clone() / div_c.clone() == Val::zero() {
                return &quotient % &make_fx(self.parameters());
            }
            //quotient = add_polyn(&quotient, &t, q);
            quotient = &quotient + &t;
            //quotient = mod_polyn(quotient, fx, q);
            //print!("/Quotient: ");
            //print_vec(&quotient.val);
            quotient = &quotient % &make_fx(self.parameters());
            //remainder = sub_polyn(&remainder, &mult_polyn(&t, &rhs.val, q), q);
            /*remainder = &remainder - &(&t * &rhs);*/
            remainder = &remainder - &(&t.old_mul(&rhs));
            //remainder = mod_polyn(remainder, fx, q);
            //print!("/remainder: ");
            //print_vec(&remainder.val);
            remainder = &remainder % &make_fx(self.parameters());
        }
        quotient
    }
}

pub fn degree<Val: Value>(a: &Vec<Val>) -> usize {
    return a.len();
}
pub fn make_fx<Val: Value>(parameters: &Params<Val>) -> Polynomial<Val> {
    let mut val = vec![Val::one()];
    for _ in 0..(parameters.borrow().degree - 1) {
        val.push(Val::zero());
    }
    val.push(Val::one());
    Polynomial {
        val,
        parameters: Rc::clone(parameters),
    }
}

pub fn single_term<Val: Value>(
    parameters: &Params<Val>,
    coefficient: Val,
    degree: usize,
) -> Polynomial<Val> {
    let mut val = Vec::new();
    let mut degree = degree as i32;
    while degree > 0 {
        degree -= 1;
        val.push(Val::zero());
    }
    val.push(coefficient);
    //val = rm_trailing_zeroes(val);
    Polynomial {
        val,
        parameters: Rc::clone(parameters),
    }
    //rm_trailing_zeroes(val)
}

// pub fn m_ladder(a: Polynomial, b: i64) -> Polynomial {
//     let mut check = false;
//     let mut x_2 = &a * &a;
//     let mut x_1 = a;
//     for x in (0..64).rev().map(|n| (b >> n) & 1) {
//         if !check {
//             // print!(":");
//         }
//         //println!("{}", x);
//         if check {
//             if (x % 2) == 0 {
//                 //println!("{}, {}", x_1, x_2);
//                 x_2 = &x_2 * &x_1;
//                 x_1 = &x_1 * &x_1;
//             } else {
//                 //println!("{}, {}", x_1, x_2);
//                 x_1 = &x_1 * &x_2;
//                 x_2 = &x_2 * &x_2;
//             }
//         }
//         if (x % 2) == 1 && !check {
//             check = true;
//         }
//         x_1 = &x_1 % &*FX;
//         x_2 = &x_2 % &*FX;
//     }
//     return x_1;
// }

pub fn rm_trailing_zeroes<Val: Value>(c: Vec<Val>) -> Vec<Val> {
    //return c;
    let first = c.iter().rev().position(|x| !x.is_zero());
    match first {
        Some(v) => {
            let mut _c = c.into_iter().rev().collect::<Vec<_>>().split_off(v);
            let __c = _c.into_iter().rev().collect::<Vec<_>>();
            __c
        }
        None => Vec::new(),
        //None -> return,
        //Some(v) -> v,
    }
}

fn mod_coeff<Val: Value>(coeff: Val, q: Val) -> Val {
    //(coeff % q.clone() + q.clone()) % q.clone()
    return coeff.rem_euclid(&q);
}

// fn rem_other(this: Polynomial) -> Polynomial {
//     let coeff_mod = Q;
//     let degree = DEGREE as usize;
//     let mut out_val = vec![0; degree];

//     // Take the polynomial mod (X^N + 1).
//     // 1. After a multiplication by X^{2N}, the polynomial is unchanged mod (X^N + 1).
//     //    Therefore, we can take the degree % 2N.
//     // 2. If degree % 2N > N, the coefficients should be negated and added to the degree % N.
//     // 3. If degree % 2N <= N, the coefficients should be added to the degree % 2N.
//     for (i, coeff) in this.val.iter().enumerate() {
//         // $ X^i == X^{i + j * 2N} mod (X^N + 1) for all j $
//         // So we can take the coeff degree mod 2N.
//         let reduced_i = i % (2 * degree);
//         if reduced_i >= degree {
//             out_val[reduced_i % degree] -= coeff;
//         } else {
//             out_val[reduced_i] += coeff;
//         }
//     }

//     // Take each coefficient % coeff_mod
//     for coeff in out_val.iter_mut() {
//         *coeff = mod_coeff(*coeff, coeff_mod)
//     }
//     Polynomial::new(out_val)
// }

#[cfg(test)]
mod tests {
    use crate::{fractions::Fraction, rlwe::RLWE};

    use super::*;

    #[test]
    fn first_step() {
        /*pub const DEGREE: i32 = 4;
        pub const Q: i32 = 7681; //1291; //9073;
        pub const T: i32 = 6;
        pub const T_RELIN: i32 = 32;
        pub const P: i32 = 4;*/
        let parameters = Parameters {
            degree: 4,
            q: 7681,
            t: vec![6],
            t_relin: 32,
            p: 4,
            log_range: 0,
        };
        let _p = &Rc::new(RefCell::new(parameters));
        let p1 = Polynomial::new(vec![1, 2, 3, 4], _p);
        let p2 = Polynomial::new(vec![5, 6, 7, 8], _p);
        let right = &p1 * &p2;
        let left = Polynomial::new(vec![5, 16, 34, 60, 61, 52, 32], _p);
        assert_eq!(left.get_mod(), right.get_mod());
    }
    #[test]
    fn trying_poly_mult() {
        let parameters = Parameters {
            degree: 4,
            q: 7681,
            t: vec![6],
            t_relin: 32,
            p: 4,
            log_range: 0,
        };
        let _p = &Rc::new(RefCell::new(parameters));

        let p1 = Polynomial::new(vec![1, 2, 3, 4], _p);
        let p2 = Polynomial::new(vec![5, 6, 7, 8], _p);
        //let left = Polynomial::new(vec![-56, -36, 2, 60], _p);
        let left = &p1 * &p2;

        //let left = left.normal_mod(_p.borrow().q);
        let left = &left.get_mod();
        //let right = &(&p1 * &p2) % &*FX;

        let right = Polynomial::new(vec![5, 16, 34, 60, 61, 52, 32], _p);
        let right = &right.get_mod();
        assert_eq!(left, right);
    }
    #[test]
    fn test_rem() {
        let parameters = Parameters {
            degree: 4,
            q: 7681,
            t: vec![6],
            t_relin: 32,
            p: 4,
            log_range: 0,
        };
        let _p = &Rc::new(RefCell::new(parameters));
        let right = Polynomial::new(vec![5, 16, 34, 60, 61, 52, 32], _p);
        assert_eq!(
            /*&right % &make_fx(_p)*/ right.get_mod().mod_q(),
            Polynomial::new(vec![7625, 7645, 2, 60], _p)
        );
    }
    //#[test]
    fn comparing_mults_big_ints() {
        let parameters: Parameters<BigInt> = Parameters {
            degree: 4,
            q: BigInt::try_from(7681).unwrap(),
            t: vec![BigInt::try_from(16).unwrap()],
            t_relin: BigInt::try_from(32).unwrap(),
            p: BigInt::try_from(4).unwrap(),
            log_range: 0,
        };
        let _p = &Rc::new(RefCell::new(parameters));
        // let p1 = Polynomial::new(
        //     vec![
        //         BigInt::try_from(1).unwrap(),
        //         BigInt::try_from(2).unwrap(),
        //         BigInt::try_from(3).unwrap(),
        //         BigInt::try_from(4).unwrap(),
        //     ],
        //     _p,
        // );
        // let p2 = Polynomial::new(
        //     vec![
        //         BigInt::try_from(5).unwrap(),
        //         BigInt::try_from(6).unwrap(),
        //         BigInt::try_from(7).unwrap(),
        //         BigInt::try_from(8).unwrap(),
        //     ],
        //     _p,
        // );
        let p1 = Polynomial::uniform_sample(_p);
        let p2 = Polynomial::uniform_sample(_p);
        let ls = p1.old_mul(&p2);
        let ls = RLWE::delta_polyn(ls.clone()).mod_q();
        println!("LS: {}", ls);
        let degree = p1.parameters().borrow().degree;
        let q = &p1.parameters().borrow().q.clone();
        let big_q = &BigInt::try_from((q - 1) * 17 + 1).unwrap();
        let psi = &Ntt::second_primitive_root(big_q, degree);
        let p1_ntt = Ntt::iter_cooley_tukey(p1, psi, big_q);
        let p2_ntt = Ntt::iter_cooley_tukey(p2, psi, big_q);
        println!("BEFORE: {}, {}", p1_ntt, p2_ntt);
        let mut p_ntt = Ntt::element_wise_mult(p1_ntt, p2_ntt, big_q);
        println!("RS: {}", p_ntt);
        let rs = Ntt::iter_gentleman_sande(p_ntt, psi, big_q);
        let rs = RLWE::delta_polyn(rs).mod_q();

        println!("AFTER RS: {}", rs);

        assert_eq!(ls.get_mod().mod_q(), rs.mod_q());
    }
    //#[test]
    fn comparing_mults() {
        let parameters = Parameters {
            degree: 4,
            q: 7681,
            t: vec![6],
            t_relin: 32,
            p: 4,
            log_range: 0,
        };
        let _p = &Rc::new(RefCell::new(parameters));
        let p1 = Polynomial::new(vec![1, 2, 3, 4], _p);
        let p2 = Polynomial::new(vec![5, 6, 7, 8], _p);
        let ls = p1.old_mul(&p2);
        let ls = RLWE::delta_polyn(ls).mod_q();
        println!("LS: {}", ls);
        let degree = p1.parameters().borrow().degree;
        let q = &p1.parameters().borrow().q.clone();
        let psi = &Ntt::second_primitive_root(q, degree);
        let p1_ntt = Ntt::iter_cooley_tukey(p1, psi, q);
        let p2_ntt = Ntt::iter_cooley_tukey(p2, psi, q);
        println!("BEFORE: {}, {}", p1_ntt, p2_ntt);
        let mut p_ntt = Ntt::element_wise_mult(p1_ntt, p2_ntt, q);
        //let mut delta = Polynomial::new(vec![p_ntt.q() / p_ntt.t(), 0, 0, 0], _p);
        //delta = Ntt::iter_cooley_tukey(delta, psi, q);
        //let p_ntt = Ntt::element_wise_div(p_ntt, delta, q);
        println!("RS: {}", p_ntt);
        let p_ntt = RLWE::delta_polyn(p_ntt).mod_q();
        println!("RS: {}", p_ntt);
        let rs = Ntt::iter_gentleman_sande(p_ntt, psi, q);

        assert_eq!(ls.get_mod().mod_q(), rs.get_mod().mod_q());
    }
    //#[test]
    fn for_extended_euclidian() {
        let parameters = Parameters {
            degree: 8,
            q: Fraction::from_i64(7681).unwrap(),
            t: vec![Fraction::from_i64(20).unwrap()],
            t_relin: Fraction::from_i64(32).unwrap(),
            p: Fraction::from_i64(4).unwrap(),
            log_range: 0,
        };
        let _p = &Rc::new(RefCell::new(parameters));
        let x = Fraction::new(1, 1);
        let mut a = Polynomial::new(
            vec![
                Fraction::new(-1, 1),
                Fraction::new(0, 1),
                Fraction::new(-1, 1),
                Fraction::new(0, 1),
                Fraction::new(2, 1),
                Fraction::new(1, 1),
            ],
            _p,
        );
        let mut b = Polynomial::new(
            vec![
                Fraction::new(1, 1),
                Fraction::new(0, 1),
                Fraction::new(0, 1),
                Fraction::new(0, 1),
                Fraction::new(1, 1),
            ],
            _p,
        );
        Polynomial::extended_gcd(&a, &b);
        println!("-------");

        let mut div = &a / &b;
        let mut rem = &a % &b;
        println!("A: {}, B: {}", a, b);
        println!("Div: {}, Rem: {}", div, rem);

        while b.degree() > 1 {
            a = b;
            b = rem;
            println!("A: {}, B: {}", a, b);
            div = &a / &b;
            rem = &a % &b;
            println!("Div: {}, Rem: {}", div, rem);
        }
    }
    #[test]
    fn trying_clpx() {
        let parameters = Parameters {
            degree: 8,
            q: Fraction::new(7681, 1),
            t: vec![Fraction::new(8, 1), Fraction::new(1, 1)],
            t_relin: Fraction::new(32, 1),
            p: Fraction::new(4, 1),
            log_range: 0,
        };
        let i_parameters = Parameters {
            degree: 8,
            q: 7681,
            t: vec![8, 1],
            t_relin: 32,
            p: 4,
            log_range: 0,
        };
        let _p = &Rc::new(RefCell::new(parameters));
        let _ip = &Rc::new(RefCell::new(i_parameters));
        let x = Fraction::new(1, 1);
        let mut a = Polynomial::new(vec![Fraction::new(8, 1), Fraction::new(1, 1)], _p);
        let mut b = Polynomial::new(
            vec![
                Fraction::new(1, 1),
                Fraction::new(0, 1),
                Fraction::new(0, 1),
                Fraction::new(0, 1),
                Fraction::new(1, 1),
            ],
            _p,
        );
        let old_s = Polynomial::extended_gcd(&a, &b);
        let mut temp = vec![];
        for x in old_s.val {
            temp.push(Fraction::new(7681, 1) * x.reduce());
        }
        let fin = Polynomial::new(temp, _p);
        println!("FINAL: {}", fin);
        let ex = Polynomial::new(vec![1, 2, 3, 4], _ip);
        let frac_ex = Fraction::frac_poly(ex, _p);
        let res = Fraction::round_multiply(&frac_ex, &fin, _ip);
        println!("RESULT: {}", res);
        let hmm = RLWE::delta_polyn(res);
        println!("HMMMM: {}", hmm);
        assert_eq!(1, 0);
    }
    // #[test]
    // fn comparing_to_other_mult() {
    //     let p1 = Polynomial::new(vec![1, 2, 3, 4]);
    //     let p2 = Polynomial::new(vec![5, 6, 7, 8]);
    //     assert_eq!(&p1 * &p2, mul_other(p1, p2));

    //     let p3 = Polynomial::new(vec![1, 6, 2, 9]);
    //     let p4 = Polynomial::new(vec![1, 13743981, 3, 12005, 3, 10000, 15]);
    //     assert_eq!((&p3 * &p4).normal_mod(Q), mul_other(p3, p4).normal_mod(Q));
    // }

    // #[test]
    // fn comparing_rems() {
    //     let p4 = Polynomial::new(vec![1, 10000, 3, 5, 3, 0, 10]);
    //     let p4 = Polynomial::uniform_sample();
    //     assert_eq!((&p4 % &*FX).normal_mod(Q), rem_other(p4))
    // }
}

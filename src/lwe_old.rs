use std::cmp::max;

use rand_distr::{Distribution, Normal, Uniform};
use std::fmt::Write as _;

pub struct LWE {
    secret: Vec<i32>,
    pub public: (Vec<i32>, Vec<i32>),
    degree: i32,
    q: i32,
    fx: Vec<i32>,
}
impl LWE {
    pub fn _new(
        secret: Vec<i32>,
        public: (Vec<i32>, Vec<i32>),
        degree: i32,
        q: i32,
        fx: Vec<i32>,
    ) -> Self {
        LWE {
            secret,
            public,
            degree,
            q,
            fx,
        }
    }
    pub fn new(degree: i32, q: i32, fx: Vec<i32>) -> Self {
        let secret = normal_sample(degree + 1);
        let e = /*vec![];*/ normal_sample(degree + 1);
        let a = uniform_sample(degree + 1, q); //mod_polyn(uniform_sample(degree + 1, q), &fx, q);
        let pk: (Vec<i32>, Vec<i32>) = (
            //mod_polyn(
            sub_polyn(
                &vec![],
                &add_polyn(&mult_polyn(&a, &secret, q), /*&vec![]*/ &e, q),
                q,
            ),
            //&fx,
            //q,
            //),
            add_polyn(&a, &vec![0], q),
        );
        /*print!("Secret RWLE: ");
        print_vec(&secret);
        print!("Noise: ");
        print_vec(&e);
        print!("Public RWLE: ");
        print_vec(&pk.0);
        print_vec(&pk.1);*/
        LWE {
            secret,
            public: pk,
            degree,
            q,
            fx,
        }
    }
    pub fn decrypt_encoding(&mut self, ct: &str, size: i32) -> String {
        let my_m = base64::decode(ct).unwrap();
        let cct: Vec<(Vec<i32>, Vec<i32>)> = bincode::deserialize(&my_m).unwrap();
        //let ct = lwe::encrypt(&a.public, &m, degree, q, &fx, size);
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
    }
    pub fn decrypt(&mut self, ct: &(Vec<i32>, Vec<i32>), size: i32) -> Vec<i32> {
        let _pt = add_polyn(&ct.0, &mult_polyn(&ct.1, &self.secret, self.q), self.q);
        let __pt = mod_polyn(_pt, &self.fx, self.q);
        let pt = add_polyn(&__pt, &vec![], self.q);
        //let pt = mult_polyn(&vec![size], &pt, self.q);
        return rm_trailing_zeroes(
            pt.into_iter()
                .map(|i| (((i * size) + (self.q / 2)) / self.q).rem_euclid(size))
                .collect(),
        );
    }

    pub fn sample(&mut self, degree: i32) -> Vec<i32> {
        //let mut rng = rand::rngs::StdRng::seed_from_u64(100);
        let normal = Normal::new(0.0, 1.0).unwrap();
        let mut sample = Vec::new();
        for _ in 0..degree {
            let x: f32 = normal.sample(&mut rand::thread_rng());
            let x = x.trunc();
            sample.push(x as i32);
        }
        sample
    }
}
pub fn degree(a: &Vec<i32>) -> usize {
    return a.len();
}
pub fn single_term(coefficient: i32, degree: usize) -> Vec<i32> {
    let mut val = Vec::new();
    let mut degree = degree as i32;
    while degree > 0 {
        degree += -1;
        val.push(0);
    }
    val.push(coefficient);
    rm_trailing_zeroes(val)
}
pub fn encrypt(
    pk: &(Vec<i32>, Vec<i32>),
    m: &Vec<i32>,
    degree: i32,
    q: i32,
    fx: &Vec<i32>,
    size: i32,
) -> (Vec<i32>, Vec<i32>) {
    let u = /*vec![1];*/ normal_sample(degree + 1);
    let e1 = normal_sample(degree + 1);
    let e2 = normal_sample(degree + 1);
    /*for i in 0..100 {
    e1 = add_polyn(&e1, &normal_sample(degree + 1), q);
    e2 = add_polyn(&e2, &normal_sample(degree + 1), q);
    }*/
    let _ct0 = add_polyn(&mult_polyn(&pk.0, &u, q), &e1, q);
    let __ct0 = add_polyn(&_ct0, &mult_polyn(&vec![q / size], m, q), q);
    let __ct1 = add_polyn(&mult_polyn(&pk.1, &u, q), &e2, q);
    //let ct0 = mod_polyn(__ct0, fx, q);
    //let ct1 = mod_polyn(__ct1, fx, q);
    (__ct0, __ct1)
}
pub fn div_polyn(dividend: Vec<i32>, divisor: Vec<i32>, fx: &Vec<i32>, q: i32) -> Vec<i32> {
    let mut quotient: Vec<i32> = Vec::new();
    let mut remainder = dividend.clone();
    while (degree(&remainder) >= degree(&divisor)) && ((&remainder).into_iter().any(|x| *x != 0)) {
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
        let rem_c = remainder.get(remainder.len() - 1).unwrap();
        let div_c = divisor.get(divisor.len() - 1).unwrap();
        let t = single_term(rem_c / div_c, remainder.len() - divisor.len());
        quotient = add_polyn(&quotient, &t, q);
        quotient = mod_polyn(quotient, fx, q);
        remainder = sub_polyn(&remainder, &mult_polyn(&t, &divisor, q), q);
        remainder = mod_polyn(remainder, fx, q);
    }
    remainder
}
pub fn sub_polyn(a: &Vec<i32>, b: &Vec<i32>, q: i32) -> Vec<i32> {
    let mut c = Vec::new();
    for i in 0..max(a.len(), b.len()) {
        let _a = a.get(i).unwrap_or(&0);
        let _b = b.get(i).unwrap_or(&0);
        c.push((_a - _b).rem_euclid(q));
    }
    rm_trailing_zeroes(c)
}
pub fn rm_trailing_zeroes(c: Vec<i32>) -> Vec<i32> {
    //return c;
    let first = c.iter().rev().position(|x| *x != 0);
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
pub fn add_polyn(a: &Vec<i32>, b: &Vec<i32>, q: i32) -> Vec<i32> {
    let mut c = Vec::new();
    for i in 0..max(a.len(), b.len()) {
        let _a = a.get(i).unwrap_or(&0);
        let _b = b.get(i).unwrap_or(&0);
        c.push((_a + _b).rem_euclid(q));
    }
    rm_trailing_zeroes(c)
}
pub fn mod_polyn(c: Vec<i32>, fx: &Vec<i32>, q: i32) -> Vec<i32> {
    div_polyn(c, fx.to_vec(), fx, q)
}
pub fn delta_polyn(c: Vec<i32>, t: i32, q: i32) -> Vec<i32> {
    return rm_trailing_zeroes(
        c.into_iter()
            .map(|i| (((i * t) + (q / 2)) / q).rem_euclid(t))
            .collect(),
    );
}
pub fn mult_polyn(a: &Vec<i32>, b: &Vec<i32>, q: i32) -> Vec<i32> {
    let mut c = Vec::new();
    for i in 1..(a.len() + b.len()) {
        let mut s = 0;
        for f in 0..i {
            let _a = a.get(f).unwrap_or(&0);
            let _b = b.get(i - f - 1).unwrap_or(&0);
            s = (s + (_a * _b)).rem_euclid(q);
        }
        c.push(s);
    }
    c = rm_trailing_zeroes(c);
    c
}

pub fn add_ct(a: &(Vec<i32>, Vec<i32>), b: &(Vec<i32>, Vec<i32>), q: i32) -> (Vec<i32>, Vec<i32>) {
    let c0 = add_polyn(&a.0, &b.0, q);
    let c1 = add_polyn(&a.1, &b.1, q);
    (c0, c1)
}

pub fn mult_ct(
    a: &(Vec<i32>, Vec<i32>),
    b: &(Vec<i32>, Vec<i32>),
    t: i32,
    q: i32,
) -> (Vec<i32>, Vec<i32>) {
    let c0 = delta_polyn(mult_polyn(&a.0, &a.0, q), t, q);
    let _c1 = delta_polyn(mult_polyn(&a.0, &b.1, q), t, 1);
    let c1 = delta_polyn(add_polyn(&_c1, &_c1, q), 1, q);
    //let c2 = delta_polyn(mult_polyn(&b.1, &b.1, q), t, q);
    (c0, c1)
}
/**
# Vector2d
Implement Vector2d
Fields:
1. degree
Return Value:
1. Vec<i32>
The sample function can be used to sample a polynomial from the Normal Distribution.
This is done by sampling the coefficients from a normal distribution with mean = 0, stdv = 0;
In this way, with a ring R/f(x) where f(x) is some x^n + 1, we can generalize the polynomial we
create as sampled from a normal distribution.

The return value represents a polynomial in summation cx^i, where c is the value at i and the index
i represents the power. So the first term is a constant, and is the first element in the vector.

We use the paramter degree in order to specify the bound of the polynomials. Our ring is R/f(x) with
f(x) = x^degree + 1, meaning that every polynomial in the ring must be under modulo f(x). As such, we
return a polynomial with degree - 1, as f(x) == 0 in the ring.
*/
pub fn normal_sample(degree: i32) -> Vec<i32> {
    //let mut rng = rand::rngs::StdRng::seed_from_u64(100);
    let normal = Normal::new(0.0, 1.0).unwrap();
    let mut sample = Vec::new();
    for _ in 0..degree {
        let x: f32 = normal.sample(&mut rand::thread_rng());
        let x = x.trunc();
        //print!("{}, ", x);
        sample.push(x as i32);
    }
    //println!("SAMPLE MEAN: {}", val / 100.0);
    sample
}
pub fn uniform_sample(degree: i32, q: i32) -> Vec<i32> {
    //let mut rng = rand::rngs::StdRng::seed_from_u64(100);
    let uniform = Uniform::try_from(-q..q + 1).unwrap();
    let mut sample = Vec::new();
    for _ in 0..degree {
        let x: i32 = uniform.sample(&mut rand::thread_rng());
        //print!("{}, ", x);
        sample.push(x as i32);
    }
    //println!("SAMPLE MEAN: {}", val / 100.0);
    sample
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_adds() {
        let v1 = vec![1, 2, 3];
        let v2 = vec![4, 5, 6];
        let v3 = add_polyn(&v1, &v2, 10);
        assert_eq!(v3, vec![5, 7, 9]);
    }
    #[test]
    fn it_adds2() {
        let v1 = vec![1, 2, 3];
        let v2 = vec![4, 5, 6];
        let v3 = add_polyn(&v1, &v2, 8);
        assert_eq!(v3, vec![5, 7, 1]);
    }
    #[test]
    fn it_mults() {
        let v1 = vec![1, 2];
        let v2 = vec![4, 5];
        let v3 = mult_polyn(&v1, &v2, 100);
        assert_eq!(v3, vec![4, 13, 10]);
    }

    #[test]
    fn it_mults2() {
        let v1 = vec![5, 1, 3];
        let v2 = vec![2, 3, 8];
        let v3 = mult_polyn(&v1, &v2, 100);
        assert_eq!(v3, vec![10, 17, 49, 17, 24]);
    }
    #[test]
    fn it_mults3() {
        let v1 = vec![5, 1, 3];
        let v2 = vec![2, 3, 8];
        let v3 = mult_polyn(&v1, &v2, 8);
        //1 + 0 x + 0 x2+ 0x3 + 0 x4+ 0 x5+ 0x6 + 0 x7+ 1 * x8
        assert_eq!(v3, vec![2, 1, 1, 1]);
    }
    #[test]
    fn it_mults4() {
        let v1 = vec![0, 0, 1];
        let v2 = vec![2, 1];
        let v3 = mult_polyn(&v1, &v2, 10);
        assert_eq!(v3, vec![0, 0, 2, 1]);
    }
    //#[test]
    fn it_divs() {
        let v1 = vec![6, 7, 4, 1];
        let v2 = vec![2, 1];
        let v3 = div_polyn(v1, v2, &vec![0, 0, 0, 0, 0, 0, 1], 100);
        assert_eq!(v3, vec![]);
    }
    //#[test]
    fn it_divs2() {
        let v1 = vec![5, 0, 9, 4, 0, 0, 0, 0, 10];
        let v2 = vec![1, 0, 0, 0, 0, 0, 0, 0, 1];
        let v3 = div_polyn(v1, v2, &vec![1, 0, 0, 0, 0, 0, 0, 0, 1], 200);
        assert_eq!(v3, vec![-5, 0, 9, 4]);
    }
    //#[test]
    fn it_divs22() {
        let v1 = vec![2, 3, 4, 5, 0, 9, 4, 0, 0, 0, 0, 10];
        let v2 = vec![3, 1, 0, 0, 0, 0, 0, 0, 0, 1];
        let v3 = mod_polyn(
            mult_polyn(&v1, &v2, 200),
            &vec![1, 0, 0, 0, 0, 0, 0, 0, 1],
            200,
        );
        assert_eq!(v3, vec![6, 9, 12, -15, 0, 27, 12]);
    }
    //#[test]
    fn it_divs3() {
        let v1 = vec![5, -3, 9, 4, 0, 0, 0, -5, 10];
        let v2 = vec![1, 0, 91, -5, 0, 0, 0, 0, 1];
        let v3 = div_polyn(v1, v2, &vec![1, 0, 0, 0, 0, 0, 0, 0, 1], 10000);
        assert_eq!(v3, vec![-5, -3, -901, 54, 0, 0, 0, -5]);
    }
}

use std::{
    fmt::Display,
    ops::{Add, Div, Mul, Rem, Sub},
    rc::Rc,
};

use num_traits::{Euclid, FromPrimitive, One, ToPrimitive, Zero};

use crate::{
    polynomial::{Parameters, Polynomial, Sampling, Value},
    rlwe::Params,
};

#[derive(Clone, Debug)]
pub struct Fraction<Val> {
    num: Val,
    den: Val,
}

impl<Val> Fraction<Val>
where
    Val: Value,
{
    pub fn new(num: Val, den: Val) -> Self {
        Self { num, den }
    }
    pub fn frac_params(parameters: &Params<Val>) -> Parameters<Fraction<Val>> {
        let p = parameters.borrow();
        let mut temp = vec![];
        for x in p.t.clone() {
            temp.push(Fraction::new(x, Val::one()));
        }
        return Parameters {
            degree: 4,
            q: Fraction::new(p.q.clone(), Val::one()),
            t: temp,
            t_relin: Fraction::new(p.t_relin.clone(), Val::one()),
            p: Fraction::new(p.p.clone(), Val::one()),
            log_range: 7681_i64.ilog(32) as usize + 1,
            root: Fraction::new(Val::zero(), Val::zero()),
        };
    }
    pub fn frac_poly(
        poly: Polynomial<Val>,
        parameters: &Params<Fraction<Val>>,
    ) -> Polynomial<Fraction<Val>> {
        let mut temp = vec![];
        for term in poly.val {
            let fraction = Fraction::new(term, Val::one());
            temp.push(fraction);
        }
        Polynomial::new(temp, parameters)
    }
    pub fn convert(&self) -> Val {
        if self.num < Val::zero() {
            (self.num.clone() - (self.den.clone() / (Val::one() + Val::one()))) / self.den.clone()
        } else {
            (self.num.clone() + (self.den.clone() / (Val::one() + Val::one()))) / self.den.clone()
        }
    }
    pub fn flatten(&mut self) {
        println!("BEFORE, {}, {}", self.num.clone(), self.den.clone());
        let times = self.convert();
        let sub = times * self.den.clone();
        self.num = self.num.clone() - sub;
        println!("HOE KAN DIT {} {}", self.num.clone(), self.den.clone())
    }
    pub fn round_multiply(
        a: &Polynomial<Self>,
        b: &Polynomial<Self>,
        parameters: &Params<Val>,
    ) -> Polynomial<Val> {
        //let c = a * b;
        let c = a.old_mul(b);
        let mut temp = vec![];
        for term in c.val {
            let x = term.convert();
            temp.push(x);
        }
        Polynomial::new(temp, parameters)
    }
    pub fn reduce(&self) -> Self {
        //https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
        let mut a = self.num.clone();
        let mut b = self.den.clone();
        if a.is_zero() {
            return Fraction::new(Val::zero(), Val::one());
        }
        while b != Val::zero() {
            let t = b.clone();
            b = a % b;
            a = t;
        }
        return Fraction::new(self.num.clone() / a.clone(), self.den.clone() / a);
    }
    fn inverse(&self, q: Val) -> Val {
        //https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
        let mut t = Val::zero();
        let mut newt = Val::one();
        let mut r = q.clone();
        let mut newr = self.den.clone();

        while newr != Val::zero() {
            let quotient = r.clone() / newr.clone();
            (t, newt) = (newt.clone(), t - quotient.clone() * newt);
            (r, newr) = (newr.clone(), r - quotient * newr);
        }
        if r > Val::one() {
            return Val::zero();
        }
        if t < Val::zero() {
            t = t + q;
        }
        return t * self.num.clone();
    }
}

impl<Val> Display for Fraction<Val>
where
    Val: Value,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}/{}", self.num, self.den)
    }
}

impl<Val> PartialOrd for Fraction<Val>
where
    Val: Value,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        return Some(
            (self.num.clone() * other.den.clone()).cmp(&(other.num.clone() * self.den.clone())),
        );
    }
}
impl<Val> PartialEq for Fraction<Val>
where
    Val: Value,
{
    fn eq(&self, other: &Self) -> bool {
        let a = self.reduce();
        let b = other.reduce();
        return (a.num == b.num) && (a.den == b.den);
    }
}
impl<Val> Sub for Fraction<Val>
where
    Val: Value,
{
    type Output = Fraction<Val>;

    fn sub(self, rhs: Self) -> Self::Output {
        if rhs.is_zero() {
            return self;
        }
        let temp = Fraction::new(
            self.num * rhs.den.clone() - rhs.num * self.den.clone(),
            self.den * rhs.den,
        );
        temp.reduce()
    }
}
impl<Val> Add for Fraction<Val>
where
    Val: Value,
{
    type Output = Fraction<Val>;

    fn add(self, rhs: Self) -> Self::Output {
        let temp = Fraction::new(
            self.num * rhs.den.clone() + rhs.num * self.den.clone(),
            self.den * rhs.den,
        );
        temp.reduce()
    }
}

impl<Val> Euclid for Fraction<Val>
where
    Val: Value,
{
    fn div_euclid(&self, v: &Self) -> Self {
        return self.clone() % v.clone();
    }

    fn rem_euclid(&self, v: &Self) -> Self {
        return self.clone() % v.clone();
    }
}

impl<Val> Div for Fraction<Val>
where
    Val: Value,
{
    type Output = Fraction<Val>;

    fn div(self, rhs: Self) -> Self::Output {
        let temp = Fraction::new(self.num * rhs.den, self.den * rhs.num);
        temp.reduce()
    }
}

impl<Val> Mul for Fraction<Val>
where
    Val: Value,
{
    type Output = Fraction<Val>;

    fn mul(self, rhs: Self) -> Self::Output {
        let temp = Fraction::new(self.num * rhs.num, self.den * rhs.den);
        temp.reduce()
    }
}

impl<Val> Rem for Fraction<Val>
where
    Val: Value,
{
    type Output = Fraction<Val>;

    fn rem(self, rhs: Self) -> Self::Output {
        let mut t = Val::zero();
        let mut newt = Val::one();
        let mut r = rhs.num.clone();
        let mut newr = self.den;

        while newr != Val::zero() {
            let quotient = r.clone() / newr.clone();
            (t, newt) = (newt.clone(), t - quotient.clone() * newt);
            (r, newr) = (newr.clone(), r - quotient * newr);
        }
        if r > Val::one() {
            return Fraction::new(Val::zero(), Val::one());
        }
        if t < Val::zero() {
            t = t + rhs.num;
        }
        return Fraction::new(t * self.num, Val::one());
    }
}

impl<Val> FromPrimitive for Fraction<Val>
where
    Val: Value,
{
    fn from_i64(n: i64) -> Option<Self> {
        return None;
        //Some(Fraction::new(n.try_into().unwrap(), Val::one()))
    }

    fn from_u64(n: u64) -> Option<Self> {
        None
        //Some(Fraction::new(n.try_into().unwrap(), 1))
    }
}

impl<Val> One for Fraction<Val>
where
    Val: Value,
{
    fn one() -> Self {
        Fraction::new(Val::one(), Val::one())
    }
}

impl<Val> Zero for Fraction<Val>
where
    Val: Value,
{
    fn zero() -> Self {
        Fraction::new(Val::zero(), Val::one())
    }

    fn is_zero(&self) -> bool {
        self.num == Val::zero()
    }
}

impl<Val> Sampling for Fraction<Val>
where
    Val: Value,
{
    type Value = Fraction<Val>;

    fn normal_sample(parameters: Self::Value, size: usize) -> Vec<Self::Value> {
        vec![Fraction::new(Val::one(), Val::one())]
    }

    fn uniform_sample(parameters: Self::Value, size: usize) -> Vec<Self::Value> {
        vec![Fraction::new(Val::one(), Val::one())]
    }
}

impl<Val> From<Fraction<Val>> for usize {
    fn from(value: Fraction<Val>) -> Self {
        1
    }
}
impl<Val> Ord for Fraction<Val>
where
    Val: Value,
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        return (self.num.clone() * other.den.clone()).cmp(&(other.num.clone() * self.den.clone()));
    }
}
impl<Val> Eq for Fraction<Val> where Val: Value {}

impl<Val> Value for Fraction<Val> where Val: Value {}

#[cfg(test)]
mod tests {
    use crate::rlwe::RLWE;

    use super::*;

    #[test]
    fn test_reduce() {
        let mut frac = Fraction::new(5, 25);
        frac = frac.reduce();
        assert_eq!(Fraction::new(1, 5), frac);

        let mut frac = Fraction::new(24, 36);
        frac = frac.reduce();
        assert_eq!(Fraction::new(2, 3), frac);
    }

    #[test]
    fn test_inverse() {
        let mut frac = Fraction::new(1, 5);
        let inverse = frac.inverse(9);
        assert_eq!(inverse, 2);

        let mut frac = Fraction::new(3, 5);
        let inverse = frac.inverse(9);
        assert_eq!(inverse, 6);

        let mut frac = Fraction::new(1, 532);
        let inverse = frac.inverse(639);
        assert_eq!(inverse, 424);
    }
    #[test]
    fn big_test_inv() {
        let mut frac = Fraction::new(1, -160001);
        let inverse = frac.inverse(7681);
        assert_eq!(inverse, 3421);
    }

    #[test]
    fn test_add_and_mult() {
        let mut a = Fraction::new(5, 9);
        let mut b = Fraction::new(4, 6);
        let mut c = Fraction::new(19, 23);
        let x = a + b;
        let y = x * c;
        println!("{}", y);
        assert_eq!(y, Fraction::new(1254, 1242));
        let last = y.reduce();
        println!("{}", last);
        assert_eq!(last, Fraction::new(209, 207))
    }
}

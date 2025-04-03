use std::{
    fmt::Display,
    ops::{Add, Div, Mul, Rem, Sub},
};

use num_traits::{Euclid, FromPrimitive, One, ToPrimitive, Zero};

use crate::polynomial::{Sampling, Value};

#[derive(Clone, Debug)]
pub struct Fraction {
    num: i64,
    den: i64,
}

impl Fraction {
    pub fn new(num: i64, den: i64) -> Self {
        Self { num, den }
    }
    fn reduce(&self) -> Self {
        //https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
        let mut a = self.num;
        let mut b = self.den;
        while b != 0 {
            let t = b;
            b = a % b;
            a = t;
        }
        return Fraction::new(self.num / a, self.den / a);
    }
    fn inverse(&self, q: i64) -> i64 {
        //https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
        let mut t = 0;
        let mut newt = 1;
        let mut r = q;
        let mut newr = self.den;

        while newr != 0 {
            let quotient = r.clone() / newr.clone();
            (t, newt) = (newt.clone(), t - quotient.clone() * newt);
            (r, newr) = (newr.clone(), r - quotient * newr);
        }
        if r > 1 {
            return 0;
        }
        if t < 0 {
            t = t + q;
        }
        return t * self.num;
    }
}

impl Display for Fraction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}/{}", self.num, self.den)
    }
}

impl PartialOrd for Fraction {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        return Some((self.num * other.den).cmp(&(other.num * self.den)));
    }
}
impl PartialEq for Fraction {
    fn eq(&self, other: &Self) -> bool {
        let a = self.reduce();
        let b = other.reduce();
        return (a.num == b.num) && (a.den == b.den);
    }
}
impl Sub for Fraction {
    type Output = Fraction;

    fn sub(self, rhs: Self) -> Self::Output {
        Fraction::new(self.num * rhs.den - rhs.num * self.den, self.den * rhs.den)
    }
}
impl Add for Fraction {
    type Output = Fraction;

    fn add(self, rhs: Self) -> Self::Output {
        Fraction::new(self.num * rhs.den + rhs.num * self.den, self.den * rhs.den)
    }
}

impl Euclid for Fraction {
    fn div_euclid(&self, v: &Self) -> Self {
        return self.clone() % v.clone();
    }

    fn rem_euclid(&self, v: &Self) -> Self {
        return self.clone() % v.clone();
    }
}

impl Div for Fraction {
    type Output = Fraction;

    fn div(self, rhs: Self) -> Self::Output {
        Fraction::new(self.num * rhs.den, self.den * rhs.num)
    }
}

impl Mul for Fraction {
    type Output = Fraction;

    fn mul(self, rhs: Self) -> Self::Output {
        Fraction::new(self.num * rhs.num, self.den * rhs.den)
    }
}

impl Rem for Fraction {
    type Output = Fraction;

    fn rem(self, rhs: Self) -> Self::Output {
        let mut t = 0;
        let mut newt = 1;
        let mut r = rhs.num;
        let mut newr = self.den;

        while newr != 0 {
            let quotient = r.clone() / newr.clone();
            (t, newt) = (newt.clone(), t - quotient.clone() * newt);
            (r, newr) = (newr.clone(), r - quotient * newr);
        }
        if r > 1 {
            return Fraction::new(0, 1);
        }
        if t < 0 {
            t = t + rhs.num;
        }
        return Fraction::new(t * self.num, 1);
    }
}

impl FromPrimitive for Fraction {
    fn from_i64(n: i64) -> Option<Self> {
        Some(Fraction::new(n, 1))
    }

    fn from_u64(n: u64) -> Option<Self> {
        Some(Fraction::new(n.to_i64().unwrap(), 1))
    }
}

impl One for Fraction {
    fn one() -> Self {
        Fraction::new(1, 1)
    }
}

impl Zero for Fraction {
    fn zero() -> Self {
        Fraction::new(0, 1)
    }

    fn is_zero(&self) -> bool {
        self.num == 0
    }
}

impl Sampling for Fraction {
    type Value = Fraction;

    fn normal_sample(parameters: Self::Value, size: usize) -> Vec<Self::Value> {
        vec![Fraction::new(1, 1)]
    }

    fn uniform_sample(parameters: Self::Value, size: usize) -> Vec<Self::Value> {
        vec![Fraction::new(1, 1)]
    }
}

impl From<Fraction> for usize {
    fn from(value: Fraction) -> Self {
        (value.num / value.den).to_usize().unwrap()
    }
}

impl Value for Fraction {}

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

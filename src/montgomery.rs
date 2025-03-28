use num_bigint::BigUint;
use num_traits::{one, FromPrimitive, One, Zero};
use std::ops::{BitAnd, Div, Rem, Shr, Sub};
use subtle::ConditionallySelectable;
use uint::construct_uint;

use crate::polynomial::Value;
construct_uint! {
    pub struct U1024(16);
}

impl ConditionallySelectable for U1024 {
    fn conditional_select(a: &Self, b: &Self, choice: subtle::Choice) -> Self {
        return U1024(<[u64; 16]>::conditional_select(&a.0, &b.0, choice));
    }
}

impl One for U1024 {
    fn one() -> Self {
        return U1024([0; 16]) + 1;
    }
}

pub fn square_multiply<M, N>(a: M, b: N, q: M) -> M
where
    M: Value,
    N: PartialOrd + Zero + Clone + One + Rem<Output = N> + Sub<Output = N> + Div<Output = N>,
{
    if b == N::zero() {
        return M::one();
    }
    let two = N::one() + N::one();
    let mut n = b.clone();
    let mut x = a.clone();
    let mut y = M::one();
    while n > N::one() {
        if (n.clone() % two.clone() == N::one()) {
            y = (x.clone() * y) % q.clone();
            n = n - N::one();
        }
        x = (x.clone() * x.clone()) % q.clone();
        n = n / two.clone();
    }
    return (x * y) % q;
}

pub fn square_multiply_no_mod<M, N>(a: M, b: N, q: M) -> M
where
    M: Value,
    N: PartialOrd + Zero + Clone + One + Rem<Output = N> + Sub<Output = N> + Div<Output = N>,
{
    if b == N::zero() {
        return M::one();
    }
    let two = N::one() + N::one();
    let mut n = b.clone();
    let mut x = a.clone();
    let mut y = M::one();
    while n > N::one() {
        if (n.clone() % two.clone() == N::one()) {
            y = (x.clone() * y); //% q.clone();
            n = n - N::one();
        }
        x = (x.clone() * x.clone()); //% q.clone();
        n = n / two.clone();
    }
    return (x * y); //% q;
}

pub fn standard_montgomery_ladder(a: i64, b: i64, q: i64) -> i64 {
    let mut check = false;
    let mut x_1 = a;
    let mut x_2 = a * a;
    if b == 0 {
        return 1;
    }
    for x in (0..64).rev().map(|n| (b >> n) & 1) {
        if check {
            if (x % 2) == 0 {
                x_2 = x_1 * x_2;
                x_1 = x_1 * x_1;
            } else {
                x_1 = x_1 * x_2;
                x_2 = x_2 * x_2;
            }
        }
        if (x % 2) == 1 && !check {
            check = true;
        }
        x_1 = x_1.rem_euclid(q);
        x_2 = x_2.rem_euclid(q);
    }
    return x_1;
}

/*pub fn secure_montgomery_ladder<M>(a: M, b: i64, q: M) -> M
where
    M: One + Rem<Output = M>,
{
    let mut registers: [M; 2] = [one(), a];

    for x in (0..64).rev().map(|n| (b >> n) & 1) {
        registers[((x + 1) % 2) as usize] = registers[0] * registers[1];
        registers[x as usize] = registers[x as usize] * registers[x as usize];
        registers[0] = registers[0] % q;
        registers[1] = registers[1] % q;
    }
    return registers[0];
}*/

pub fn more_secure_mont_ladder<M>(a: M, b: M, q: M) -> M
where
    M: Shr<Output = M>
        + One
        + Rem<Output = M>
        + ConditionallySelectable
        + FromPrimitive
        + BitAnd<Output = M>
        + Zero
        + PartialEq,
{
    let (mut l, mut h) = (one(), a);
    for x in (0..32)
        .rev()
        .map(|n| (b >> M::from_i32(n).unwrap()) & one())
    {
        let (_l, _h) = (l, h);
        let res = x % (M::one() + one()) == M::zero();
        l = M::conditional_select(&(_l * _l), &(_l * _h), (!res as u8).into()) % q;
        h = M::conditional_select(&(_l * _h), &(_h * _h), (!res as u8).into()) % q;
    }
    return l;
}

pub fn secure_exponentiation(a: i64, b: i64, q: i64) -> i64 {
    let mut registers = [1, a];

    for x in (0..64).rev().map(|n| (b >> n) & 1) {
        registers[((x + 1) % 2) as usize] = registers[0] * registers[1];
        registers[x as usize] = registers[x as usize] * registers[x as usize];
        registers[0] = registers[0].rem_euclid(q);
        registers[1] = registers[1].rem_euclid(q);
    }
    return registers[0];
}

#[cfg(test)]
mod tests {

    use num_traits::ToPrimitive;

    use super::*;
    #[test]
    fn test_square_multiply() {
        let ls = square_multiply(2, 8, 1000);
        assert_eq!(ls, 256);
    }
    /*#[test]
    fn test_secure_ladder() {
        let ls = secure_montgomery_ladder(2, 10, 10000);
        assert_eq!(ls, 1024);
    }*/
    #[test]
    fn subtle_test() {
        let ls = more_secure_mont_ladder(2, 5, 10000);
        assert_eq!(ls, 32);
    }
    /*#[test]
    fn with_subtle() {
        for i in 0..10 {
            for j in 0..10 {
                println!("i{} : j{}", i, j);
                let ls = secure_montgomery_ladder(i, j, 10000);
                let rs = more_secure_mont_ladder(i, j, 10000);
                assert_eq!(ls, rs);
            }
        }
    }*/
    #[test]
    fn num_uint_test() {
        let mut a = num_bigint::BigInt::new(num_bigint::Sign::Plus, vec![123]);
        println!("{}", a);
        let mut b = num_bigint::BigInt::new(num_bigint::Sign::Plus, vec![2]);
        println!("{}", b);
        let mut q = num_bigint::BigInt::new(num_bigint::Sign::Plus, vec![9, 9, 9, 9, 9]);
        let ls = square_multiply(a, b, q);
        let rs = square_multiply(123, 2, 99999);
        println!("{}, {}", ls, rs);
        assert_eq!(ls.to_i64().unwrap(), rs);
    }
}

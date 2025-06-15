use montgomery::square_multiply as m_ladder;
use std::rc::Rc;

use crate::{
    montgomery,
    polynomial::{Polynomial, Value},
};

//trait Val: Rem + PartialOrd + Zero + Clone + ConditionallySelectable + One {}
//impl<T> Val for T where T: Rem + PartialOrd + Zero + Clone + ConditionallySelectable + One {}

pub struct Ntt<Val> {
    q: Val,
}
impl<Val> Ntt<Val>
where
    Val: Value,
{
    pub fn ntt_mult(a: Polynomial<Val>, b: Polynomial<Val>, q: &Val, psi: &Val) -> Polynomial<Val> {
        let a = a.get_mod().mod_q();
        let b = b.get_mod().mod_q();
        let a_ntt = Ntt::iter_cooley_tukey(a, psi, q);
        let b_ntt = Ntt::iter_cooley_tukey(b, psi, q);
        let a_b_ntt = Ntt::element_wise_mult(a_ntt, b_ntt, q);
        return Ntt::iter_gentleman_sande(a_b_ntt, psi, q);
    }
    fn primitive_root(q: &Val, degree: i64) -> Val {
        let mut i: Val = Val::zero();
        while &i < q {
            if m_ladder(i.clone(), Val::from_i64(degree).unwrap(), q.clone()).rem_euclid(&q)
                == Val::one()
            {
                let mut is_root = true;
                for j in 1..(degree - 1) {
                    if m_ladder(i.clone(), Val::from_i64(j).unwrap(), q.clone()).rem_euclid(&q)
                        == Val::one()
                    {
                        is_root = false;
                        //break;
                    }
                }
                if is_root {
                    return i;
                }
            }
            i = i + Val::one();
        }
        return Val::zero();
    }

    pub fn second_primitive_root(q: &Val, degree: i64) -> Val {
        let omega = Ntt::primitive_root(q, degree);
        let mut i: Val = Val::zero();
        while &i < q {
            if m_ladder(i.clone(), Val::one() + Val::one(), q.clone()).rem_euclid(&q) == omega {
                if m_ladder(i.clone(), Val::from_i64(degree).unwrap(), q.clone()).rem_euclid(&q)
                    == q.clone() - Val::one()
                {
                    return i;
                }
            }
            i = i + Val::one();
        }
        return Val::zero();
    }
    fn ntt(x: Polynomial<Val>, omega: &Val, q: &Val) -> Polynomial<Val> {
        let mut res = vec![];
        for i in 0..x.degree() {
            let mut res_val = Val::zero();
            for j in 0..x.degree() {
                let n = (i * j).rem_euclid(x.degree());
                //let val = omega.pow(n as u32).rem_euclid(q);
                let val = m_ladder(
                    omega.clone(),
                    Val::from_i64(n as i64).unwrap().clone(),
                    q.clone(),
                );
                res_val = res_val + x.val[j].clone() * val;
                res_val = res_val.rem_euclid(&q);
            }
            res.push(res_val);
        }
        return Polynomial::new(res, x.parameters());
    }
    pub fn second_ntt(x: Polynomial<Val>, psi: &Val, q: &Val) -> Polynomial<Val> {
        let mut res = vec![];
        for i in 0..x.degree() {
            let mut res_val = Val::zero();
            for j in 0..x.degree() {
                let n = ((2 * i * j) + j).rem_euclid(2 * x.degree());
                let val = m_ladder(psi.clone(), Val::from_i64(n as i64).unwrap(), q.clone());
                //print!("{} ", n);
                res_val = res_val + x.val[j].clone() * val;
                res_val = res_val.rem_euclid(&q);
            }
            //println!("");
            res.push(res_val);
        }
        return Polynomial::new(res, x.parameters());
    }
    fn transpose_second_ntt(x: Polynomial<Val>, psi: &Val, q: &Val) -> Polynomial<Val> {
        let mut res = vec![];
        for i in 0..x.degree() {
            let mut res_val = Val::zero();
            for j in 0..x.degree() {
                let n = ((2 * i * j) + i).rem_euclid(2 * x.degree());
                //let val = omega.pow(n as u32).rem_euclid(q);
                let val = m_ladder(psi.clone(), Val::from_i64(n as i64).unwrap(), q.clone());
                res_val = res_val + x.val[j].clone() * val.clone();
                res_val = res_val.rem_euclid(&q);
            }
            res.push(res_val);
        }
        return Polynomial::new(res, x.parameters());
    }

    fn inv_ntt(x: Polynomial<Val>, omega: &Val, q: &Val) -> Polynomial<Val> {
        let x = Ntt::ntt(x, &Ntt::inverse(omega, q), q);
        let scaling_factor = Ntt::inverse(&Val::from_i32(x.degree() as i32).unwrap(), q);
        let mut res = vec![];
        for i in 0..x.degree() {
            res.push((scaling_factor.clone() * x.val[i].clone()).rem_euclid(&q));
        }
        Polynomial::new(res, x.parameters())
    }

    pub fn second_inv_ntt(x: Polynomial<Val>, psi: &Val, q: &Val) -> Polynomial<Val> {
        let x = Ntt::transpose_second_ntt(x, &Ntt::inverse(psi, q), q);
        let scaling_factor = Ntt::inverse(&Val::from_i32(x.degree() as i32).unwrap(), q);
        let mut res = vec![];
        for i in 0..x.degree() {
            res.push((scaling_factor.clone() * x.val[i].clone()).rem_euclid(&q));
        }
        Polynomial::new(res, x.parameters())
    }

    pub fn element_wise_mult(a: Polynomial<Val>, b: Polynomial<Val>, q: &Val) -> Polynomial<Val> {
        let mut res = vec![];
        for i in 0..(a.parameters().borrow().degree as usize) {
            res.push(a.val[i].clone() * b.val[i].clone());
        }
        Polynomial::new(res, a.parameters())
    }
    pub fn element_wise_div(a: Polynomial<Val>, b: Polynomial<Val>, q: &Val) -> Polynomial<Val> {
        let mut res = vec![];
        for i in 0..(a.parameters().borrow().degree as usize) {
            if b.val[i] != Val::zero() {
                res.push(a.val[i].clone() / b.val[i].clone());
            }
        }
        Polynomial::new(res, a.parameters())
    }

    fn inverse(a: &Val, q: &Val) -> Val {
        //https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
        let mut t = Val::zero();
        let mut newt = Val::one();
        let mut r = q.clone();
        let mut newr = a.clone();

        while newr != Val::zero() {
            let quotient = r.clone() / newr.clone();
            (t, newt) = (newt.clone(), t - quotient.clone() * newt);
            (r, newr) = (newr.clone(), r - quotient * newr);
        }
        if r > Val::one() {
            return Val::zero();
        }
        if t < Val::zero() {
            t = t + q.clone();
        }
        return t;
    }
    fn poly_step_by(x: Polynomial<Val>, offset: usize, n: i32) -> Polynomial<Val> {
        let parameters = Rc::clone(x.parameters());
        let Val = x.val;
        let res = Val
            .iter()
            .skip(offset)
            .step_by(2)
            .cloned()
            .collect::<Vec<_>>();
        Polynomial::new(res, &parameters)
    }

    fn cooley_tukey(x: Polynomial<Val>, psi: &Val, q: &Val) -> Polynomial<Val> {
        println!("degree: {}", x.degree());
        let psi = &Ntt::second_primitive_root(q, x.degree() as i64);
        if x.degree() <= 4 {
            return Ntt::second_ntt(x, psi, q);
        }
        let x_1 = Ntt::poly_step_by(x.clone(), 0, 2);
        let x_2 = Ntt::poly_step_by(x.clone(), 1, 2);
        let mut x_1 = Ntt::cooley_tukey(
            x_1,
            &Ntt::second_primitive_root(q, x.degree() as i64 / 2),
            q,
        );
        let mut x_2 = Ntt::cooley_tukey(
            x_2,
            &Ntt::second_primitive_root(q, x.degree() as i64 / 2),
            q,
        );
        //let mut x_1 = cooley_tukey(x_1, psi, q);
        //let mut x_2 = cooley_tukey(x_2, psi, q);
        //let mut x_res: Vec<i32> = interleave(x_1.val, x_2.val).collect::<Vec<_>>();
        let mut x_res = x_1.val.clone();
        x_res.append(&mut x_2.val);
        let psi = Ntt::second_primitive_root(q, x.degree() as i64 / 2);
        for i in 0..x.degree() / 2 {
            let a = x_res[i].clone();
            let factor = m_ladder(psi.clone(), Val::from_i64(i as i64).unwrap(), q.clone());
            let b = factor * x_res[i + x.degree() / 2].clone();
            x_res[i] = (a.clone() + b.clone()).rem_euclid(&q);
            x_res[i + x.degree() / 2] = (a - b).rem_euclid(&q);
        }
        return Polynomial::new(x_res, x.parameters());
    }
    fn bit_reverse(num: i64, degree: i64, q: Val) -> usize {
        let mut counter = Val::zero();
        let mut result = Val::zero();
        for x in (0..degree.checked_ilog2().unwrap())
            .rev()
            .map(|n| (num >> n) & 1)
        {
            result = result
                + Val::from_i64(x).unwrap()
                    * m_ladder(Val::one() + Val::one(), counter.clone(), q.clone());
            counter = counter + Val::one();
        }
        return result.try_into().unwrap_or(1);
    }
    fn stored_roots(psi: &Val, q: &Val, degree: i64) -> Vec<Val> {
        let mut temp = vec![Val::zero(); degree as usize];
        for i in 0..degree {
            temp[i as usize] = m_ladder(
                psi.clone(),
                Ntt::bit_reverse(i, degree, q.clone()),
                q.clone(),
            );
        }
        return temp;
    }

    fn bit_reverse_Polynomial(a: Polynomial<Val>, q: Val) -> Polynomial<Val> {
        let mut res = vec![];
        for i in 0..a.degree() {
            let index: usize = Ntt::bit_reverse(i as i64, a.degree() as i64, q.clone());
            res.push(a.val[index].clone());
        }
        return Polynomial::new(res, a.parameters());
    }

    pub fn iter_cooley_tukey(x: Polynomial<Val>, psi: &Val, q: &Val) -> Polynomial<Val> {
        let degree = x.parameters().borrow().degree as usize;
        let roots = Ntt::stored_roots(psi, q, degree as i64);
        let mut a = x.val.clone();
        let mut t = degree;
        let mut i = 1;
        while i < degree {
            t = t / 2;
            for j in 0..i {
                let j1 = 2 * j * t;
                let j2 = j1 + t;
                let root = roots[i + j].clone();
                //Also works:
                //let root = m_ladder(psi, bit_reverse((i + j) as i64, x.degree() as i64, q), q);
                for k in j1..j2 {
                    let u = a[k].clone();
                    let v = a[k + t].clone() * root.clone();
                    a[k] = (u.clone() + v.clone()).rem_euclid(&q);
                    a[k + t] = (u - v).rem_euclid(&q);
                }
            }
            i = i * 2;
        }
        return Polynomial::new(a, x.parameters());
    }
    pub fn iter_gentleman_sande(x: Polynomial<Val>, psi: &Val, q: &Val) -> Polynomial<Val> {
        let degree = x.parameters().borrow().degree as usize;
        let roots = Ntt::stored_roots(&Ntt::inverse(psi, q), q, degree as i64);
        let mut a = x.val.clone();
        let mut t = 1;
        let mut i = degree;
        while i > 1 {
            i = i / 2;
            let mut j1 = 0;
            for j in 0..i {
                let j2 = j1 + t;
                let root = roots[i + j].clone();
                //Also works:
                //let root = m_ladder(inverse(psi, q), bit_reverse((i + j) as i64, x.degree() as i64, q), q);
                for k in j1..j2 {
                    let u = a[k].clone();
                    let v = a[k + t].clone();
                    a[k] = (u.clone() + v.clone()).rem_euclid(&q);
                    a[k + t] = ((u - v) * root.clone()).rem_euclid(&q);
                }
                j1 = j1 + 2 * t;
            }
            t = t * 2;
        }
        let scaling_factor = Ntt::inverse(&Val::from_i32(x.degree() as i32).unwrap(), q);
        for i in 0..x.degree() {
            a[i] = (scaling_factor.clone() * a[i].clone()).rem_euclid(&q);
        }
        return Polynomial::new(a, x.parameters());
    }
}

#[cfg(test)]
mod tests {
    use std::cell::RefCell;

    use crate::{polynomial::Parameters, print_vec, rlwe::get_def_params};

    use super::*;
    #[test]

    fn test_cooley_tukey() {
        let _p = &Rc::new(RefCell::new(get_def_params()));
        let q = &_p.borrow().q;
        let degree = _p.borrow().degree;
        let psi = &Ntt::second_primitive_root(q, degree);
        println!("PSI FOR DEGREE 2: {}", Ntt::second_primitive_root(q, 1));
        let ls = Ntt::cooley_tukey(Polynomial::new(vec![1, 2, 3, 4], _p), psi, q);
        print_vec(&ls.val);
        assert_eq!(ls, Polynomial::new(vec![1467, 2807, 3471, 7621], _p));
    }
    //#[test]
    fn test_cooley_tukey_2() {
        let _p = &Rc::new(RefCell::new(get_def_params()));
        let q = &_p.borrow().q;
        let psi = &Ntt::second_primitive_root(q, 8);
        let ls = Ntt::cooley_tukey(Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7, 8], _p), psi, q);
        println!("Cooley-Tukey result: ");
        print_vec(&ls.val);
        let test_val = Ntt::second_ntt(Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7, 8], _p), psi, q);
        println!("naive result: ");
        print_vec(&test_val.val);
        let ls = Ntt::second_inv_ntt(ls, psi, q);
        assert_eq!(ls, Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7, 8], _p));
    }
    #[test]
    fn test_iter_cooley_tukey() {
        let _p = &Rc::new(RefCell::new(get_def_params()));

        let q = &_p.borrow().q;
        let psi = &Ntt::second_primitive_root(q, 4);
        let mut ls = Ntt::iter_cooley_tukey(Polynomial::new(vec![1, 2, 3, 4], _p), psi, q);
        println!("Cooley-Tukey result: ");
        ls = Ntt::bit_reverse_Polynomial(ls, q.clone());
        print_vec(&ls.val);
        let test_val = Ntt::second_ntt(Polynomial::new(vec![1, 2, 3, 4], _p), psi, q);
        println!("naive result: ");
        print_vec(&test_val.val);
        let ls = Ntt::second_inv_ntt(ls, psi, q);
        assert_eq!(ls, Polynomial::new(vec![1, 2, 3, 4], _p));
    }
    #[test]
    fn test_iter_cooley_tukey_2() {
        let params = Parameters {
            degree: 8,
            q: 7681,
            t: vec![6],
            t_relin: 32,
            p: 4,
            log_range: 0,
            root: 0,
        };
        let _p = &Rc::new(RefCell::new(params));

        let q = &_p.borrow().q;
        let psi = &Ntt::second_primitive_root(q, 8);
        let mut ls =
            Ntt::iter_cooley_tukey(Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7, 8], _p), psi, q);
        ls = Ntt::bit_reverse_Polynomial(ls, q.clone());
        println!("Cooley-Tukey result: ");
        print_vec(&ls.val);
        let test_val = Ntt::second_ntt(Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7, 8], _p), psi, q);
        println!("naive result: ");
        print_vec(&test_val.val);
        let ls = Ntt::second_inv_ntt(ls, psi, q);
        assert_eq!(ls, Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7, 8], _p));
    }
    #[test]
    fn test_iter_gentleman_sande() {
        let _p = &Rc::new(RefCell::new(get_def_params()));

        let q = &7681;
        let psi = &Ntt::second_primitive_root(q, 4);
        let mut ls = Ntt::iter_cooley_tukey(Polynomial::new(vec![1, 2, 3, 4], _p), psi, q);
        println!("Cooley-Tukey result: ");
        //ls = bit_reverse_Polynomial(ls, q);
        print_vec(&ls.val);
        let test_val = Ntt::second_ntt(Polynomial::new(vec![1, 2, 3, 4], _p), psi, q);
        println!("naive result: ");
        print_vec(&test_val.val);
        let ls = Ntt::iter_gentleman_sande(ls, psi, q);
        assert_eq!(ls, Polynomial::new(vec![1, 2, 3, 4], _p));
    }
    #[test]
    fn m_ladder_test_1() {
        assert_eq!(3383_i64.pow(3).rem_euclid(5), m_ladder(3383, 3, 5));
        assert_eq!(3383_i64.pow(3).rem_euclid(7681), m_ladder(3383, 3, 7681));
    }
    #[test]
    fn m_ladder_test_2() {
        let omega: i64 = 3383;
        //let n = 3;
        let q = 7681;
        for i in 0_i64..16 {
            let n = i.rem_euclid(4);
            let val_1 = omega.pow(n as u32).rem_euclid(q);
            let val_2 = m_ladder(omega, n as i64, q as i64);
            assert_eq!(val_1, val_2);
        }
    }
    #[test]
    fn roots() {
        let ls = Ntt::primitive_root(&7681, 4);
        assert_eq!(ls, 3383);
        let ls = Ntt::primitive_root(&7681, 9);
        assert_eq!(ls, 0);
    }

    #[test]
    fn ntt_result() {
        let _p = &Rc::new(RefCell::new(get_def_params()));

        let omega = Ntt::primitive_root(&7681, 4);
        let ls = Ntt::ntt(Polynomial::new(vec![1, 2, 3, 4], _p), &omega, &7681);
        print_vec(&ls.val);
        assert_eq!(ls, Polynomial::new(vec![10, 913, 7679, 6764], _p));
    }
    #[test]
    fn second_ntt_result() {
        let _p = &Rc::new(RefCell::new(get_def_params()));

        let psi = &Ntt::second_primitive_root(&7681, 4);
        let ls = Ntt::second_ntt(Polynomial::new(vec![1, 2, 3, 4], _p), psi, &7681);
        print_vec(&ls.val);
        assert_eq!(ls, Polynomial::new(vec![1467, 2807, 3471, 7621], _p));
    }
    #[test]
    fn t_inverse() {
        let omega = &Ntt::primitive_root(&7681, 4);
        assert_eq!(Ntt::inverse(omega, &7681), 4298);
        assert_eq!(Ntt::inverse(&4, &7681), 5761);
    }
    #[test]
    fn test_inv_ntt() {
        let _p = &Rc::new(RefCell::new(get_def_params()));

        let initial_vector = Polynomial::new(vec![1, 2, 3, 4], _p);
        let q = &7681;
        let omega = &Ntt::primitive_root(q, 4);
        let step_1 = Ntt::ntt(initial_vector.clone(), omega, q);
        let ls = Ntt::inv_ntt(step_1, omega, q);
        assert_eq!(ls, initial_vector);
    }
    #[test]
    fn test_ntt_mult() {
        let _p = &Rc::new(RefCell::new(get_def_params()));

        let poly_1 = Polynomial::new(vec![1, 2, 3, 4], _p);
        let poly_2 = Polynomial::new(vec![5, 6, 7, 8], _p);
        let q = &_p.borrow().q;
        //let q =
        let omega = &Ntt::primitive_root(q, 4);
        let dft_1 = Ntt::ntt(poly_1.clone(), omega, q);
        let dft_2 = Ntt::ntt(poly_2.clone(), omega, q);
        println!("DFTS:");
        print_vec(&dft_1.val);
        print_vec(&dft_2.val);
        println!("----");
        let mult_dft = Ntt::element_wise_mult(dft_1, dft_2, q);
        println!("mult_DFT:");
        print_vec(&mult_dft.val);
        let ls = Ntt::inv_ntt(mult_dft.mod_q(), omega, q);
        print_vec(&ls.val);
        assert_eq!(ls, Polynomial::new(vec![66, 68, 66, 60], _p));
    }
    #[test]
    fn test_second_prim_root() {
        let ls = Ntt::second_primitive_root(&7681, 4);
        assert_eq!(ls, 1925);
    }
    //#[test]
    fn test_second_inv_ntt() {
        let _p = &Rc::new(RefCell::new(get_def_params()));

        let initial_vector = Polynomial::new(vec![1, 2, 3, 4, 5, 6, 7], _p);
        let q = &7681;
        let psi = &Ntt::second_primitive_root(q, 7);
        let step_1 = Ntt::second_ntt(initial_vector.clone(), psi, q);
        print_vec(&step_1.val);
        let ls = Ntt::second_inv_ntt(step_1, psi, q);
        assert_eq!(ls, initial_vector);
    }
    #[test]
    fn test_second_ntt_mult() {
        let _p = &Rc::new(RefCell::new(get_def_params()));

        let poly_1 = Polynomial::new(vec![1, 2, 3, 4], _p);
        let poly_2 = Polynomial::new(vec![5, 6, 7, 8], _p);
        let q = &_p.borrow().q;
        let psi = &Ntt::second_primitive_root(q, 4);
        let dft_1 = Ntt::second_ntt(poly_1.clone(), psi, q);
        let dft_2 = Ntt::second_ntt(poly_2.clone(), psi, q);
        println!("DFTS:");
        print_vec(&dft_1.val);
        print_vec(&dft_2.val);
        println!("----");
        let mult_dft = Ntt::element_wise_mult(dft_1, dft_2, q);
        println!("mult_DFT:");
        print_vec(&mult_dft.val);
        let ls = Ntt::second_inv_ntt(mult_dft.mod_q(), psi, q);
        let rs = (&poly_1 * &poly_2).get_mod().mod_q();
        print_vec(&ls.val);
        assert_eq!(ls, rs);
    }
    #[test]
    fn test_bit_reversal() {
        let x = Ntt::bit_reverse(7, 8, 7681);
        assert_eq!(x, 7);
    }
    #[test]
    fn test_storing_roots() {
        let q = &7681;
        let degree = 8;
        let psi = &Ntt::second_primitive_root(q, degree);
        let res = Ntt::stored_roots(psi, q, degree);
        for i in 0..degree {
            print!(", {}", m_ladder(psi.clone(), i, q.clone()));
        }
        println!();
        for i in 0..degree {
            print!(", {}", res[i as usize]);
        }
        assert_eq!(res, [1, 4298, 1213, 5756, 527, 6832, 1728, 7098]);
    }
    #[test]
    fn final_test_ntt_mult() {
        let _p = &Rc::new(RefCell::new(get_def_params()));

        let q = &_p.borrow().q;
        let degree = _p.borrow().degree;
        let psi = &Ntt::second_primitive_root(q, degree);
        let a = Polynomial::new(vec![1, 2, 3, 4], _p);
        let b = Polynomial::new(vec![5, 6, 7, 8], _p);
        let ls = Ntt::ntt_mult(a.clone(), b.clone(), q, psi);
        print_vec(&(&a * &b).val);
        let rs = (&a.old_mul(&b)).get_mod().mod_q();
        print_vec(&ls.val);
        print_vec(&rs.val);
        assert_eq!(ls, rs);
        //assert_eq!(1, 2);
    }
    #[test]
    fn generic_testing_ntt() {
        let ls = Ntt::primitive_root(&7681, 4);
        assert_eq!(ls, 3383);

        let ls = Ntt::primitive_root(&7681, 9);
        assert_eq!(ls, 0);

        let ls = Ntt::second_primitive_root(&7681, 4);
        assert_eq!(ls, 1925);
    }
}

use crate::print_vec;

pub fn rec_dft<'a>(a: Vec<i32>, n: i32) -> Vec<i32> {
    let mut _a = vec![0; n as usize];
    if n == 1 {
        return vec![a[0].clone()];
        //return a;
    }
    let w_n = n;
    let mut w = 1;
    let even = a.iter().cloned().step_by(2).collect::<Vec<_>>();
    let odd = a.iter().cloned().skip(1).step_by(2).collect::<Vec<_>>();
    let _a1 = rec_dft(even, n / 2);
    let _a2 = rec_dft(odd, n / 2);
    for k in 0..(n / 2) {
        _a[k as usize] = _a1[k as usize] + w * _a2[k as usize];
        _a[(k + n / 2) as usize] = _a1[k as usize] - w * _a2[k as usize];
        w *= w_n;
    }
    print_vec(&_a);
    _a
}

pub fn m_ladder(a: i64, b: i64) -> i64 {
    let mut check = false;
    let mut x_1 = a;
    let mut x_2 = a * a;
    for x in (0..64).rev().map(|n| (b >> n) & 1) {
        if !check {
            print!(":");
        }
        println!("{}", x);
        if check {
            if (x % 2) == 0 {
                println!("{}, {}", x_1, x_2);
                x_2 = x_1 * x_2;
                x_1 = x_1 * x_1;
            } else {
                println!("{}, {}", x_1, x_2);
                x_1 = x_1 * x_2;
                x_2 = x_2 * x_2;
            }
        }
        if (x % 2) == 1 && !check {
            check = true;
        }
    }
    return x_1;
}

#[cfg(test)]
mod tests {
    use super::*;

    //#[test]
    fn testing() {
        let result = rec_dft(vec![1, 4, 7, 4], 4);
        assert_eq!(vec![1, 3], result);
    }

    #[test]
    fn ladder_mult() {
        assert_eq!(m_ladder(2, 2), 4);
        assert_eq!(m_ladder(4, 2), 16);
        assert_eq!(m_ladder(6, 3), 216);
    }
}

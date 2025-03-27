use crate::lwe_old as lwe;
use std::fmt::Display;
use std::ops::AddAssign;
use std::{collections::HashMap, fmt::Write};
use std::{fs, thread};

#[derive(PartialEq, Eq, Hash)]
pub struct Matrix {
    q1: Vec<i32>,
    q2: Vec<i32>,
    q3: Vec<i32>,
    q4: Vec<i32>,
}
impl Matrix {
    pub fn new(q1: Vec<i32>, q2: Vec<i32>, q3: Vec<i32>, q4: Vec<i32>) -> Self {
        Self { q1, q2, q3, q4 }
    }
}

#[derive(PartialEq, Eq, Hash)]
pub struct Matrix2 {
    data: Vec<Vec<Vec<i32>>>,
}
impl Matrix2 {
    pub fn new(data: Vec<Vec<Vec<i32>>>) -> Self {
        Self { data }
    }
}

pub struct Poly {
    val: Vec<i32>,
    degree: i32,
    q: i32,
    cur: usize,
}
impl Poly {
    pub fn new(_degree: i32, _q: i32) -> Self {
        Self {
            val: Vec::new(),
            degree,
            q,
            cur: 0,
        }
    }
}
impl Iterator for Poly {
    type Item = Vec<i32>;
    fn next(&mut self) -> Option<Self::Item> {
        //print!("EA");
        //print!("EA");
        while self.val.len() < self.degree as usize {
            self.val.push(0);
            //print!("EAFA");
        }
        if self.val.len() == self.degree as usize && self.val.last().unwrap() < &self.q {
            self.val.last_mut().unwrap().add_assign(1);
            return Some(self.val.clone());
        }
        while self.val.len() > 1 {
            self.val.pop();
            if self.val.last().unwrap() < &self.q {
                self.val.last_mut().unwrap().add_assign(1);
                //print!("HMM: ");
                //print_vec(&self.val);
                return self.next();
            }
        }

        /*if self.val.last().unwrap() < &self.q {
            self.val.last_mut().unwrap().add_assign(1);
            return Some(self.val.clone());
        } else if self.val.len() < self.degree as usize {
            self.val.push(-self.q);
            return Some(self.val.clone());
            }*/
        None
    }
}

impl Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        //write!(f, "'");
        write_vec(f, &self.q1);
        //write!(f, "\n");
        write_vec(f, &self.q2);
        //write!(f, "\n");
        write_vec(f, &self.q3);
        //write!(f, "\n");
        write_vec(f, &self.q4);
        write!(f, "'")
    }
}
impl Display for Matrix2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        //write!(f, "'");
        for i in &self.data {
            for j in i {
                write_vec(f, &j);
            }
        }
        write!(f, "'")
    }
}
pub fn write_vec(f: &mut std::fmt::Formatter<'_>, some_vec: &Vec<impl Display>) {
    write!(f, "[");
    for i in 0..some_vec.len() {
        write!(f, "{}", some_vec.get(i).unwrap());
        if i < some_vec.len() - 1 {
            write!(f, ", ");
        }
    }
    write!(f, "]");
}
const p: i32 = 10;
const q: i32 = 2;
const degree: i32 = 2;

pub fn make_matrix2(q1: i32, degree1: i32, fx1: &Vec<i32>, e: bool) -> Matrix2 {
    let mut results = Vec::new();
    let mut au = Vec::new();
    let mut eu = Vec::new();
    let mut u = Vec::new();
    let mut ah = Vec::new();
    let mut eh = Vec::new();
    for _ in 0..2 {
        u.push(lwe::uniform_sample(degree1, q1));
        au.push(lwe::uniform_sample(degree1, q1));
        eu.push(vec![]); //lwe::uniform_sample(degree1, q1));
    }
    for _ in 0..2 {
        ah.push(lwe::uniform_sample(degree1, q1));
        //eh.push(lwe::uniform_sample(degree1, q1));
        eh.push(vec![]); //lwe::uniform_sample(degree1, q1));
    }
    for i in 0..2 {
        //let ah = lwe::uniform_sample(degree1, q1);
        //let eh = vec![]; //lwe::uniform_sample(degree1, q1);
        let mut temp = Vec::new();
        for j in 0..2 {
            /*let a = lwe::add_polyn(&ah[i], &au[j], q1);
            let e = lwe::add_polyn(&eh[i], &eu[j], q1);
            let result = lwe::add_polyn(&lwe::mult_polyn(&a, &u[j], q1), &e, q1);
            let result = lwe::mod_polyn(result, &fx1, q1);
            temp.push(result);*/
            print!("i : {}", i);
            print!(" and j : {}", j);
            let a = lwe::add_polyn(&ah[i], &au[j], q1);
            let e = lwe::add_polyn(&eh[i], &eu[j], q1);
            let result = lwe::add_polyn(&lwe::mult_polyn(&a, &u[j], q1), &e, q1);
            let result = lwe::mod_polyn(result, &fx1, q1);
            temp.push(result);
        }
        println!();
        results.push(temp);
    }
    println!();
    Matrix2::new(results)
}

pub fn make_matrix(q1: i32, degree1: i32, fx1: &Vec<i32>, e: bool) -> Matrix2 {
    let mut au = Vec::new();
    let mut eu = Vec::new();
    let mut u = Vec::new();
    let mut ah = Vec::new();
    let mut eh = Vec::new();
    for _ in 0..2 {
        u.push(lwe::uniform_sample(degree1, q1));
        au.push(lwe::uniform_sample(degree1, q1));
        eu.push(vec![]); //lwe::uniform_sample(degree1, q1));
    }
    for _ in 0..2 {
        ah.push(lwe::uniform_sample(degree1, q1));
        //eh.push(lwe::uniform_sample(degree1, q1));
        eh.push(vec![]); //lwe::uniform_sample(degree1, q1));
    }

    let a = lwe::add_polyn(&ah[0], &au[0], q1);
    let e = lwe::add_polyn(&eh[0], &eu[0], q1);
    let topleft = lwe::add_polyn(&lwe::mult_polyn(&a, &u[0], q1), &e, q1);
    let topleft = lwe::mod_polyn(topleft, &fx1, q1);
    //print_vec(&topleft);
    //topright
    let a = lwe::add_polyn(&ah[0], &au[1], q1);
    let e = lwe::add_polyn(&eh[0], &au[1], q1);
    let topright = lwe::add_polyn(&lwe::mult_polyn(&a, &u[1], q1), &e, q1);
    let topright = lwe::mod_polyn(topright, &fx1, q1);
    //botleft
    let a = lwe::add_polyn(&ah[1], &au[0], q1);
    let e = lwe::add_polyn(&eh[1], &eu[0], q1);
    let botleft = lwe::add_polyn(&lwe::mult_polyn(&a, &u[0], q1), &e, q1);
    let botleft = lwe::mod_polyn(botleft, &fx1, q1);
    //botright
    let a = lwe::add_polyn(&ah[1], &au[1], q1);
    let e = lwe::add_polyn(&eh[1], &eu[1], q1);
    let botright = lwe::add_polyn(&lwe::mult_polyn(&a, &u[1], q1), &e, q1);
    let botright = lwe::mod_polyn(botright, &fx1, q1);

    Matrix2::new(vec![vec![topleft, topright], vec![botleft, botright]])

    /*let u1 = lwe::uniform_sample(degree1, q1);
    let u2 = lwe::uniform_sample(degree1, q1);
    let au1 = lwe::uniform_sample(degree1, q1);
    let au2 = lwe::uniform_sample(degree1, q1);
    let mut eu1 = vec![];
    let mut eu2 = vec![];
    if (e) {
        eu1 = /*vec![]; */lwe::uniform_sample(degree1, q1);
        eu2 = /*vec![]; */lwe::uniform_sample(degree1, q1);
    }

    let ah1 = lwe::uniform_sample(degree1, q1);
    let ah2 = lwe::uniform_sample(degree1, q1);
    let eh1 = vec![]; //lwe::uniform_sample(degree1, q1);
    let eh2 = vec![]; //lwe::uniform_sample(degree1, q1);
                      //print_vec(&eh1);
                      //topleft
    let a = lwe::add_polyn(&ah1, &au1, q1);
    let e = lwe::add_polyn(&eh1, &eu1, q1);
    let topleft = lwe::add_polyn(&lwe::mult_polyn(&a, &u1, q1), &e, q1);
    let topleft = lwe::mod_polyn(topleft, &fx1, q1);
    //print_vec(&topleft);
    //topright
    let a = lwe::add_polyn(&ah1, &au2, q1);
    let e = lwe::add_polyn(&eh1, &au2, q1);
    let topright = lwe::add_polyn(&lwe::mult_polyn(&a, &u2, q1), &e, q1);
    let topright = lwe::mod_polyn(topright, &fx1, q1);
    //botleft
    let a = lwe::add_polyn(&ah2, &au1, q1);
    let e = lwe::add_polyn(&eh2, &eu1, q1);
    let botleft = lwe::add_polyn(&lwe::mult_polyn(&a, &u1, q1), &e, q1);
    let botleft = lwe::mod_polyn(botleft, &fx1, q1);
    //botright
    let a = lwe::add_polyn(&ah2, &au2, q1);
    let e = lwe::add_polyn(&eh2, &eu2, q1);
    let botright = lwe::add_polyn(&lwe::mult_polyn(&a, &u2, q1), &e, q1);
    let botright = lwe::mod_polyn(botright, &fx1, q1);

    Matrix2::new(vec![vec![topleft, topright], vec![botleft, botright]])*/
}
pub fn parallel_test(q1: i32, degree1: i32, fx1: Vec<i32>, e: bool, path1: &str, path2: &str) {
    let worker_count = 1;
    let n_elements = 100000; //100000;

    let handles: Vec<_> = (0..worker_count)
        .map(|_| {
            thread::spawn(move || {
                let mut local_hashmap = HashMap::new();

                for _ in 0..n_elements {
                    let ffx = vec![1, 0, 1];
                    let item = make_matrix2(q1, degree1, &ffx, e);

                    ////

                    *local_hashmap.entry(item).or_insert(0) += 1;
                }

                local_hashmap
            })
        })
        .collect();

    let mut freq_map = HashMap::new();

    for handle in handles {
        let local_hashmap = handle.join().unwrap();
        for (c, count) in local_hashmap {
            *freq_map.entry(c).or_insert(0) += count;
        }
    }
    let mut it = freq_map.iter().peekable();
    let mut buf = String::new();
    let mut lab = String::new();
    //write!(&mut buf, "[");
    while let Some((key, val)) = it.next() {
        write!(&mut buf, "{}", val);
        write!(&mut lab, "{}", key);
        if !it.peek().is_none() {
            write!(&mut buf, " ");
            write!(&mut lab, " ");
        }
    }

    //write!(&mut buf, "]");
    //let path = "/Users/rensweerman/PycharmProjects/pythonProject/valuesp2.txt";
    fs::write(path1, buf).expect("Unable to write file");
    //let path = "/Users/rensweerman/PycharmProjects/pythonProject/labelsp2.txt";
    fs::write(path2, lab).expect("Unable to write file");
}

// pub fn test() {
//     let mut map = HashMap::new();
//     let fx = vec![1, 0, 0, 0, 0, 1];

//     //8924237
//     //let u1 = lwe::uniform_sample(degree, q); //vec![];
//     //let u2 = lwe::uniform_sample(degree, q);
//     for i in 0..10000000 {
//         //let item = make_matrix(true);
//         let count = map.entry(item).or_insert(0);
//         *count += 1;
//     }
//     let mut it = map.iter().peekable();
//     let mut buf = String::new();
//     let mut lab = String::new();
//     //write!(&mut buf, "[");
//     while let Some((key, val)) = it.next() {
//         write!(&mut buf, "{}", val);
//         write!(&mut lab, "{}", key);
//         if !it.peek().is_none() {
//             write!(&mut buf, " ");
//             write!(&mut lab, " ");
//         }
//     }

//     //write!(&mut buf, "]");
//     let path = "/Users/rensweerman/PycharmProjects/pythonProject/valuestest.txt";
//     fs::write(path, buf).expect("Unable to write file");
//     let path = "/Users/rensweerman/PycharmProjects/pythonProject/labelstest.txt";
//     fs::write(path, lab).expect("Unable to write file");
// }

pub fn test_no_e() {
    let fx = vec![1, 0, 1];
    let mut map = HashMap::new();
    //8924237
    //let u1 = lwe::uniform_sample(degree, q); //vec![];
    //let u2 = lwe::uniform_sample(degree, q);
    for _ in 0..1000000 {
        let u1 = lwe::uniform_sample(degree, q);
        let u2 = lwe::uniform_sample(degree, q);
        let au1 = lwe::uniform_sample(degree, q);
        let au2 = lwe::uniform_sample(degree, q);
        let eu1 = vec![]; //lwe::uniform_sample(degree, q);
        let eu2 = vec![]; //lwe::uniform_sample(degree, q);

        let ah1 = lwe::uniform_sample(degree, q);
        let ah2 = lwe::uniform_sample(degree, q);
        let eh1 = vec![]; //lwe::uniform_sample(degree, q);
        let eh2 = vec![]; //lwe::uniform_sample(degree, q);
                          //print_vec(&eh1);
                          //topleft
        let a = lwe::add_polyn(&ah1, &au1, q);
        let e = lwe::add_polyn(&eh1, &eu1, q);
        let topleft = lwe::add_polyn(&lwe::mult_polyn(&a, &u1, q), &e, q);
        let topleft = lwe::mod_polyn(topleft, &fx, q);
        //print_vec(&topleft);
        //topright
        let a = lwe::add_polyn(&ah1, &au2, q);
        let e = lwe::add_polyn(&eh1, &au2, q);
        let topright = lwe::add_polyn(&lwe::mult_polyn(&a, &u2, q), &e, q);
        let topright = lwe::mod_polyn(topright, &fx, q);
        //botleft
        let a = lwe::add_polyn(&ah2, &au1, q);
        let e = lwe::add_polyn(&eh2, &eu1, q);
        let botleft = lwe::add_polyn(&lwe::mult_polyn(&a, &u1, q), &e, q);
        let botleft = lwe::mod_polyn(botleft, &fx, q);
        //botright
        let a = lwe::add_polyn(&ah2, &au2, q);
        let e = lwe::add_polyn(&eh2, &eu2, q);
        let botright = lwe::add_polyn(&lwe::mult_polyn(&a, &u2, q), &e, q);
        let botright = lwe::mod_polyn(botright, &fx, q);

        let item = Matrix::new(topleft, topright, botleft, botright);
        let count = map.entry(item).or_insert(0);
        *count += 1;
    }
    // for (key, val) in map.iter() {
    //     print!("{}", val);
    //     if val < map {
    //         print!(", ");
    //     }
    // }
    let mut it = map.iter().peekable();
    let mut buf = String::new();
    let mut lab = String::new();
    //write!(&mut buf, "[");
    while let Some((key, val)) = it.next() {
        write!(&mut buf, "{}", val);
        write!(&mut lab, "{}", key);
        if !it.peek().is_none() {
            write!(&mut buf, " ");
            write!(&mut lab, " ");
        }
    }

    //write!(&mut buf, "]");
    let path = "/Users/rensweerman/PycharmProjects/pythonProject/valuesnoe.txt";
    fs::write(path, buf).expect("Unable to write file");
    let path = "/Users/rensweerman/PycharmProjects/pythonProject/labelsnoe.txt";
    fs::write(path, lab).expect("Unable to write file");
}

pub fn test2() {
    let fx = vec![1, 0, 1];
    let mut map = HashMap::new();
    //8924237
    Poly::new(degree, q).for_each(|u1| {
        Poly::new(degree, q).for_each(|u2| {
            Poly::new(degree, q).for_each(|au1| {
                Poly::new(degree, q).for_each(|au2| {
                    Poly::new(degree, q).for_each(|eu1| {
                        Poly::new(degree, q).for_each(|eu2| {
                            Poly::new(degree, q).for_each(|ah1| {
                                Poly::new(degree, q).for_each(|ah2| {
                                    Poly::new(degree, q).for_each(|eh1| {
                                        Poly::new(degree, q).for_each(|eh2| {
                                            //topleft
                                            let a = lwe::add_polyn(&ah1, &au1, q);
                                            let e = lwe::add_polyn(&eh1, &eu1, q);
                                            let topleft =
                                                lwe::add_polyn(&lwe::mult_polyn(&a, &u1, q), &e, q);
                                            let topleft = lwe::mod_polyn(topleft, &fx, q);
                                            //print_vec(&topleft);
                                            //topright
                                            let a = lwe::add_polyn(&ah1, &au2, q);
                                            let e = lwe::add_polyn(&eh1, &au2, q);
                                            let topright =
                                                lwe::add_polyn(&lwe::mult_polyn(&a, &u2, q), &e, q);
                                            let topright = lwe::mod_polyn(topright, &fx, q);
                                            //botleft
                                            let a = lwe::add_polyn(&ah2, &au1, q);
                                            let e = lwe::add_polyn(&eh2, &eu1, q);
                                            let botleft =
                                                lwe::add_polyn(&lwe::mult_polyn(&a, &u1, q), &e, q);
                                            let botleft = lwe::mod_polyn(botleft, &fx, q);
                                            //botright
                                            let a = lwe::add_polyn(&ah2, &au2, q);
                                            let e = lwe::add_polyn(&eh2, &eu2, q);
                                            let botright =
                                                lwe::add_polyn(&lwe::mult_polyn(&a, &u2, q), &e, q);
                                            let botright = lwe::mod_polyn(botright, &fx, q);

                                            let item =
                                                Matrix::new(topleft, topright, botleft, botright);
                                            let count = map.entry(item).or_insert(0);
                                            *count += 1;
                                        });
                                    });
                                });
                            });
                        });
                    });
                });
            });
        });
    });

    // for (key, val) in map.iter() {
    //     print!("{}", val);
    //     if val < map {
    //         print!(", ");
    //     }
    // }
    let mut it = map.iter().peekable();
    let mut buf = String::new();
    let mut lab = String::new();
    //write!(&mut buf, "[");
    while let Some((key, val)) = it.next() {
        write!(&mut buf, "{}", val);
        write!(&mut lab, "{}", key);
        if !it.peek().is_none() {
            write!(&mut buf, " ");
            write!(&mut lab, " ");
        }
    }

    //write!(&mut buf, "]");
    let path = "/Users/rensweerman/PycharmProjects/pythonProject/values.txt";
    fs::write(path, buf).expect("Unable to write file");
    let path = "/Users/rensweerman/PycharmProjects/pythonProject/labels.txt";
    fs::write(path, lab).expect("Unable to write file");
}

pub fn test3() {
    let mut count = 0;
    let fx = vec![1, 0, 1];
    let mut map = HashMap::new();
    //8924237
    let u1 = lwe::uniform_sample(degree, q);
    let u2 = lwe::uniform_sample(degree, q);
    let au1 = lwe::uniform_sample(degree, q);
    let au2 = lwe::uniform_sample(degree, q);
    let eu1 = vec![]; //lwe::uniform_sample(degree, q);
    let eu2 = vec![]; //lwe::uniform_sample(degree, q);
    Poly::new(degree, q).for_each(|u1| {
        println!("1: {}", count);
        Poly::new(degree, q).for_each(|u2| {
            println!("2: {}", count);
            Poly::new(degree, q).for_each(|au1| {
                //println!("3: {}", count);
                Poly::new(degree, q).for_each(|au2| {
                    Poly::new(degree, q).for_each(|ah1| {
                        Poly::new(degree, q).for_each(|ah2| {
                            let eh1 = vec![];
                            let eh2 = vec![];
                            //Poly::new(degree, q).for_each(|eh1| {
                            //Poly::new(degree, q).for_each(|eh2| {
                            count += 1;
                            //topleft
                            let a = lwe::add_polyn(&ah1, &au1, q);
                            let e = lwe::add_polyn(&eh1, &eu1, q);
                            let topleft = lwe::add_polyn(&lwe::mult_polyn(&a, &u1, q), &e, q);
                            let topleft = lwe::mod_polyn(topleft, &fx, q);
                            //print_vec(&topleft);
                            //topright
                            let a = lwe::add_polyn(&ah1, &au2, q);
                            let e = lwe::add_polyn(&eh1, &au2, q);
                            let topright = lwe::add_polyn(&lwe::mult_polyn(&a, &u2, q), &e, q);
                            let topright = lwe::mod_polyn(topright, &fx, q);
                            //botleft
                            let a = lwe::add_polyn(&ah2, &au1, q);
                            let e = lwe::add_polyn(&eh2, &eu1, q);
                            let botleft = lwe::add_polyn(&lwe::mult_polyn(&a, &u1, q), &e, q);
                            let botleft = lwe::mod_polyn(botleft, &fx, q);
                            //botright
                            let a = lwe::add_polyn(&ah2, &au2, q);
                            let e = lwe::add_polyn(&eh2, &eu2, q);
                            let botright = lwe::add_polyn(&lwe::mult_polyn(&a, &u2, q), &e, q);
                            let botright = lwe::mod_polyn(botright, &fx, q);

                            let item = Matrix::new(topleft, topright, botleft, botright);
                            let count = map.entry(item).or_insert(0);
                            *count += 1;
                        });
                    });
                });
            });
        });
    });
    //});
    //});

    // for (key, val) in map.iter() {
    //     print!("{}", val);
    //     if val < map {
    //         print!(", ");
    //     }
    // }
    let mut it = map.iter().peekable();
    let mut buf = String::new();
    let mut lab = String::new();
    //write!(&mut buf, "[");
    while let Some((key, val)) = it.next() {
        write!(&mut buf, "{}", val);
        write!(&mut lab, "{}", key);
        if !it.peek().is_none() {
            write!(&mut buf, " ");
            write!(&mut lab, " ");
        }
    }

    //write!(&mut buf, "]");
    let path = "/Users/rensweerman/PycharmProjects/pythonProject/values.txt";
    fs::write(path, buf).expect("Unable to write file");
    let path = "/Users/rensweerman/PycharmProjects/pythonProject/labels.txt";
    fs::write(path, lab).expect("Unable to write file");

    println!("{}", count);
}

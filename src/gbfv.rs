use std::{
    cell::RefCell,
    fmt::Display,
    ops::{BitAnd, Div, Sub},
    rc::Rc,
};

use crate::{
    polynomial::{Polynomial, Value},
    rlwe::Params,
};

fn encode_m<Val: Value>(m: &Polynomial<Val>, p: &Params<Val>) -> Polynomial<Val> {
    //naive:
    let m_0 = m.val.get(0).unwrap_or(&Val::zero());
    let p_0 = m.val.get(0).unwrap_or(&Val::zero());
    return Polynomial::new(vec![Val::one()], p);
}

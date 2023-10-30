use std::ops::{Add, Div, Mul, Sub};

#[inline]
pub fn add_scalar<T: Add<Output = T> + Copy>(dx: &[T], da: T) -> Vec<T> {
    let mut r = Vec::<T>::new();
    let n = dx.len();
    let m = n % 4;
    if m != 0 {
        for i in 0..m {
            r.push(dx[i] + da);
        }
        if n < 4 {
            return r;
        }
    }
    for i in (m..n).step_by(4) {
        r.push(dx[i] + da);
        r.push(dx[i + 1] + da);
        r.push(dx[i + 2] + da);
        r.push(dx[i + 3] + da);
    }
    r
}

#[inline]
pub fn mul_scalar<T: Mul<Output = T> + Copy>(dx: &[T], da: T) -> Vec<T> {
    let mut r = Vec::<T>::new();
    let n = dx.len();
    let m = n % 4;
    if m != 0 {
        for i in 0..m {
            r.push(dx[i] * da);
        }
        if n < 4 {
            return r;
        }
    }
    for i in (m..n).step_by(4) {
        r.push(dx[i] * da);
        r.push(dx[i + 1] * da);
        r.push(dx[i + 2] * da);
        r.push(dx[i + 3] * da);
    }
    r
}

#[inline]
pub fn sub_scalar<T: Sub<Output = T> + Copy>(dx: &[T], da: T) -> Vec<T> {
    let mut r = Vec::<T>::new();
    let n = dx.len();
    let m = n % 4;
    if m != 0 {
        for i in 0..m {
            r.push(dx[i] - da);
        }
        if n < 4 {
            return r;
        }
    }
    for i in (m..n).step_by(4) {
        r.push(dx[i] - da);
        r.push(dx[i + 1] - da);
        r.push(dx[i + 2] - da);
        r.push(dx[i + 3] - da);
    }
    r
}

#[inline]
pub fn div_scalar<T: Div<Output = T> + Copy>(dx: &[T], da: T) -> Vec<T> {
    let mut r = Vec::<T>::new();
    let n = dx.len();
    let m = n % 4;
    if m != 0 {
        for i in 0..m {
            r.push(dx[i] / da);
        }
        if n < 4 {
            return r;
        }
    }
    for i in (m..n).step_by(4) {
        r.push(dx[i] / da);
        r.push(dx[i + 1] / da);
        r.push(dx[i + 2] / da);
        r.push(dx[i + 3] / da);
    }
    r
}

#[inline]
pub fn add<T: Add<Output = T> + Copy>(dx: &[T], dy: &[T]) -> Vec<T> {
    assert!(dx.len() == dy.len(), "Can't broadcast between dx and dy");
    let mut r = Vec::<T>::new();
    let n = dx.len();
    let m = n % 4;
    if m != 0 {
        for i in 0..m {
            r.push(dx[i] + dy[i]);
        }
        if n < 4 {
            return r;
        }
    }
    for i in (m..n).step_by(4) {
        r.push(dx[i] + dy[i]);
        r.push(dx[i + 1] + dy[i + 1]);
        r.push(dx[i + 2] + dy[i + 2]);
        r.push(dx[i + 3] + dy[i + 2]);
    }
    r
}

#[inline]
pub fn sub<T: Sub<Output = T> + Copy>(dx: &[T], dy: &[T]) -> Vec<T> {
    assert!(dx.len() == dy.len(), "Can't broadcast between dx and dy");
    let mut r = Vec::<T>::new();
    let n = dx.len();
    let m = n % 4;
    if m != 0 {
        for i in 0..m {
            r.push(dx[i] - dy[i]);
        }
        if n < 4 {
            return r;
        }
    }
    for i in (m..n).step_by(4) {
        r.push(dx[i] - dy[i]);
        r.push(dx[i + 1] - dy[i + 1]);
        r.push(dx[i + 2] - dy[i + 2]);
        r.push(dx[i + 3] - dy[i + 2]);
    }
    r
}

#[inline]
pub fn mul<T: Mul<Output = T> + Copy>(dx: &[T], dy: &[T]) -> Vec<T> {
    assert!(dx.len() == dy.len(), "Can't broadcast between dx and dy");
    let mut r = Vec::<T>::new();
    let n = dx.len();
    let m = n % 4;
    if m != 0 {
        for i in 0..m {
            r.push(dx[i] * dy[i]);
        }
        if n < 4 {
            return r;
        }
    }
    for i in (m..n).step_by(4) {
        r.push(dx[i] * dy[i]);
        r.push(dx[i + 1] * dy[i + 1]);
        r.push(dx[i + 2] * dy[i + 2]);
        r.push(dx[i + 3] * dy[i + 2]);
    }
    r
}

#[inline]
pub fn div<T: Div<Output = T> + Copy>(dx: &[T], dy: &[T]) -> Vec<T> {
    assert!(dx.len() == dy.len(), "Can't broadcast between dx and dy");
    let mut r = Vec::<T>::new();
    let n = dx.len();
    let m = n % 4;
    if m != 0 {
        for i in 0..m {
            r.push(dx[i] / dy[i]);
        }
        if n < 4 {
            return r;
        }
    }
    for i in (m..n).step_by(4) {
        r.push(dx[i] / dy[i]);
        r.push(dx[i + 1] / dy[i + 1]);
        r.push(dx[i + 2] / dy[i + 2]);
        r.push(dx[i + 3] / dy[i + 2]);
    }
    r
}

#[inline]
pub fn absf(dx: &[f64]) -> Vec<f64> {
    let mut r = Vec::<f64>::new();
    let n = dx.len();
    let m = n % 4;
    if m != 0 {
        for i in 0..m {
            r.push(dx[i].abs());
        }
        if n < 4 {
            return r;
        }
    }
    for i in (m..n).step_by(4) {
        r.push(dx[i].abs());
        r.push(dx[i + 1].abs());
        r.push(dx[i + 2].abs());
        r.push(dx[i + 3].abs());
    }
    r
}

#[inline]
pub fn absi(dx: &[i32]) -> Vec<i32> {
    let mut r = Vec::<i32>::new();
    let n = dx.len();
    let m = n % 4;
    if m != 0 {
        for i in 0..m {
            r.push(dx[i].abs());
        }
        if n < 4 {
            return r;
        }
    }
    for i in (m..n).step_by(4) {
        r.push(dx[i].abs());
        r.push(dx[i + 1].abs());
        r.push(dx[i + 2].abs());
        r.push(dx[i + 3].abs());
    }
    r
}

#[inline]
pub fn sumf(dx: &[f64]) -> f64 {
    let mut r = 0f64;
    let n = dx.len();
    let m = n % 4;
    if m != 0 {
        for i in 0..m {
            r += dx[i];
        }
        if n < 4 {
            return r;
        }
    }
    for i in (m..n).step_by(4) {
        r += dx[i] + dx[i + 1] + dx[i + 2] + dx[i + 3];
    }
    r
}

pub fn sumi(dx: &[i32]) -> i32 {
    let mut r = 0i32;
    let n = dx.len();
    let m = n % 4;
    if m != 0 {
        for i in 0..m {
            r += dx[i];
        }
        if n < 4 {
            return r;
        }
    }
    for i in (m..n).step_by(4) {
        r += dx[i] + dx[i + 1] + dx[i + 2] + dx[i + 3];
    }
    r
}

#[inline]
pub fn powi(dx: &[f64], p: i32) -> Vec<f64> {
    let mut r = Vec::<f64>::new();
    let n = dx.len();
    let m = n % 4;
    if m != 0 {
        for i in 0..m {
            r.push(dx[i].powi(p));
        }
        if n < 4 {
            return r;
        }
    }
    for i in (m..n).step_by(4) {
        r.push(dx[i].powi(p));
        r.push(dx[i + 1].powi(p));
        r.push(dx[i + 2].powi(p));
        r.push(dx[i + 3].powi(p));
    }
    r
}

#[inline]
pub fn powf(dx: &[f64], p: f64) -> Vec<f64> {
    let mut r = Vec::<f64>::new();
    let n = dx.len();
    let m = n % 4;
    if m != 0 {
        for i in 0..m {
            r.push(dx[i].powf(p));
        }
        if n < 4 {
            return r;
        }
    }
    for i in (m..n).step_by(4) {
        r.push(dx[i].powf(p));
        r.push(dx[i + 1].powf(p));
        r.push(dx[i + 2].powf(p));
        r.push(dx[i + 3].powf(p));
    }
    r
}

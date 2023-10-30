use std::iter::zip;
#[no_mangle]
#[inline(always)]
pub fn ddot(n: usize, dx: &[f64], dy: &[f64], incx: usize, incy: usize) -> f64 {
    // assert!(n + 1 <= dx.len() && n + 1 <= dy.len(), "index out of bound");
    let mut r = 0f64;
    let m: usize;
    let mut ix: usize;
    let mut iy: usize;

    if incx == 0 && incy == 0 {
        m = n % 5;
        if m != 0 {
            for i in 1..m + 1 {
                r += dx[i] * dy[i];
            }
            if n < 5 {
                return r;
            }
        }

        for i in (m + 1..n + 1).step_by(5) {
            r += dx[i] * dy[i]
                + dx[i + 1] * dy[i + 1]
                + dx[i + 2] * dy[i + 2]
                + dx[i + 3] * dy[i + 3]
                + dx[i + 4] * dy[i + 4]
        }
    } else {
        ix = 1;
        iy = 1;
        for _ in 1..n + 1 {
            r += dx[ix] * dy[iy];
            ix += incx;
            iy += incy;
        }
    }
    r
}

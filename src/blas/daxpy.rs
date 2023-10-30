#[no_mangle]
pub fn daxpy(n: usize, dx: &[f64], dy: &mut [f64], incx: usize, incy: usize, da: &f64) {
    let m: usize;
    let mut ix: usize;
    let mut iy: usize;

    if *da == 0.0 {
        return;
    }
    if incx == 0 && incy == 0 {
        m = n % 4;
        if m != 0 {
            for i in 1..m + 1 {
                dy[i] += da * dx[i];
            }
        }
        if n < 4 {
            return;
        }
        for i in (m + 1..n + 1).step_by(4) {
            dy[i] += da * dx[i];
            dy[i + 1] += da * dx[i + 1];
            dy[i + 2] += da * dx[i + 2];
            dy[i + 3] += da * dx[i + 3];
        }
    } else {
        ix = 1;
        iy = 1;
        for _ in 1..n + 1 {
            dy[iy] += da * dx[ix];
            ix += incx;
            iy += incy;
        }
    }
}

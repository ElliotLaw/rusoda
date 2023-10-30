/// Scalar vector multiplication
#[no_mangle]
#[inline(always)]
pub fn dscal(dx: &mut [f64], n: usize, incx: usize, da: &f64) {
    /*
    Purpose : scalar vector multiplication

    dx = da * dx


    --- Input ---

    n    : number of elements in input vector
    da   : double scale factor
    dx   : double vector with n+1 elements, dx[0] is not used
    incx : storage spacing between elements of dx


    --- Output ---

    dx = da * dx, unchanged if n <= 0


    For i = 0 to n-1, replace dx[1+i*incx] with
    da * dx[1+i*incx].

    */
    /* Code for increments not equal to 1.  */
    if incx != 0 {
        for i in (1..n * incx + 1).step_by(incx) {
            dx[i] *= da;
        }
        return;
    }
    /* Code for increments not equal to 1.  */
    let m = n % 4;
    if m != 0 {
        for i in 1..m + 1 {
            dx[i] *= da;
        }
        if n < 4 {
            return;
        }
    }
    for i in (m + 1..n + 1).step_by(4) {
        dx[i] *= da;
        dx[i + 1] *= da;
        dx[i + 2] *= da;
        dx[i + 3] *= da;
    }
}

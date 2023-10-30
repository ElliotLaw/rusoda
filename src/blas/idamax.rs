#[no_mangle]
#[inline(always)]
pub fn idamax(dx: &[f64], n: usize, incx: usize) -> usize {
    /*
    [2023-10-17]

    Impl of linpack idamax.f by Elliot Law.

    Purpose : Find largest component of double vector dx


    --- Input ---

    n    : number of elements in input vector
    dx   : double vector with n+1 elements, dx[0] is not used
    incx : storage spacing between elements of dx


    --- Output ---

    idamax : smallest index, 0 if n <= 0


    Find smallest index of maximum magnitude of dx.
    idamax = first i, i=1 to n, to minimize fabs( dx[1-incx+i*incx] ).

    [Author]
    Univ. of Tennessee
    Univ. of California Berkeley
    Univ. of Colorado Denver
    NAG Ltd.

    */
    let mut idmax = 1usize;
    let mut dmax: f64;
    let mut ix: usize;

    if n == 1 {
        return idmax;
    }

    if incx == 0 {
        dmax = dx[1].abs();
        for i in 2..n + 1 {
            if dx[i].abs() > dmax {
                idmax = i;
                dmax = dx[i].abs();
            }
        }
    } else {
        ix = 1;
        dmax = dx[1].abs();
        ix += incx;
        for i in 2..n + 1 {
            if dx[ix].abs() > dmax {
                idmax = i;
                dmax = dx[ix].abs();
            }
            ix += incx;
        }
    }
    idmax
}

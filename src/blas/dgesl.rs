use super::{daxpy::daxpy, ddot::ddot};

/// Solves the double precision system A * x = b  or  trans(A) * x = b using the factors computed by dgeco or dgefa.
#[no_mangle]
#[inline(always)]
pub fn dgesl(a: &mut Vec<Vec<f64>>, n: usize, ipvt: &[usize], b: &mut [f64], job: usize) {
    /*
    [2023-10-17]

    Impl of linpack dgesl.f by Elliot Law.

    [Purpose]
    DGESL solves the double precision system
    A * x = b  or  trans(A) * x = b
    using the factors computed by dgeco or dgefa.

    [Author]
    Cleve Moler,
    University of New Mexico,
    Argonne National Lab.

    */
    let mut t: f64;
    let mut j: usize;
    /*
       Job = 0, solve a * x = b.
    */
    if job == 0 {
        /*
           First solve L * y = b.
        */
        for k in 1..n + 1 {
            t = ddot(k - 1, &a[k], b, 1, 1);
            b[k] = (b[k] - t) / a[k][k];
        }
        /*
           Now solve U * x = y.
        */
        for k in (1..n).rev() {
            (*b)[k] += ddot(n - k, &a[k][k..], &b[k..], 1, 1);
            j = ipvt[k];
            if j != k {
                t = b[j];
                b[j] = b[k];
                b[k] = t;
            }
        }
        return;
    }
    /*
       Job = nonzero, solve Transpose(a) * x = b.

       First solve Transpose(U) * y = b.
    */
    for k in 1..n {
        j = ipvt[k];
        t = b[j];
        if j != k {
            b[j] = b[k];
            b[k] = t;
        }
        daxpy(n - k, &mut a[k][k..], &mut b[k..], 1, 1, &t);
    }
    /*
       Now solve Transpose(L) * x = y.
    */
    for k in (1..n + 1).rev() {
        b[k] /= a[k][k];
        t = -b[k];
        daxpy(k - 1, &mut a[k], b, 1, 1, &t);
    }
}

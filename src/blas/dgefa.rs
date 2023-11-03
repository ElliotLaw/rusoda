use super::{dscal::dscal, idamax::idamax};

/// Factors a double matrix by Gaussian elimination.
#[no_mangle]
#[inline(always)]
pub fn dgefa(a: &mut Vec<Vec<f64>>, n: usize, ipvt: &mut [usize], info: &mut usize) {
    /*
       Purpose : dgefa factors a double matrix by Gaussian elimination.

       dgefa is usually called by dgeco, but it can be called directly
       with a saving in time if rcond is not needed.
       (Time for dgeco) = (1+9/n)*(time for dgefa).

       This c version uses algorithm kji rather than the kij in dgefa.f.
       Note that the fortran version input variable lda is not needed.


       On Entry :

          a   : double matrix of dimension ( n+1, n+1 ),
                the 0-th row and column are not used.
                a is created using NewDoubleMatrix, hence
                lda is unnecessary.
          n   : the row dimension of a.

       On Return :

          a     : a lower triangular matrix and the multipliers
                  which were used to obtain it.  The factorization
                  can be written a = L * U where U is a product of
                  permutation and unit upper triangular matrices
                  and L is lower triangular.
          ipvt  : an n+1 integer vector of pivot indices.
          *info : = 0 normal value,
                  = k if U[k][k] == 0.  This is not an error
                    condition for this subroutine, but it does
                    indicate that dgesl or dgedi will divide by
                    zero if called.  Use rcond in dgeco for
                    a reliable indication of singularity.

                    Notice that the calling program must use &info.

       BLAS : daxpy, dscal, idamax
    */
    let mut t: f64;
    let mut j: usize;
    *info = 0;
    for k in 1..n {
        /*
           Find j = pivot index.  Note that a[k]+k-1 is the address of
           the 0-th element of the row vector whose 1st element is a[k][k].
        */
        j = idamax(&a[k][k - 1..], n - k + 1, 1) + k - 1;
        ipvt[k] = j;
        /*
           Zero pivot implies this row already triangularized.
        */
        if a[k][j] == 0. {
            *info = k;
            continue;
        }
        /*
           Interchange if necessary.
        */
        if j != k {
            t = a[k][j];
            a[k][j] = a[k][k];
            a[k][k] = t;
        }
        /*
           Compute multipliers.
        */
        t = -1f64 / a[k][k];
        dscal(&mut a[k][k..], n - k, 1, &t);
        /*
           Column elimination with row indexing.
        */
        for i in k + 1..n + 1 {
            t = a[i][j];
            if j != k {
                a[i][j] = a[i][k];
                a[i][k] = t;
            }
            daxpy_inner(n - k, a, k, i, k, &t);
        }
    }
    ipvt[n] = n;
    if a[n][n] == 0f64 {
        *info = n;
    }
}

fn daxpy_inner(n: usize, a: &mut Vec<Vec<f64>>, idx: usize, idy: usize, k_: usize, da: &f64) {
    /*daxpy for dgefa */
    let m: usize;

    if *da == 0.0 {
        return;
    }

    m = n % 4;
    if m != 0 {
        for i in 1 + k_..m + 1 + k_ {
            a[idy][i] += da * a[idx][i];
        }
    }
    if n < 4 {
        return;
    }
    for i in (m + 1 + k_..n + 1 + k_).step_by(4) {
        a[idy][i] += da * a[idx][i];
        a[idy][i + 1] += da * a[idx][i + 1];
        a[idy][i + 2] += da * a[idx][i + 2];
        a[idy][i + 3] += da * a[idx][i + 3];
    }
}

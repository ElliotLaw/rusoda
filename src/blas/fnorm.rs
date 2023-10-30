#[no_mangle]
pub fn fnorm(n: usize, a: &Vec<Vec<f64>>, w: &[f64]) -> f64 {
    /*
       This subroutine computes the norm of a full n by n matrix,
       stored in the array a, that is consistent with the weighted max-norm
       on vectors, with weights stored in the array w.

        fnorm = max(i=1,...,n) ( w[i] * sum(j=1,...,n) fabs( a[i][j] ) / w[j] )
    */
    let mut an = 0.;
    let mut sum: f64;
    let mut ap1: &[f64];
    for i in 1..n + 1 {
        sum = 0.;
        ap1 = &a[i];
        for j in 1..n + 1 {
            sum += ap1[j].abs() / w[j];
        }
        if sum * w[i] > an {
            an = sum * w[i];
        }
    }
    an
}



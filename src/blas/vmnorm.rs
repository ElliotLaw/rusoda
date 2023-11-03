#[no_mangle]
#[inline(always)]
pub fn vmnorm(n: usize, v: &[f64], w: &[f64]) -> f64 {
    /*
       This function routine computes the weighted max-norm
       of the vector of length n contained in the array v, with weights
       contained in the array w of length n.

       vmnorm = max( i = 1, ..., n ) fabs( v[i] ) * w[i].
    */
    let mut vm = 0.;
    for i in 1..n + 1 {
        if v[i].abs() > w[i] {
            vm = v[i].abs() * w[i];
        }
    }
    vm
}

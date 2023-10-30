use rand::Rng;

pub fn ddot_test() {
    let mut rng = rand::thread_rng();
    let n = 10000;
    let mut dx = Vec::<f64>::new();
    let mut dy = Vec::<f64>::new();
    for _ in 0..n + 1 {
        dx.push(rng.gen::<f64>());
        dy.push(rng.gen::<f64>());
    }

    use crate::blas::ddot::ddot;
    use std::time::Instant;

    let tt = Instant::now();
    let rr = ddot(n, &dx, &dy, 0, 0);
    println!("{}, {}ns", rr, tt.elapsed().as_nanos());
}

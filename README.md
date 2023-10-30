# rusoda

#### Description

rust implementation of odepack dlsoda

#### Instructions

1.  Build a struct which contains your ode system coeffs and number of eqs.

```
struct Oral1Cpt {
    incalc_par: Vec<f64>,
    neq: usize,
}
```

2.  Implement OdeSystem trait for your struct(define the ode function).

```
use rusoda::OdeSystem;

impl OdeSystem for Oral1Cpt {
    fn func(&self, _t: f64, _y: &mut [f64], _dy: &mut [f64]) {
        let ka = self.incalc_par[0];
        let cl = self.incalc_par[1];
        let v = self.incalc_par[2];
        (*_dy)[0] = -ka * (*_y)[0];
        (*_dy)[1] = ka * (*_y)[0] - cl / v * (*_y)[1];
    }
}
```

3.  Initialize LSODA solver and call "solve" method.

```
use rusoda::IStateInput;
use rusoda::LSODA;

use std::time::Instant;
    env_logger::init();
    let mut t = 0.;
    let tout = 12.;
    let y = [4.0, 0.0];
    let mut lsoda = LSODA::init();
    let sys = Oral1Cpt::init(vec![3.09, 32., 648.]);

    let mut state = IStateInput::InitialCall;
    let tt = Instant::now();

    let res = lsoda.solve(
        &sys, sys.neq, &y, &mut t, tout, &mut state, 1e-3, 1e-6, false,
    );

    println!("{:?},TIME:{}MS", res, tt.elapsed().as_millis())
```

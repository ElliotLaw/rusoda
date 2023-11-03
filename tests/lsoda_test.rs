extern crate rusoda;
use rusoda::IStateInput;
use rusoda::{OdeSystem, LSODA};

struct Oral1Cpt {
    incalc_par: Vec<f64>,
    neq: usize,
}

impl Oral1Cpt {
    fn init(par: Vec<f64>) -> Self {
        Self {
            incalc_par: par.clone(),
            neq: 2,
        }
    }
}

impl OdeSystem for Oral1Cpt {
    fn func(&self, _t: f64, _y: &mut [f64], _dy: &mut [f64]) {
        let ka = self.incalc_par[0];
        let cl = self.incalc_par[1];
        let v = self.incalc_par[2];
        (*_dy)[0] = -ka * (*_y)[0];
        (*_dy)[1] = ka * (*_y)[0] - cl / v * (*_y)[1];
    }
    fn reverse_func(&self, _t: f64, _y: &mut [f64], _dy: &mut [f64]) {
        let ka = self.incalc_par[0];
        let cl = self.incalc_par[1];
        let v = self.incalc_par[2];
        (*_dy)[0] = ka * (*_y)[0];
        (*_dy)[1] = -ka * (*_y)[0] + cl / v * (*_y)[1];
    }
    fn get_neq(&self) -> usize {
        self.neq
    }
}

pub fn oral_1cpt_test() {
    use std::time::Instant;
    env_logger::init();
    let mut t = 0.;
    let tout = 12.;
    let y = [4.0, 0.0];
    let mut lsoda = LSODA::init();
    let sys = Oral1Cpt::init(vec![3.09, 32., 648.]);

    let mut state = IStateInput::InitialCall;
    let tt = Instant::now();

    let res = lsoda.solve(&sys, &y, &mut t, tout, &mut state, 1e-3, 1e-6, false, false);

    println!("{:?},TIME:{}MS", res, tt.elapsed().as_millis())
}

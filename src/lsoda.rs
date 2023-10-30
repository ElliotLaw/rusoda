use crate::blas::prelude::*;
use crate::common::*;

use log::{error, info, warn};

// use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq)]
enum OrderFlag {
    NoChangeinHNq,
    HChanged,
    BothHNqChanged,
}
impl Default for OrderFlag {
    fn default() -> Self {
        Self::NoChangeinHNq
    }
}

#[derive(Debug, Default)]
pub struct LSODA {
    lastoutcome: IStateOutput,
    ml: usize,
    mu: usize,
    imxer: usize,
    sqrteta: f64,
    mord: [usize; 3],
    sm1: [f64; 13],
    el: [f64; 14],
    cm1: [f64; 13],
    cm2: [f64; 6],
    elco: [[f64; 14]; 13],
    tesco: [[f64; 4]; 13],
    illin: usize,
    init: InitializationState,
    ierpj: usize,
    iersl: usize,
    jcur: usize,
    l: usize,
    miter: CorrectionIterMethod,
    maxord: usize,
    maxcor: usize,
    msbp: usize,
    mxncf: usize,
    kflag: KFlag,
    jstart: JStart,
    ixpr: Ixpr,
    jtyp: CorrectionIterMethod,
    mused: OdeMethod,
    mxordn: usize,
    mxords: usize,
    meth_: OdeMethod,
    n: usize,
    nq: usize,
    nst: usize,
    nfe: usize,
    nje: usize,
    nqu: usize,
    mxstep: usize,
    mxhnil: usize,
    nslast: usize,
    nhnil: usize,
    ntrep: usize,
    nyh: usize,
    ccmax: f64,
    el0: f64,
    h_: f64,
    hmin: f64,
    hmxi: f64,
    hu: f64,
    rc: f64,
    tn_: f64,
    tsw: f64,
    pdnorm: f64,
    conit: f64,
    crate_: f64,
    hold: f64,
    rmax: f64,
    ialth: usize,
    ipup: CorrectionIterMethod,
    lmax: usize,
    nslp: usize,
    pdest: f64,
    pdlast: f64,
    ratio: f64,
    icount: i32,
    irflag: i32,
    ewt: Vec<f64>,
    savf: Vec<f64>,
    acor: Vec<f64>,
    yh_: Vec<Vec<f64>>,
    wm_: Vec<Vec<f64>>,
    ipvt: Vec<usize>,
    itol_: ITol,
    rtol_: Vec<f64>,
    atol_: Vec<f64>,
}

impl LSODA {
    pub fn init() -> Self {
        Self {
            mxords: 12,
            h_: 0.0,
            tn_: 0.0,
            mord: [12, 5, 0],
            sm1: [
                0., 0.5, 0.575, 0.55, 0.45, 0.35, 0.25, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025,
            ],
            el: [0.0; 14],
            cm1: [0.0; 13],
            cm2: [0.0; 6],
            ccmax: 0.3,
            maxcor: 3,
            msbp: 20,
            mxncf: 10,
            ..Self::default()
        }
    }

    fn successreturn(&mut self, y: &mut [f64], t: &mut f64, itask: ITask, ihit: bool, tcrit: f64) {
        for i in 1..self.n + 1 {
            y[i] = self.yh_[1][i];
        }
        *t = self.tn_;
        match itask {
            ITask::NormalComputationWithoutTCrit | ITask::SingleStepWithoutTCrit => {
                if ihit {
                    *t = tcrit;
                }
            }
            _ => (),
        }
        self.lastoutcome = IStateOutput::IntegrationSuccess;
        self.illin = 0;
    }

    pub fn terminate(&mut self) {
        if self.illin == 5 {
            error!("[lsoda] repeated occurrence of illegal input. run aborted.. \n apparent infinite loop.")
        } else {
            self.illin += 1;
            self.lastoutcome = IStateOutput::IllegalInput;
        }
    }

    fn terminate_inner(&mut self, y: &mut Vec<f64>, t: &mut f64) {
        for i in 1..self.n + 1 {
            y[i] = self.yh_[1][i];
        }
        *t = self.tn_;
        self.illin = 0;
    }

    fn ewset(&mut self, ycur: &Vec<f64>) {
        match self.itol_ {
            ITol::SingleRtolSingleAtol => {
                for i in 1..self.n + 1 {
                    self.ewt[i] = self.rtol_[1] * ycur[i].abs() + self.atol_[1];
                }
            }
            ITol::SingleRtolMultiAtol => {
                for i in 1..self.n + 1 {
                    self.ewt[i] = self.rtol_[1] * ycur[i].abs() + self.atol_[i];
                }
            }
            ITol::MultiRtolSingleAtol => {
                for i in 1..self.n + 1 {
                    self.ewt[i] = self.rtol_[i] * ycur[i].abs() + self.atol_[1];
                }
            }
            ITol::MultiRtolMultiAtol => {
                for i in 1..self.n + 1 {
                    self.ewt[i] = self.rtol_[i] * ycur[i].abs() + self.atol_[i];
                }
            }
        }
    }

    fn resetcoef(&mut self) {
        /*
           The el vector and related constants are reset
           whenever the order nq is changed, or at the start of the problem.
        */
        for i in 1..self.l + 1 {
            self.el[i] = self.elco[self.nq][i];
        }
        self.rc *= self.el[1] / self.el0;
        self.el0 = self.el[1];
        self.conit = 0.5 / (self.nq + 2) as f64;
    }

    fn scaleh(&mut self, rh: &mut f64, pdh: &mut f64) {
        let mut r: f64;
        /*
           If h_ is being changed, the h_ ratio rh is checked against rmax, hmin,
           and hmxi, and the yh_ array is rescaled.  ialth is set to l = nq + 1
           to prevent a change of h_ for that many steps, unless forced by a
           convergence or error test failure.
        */
        *rh = rh.min(self.rmax);
        *rh /= 1f64.max(self.h_.abs() * self.hmxi * (*rh));
        /*
           If meth_ = 1, also restrict the new step size by the stability region.
           If this reduces h_, set irflag to 1 so that if there are roundoff
           problems later, we can assume that is the cause of the trouble.
        */
        if self.meth_ == OdeMethod::Adam {
            self.irflag = 0;
            *pdh = 0.000001_f64.max(self.h_.abs() * self.pdlast);
            if (*rh * *pdh * 1.00001) >= self.sm1[self.nq] {
                *rh = self.sm1[self.nq] / *pdh;
                self.irflag = 1;
            }
        }
        r = 1.;
        for j in 2..self.l + 1 {
            r *= *rh;
            for i in 1..self.n + 1 {
                self.yh_[j][i] *= r;
            }
        }
        self.h_ *= *rh;
        self.rc *= *rh;
        self.ialth = self.l;
    }

    fn corfailure(&mut self, told: &mut f64, rh: &mut f64, ncf: &mut usize, corflag: &mut usize) {
        *ncf += 1;
        self.rmax = 2.;
        self.tn_ = *told;
        for j in (1..self.nq + 1).rev() {
            for i1 in j..self.nq + 1 {
                for i in 1..self.n + 1 {
                    self.yh_[i1][i] -= self.yh_[i1 + 1][i];
                }
            }
        }
        if self.h_.abs() <= self.hmin * 1.00001 || *ncf == self.mxncf {
            *corflag = 2;
            return;
        }
        *corflag = 1;
        *rh = 0.25;
        self.ipup = self.miter;
    }

    #[allow(unused_assignments)]
    fn correction(
        &mut self,
        neq: &usize,
        y: &mut [f64],
        system: &dyn OdeSystem,
        corflag: &mut usize,
        pnorm: f64,
        del: &mut f64,
        delp: &mut f64,
        told: &mut f64,
        ncf: &mut usize,
        rh: &mut f64,
        m: &mut usize,
    ) {
        let mut rm: f64 = 0.;
        let mut rate: f64 = 0.;
        let mut dcon: f64 = 0.;
        /*
           Up to maxcor corrector iterations are taken.  A convergence test is
           made on the r.m.s. norm of each correction, weighted by the error
           weight vector ewt.  The sum of the corrections is accumulated in the
           vector acor[i].  The yh_ array is not altered in the corrector loop.
        */
        *m = 0;
        *corflag = 0;
        *del = 0.;
        for i in 1..self.n + 1 {
            y[i] = self.yh_[1][i];
        }
        system.func(self.tn_, &mut y[1..], &mut self.savf[1..]);
        self.nfe += 1;
        /*
           If indicated, the matrix P = I - h_ * el[1] * J is reevaluated and
           preprocessed before starting the corrector iteration.  ipup is set
           to 0 as an indicator that this has been done.
        */
        loop {
            if *m == 0 {
                match self.ipup {
                    CorrectionIterMethod::NoJac => (),
                    _ => {
                        self.prja(*neq, y, system);
                        self.ipup = CorrectionIterMethod::NoJac;
                        self.rc = 1.;
                        self.nslp = self.nst;
                        self.crate_ = 0.7;
                        if self.ierpj != 0 {
                            self.corfailure(told, rh, ncf, corflag);
                            return;
                        }
                    }
                }

                for i in 1..self.n + 1 {
                    self.acor[i] = 0.;
                }
            }
            if self.miter == CorrectionIterMethod::NoJac {
                /*
                   In case of functional iteration, update y directly from
                   the result of the last function evaluation.
                */
                for i in 1..self.n + 1 {
                    self.savf[i] = self.h_ * self.savf[i] - self.yh_[2][i];
                    y[i] = self.savf[i] - self.acor[i];
                }
                *del = vmnorm(self.n, y, &self.ewt);
                for i in 1..self.n + 1 {
                    y[i] = self.yh_[1][i] + self.el[1] * self.savf[i];
                    self.acor[i] = self.savf[i];
                }
            } else {
                for i in 1..self.n + 1 {
                    y[i] = self.h_ * self.savf[i] - self.yh_[2][i] - self.acor[i];
                }
                self.solsy(y);
                *del = vmnorm(self.n, y, &self.ewt);
                for i in 1..self.n + 1 {
                    self.acor[i] += y[i];
                    y[i] = self.yh_[1][i] + self.el[1] * self.acor[i];
                }
            }
            /*
               Test for convergence.  If *m > 0, an estimate of the convergence
               rate constant is stored in crate, and this is used in the test.

               We first check for a change of iterates that is the size of
               roundoff error.  If this occurs, the iteration has converged, and a
               new rate estimate is not formed.
               In all other cases, force at least two iterations to estimate a
               local Lipschitz constant estimate for Adams method.
               On convergence, form pdest = local maximum Lipschitz constant
               estimate.  pdlast is the most recent nonzero estimate.
            */
            if *del <= 100. * pnorm * ETA {
                break;
            }
            match self.meth_ {
                OdeMethod::BDF => {
                    if *m != 0 {
                        rm = 1024.;
                        if *del <= *delp * 1024. {
                            rm = *del / (*delp);
                        }
                        rate = rate.max(rm);
                        self.crate_ = rm.max(0.2 * self.crate_);
                    }
                    dcon =
                        *del * 1f64.min(1.5 * self.crate_) / (self.tesco[self.nq][2] * self.conit);
                    if dcon <= 1. {
                        self.pdest = self.pdest.max(rate / (self.h_ * self.el[1]).abs());
                        if self.pdest != 0. {
                            self.pdlast = self.pdest;
                        }
                        break;
                    }
                }
                _ => (),
            }
            /*
               The corrector iteration failed to converge.
               If miter != 0 and the Jacobian is out of date, prja is called for
               the next try.   Otherwise the yh_ array is retracted to its values
               before prediction, and h_ is reduced, if possible.  If h_ cannot be
               reduced or mxncf failures have occured, exit with corflag = 2.
            */
            *m += 1;
            if *m == self.maxcor || (*m >= 2 && *del > 2. * *delp) {
                if self.miter == CorrectionIterMethod::NoJac || self.jcur == 1 {
                    self.corfailure(told, rh, ncf, corflag);
                    return;
                }
                self.ipup = self.miter;
                /*
                   Restart corrector if Jacobian is recomputed.
                */
                *m = 0;
                rate = 0.;
                *del = 0.;
                for i in 1..self.n + 1 {
                    y[i] = self.yh_[1][i];
                }
                system.func(self.tn_, &mut y[1..], &mut self.savf[1..]);
                self.nfe += 1;
            } else {
                *delp = *del;
                system.func(self.tn_, &mut y[1..], &mut self.savf[1..]);
                self.nfe += 1;
            }
        }
    }

    #[allow(unused_assignments)]
    fn methodswitch(&mut self, dsm: f64, pnorm: f64, pdh: &mut f64, rh: &mut f64) {
        let lm1: usize;
        let lm1p1: usize;
        let lm2: usize;
        let lm2p1: usize;
        let nqm1: usize;
        let nqm2: usize;
        let mut rh1: f64;
        let rh2: f64;
        let mut rh1it: f64;
        let exm2: f64;
        let dm2: f64;
        let exm1: f64;
        let mut dm1: f64;
        let alpha: f64;
        let exsm: f64;
        /*
           We are current using an Adams method.  Consider switching to bdf.
           If the current order is greater than 5, assume the problem is
           not stiff, and skip this section.
           If the Lipschitz constant and error estimate are not polluted
           by roundoff, perform the usual test.
           Otherwise, switch to the bdf methods if the last step was
           restricted to insure stability ( irflag = 1 ), and stay with Adams
           method if not.  When switching to bdf with polluted error estimates,
           in the absence of other information, double the step size.

           When the estimates are ok, we make the usual test by computing
           the step size we could have (ideally) used on this step,
           with the current (Adams) method, and also that for the bdf.
           If nq > mxords, we consider changing to order mxords on switching.
           Compare the two step sizes to decide whether to switch.
           The step size advantage must be at least ratio = 5 to switch.
        */
        if self.meth_ == OdeMethod::Adam {
            if self.nq > 5 {
                return;
            }
            if dsm <= 100. * pnorm * ETA || self.pdest == 0. {
                match self.irflag {
                    0 => return,
                    _ => (),
                }
                rh2 = 2.;
                nqm2 = self.nq.min(self.mxords);
            } else {
                exsm = 1. / self.l as f64;
                rh1 = 1. / (1.2 * dsm.powf(exsm) + 0.0000012);
                rh1it = 2. * rh1;
                *pdh = self.pdlast * self.h_.abs();
                if *pdh * rh1 > 0.00001 {
                    rh1it = self.sm1[self.nq] / (*pdh);
                }
                rh1 = rh1.min(rh1it);
                if self.nq > self.mxords {
                    nqm2 = self.mxords;
                    lm2 = self.mxords + 1;
                    exm2 = 1. / lm2 as f64;
                    lm2p1 = lm2 + 1;
                    dm2 = vmnorm(self.n, &self.yh_[lm2p1], &self.ewt);
                    rh2 = 1. / (1.2 * dm2.powf(exm2) + 0.0000012)
                } else {
                    dm2 = dsm * (self.cm1[self.nq] / self.cm2[self.nq]);
                    rh2 = 1. / (1.2 * dm2.powf(exsm) + 0.0000012);
                    nqm2 = self.nq;
                }
                if rh2 < self.ratio * rh1 {
                    return;
                }
            }
            /*
               The method switch test passed.  Reset relevant quantities for bdf.
            */
            *rh = rh2;
            self.icount = 20;
            self.meth_ = OdeMethod::BDF;
            self.miter = self.jtyp;
            self.pdlast = 0.;
            self.nq = nqm2;
            self.l = self.nq + 1;
            return;
        }
        /*
           We are currently using a bdf method, considering switching to Adams.
           Compute the step size we could have (ideally) used on this step,
           with the current (bdf) method, and also that for the Adams.
           If nq > mxordn, we consider changing to order mxordn on switching.
           Compare the two step sizes to decide whether to switch.
           The step size advantage must be at least 5/ratio = 1 to switch.
           If the step size for Adams would be so small as to cause
           roundoff pollution, we stay with bdf.
        */
        exsm = 1. / self.l as f64;
        if self.mxordn < self.nq {
            nqm1 = self.mxordn;
            lm1 = self.mxordn + 1;
            exm1 = 1. / lm1 as f64;
            lm1p1 = lm1 + 1;
            dm1 = vmnorm(self.n, &self.yh_[lm1p1], &self.ewt) / self.cm1[self.mxordn];
            rh1 = 1. / (1.2 * dm1.powf(exm1) + 0.0000012);
        } else {
            dm1 = dsm * (self.cm1[self.nq] / self.cm2[self.nq]);
            rh1 = 1. / (1.2 * dm1.powf(exsm) + 0.0000012);
            nqm1 = self.nq;
            exm1 = exsm;
        }
        rh1it = 2. * rh1;
        *pdh = self.pdnorm * self.h_.abs();
        if *pdh * rh1 > 0.00001 {
            rh1it = self.sm1[nqm1] / (*pdh);
        }
        rh1 = rh1.min(rh1it);
        rh2 = 1. / (1.2 * dsm.powf(exsm) + 0.0000012);
        if rh1 * self.ratio < 5. * rh2 {
            return;
        }
        alpha = rh1.max(0.001);
        dm1 *= alpha.powf(rh1);
        if dm1 <= 1000. * ETA * pnorm {
            return;
        }
        /*
           The switch test passed.  Reset relevant quantities for Adams.
        */
        *rh = rh1;
        self.icount = 20;
        self.meth_ = OdeMethod::Adam;
        self.miter = CorrectionIterMethod::NoJac;
        self.pdlast = 0.;
        self.nq = nqm1;
        self.l = self.nq + 1;
    }

    fn endstoda(&mut self) {
        let r: f64 = 1. / self.tesco[self.nqu][2];
        for i in 1..self.n + 1 {
            self.acor[i] *= r;
        }
        self.hold = self.h_;
        self.jstart = JStart::ContinueStep;
    }

    #[allow(unused_assignments)]
    fn orderswitch(
        &mut self,
        rhup: &mut f64,
        dsm: f64,
        pdh: &mut f64,
        rh: &mut f64,
        orderflag: &mut OrderFlag,
    ) {
        let mut newq = 0usize;
        let exsm: f64;
        let mut rhdn: f64;
        let mut rhsm: f64;
        let ddn: f64;
        let exdn: f64;
        let r: f64;
        *orderflag = OrderFlag::NoChangeinHNq;
        exsm = 1. / self.l as f64;
        rhsm = 1. / (1.2 * dsm.powf(exsm) + 0.0000012);
        rhdn = 0.;
        if self.nq != 1 {
            ddn = vmnorm(self.n, &self.yh_[self.l], &self.ewt);
            exdn = 1. / self.nq as f64;
            rhdn = 1. / (1.3 * ddn.powf(exdn) + 0.0000013);
        }
        /*
           If meth_ is Adam, limit rh accordinfg to the stability region also.
        */
        match self.meth_ {
            OdeMethod::Adam => {
                *pdh = (self.h_.abs() * self.pdlast).max(0.000001);
                if self.l < self.lmax {
                    *rhup = rhup.min(self.sm1[self.l] / (*pdh));
                }
                rhsm = rhsm.min(self.sm1[self.nq] / (*pdh));
                if self.nq > 1 {
                    rhdn = rhdn.min(self.sm1[self.nq - 1] / (*pdh));
                }
                self.pdest = 0.;
            }
            _ => (),
        }
        if rhsm >= *rhup {
            if rhsm >= rhdn {
                newq = self.nq;
                *rh = rhsm;
            } else {
                newq = self.nq - 1;
                *rh = rhdn;
                match self.kflag {
                    KFlag::StepSuccess => (),
                    _ => {
                        if *rh > 1. {
                            *rh = 1.;
                        }
                    }
                }
            }
        } else {
            if *rhup <= rhdn {
                newq = self.nq - 1;
                *rh = rhdn;
                match self.kflag {
                    KFlag::StepSuccess => (),
                    _ => {
                        if *rh > 1. {
                            *rh = 1.;
                        }
                    }
                }
            } else {
                *rh = *rhup;
                if *rh >= 1.1 {
                    r = self.el[self.l] / self.l as f64;
                    self.nq = self.l;
                    self.l = self.nq + 1;
                    for i in 1..self.n + 1 {
                        self.yh_[self.l][i] = self.acor[i] * r;
                    }
                    *orderflag = OrderFlag::BothHNqChanged;
                    return;
                } else {
                    self.ialth = 3;
                    return;
                }
            }
        }
        /*
           If meth_ = Adam and h_ is restricted by stability, bypass 10 percent test.
        */
        match self.meth_ {
            OdeMethod::Adam => {
                if *rh * *pdh * 1.00001 < self.sm1[newq] {
                    match self.kflag {
                        KFlag::StepSuccess => {
                            if *rh < 1.1 {
                                self.ialth = 3;
                                return;
                            }
                        }
                        _ => (),
                    }
                }
            }
            _ => match self.kflag {
                KFlag::StepSuccess => {
                    if *rh < 1.1 {
                        self.ialth = 3;
                        return;
                    }
                }
                _ => (),
            },
        }
        match self.kflag {
            KFlag::CouldNotAchieveConvergence | KFlag::FatalinPJACoSLVS => *rh = rh.min(0.2),
            _ => (),
        }
        /*
           If there is a change of order, reset nq, l, and the coefficients.
           In any case h_ is reset according to rh and the yh_ array is rescaled.
           Then exit or redo the step.
        */
        if newq == self.nq {
            *orderflag = OrderFlag::HChanged;
            return;
        }
        self.nq = newq;
        self.l = self.nq + 1;
        *orderflag = OrderFlag::BothHNqChanged;
    }

    fn solsy(&mut self, y: &mut [f64]) {
        self.iersl = 0;
        if self.miter != CorrectionIterMethod::InternalFullJac {
            warn!("[LSODA] solsy -- miter != 2\n");
            return;
        } else {
            dgesl(&mut self.wm_, self.n, &self.ipvt, y, 0);
        }
    }

    #[allow(unused_assignments)]
    fn intdy(&self, t: f64, k: usize, dky: &mut [f64], iflag: &mut IFlag) {
        let mut ic = 0usize;
        let mut jp1 = 0usize;
        let r: f64;
        let s: f64;
        let tp: f64;
        *iflag = IFlag::KTLegal;
        if k > self.nq {
            error!("[intdy] k = {} illegal\n", k);
            *iflag = IFlag::IllegalK;
            return;
        }
        tp = self.tn_ - self.hu - 100. * ETA * (self.tn_ + self.hu);
        if (t - tp) * (t - self.tn_) > 0. {
            error!(
                "[intdy] t = {} illegal.t not in interval tcur - hu to tcur\n",
                t
            );
            *iflag = IFlag::IllegalT;
            return;
        }
        s = (t - self.tn_) / self.h_;
        ic = 1;
        for jj in (self.l - k)..self.nq + 1 {
            ic *= jj;
        }
        for i in 1..self.n + 1 {
            dky[i] = ic as f64 * self.yh_[self.l][i];
        }
        for j in (k..self.nq).rev() {
            jp1 = j + 1;
            ic = 1;
            for jj in (jp1 - k)..j + 1 {
                ic *= jj;
            }
            for i in 1..self.n + 1 {
                dky[i] = ic as f64 * self.yh_[jp1][i] + s * dky[i];
            }
        }
        if k == 0 {
            return;
        }
        r = self.h_.powi(-(k as i32));
        for i in 1..self.n + 1 {
            dky[i] *= r;
        }
    }

    #[allow(unused_assignments)]
    fn prja(&mut self, _neq: usize, y: &mut [f64], system: &dyn OdeSystem) {
        // let mut i = 0usize;
        let mut ier = 0usize;
        // let mut j = 0usize;
        let mut fac = 0.0;
        let mut hl0 = 0.0;
        let mut r = 0.0;
        let mut r0 = 0.0;
        let mut yj = 0.0;
        /*
           prja is called by stoda to compute and process the matrix
           P = I - h_ * el[1] * J, where J is an approximation to the Jacobian.
           Here J is computed by finite differencing.
           J, scaled by -h_ * el[1], is stored in wm_.  Then the norm of J ( the
           matrix norm consistent with the weighted max-norm on vectors given
           by vmnorm ) is computed, and J is overwritten by P.  P is then
           subjected to LU decomposition in preparation for later solution
           of linear systems with p as coefficient matrix.  This is done
           by dgefa if miter = 2, and by dgbfa if miter = 5.
        */
        self.nje += 1;
        self.ierpj = 0;
        self.jcur = 1;
        hl0 = self.h_ * self.el0;
        /*
           If miter = 2, make n calls to f to approximate J.
        */
        if self.miter == CorrectionIterMethod::InternalFullJac {
            fac = vmnorm(self.n, &self.savf, &self.ewt);
            r0 = 1000.0 * self.h_.abs() * ETA * self.n as f64 * fac;
            if r0 == 0.0 {
                r0 = 1.0;
            }
            for j in 1..self.n + 1 {
                yj = y[j];
                r = (self.sqrteta * yj.abs()).max(r0 / self.ewt[j]);
                y[j] += r;
                system.func(self.tn_, &mut y[1..], &mut self.acor[1..]);
                for i in 1..self.n + 1 {
                    self.wm_[i][j] = (self.acor[i] - self.savf[i]) * fac;
                }
                y[j] = yj;
            }
            self.nfe += self.n;
            /*
               Compute norm of Jacobian.
            */
            self.pdnorm = fnorm(self.n, &self.wm_, &self.ewt) / hl0.abs();
            /*
               Add identity matrix.
            */
            for i in 1..self.n + 1 {
                self.wm_[i][i] += 1.;
            }
            /*
               Do LU decomposition on P.
            */
            dgefa(&mut self.wm_, self.n, &mut self.ipvt, &mut ier);
        }
    }

    fn cfode(&mut self, method: OdeMethod) {
        /*
           cfode is called by the integrator routine to set coefficients
           needed there.  The coefficients for the current method, as
           given by the value of meth, are set for all orders and saved.
           The maximum order assumed here is 12 if meth = 1 and 5 if meth = 2.
           ( A smaller value of the maximum order is also allowed. )
           cfode is called once at the beginning of the problem, and
           is not called again unless and until meth is changed.

           The _C(elco) array contains the basic method coefficients.
           The coefficients el[i], 1 < i < nq+1, for the method of
           order nq are stored in _C(elco)[nq][i].  They are given by a generating
           polynomial, i.e.,

              l(x) = el[1] + el[2]*x + ... + el[nq+1]*x^nq.

           For the implicit Adams method, l(x) is given by

              dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),   l(-1) = 0.

           For the bdf methods, l(x) is given by

              l(x) = (x+1)*(x+2)*...*(x+nq)/k,

           where   k = factorial(nq)*(1+1/2+...+1/nq).

           The _C(tesco) array contains test constants used for the
           local error test and the selection of step size and/or order.
           At order nq, _C(tesco)[nq][k] is used for the selection of step
           size at order nq-1 if k = 1, at order nq if k = 2, and at order
           nq+1 if k = 3.
        */
        let mut nqm1: usize;
        let mut nqp1: usize;
        let mut agamq: f64;
        let mut fnq: f64;
        let mut fnqm1: f64;
        let mut pc: [f64; 13] = [0.; 13];
        let mut pint: f64;
        let mut ragq: f64;
        let mut rqfac: f64;
        let mut rq1fac: f64;
        let mut tsign: f64;
        let mut xpin: f64;
        if method == OdeMethod::Adam {
            self.elco[1][1] = 1.;
            self.elco[1][2] = 1.;
            self.tesco[1][1] = 0.;
            self.tesco[1][2] = 2.;
            self.tesco[2][1] = 1.;
            self.tesco[12][3] = 0.;
            pc[1] = 1.;
            rqfac = 1.;
            for nq in 2..13 {
                /*
                   The pc array will contain the coefficients of the polynomial

                      p(x) = (x+1)*(x+2)*...*(x+nq-1).

                   Initially, p(x) = 1.
                */
                rq1fac = rqfac;
                rqfac /= nq as f64;
                nqm1 = nq - 1;
                fnqm1 = nqm1 as f64;
                nqp1 = nq + 1;
                /*
                   Form coefficients of p(x)*(x+nq-1).
                */
                pc[nq] = 0.;
                for i in (2..nq + 1).rev() {
                    pc[i] = pc[i - 1] + fnqm1 * pc[i];
                }
                pc[1] *= fnqm1;
                /*
                   Compute integral, -1 to 0, of p(x) and x*p(x).
                */
                pint = pc[1];
                xpin = pc[1] / 2.;
                tsign = 1.;
                for i in 2..nq + 1 {
                    tsign = -tsign;
                    pint += tsign * pc[i] / i as f64;
                    xpin += tsign * pc[i] / (i + 1) as f64;
                }
                /*
                   Store coefficients in elco and tesco.
                */
                self.elco[nq][1] = pint * rq1fac;
                self.elco[nq][2] = 1.;
                for i in 2..nq + 1 {
                    self.elco[nq][i + 1] = rq1fac * pc[i] / i as f64;
                }
                agamq = rqfac * xpin;
                ragq = 1. / agamq;
                self.tesco[nq][2] = ragq;
                if nq < 12 {
                    self.tesco[nqp1][1] = ragq * rqfac / nqp1 as f64;
                }
                self.tesco[nqm1][3] = ragq;
            }
            return;
        }
        /* meth_ = BDF */
        pc[1] = 1.;
        rq1fac = 1.;
        /*
           The pc array will contain the coefficients of the polynomial
              p(x) = (x+1)*(x+2)*...*(x+nq).
           Initially, p(x) = 1.
        */
        for nq in 1..6 {
            fnq = nq as f64;
            nqp1 = nq + 1;
            /*
               Form coefficients of p(x)*(x+nq).
            */
            pc[nqp1] = 0.;
            for i in (2..nq + 2).rev() {
                pc[i] = pc[i - 1] + fnq * pc[i];
            }
            pc[1] *= fnq;
            /*
               Store coefficients in elco and tesco.
            */
            for i in 1..nqp1 + 1 {
                self.elco[nq][i] = pc[i] / pc[2];
            }
            self.elco[nq][2] = 1.;
            self.tesco[nq][1] = rq1fac;
            self.tesco[nq][2] = nqp1 as f64 / self.elco[nq][1];
            self.tesco[nq][3] = (nq + 2) as f64 / self.elco[nq][1];
            rq1fac /= fnq;
        }
        return;
    }

    #[allow(unused_assignments)]
    fn stoda(&mut self, neq: usize, y: &mut [f64], system: &dyn OdeSystem) {
        assert!(neq + 1 == y.len());
        let mut corflag = 0usize;
        let mut orderflag = OrderFlag::NoChangeinHNq;
        let mut m = 0usize;
        let mut ncf = 0usize;
        let mut del = 0.0;
        let mut delp = 0.0;
        let mut dsm = 0.0;
        let mut dup = 0.0;
        let mut exup = 0.0;
        let mut r = 0.0;
        let mut rh = 0.0;
        let mut rhup = 0.0;
        let mut told = 0.0;
        let mut pdh = 0.0;
        let mut pnorm = 0.0;
        /*
           stoda performs one step of the integration of an initial value
           problem for a system of ordinary differential equations.
           Note.. stoda is independent of the value of the iteration method
           indicator miter, when this is != 0, and hence is independent
           of the type of chord method used, or the Jacobian structure.
           Communication with stoda is done with the following variables:

           jstart = enum for input only, with the following
                    values and meanings:

                    FirstStep:         perform the first step,
                    ContinueStep:      take a new step continuing from the last,
                    NewHNextStep:      take the next step with a new value of h_,
                                       n, meth_, miter, and/or matrix parameters.
                    NewParamNextStep:  take the next step with a new value of h_,
                                       but with other inputs unchanged.

           kflag = a completion enum with the following meanings:

                    StepSuccess                 the step was successful,
                    CouldNotAchieveReqErr       the requested error could not
                                                be achieved,
                    CouldNotAchieveConvergence  corrector convergence could not
                                                be achieved,
                    FatalinPJACoSLVS            fatal error in prja or solsy.

           miter = corrector iteration method:

                     0  functional iteration,
                    >0  a chord method corresponding to jacobian type jt.

        */
        told = self.tn_;
        self.ierpj = 0;
        self.iersl = 0;
        self.jcur = 0;
        /*
           On the first call, the order is set to 1, and other variables are
           initialized.  rmax is the maximum ratio by which h_ can be increased
           in a single step.  It is initially 1.e4 to compensate for the small
           initial h_, but then is normally equal to 10.  If a filure occurs
           (in corrector convergence or error test), rmax is set at 2 for
           the next increase.
           cfode is called to get the needed coefficients for both methods.
        */
        match self.jstart {
            JStart::FirstStep => {
                self.lmax = self.maxord + 1;
                self.nq = 1;
                self.l = 2;
                self.ialth = 2;
                self.rmax = 10000.;
                self.rc = 0.;
                self.el0 = 1.;
                self.crate_ = 0.7;
                self.hold = self.h_;
                self.nslp = 0;
                self.ipup = self.miter;
                /*
                   Initialize switching parameters.  meth_ = 1 is assumed initially.
                */
                self.icount = 20;
                self.irflag = 0;
                self.pdest = 0.;
                self.pdlast = 0.;
                self.ratio = 5.;
                self.cfode(OdeMethod::BDF);
                for i in 1..6 {
                    self.cm2[i] = self.tesco[i][2] * self.elco[i][i + 1];
                }
                self.cfode(OdeMethod::Adam);
                for i in 1..13 {
                    self.cm1[i] = self.tesco[i][2] * self.elco[i][i + 1];
                }
                self.resetcoef();
            }
            /*
               The following block handles preliminaries needed when jstart = -1.
               ipup is set to miter to force a matrix update.
               If an order increase is about to be considered ( ialth = 1 ),
               ialth is reset to 2 to postpone consideration one more step.
               If the caller has changed meth_, cfode is called to reset
               the coefficients of the method.
               If h_ is to be changed, yh_ must be rescaled.
               If h_ or meth_ is being changed, ialth is reset to l = nq + 1
               to prevent further changes in h_ for that many steps.
            */
            JStart::ContinueStep => {
                self.ipup = self.miter;
                self.lmax = self.maxord + 1;
                if self.ialth == 1 {
                    self.ialth = 2;
                }
                if self.meth_ != self.mused {
                    self.cfode(self.meth_);
                    self.ialth = self.l;
                    self.resetcoef();
                }
                if self.h_ != self.hold {
                    rh = self.h_ / self.hold;
                    self.h_ = self.hold;
                    self.scaleh(&mut rh, &mut pdh);
                }
            }
            JStart::NewHNextStep => {
                if self.h_ != self.hold {
                    rh = self.h_ / self.hold;
                    self.h_ = self.hold;
                    self.scaleh(&mut rh, &mut pdh);
                }
            }
            JStart::NewParamNextStep => (),
        }
        /*
           Prediction.
           This section computes the predicted values by effectively
           multiplying the yh_ array by the pascal triangle matrix.
           rc is the ratio of new to old values of the coefficient h_ * el[1].
           When rc differs from 1 by more than ccmax, ipup is set to miter
           to force pjac to be called, if a jacobian is involved.
           In any case, prja is called at least every msbp steps.
        */
        loop {
            loop {
                if self.rc - 1. > self.ccmax || self.nst >= self.nslp + self.msbp {
                    self.ipup = self.miter;
                }
                self.tn_ += self.h_;
                for j in (1..self.nq + 1).rev() {
                    for i1 in j..self.nq + 1 {
                        for i in 1..self.n + 1 {
                            self.yh_[i1][i] += self.yh_[i1 + 1][i];
                        }
                    }
                }
                pnorm = vmnorm(self.n, &self.yh_[1], &self.ewt);
                self.correction(
                    &neq,
                    y,
                    system,
                    &mut corflag,
                    pnorm,
                    &mut del,
                    &mut delp,
                    &mut told,
                    &mut ncf,
                    &mut rh,
                    &mut m,
                );
                if corflag == 0 {
                    break;
                } else if corflag == 1 {
                    rh = rh.max(self.hmin / self.h_.abs());
                    self.scaleh(&mut rh, &mut pdh);
                    continue;
                } else {
                    self.kflag = KFlag::CouldNotAchieveConvergence;
                    self.hold = self.h_;
                    self.jstart = JStart::ContinueStep;
                    return;
                }
            }
            /*
               The corrector has converged.  jcur is set to 0
               to signal that the Jacobian involved may need updating later.
               The local error test is done now.
            */
            self.jcur = 0;
            if m == 0 {
                dsm = del / self.tesco[self.nq][2];
            } else if m > 0 {
                dsm = vmnorm(self.n, &self.acor, &self.ewt);
            }
            if dsm <= 1. {
                /*
                   After a successful step, update the yh_ array.
                   Decrease icount by 1, and if it is -1, consider switching methods.
                   If a method switch is made, reset various parameters,
                   rescale the yh_ array, and exit.  If there is no switch,
                   consider changing h_ if ialth = 1.  Otherwise decrease ialth by 1.
                   If ialth is then 1 and nq < maxord, then acor is saved for
                   use in a possible order increase on the next step.
                   If a change in h_ is considered, an increase or decrease in order
                   by one is considered also.  A change in h_ is made only if it is by
                   a factor of at least 1.1.  If not, ialth is set to 3 to prevent
                   testing for that many steps.
                */
                self.kflag = KFlag::StepSuccess;
                self.nst += 1;
                self.hu = self.h_;
                self.nqu = self.nq;
                self.mused = self.meth_;
                for j in 1..self.l + 1 {
                    r = self.el[j];
                    for i in 1..self.n + 1 {
                        self.yh_[j][i] += r * self.acor[i];
                    }
                }
                self.icount -= 1;
                if self.icount < 0 {
                    self.methodswitch(dsm, pnorm, &mut pdh, &mut rh);
                    if self.meth_ != self.mused {
                        rh = rh.max(self.hmin / self.h_.abs());
                        self.scaleh(&mut rh, &mut pdh);
                        self.rmax = 10.;
                        self.endstoda();
                        break;
                    }
                }
                /*
                   No method switch is being made.  Do the usual step/order selection.
                */
                self.ialth -= 1;
                if self.ialth == 0 {
                    rhup = 0.;
                    if self.l != self.lmax {
                        for i in 1..self.n + 1 {
                            self.savf[i] = self.acor[i] - self.yh_[self.lmax][i];
                        }
                        dup = vmnorm(self.n, &self.savf, &self.ewt);
                        exup = 1. / (self.l + 1) as f64;
                        rhup = 1. / (1.4 * dup.powf(exup) + 0.0000014);
                    }
                    self.orderswitch(&mut rhup, dsm, &mut pdh, &mut rh, &mut orderflag);
                    match orderflag {
                        OrderFlag::NoChangeinHNq => {
                            /*
                               No change in h_ or nq.
                            */
                            self.endstoda();
                            break;
                        }
                        OrderFlag::HChanged => {
                            /*
                               h_ is changed, but not nq.
                            */
                            rh = rh.max(self.hmin / self.h_.abs());
                            self.scaleh(&mut rh, &mut pdh);
                            self.rmax = 10.;
                            self.endstoda();
                            break;
                        }
                        OrderFlag::BothHNqChanged => {
                            /*
                               both nq and h_ are changed.
                            */
                            self.resetcoef();
                            rh = rh.max(self.hmin / self.h_.abs());
                            self.scaleh(&mut rh, &mut pdh);
                            self.rmax = 10.;
                            self.endstoda();
                            break;
                        }
                    }
                }
                if self.ialth > 1 || self.l == self.lmax {
                    self.endstoda();
                    break;
                }
                for i in 1..self.n + 1 {
                    self.yh_[self.lmax][i] = self.acor[i];
                }
                self.endstoda();
                break;
            }
            /*
              The error test failed.  kflag keeps track of multiple failures.
              Restore tn_ and the yh_ array to their previous values, and prepare
              to try the step again.  Compute the optimum step size for this or
              one lower.  After 2 or more failures, h_ is forced to decrease
              by a factor of 0.2 or less.
            */
            else {
                match self.kflag {
                    KFlag::StepSuccess => self.kflag = KFlag::CouldNotAchieveReqError,
                    KFlag::CouldNotAchieveReqError => {
                        self.kflag = KFlag::CouldNotAchieveConvergence
                    }
                    KFlag::CouldNotAchieveConvergence => self.kflag = KFlag::FatalinPJACoSLVS,
                    KFlag::FatalinPJACoSLVS => self.kflag = KFlag::MultiFails(4),
                    KFlag::MultiFails(x) => self.kflag = KFlag::MultiFails(x + 1),
                }
                self.tn_ = told;
                for j in (1..self.nq + 1).rev() {
                    for i1 in j..self.nq + 1 {
                        for i in 1..self.n + 1 {
                            self.yh_[i1][i] -= self.yh_[i1 + 1][i];
                        }
                    }
                }
                self.rmax = 2.;
                if self.h_.abs() <= self.hmin * 1.00001 {
                    self.kflag = KFlag::CouldNotAchieveReqError;
                    self.hold = self.h_;
                    self.jstart = JStart::ContinueStep;
                    break;
                }
                match self.kflag {
                    KFlag::StepSuccess
                    | KFlag::CouldNotAchieveReqError
                    | KFlag::CouldNotAchieveConvergence => {
                        rhup = 0.;
                        self.orderswitch(&mut rhup, dsm, &mut pdh, &mut rh, &mut orderflag);
                        match orderflag {
                            OrderFlag::NoChangeinHNq | OrderFlag::HChanged => {
                                if orderflag == OrderFlag::NoChangeinHNq {
                                    rh = rh.min(0.2);
                                }
                                rh = rh.max(self.hmin / self.h_.abs());
                                self.scaleh(&mut rh, &mut pdh);
                            }
                            _ => {
                                self.resetcoef();
                                rh = rh.max(self.hmin / self.h_.abs());
                                self.scaleh(&mut rh, &mut pdh);
                            }
                        }
                        continue;
                    }
                    KFlag::FatalinPJACoSLVS => {
                        rh = 0.1;
                        rh = rh.min(self.hmin / self.h_.abs());
                        self.h_ *= rh;
                        for i in 1..self.n + 1 {
                            y[i] = self.yh_[1][i];
                        }
                        system.func(self.tn_, &mut y[1..], &mut self.savf[1..]);
                        self.nfe += 1;
                        for i in 1..self.n + 1 {
                            self.yh_[2][i] = self.h_ * self.savf[i];
                        }
                        self.ipup = self.miter;
                        self.ialth = 5;
                        if self.nq == 1 {
                            continue;
                        }
                        self.nq = 1;
                        self.l = 2;
                        self.resetcoef();
                        continue;
                    }
                    KFlag::MultiFails(x) => {
                        if x == 10 {
                            self.kflag = KFlag::CouldNotAchieveReqError;
                            self.hold = self.h_;
                            self.jstart = JStart::ContinueStep;
                            break;
                        } else {
                            rh = 0.1;
                            rh = rh.min(self.hmin / self.h_.abs());
                            self.h_ *= rh;
                            for i in 1..self.n + 1 {
                                y[i] = self.yh_[1][i];
                            }
                            system.func(self.tn_, &mut y[1..], &mut self.savf[1..]);
                            self.nfe += 1;
                            for i in 1..self.n + 1 {
                                self.yh_[2][i] = self.h_ * self.savf[i];
                            }
                            self.ipup = self.miter;
                            self.ialth = 5;
                            if self.nq == 1 {
                                continue;
                            }
                            self.nq = 1;
                            self.l = 2;
                            self.resetcoef();
                            continue;
                        }
                    }
                }
                /*
                  Control reaches this section if 3 or more failures have occurred.
                  If 10 failures have occurred, exit with kflag = -1.
                  It is assumed that the derivatives that have accumulated in the
                  yh_ array have errors of the wrong order.  Hence the first
                  derivative is recomputed, and the order is set to 1.  Then
                  h_ is reduced by a factor of 10, and the step is retried,
                  until it succeeds or h_ reaches hmin.
                */
            }
        }
    }

    #[allow(unused_assignments)]
    pub fn step(
        &mut self,
        system: &dyn OdeSystem,
        neq: usize,
        y: &mut Vec<f64>,
        t: &mut f64,
        tout: f64,
        itask: ITask,
        istatein: &mut IStateInput,
        iopt: IOpt,
        jt: CorrectionIterMethod,
        iworks: IWork,
        rworks: Rworks,
    ) {
        // println!("{:?}", y);
        if tout > *t {
        } else {
            self.terminate();
            error!("[LSODA] tout = {} < t = {}", tout, t);
            return;
        }

        let mut iflag = IFlag::KTLegal;
        let mut lenyh = 0;
        let mut ihit = false;
        let mut atoli = 0.0;
        let mut ayi = 0.0;
        let mut big = 0.0;
        let mut h0 = 0.0;
        let mut hmax = 0.0;
        let mut hmx = 0.0;
        let mut rh = 0.0;
        let mut rtoli = 0.0;
        let mut tcrit = 0.0;
        let mut tdist = 0.0;
        let mut tnext = 0.0;
        let mut tol = 0.0;
        let mut tolsf = 0.0;
        let mut tp = 0.0;
        let mut size = 0.0;
        let mut sum = 0.0;
        let mut w0 = 0.0;
        /*
           Block a.
           This code block is executed on every call.
           If this is the initial call but the flag init shows that initialization has not
           yet been done, an error return occurs.
           If this is the initial call and tout = t, return immediately(Do Nothing).
        */
        match istatein {
            IStateInput::ParamChange | IStateInput::ParamNotChange => match self.init {
                InitializationState::NotYet => {
                    self.terminate();
                    error!("[lsoda] calling solve or recall but lsoda not initialized\n");
                    return;
                }
                _ => (),
            },
            _ => (),
        };

        /*
           Block b.
           The next code block is executed for the initial call,
           or for a continuation call with parameter changes.
           It contains checking of all inputs and various initializations.

           First check legality of the non-optional inputs neq, itol, iopt,
           jt, ml, and mu.
        */
        match istatein {
            IStateInput::InitialCall | IStateInput::ParamChange => {
                self.ntrep = 0;

                if neq <= 0 {
                    self.terminate();
                    error!("[lsoda] neq is less than 1, giving {}", neq);
                    return;
                }

                match istatein {
                    IStateInput::ParamChange => {
                        if neq > self.n {
                            self.terminate();
                            error!("[lsoda] param changed and neq increased");
                            return;
                        }
                    }
                    _ => (),
                };
                self.n = neq;
                self.jtyp = jt;
                match jt {
                    CorrectionIterMethod::UserBandJac | CorrectionIterMethod::InternalBandJac => {
                        if iworks.0 >= self.n {
                            self.terminate();
                            error!("[lsoda] ml not between 1 and neq");
                            return;
                        }

                        if iworks.1 >= self.n {
                            self.terminate();
                            error!("[lsoda] mu not between 1 and neq");
                            return;
                        }
                        self.ml = iworks.0;
                        self.mu = iworks.1;
                    }
                    _ => (),
                }
                /* Next process and check the optional inpus.   */
                match iopt {
                    /* Default options.*/
                    IOpt::NoOptionalInputs => {
                        self.ixpr = Ixpr::NoExtraPrint;
                        self.mxstep = 500usize;
                        self.mxhnil = 10usize;
                        self.hmxi = 0.0;
                        self.hmin = 0.0;
                        match istatein {
                            IStateInput::InitialCall => {
                                h0 = 0.0;
                                self.mxordn = self.mord[0];
                                self.mxords = self.mord[1];
                            }
                            _ => (),
                        }
                    }
                    /* Optional inputs.*/
                    _ => {
                        self.ixpr = iworks.2;
                        self.mxstep = iworks.3;
                        self.mxhnil = iworks.4;
                        match istatein {
                            IStateInput::InitialCall => {
                                h0 = rworks.h0;
                                self.mxordn = iworks.5;
                                if self.mxordn == 0 {
                                    self.mxordn = 100usize;
                                }
                                self.mxordn = self.mxordn.min(self.mord[0]);
                                self.mxords = iworks.6;
                                if self.mxords == 0 {
                                    self.mxords = 100usize;
                                }
                                self.mxords = self.mxords.min(self.mord[1]);
                                if (tout - *t) * h0 < 0.0 {
                                    self.terminate();
                                    error!("[lsoda] tout = {} behind t = {}. integration direction is given by {}", tout, t, h0);
                                    return;
                                }
                            }
                            _ => (),
                        }
                        hmax = rworks.hmax;
                        if hmax < 0.0 {
                            self.terminate();
                            error!("[lsoda] rworks.hmax < 0.");
                            return;
                        }

                        self.hmxi = 0.0;
                        if hmax != 0.0 {
                            self.hmxi = 1.0 / hmax;
                        }
                        self.hmin = rworks.hmin;
                        if self.hmin < 0.0 {
                            self.terminate();
                            error!("[lsoda] rworks.hmin < 0.");
                            return;
                        }
                    }
                }
            }
            _ => (),
        };

        match istatein {
            /*
                for the initial call only
                It contains all remaining initializations, the initial call to f,
                and the calculation of the initial step size.
                The error weights in ewt are inverted after being loaded.
            */
            IStateInput::InitialCall => {
                /*
                    Allocate Memory.
                */
                self.sqrteta = SQRTETA;
                self.meth_ = OdeMethod::Adam;
                self.nyh = self.n;
                lenyh = 1 + self.mxordn.max(self.mxords);
                self.yh_.resize(lenyh + 1, vec![0.0; self.nyh + 1]);
                self.wm_.resize(self.nyh + 1, vec![0.0; self.nyh + 1]);
                self.ewt.resize(self.nyh + 1, 0.0);
                self.savf.resize(self.nyh + 1, 0.0);
                self.acor.resize(self.nyh + 1, 0.0);
                self.ipvt.resize(self.nyh + 1, 0);
            }
            _ => (),
        };

        match istatein {
            IStateInput::InitialCall | IStateInput::ParamChange => {
                rtoli = self.rtol_[1];
                atoli = self.atol_[1];
                for i in 1..self.n + 1 {
                    match self.itol_ {
                        ITol::MultiRtolMultiAtol => {
                            rtoli = self.rtol_[i];
                            atoli = self.atol_[i];
                        }
                        ITol::MultiRtolSingleAtol => {
                            rtoli = self.rtol_[i];
                        }
                        ITol::SingleRtolMultiAtol => {
                            atoli = self.atol_[i];
                        }
                        _ => (),
                    }
                }
                if rtoli < 0.0 {
                    self.terminate();
                    error!("[lsoda] rtol = {} is less than 0.\n", rtoli);
                    return;
                }
            }
            _ => (),
        };
        if *istatein == IStateInput::ParamChange {
            self.jstart = JStart::NewHNextStep;
        }

        match istatein {
            IStateInput::InitialCall => {
                self.tn_ = *t;
                self.tsw = *t;
                self.maxord = self.mxordn;
                match itask {
                    ITask::NormalComputationWithoutTCrit | ITask::SingleStepWithoutTCrit => {
                        tcrit = rworks.tcrit;
                        if (tcrit - tout) * (tout - *t) < 0.0 {
                            self.terminate();
                            error!("[lsoda] itask without tcrit and tcrit behind tout\n");
                            return;
                        }
                        if h0 != 0. && (*t + h0 - tcrit) * h0 > 0.0 {
                            h0 = tcrit - *t;
                        }
                    }
                    _ => (),
                }
                self.jstart = JStart::FirstStep;
                self.nhnil = 0;
                self.nst = 0;
                self.nje = 0;
                self.nslast = 0;
                self.hu = 0.;
                self.nqu = 0;
                self.mused = OdeMethod::Adam;
                self.miter = CorrectionIterMethod::NoJac;
                self.ccmax = 0.3;
                self.maxcor = 3;
                self.msbp = 20;
                self.mxncf = 10;
                /* Initial call to f.  */
                if self.yh_.len() != lenyh + 1 {
                    error!("length of yh_ is not right");
                    return;
                }
                if self.yh_[0].len() != self.nyh + 1 {
                    error!("length of yh_[0] is not right");
                    return;
                }
                system.func(*t, &mut y[1..], &mut self.yh_[2][1..]);

                self.nfe = 1;

                /* Load the initial value vector in yh_.  */
                for i in 1..self.n + 1 {
                    self.yh_[1][i] = y[i];
                }

                /* Load and invert the ewt array.  ( h_ is temporarily set to 1. ) */
                self.nq = 1;
                self.h_ = 1.0;
                self.ewset(y);
                for i in 1..self.n + 1 {
                    if self.ewt[i] <= 0.0 {
                        self.terminate_inner(y, t);
                        error!("[LSODA] ewt[{}] = {} <= 0.", i, self.ewt[i]);
                    }
                    self.ewt[i] = 1.0 / self.ewt[i];
                }
                /*
                   The coding below computes the step size, h0, to be attempted on the
                   first step, unless the user has supplied a value for this.
                   First check that tout - *t differs significantly from zero.
                   A scalar tolerance quantity tol is computed, as max(rtol[i])
                   if this is positive, or max(atol[i]/fabs(y[i])) otherwise, adjusted
                   so as to be between 100*ETA and 0.001.
                   Then the computed value h0 is given by

                      h0^(-2) = 1. / ( tol * w0^2 ) + tol * ( norm(f) )^2

                   where   w0     = max( abs(*t), abs(tout) ),
                           f      = the initial value of the vector f(t,y), and
                           norm() = the weighted vector norm used throughout, given by
                                    the vmnorm function routine, and weighted by the
                                    tolerances initially loaded into the ewt array.

                   The sign of h0 is inferred from the initial values of tout and *t.
                   abs(h0) is made < fabs(tout-*t) in any case.
                */
                if h0 == 0.0 {
                    tdist = (tout - *t).abs();
                    if t.abs() > tout.abs() {
                        w0 = *t;
                    } else {
                        w0 = tout;
                    }
                    if tdist < 2.0 * ETA * w0 {
                        self.terminate();
                        error!("[LSODA] tout too close to t to start integration\n ");
                        return;
                    }
                    tol = self.rtol_[1];
                    match self.itol_ {
                        ITol::MultiRtolSingleAtol | ITol::MultiRtolMultiAtol => {
                            for i in 2..self.n + 1 {
                                if tol < self.rtol_[i] {
                                    tol = self.rtol_[i];
                                }
                            }
                        }
                        _ => (),
                    }
                    if tol <= 0.0 {
                        atoli = self.atol_[1];
                        for i in 1..self.n + 1 {
                            match self.itol_ {
                                ITol::SingleRtolMultiAtol | ITol::MultiRtolMultiAtol => {
                                    atoli = self.atol_[i];
                                }
                                _ => (),
                            }
                            ayi = y[i].abs();
                            if ayi != 0.0 {
                                if tol < atoli / ayi {
                                    tol = atoli / ayi;
                                }
                            }
                        }
                    }
                    if tol < 100.0 * ETA {
                        tol = 100.0 * ETA;
                    }
                    if tol > 0.001 {
                        tol = 0.001;
                    }
                    sum = vmnorm(self.n, &self.yh_[2], &self.ewt);
                    sum = 1.0 / (tol * w0 * w0) + tol * sum * sum;
                    h0 = 1.0 / sum.sqrt();
                    h0 = h0.min(tdist);
                    if (tout - *t) < 0.0 {
                        h0 = -h0;
                    }
                } /* end if ( h0 == 0. )   */
                /*
                   Adjust h0 if necessary to meet hmax bound.
                */
                rh = h0.abs() * self.hmxi;
                if rh > 1.0 {
                    h0 /= rh;
                }
                /*
                   Load h_ with h0 and scale yh_[2] by h0.
                */

                self.h_ = h0;
                for i in 1..self.n + 1 {
                    self.yh_[2][i] *= h0;
                }
            } /* end initial call configuring  */
            _ => (),
        }

        match istatein {
            /*
               Block d.
               The next code block is for continuation calls only ( *istate = 2 or 3 )
               and is to check stop conditions before taking a step.
            */
            IStateInput::ParamChange | IStateInput::ParamNotChange => {
                self.nslast = self.nst;
                match itask {
                    ITask::NormalComputation => {
                        if (self.tn_ - tout) * self.h_ >= 0. {
                            self.intdy(tout, 0, y, &mut iflag);
                            match iflag {
                                IFlag::KTLegal => (),
                                _ => {
                                    error!("[lsoda] trouble from intdy, itask = NormalComputation, tout = {}\n",
                                tout);
                                    self.terminate();
                                    return;
                                }
                            }
                            *t = tout;
                            self.lastoutcome = IStateOutput::IntegrationSuccess;
                            self.illin = 0;
                            return;
                        }
                    }
                    ITask::SingleStep => (),
                    ITask::StopAtFirstMeshPoint => {
                        tp = self.tn_ - self.hu * (1. + 100. * ETA);
                        if (tp - tout) * self.h_ > 0. {
                            error!(
                                "[lsoda] itask = StopAtFirstMeshPoint and tout behind tcur - hu\n"
                            );
                            self.terminate();
                            return;
                        }
                        if (self.tn_ - tout) * self.h_ < 0. {
                            ()
                        } else {
                            self.successreturn(y, t, itask, ihit, tcrit);
                            return;
                        }
                    }
                    ITask::NormalComputationWithoutTCrit => {
                        tcrit = rworks.tcrit;
                        if (self.tn_ - tcrit) * self.h_ > 0. {
                            error!("[lsoda] itask = NormalComputationWithoutTCrit and tcrit behind tcur\n");
                            self.terminate();
                            return;
                        }
                        if (tcrit - tout) * self.h_ < 0. {
                            error!("[lsoda] itask = NormalComputationWithoutTCrit and tout behind tcrit\n");
                            self.terminate();
                            return;
                        }
                        if (self.tn_ - tout) * self.h_ >= 0. {
                            self.intdy(tout, 0, y, &mut iflag);
                            match iflag {
                                IFlag::KTLegal => (),
                                _ => {
                                    error!("[lsoda] trouble from intdy, itask = NormalComputationWithoutTCrit, tout = {}\n", tout);
                                    self.terminate();
                                    return;
                                }
                            }
                            *t = tout;
                            self.lastoutcome = IStateOutput::IntegrationSuccess;
                            self.illin = 0;
                            return;
                        }
                    }
                    ITask::SingleStepWithoutTCrit => {
                        tcrit = rworks.tcrit;
                        if (self.tn_ - tcrit) * self.h_ > 0. {
                            error!(
                                "[lsoda] itask = SingleStepWithoutTCrit and tcrit behind tcur\n"
                            );
                            self.terminate();
                            return;
                        }
                        hmx = self.tn_.abs() + self.h_.abs();
                        ihit = (self.tn_ - tcrit).abs() <= (100. * ETA * hmx);
                        if ihit {
                            *t = tcrit;
                            self.successreturn(y, t, itask, ihit, tcrit);
                            return;
                        }
                        tnext = self.tn_ + self.h_ * (1. + 4. * ETA);
                        if (tnext - tcrit) * self.h_ <= 0. {
                            ()
                        } else {
                            self.h_ = (tcrit - self.tn_) * (1. + 4. * ETA);
                            match self.lastoutcome {
                                IStateOutput::IntegrationSuccess => {
                                    self.jstart = JStart::NewHNextStep
                                }
                                _ => (),
                            }
                        }
                    }
                } /* end match task   */
            }
            _ => (),
        }
        /*
           Block e.
           The next block is normally executed for all calls and contains
           the call to the one-step core integrator stoda.

           This is a looping point for the integration steps.

           First check for too many steps being taken, update ewt ( if not at
           start of problem).  Check for too much accuracy being requested, and
           check for h_ below the roundoff level in *t.
        */
        loop {
            if *istatein != IStateInput::InitialCall || self.nst != 0 {
                if self.nst - self.nslast >= self.mxstep {
                    error!("[LSODA] {}  steps taken before reaching tout", self.mxstep);
                    self.lastoutcome = IStateOutput::ExcessiveAmt;
                    self.terminate_inner(y, t);
                    return;
                }
                self.ewset(&self.yh_[1].clone());
                for i in 1..self.n + 1 {
                    if self.ewt[i] <= 0. {
                        error!("[LSODA] ewt[{}] <= 0.", i);
                        self.lastoutcome = IStateOutput::EWTiZero;
                        self.terminate_inner(y, t);
                        return;
                    }
                    self.ewt[i] = 1. / self.ewt[i];
                }
            }
            tolsf = ETA * vmnorm(self.n, &self.yh_[1], &self.ewt);
            if tolsf > 0.01 {
                tolsf *= 200.;
                if self.nst == 0 {
                    error!(
                        "[LSODA] at start of problem, too much accuracy\n
                    requested for precision of machine,\n
                    suggested scaling factor = {}\n",
                        tolsf
                    );
                    self.terminate();
                    return;
                }
                error!(
                    "[LSODA] at t = {}, too much accuracy requested\n
                for precision of machine, suggested\n
                scaling factor = {}\n",
                    t, tolsf
                );
                self.lastoutcome = IStateOutput::TooMuchAccuracy;
                self.terminate_inner(y, t);
                return;
            }

            if self.tn_ + self.h_ == self.tn_ {
                self.nhnil += 1;
                if self.nhnil <= self.mxhnil {
                    warn!(
                        "[LSODA] warning..internal t = {} and h_ = {} are\n
                    such that in the machine, t + h_ = t on the next step\n
                    solver will continue anyway.\n",
                        self.tn_, self.h_
                    );
                }
                if self.nhnil == self.mxhnil {
                    error!("[LSODA] above warning has been issued {} times, it will not be issued again for this problem", self.nhnil);
                }
            }
            /* Call stoda */
            self.stoda(neq, y, system);
            match self.kflag {
                KFlag::StepSuccess => {
                    /*
                       Block f.
                       The following block handles the case of a successful return from the
                       core integrator ( kflag = 0 ).
                       If a method switch was just made, record tsw, reset maxord,
                       set jstart to -1 to signal stoda to complete the switch,
                       and do extra printing of data if ixpr = 1.
                       Then, in any case, check for stop conditions.
                    */
                    self.init = InitializationState::Done;
                    if self.meth_ != self.mused {
                        self.tsw = self.tn_;
                        self.maxord = self.mxordn;
                        match self.meth_ {
                            OdeMethod::BDF => self.maxord = self.mxords,
                            OdeMethod::Adam => (),
                        }
                        self.jstart = JStart::NewParamNextStep;
                        match self.ixpr {
                            Ixpr::PrintOnSwitch => match self.meth_ {
                                OdeMethod::BDF => {
                                    info!("[LSODA] a switch to the stiff method has occurred")
                                }
                                OdeMethod::Adam => {
                                    info!("[LSODA] a switch to the nonstiff method has occurred")
                                }
                            },
                            _ => (),
                        }
                    }
                    match itask {
                        ITask::NormalComputation => {
                            if (self.tn_ * tout) * self.h_ < 0. {
                                continue;
                            }
                            self.intdy(tout, 0, y, &mut iflag);
                            *t = tout;
                            self.lastoutcome = IStateOutput::IntegrationSuccess;
                            self.illin = 0;
                            return;
                        }
                        ITask::SingleStep => {
                            self.successreturn(y, t, itask, ihit, tcrit);
                            return;
                        }
                        ITask::StopAtFirstMeshPoint => {
                            if (self.tn_ * tout) * self.h_ >= 0. {
                                self.successreturn(y, t, itask, ihit, tcrit);
                                return;
                            }
                            continue;
                        }
                        ITask::NormalComputationWithoutTCrit => {
                            if (self.tn_ * tout) * self.h_ >= 0. {
                                self.intdy(tout, 0, y, &mut iflag);
                                *t = tout;
                                self.lastoutcome = IStateOutput::IntegrationSuccess;
                                self.illin = 0;
                                return;
                            } else {
                                hmx = self.tn_.abs() + self.h_.abs();
                                ihit = (self.tn_ - tcrit).abs() <= 100. * ETA * hmx;
                                if ihit {
                                    self.successreturn(y, t, itask, ihit, tcrit);
                                    return;
                                }
                                tnext = self.tn_ + self.h_ * (1. + 4. * ETA);
                                if (tnext - tcrit) * self.h_ <= 0. {
                                    continue;
                                }
                                self.h_ = (tcrit - self.tn_) * (1. + 4. * ETA);
                                self.jstart = JStart::NewHNextStep;
                                continue;
                            }
                        }
                        ITask::SingleStepWithoutTCrit => {
                            hmx = self.tn_.abs() + self.h_.abs();
                            ihit = (self.tn_ - tcrit).abs() <= 100. * ETA * hmx;
                            if ihit {
                                self.successreturn(y, t, itask, ihit, tcrit);
                                return;
                            }
                        }
                    }
                }
                KFlag::CouldNotAchieveReqError => {
                    error!("[LSODA] at t = {} and step size h = {}, the error test failed repeatedly or
                    with abs(h) = hmin\n", self.tn_, self.h_);
                    self.lastoutcome = IStateOutput::RepTestFailure;
                    big = 0.;
                    self.imxer = 1;
                    for i in 1..self.n + 1 {
                        size = self.acor[i].abs() * self.ewt[i];
                        if big < size {
                            big = size;
                            self.imxer = i;
                        }
                    }
                    self.terminate_inner(y, t);
                    return;
                }
                KFlag::CouldNotAchieveConvergence => {
                    error!("[LSODA] at t = {} and step size h = {}, corrector convergence failed repeatedly or
                    with abs(h) = hmin\n", self.tn_, self.h_);
                    self.lastoutcome = IStateOutput::RepConvFailure;
                    big = 0.;
                    self.imxer = 1;
                    for i in 1..self.n + 1 {
                        size = self.acor[i].abs() * self.ewt[i];
                        if big < size {
                            big = size;
                            self.imxer = i;
                        }
                    }
                    self.terminate_inner(y, t);
                    return;
                }
                _ => (),
            }
        }
    }

    pub fn solve(
        &mut self,
        system: &dyn OdeSystem,
        neq: usize,
        y: &[f64],
        tin: &mut f64,
        tout: f64,
        istate: &mut IStateInput,
        rtol: f64,
        atol: f64,
        dense: bool,
    ) -> (Vec<f64>, Vec<Vec<f64>>) {
        let limit = (tout * 10.) as usize;
        let mut y_i = vec![0.; neq + 1];
        let mut all_t = Vec::<f64>::new();
        let mut all_y = Vec::<Vec<f64>>::new();

        let iworks = IWork::default();
        let rworks = Rworks::default();
        let itask: ITask = ITask::NormalComputation;
        let iopt: IOpt = IOpt::NoOptionalInputs;
        let jt: CorrectionIterMethod = CorrectionIterMethod::InternalFullJac;
        // Set the tolerance. We should do it only once.
        self.rtol_.resize(neq + 1, rtol);
        self.atol_.resize(neq + 1, atol);
        self.rtol_[0] = 0.0;
        self.atol_[0] = 0.0;
        // Fill-in values.
        for i in 1..neq + 1 {
            y_i[i] = y[i - 1];
        }

        let f_step = 0.01;

        self.step(
            system,
            neq,
            &mut y_i,
            tin,
            *tin + f_step,
            itask,
            istate,
            iopt,
            jt,
            iworks,
            rworks,
        );
        all_t.push(*tin);
        all_y.push(y_i.to_vec());

        let mut count = 0usize;
        while *tin + 2. * self.h_ < tout {
            self.step(
                system,
                neq,
                &mut y_i,
                tin,
                *tin + self.h_,
                itask,
                istate,
                iopt,
                jt,
                iworks,
                rworks,
            );
            *tin += self.h_;
            count += 1;
            if !dense && count == limit {
                all_t.push(*tin);
                all_y.push(y_i.to_vec());
                count = 0;
            }
        }

        self.step(
            system, neq, &mut y_i, tin, tout, itask, istate, iopt, jt, iworks, rworks,
        );
        *tin = tout;
        all_t.push(*tin);
        all_y.push(y_i.to_vec());

        (all_t, all_y)
    }
}

pub trait OdeSystem {
    fn func(&self, _t: f64, _y: &mut [f64], _dy: &mut [f64]);
}

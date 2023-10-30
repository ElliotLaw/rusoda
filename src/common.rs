pub static SQRTETA: f64 = 1.4901161193847656e-08;
pub static ETA: f64 = 2.2204460492503131e-16;

#[derive(Debug, Default, Clone)]
pub struct LSODACommon {
    pub yh: Vec<Vec<f64>>,
    pub wm: Vec<Vec<f64>>,
    pub ewt: Vec<f64>,
    pub savf: Vec<f64>,
    pub acor: Vec<f64>,
    pub ipvt: Vec<f64>,
    /* static variables for LSODA */
    pub h: f64,
    /*
        [Ref]
        odepack/src/M_main/dlsoda.inc

       the step size in t last used (successfully).
    */
    pub hu: f64,
    pub rc: f64,
    pub tn: f64,
    pub tsw: f64,
    pub pdnorm: f64,
    /* no static variable for prja(), solsy() */
    /* static variables for stoda() */
    pub crate_: f64,
    pub el: [f64; 14],
    pub elco: [[f64; 14]; 13],
    pub tesco: [[f64; 4]; 13],
    pub hold: f64,
    pub rmax: f64,
    pub pdest: f64,
    pub pdlast: f64,
    /* static variables for various vectors and the Jacobian. */
    pub ialth: usize,
    pub ipup: usize,
    pub nslp: usize,
    pub icount: usize,
    pub irflag: usize,
    pub imxer: usize,
    pub illin: usize,
    pub nhnil: usize,
    pub nslast: usize,
    pub jcur: usize,
    pub meth: usize,
    pub mused: usize,
    pub nq: usize,
    pub nst: usize,
    pub ncf: usize,
    pub nfe: usize,
    pub nje: usize,
    pub nqu: usize,
    pub miter: usize,
}

/*
    LSODA Switch between AdamBash. and BDF methods to solve
    nonstiff or stiff problem.
*/
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OdeMethod {
    Adam,
    BDF,
}
impl Default for OdeMethod {
    fn default() -> Self {
        Self::Adam
    }
}

/*
    LSODA call state.
*/
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IStateInput {
    /*istate code 1 */
    InitialCall,
    /*istate code 2 */
    ParamChange,
    /*istate code 3 */
    ParamNotChange,
}
impl Default for IStateInput {
    fn default() -> Self {
        Self::InitialCall
    }
}

/*
    LSODA output state.
*/
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IStateOutput {
    /*istate code 1 */
    DoNothing,
    /*istate code 2 */
    IntegrationSuccess,
    /*istate code -1 */
    ExcessiveAmt,
    /*istate code -2 */
    TooMuchAccuracy,
    /*istate code -3 */
    IllegalInput,
    /*istate code -4 */
    RepTestFailure,
    /*istate code -5 */
    RepConvFailure,
    /*istate code -6 */
    EWTiZero,
    /*istate code -7 */
    IRWorkTooShort,
}
impl Default for IStateOutput {
    fn default() -> Self {
        Self::DoNothing
    }
}

/*
    Task type in LSODA for input.
*/
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ITask {
    /*itask code 1 */
    NormalComputation,
    /*itask code 2 */
    SingleStep,
    /*itask code 3 */
    StopAtFirstMeshPoint,
    /*itask code 4 */
    NormalComputationWithoutTCrit,
    /*itask code 5 */
    SingleStepWithoutTCrit,
}
impl Default for ITask {
    fn default() -> Self {
        Self::NormalComputation
    }
}

/*
    Tolerance type in LSODA for input.
*/
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ITol {
    /*itol code 1 */
    SingleRtolSingleAtol,
    /*itol code 2 */
    SingleRtolMultiAtol,
    /*itol code 3 */
    MultiRtolSingleAtol,
    /*itol code 4 */
    MultiRtolMultiAtol,
}
impl Default for ITol {
    fn default() -> Self {
        Self::SingleRtolMultiAtol
    }
}

/*
    Whether optional input is passed into LSODA.
*/
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IOpt {
    /*iopt code 0 */
    OptionalInputs,
    /*iopt code 1 */
    NoOptionalInputs,
}
impl Default for IOpt {
    fn default() -> Self {
        Self::NoOptionalInputs
    }
}

/*
    Jacobian type in LSODA
*/
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CorrectionIterMethod {
    /*code 0 */
    NoJac,
    /*code 1 */
    UserFullJac,
    /*code 2 */
    InternalFullJac,
    /*code 4 */
    UserBandJac,
    /*code 5 */
    InternalBandJac,
}
impl Default for CorrectionIterMethod {
    fn default() -> Self {
        Self::NoJac
    }
}

/*
    flag to generate extra printing at method switches.
    IXPR = 0 means no extra printing (the default).
    IXPR = 1 means print data on each switch.
    T, H, and NST will be printed on the same logical
    unit as used for error messages.
*/
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Ixpr {
    /*ixpr code 0 */
    NoExtraPrint,
    /*ixpr code 1 */
    PrintOnSwitch,
}
impl Default for Ixpr {
    fn default() -> Self {
        Self::NoExtraPrint
    }
}

/*
    Subroutine STODA in LSODA step type.
*/
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum JStart {
    /*jstart code 0 */
    FirstStep,
    /*jstart code > 0 */
    ContinueStep,
    /*jstart code -1 */
    NewParamNextStep,
    /*jstart code -2 */
    NewHNextStep,
}
impl Default for JStart {
    fn default() -> Self {
        Self::FirstStep
    }
}

/*
    Whether LSODA is initialized.
*/
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum InitializationState {
    Done,
    NotYet,
}
impl Default for InitializationState {
    fn default() -> Self {
        InitializationState::Done
    }
}

/*
    If K and T were legal
*/
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IFlag {
    /*iflag code 0 */
    KTLegal,
    /*iflag code -1 */
    IllegalK,
    /*iflag code -2 */
    IllegalT,
}
impl Default for IFlag {
    fn default() -> Self {
        Self::KTLegal
    }
}

/*
    completion state of STODA
*/
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KFlag {
    /*kflag code 0 */
    StepSuccess,
    /*kflag code -1 */
    CouldNotAchieveReqError,
    /*kflag code -2 */
    CouldNotAchieveConvergence,
    /*kflag code -3 */
    FatalinPJACoSLVS,
    /*kflag code <-3 */
    MultiFails(u8),
}
impl Default for KFlag {
    fn default() -> Self {
        Self::StepSuccess
    }
}

/*
    Used for work space.(usize or enum)
    IWork.0 = ml: usize
    IWork.1 = mu: usize
        These are the lower and upper half-bandwidths, respectively, of the
        banded Jacobian, excluding the main diagonal.  The band is defined by
        the matrix locations (i,j) with i-ML .le. j .le. i+MU.  ML and MU must
        satisfy  0 .le.  ML,MU  .le. NEQ-1.  These are required if JT is 4 or 5,
        and ignored otherwise.  ML and MU may in fact be the band parameters
        for a matrix to which df/dy is only approximately equal.

    IWork.2 = ixpr: Ixpr
        See enum Ixpr.

    IWork.3 = mxstep: usize
        maximum number of (internally defined) steps
        allowed during one call to the solver.
        The default value is 500.

    IWork.4 = mxhnil: usize
        maximum number of messages printed (per problem)
        warning that T + H = T on a step (H = step size).
        This must be positive to result in a non-default
        value.
        The default value is 10.

    IWork.5 = mxordn: usize
        the maximum order to be allowed for the nonstiff
        (Adams) method.  the default value is 12.
        if MXORDN exceeds the default value, it will
        be reduced to the default value.
        MXORDN is held constant during the problem.

    IWork.6 = mxords: usize
        the maximum order to be allowed for the stiff
        (BDF) method.  The default value is 5.
        If MXORDS exceeds the default value, it will
        be reduced to the default value.
        MXORDS is held constant during the problem.

*/
#[derive(Debug, Default, Clone, Copy)]
pub struct IWork(
    pub usize,
    pub usize,
    pub Ixpr,
    pub usize,
    pub usize,
    pub usize,
    pub usize,
);

/*
    Used for work space.(f64)
*/
#[derive(Debug, Default, Clone, Copy)]
pub struct Rworks {
    /*
    Critical point to return.
    TCRIT may be equal to or beyond TOUT, but not behind it
    in the direction of integration.
    This option is useful if the problem has a singularity
    at or beyond t = TCRIT.
    */
    pub tcrit: f64,
    /*
    the step size to be attempted on the first step.
    The default value is determined by the solver.
     */
    pub h0: f64,
    /*
    the maximum absolute step size allowed.
    The default value is infinite.
     */
    pub hmax: f64,
    /*
    the minimum absolute step size allowed.
    The default value is 0.  (This lower bound is not
    enforced on the final step before reaching TCRIT
    when ITASK code is 4 or 5.)
     */
    pub hmin: f64,
}

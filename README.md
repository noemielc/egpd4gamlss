# egpd4gamlss
An add-on to the GAMLSS R-package [1] allowing to fit non-stationary marginal distributions following the Extended Generalised Pareto model (EGPD) in the sense of [2].
The add-on is made generic by means of symbolic derivation so that any parametric form satisfying the requirements of the EGPD family can be used and its parameter fitted and regressed on covariates by means of Generalized Additive Modelling forms.

The repository contains a tutorial with background theory about EGPD and main functions of GAMLSS for model fitting, analysis and selection, as well as a short version of the associated tutorial code. The core code of the add-on is GenericEGPDg.R.


References:
[1] Stasinopoulos, D. Mikis, and Robert A. Rigby. "Generalized additive models for location scale and shape (GAMLSS) in R." Journal of Statistical Software 23 (2008): 1-46.
[2] Naveau, Philippe, et al. "Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection." Water Resources Research 52.4 (2016): 2753-2769.


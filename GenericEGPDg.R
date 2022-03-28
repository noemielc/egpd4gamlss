##                                                     ##
##    Extended GPD Family implemented in GAMLSS        ##
##                                                     ##
## Ext-GPD in the sense of                             ##
## Naveau, Philippe, et al. "Modeling jointly low,     ##
## moderate, and heavy rainfall intensities without a  ##
## threshold selection." Water Resources Research 52.4 ##
## (2016): 2753-2769.                                  ##
##                                                     ##
## Version 14.01.2022                                  ##
## Noémie Le Carrer & Nicolas Berthier                 ##
## Correspondance: noemie.le-carrer@imt-atlantique.fr  ##
##                                                     ##

library(Deriv)

EGPD_G <- function (G) {
  ## G can be a function (of z), or simply an expression
  
  ## CDF of the GPD (mu)
  Hfun <- function (z, mu) {
      if (mu[1] != 0) {
          1 - (1 + mu * z) ^ (- 1 / mu)
      } else {
          1 - exp(-z)
      }
  }

  ## PDF of the GPD (mu)
  hfun <- function (z, mu) {
      if (mu[1] != 0) {
          (1 + mu * z) ^ (- (mu + 1) / mu)
      } else {
          exp (-z)
      }
  }
  .Gfun <- if (is.function (G)) body (G) else G
  .gfun <- Deriv (.Gfun, "z", combine = "cbind")
  .Hfunys <- do.call ('substitute', list (body (Hfun),
                                          list (z = quote (y / sigma))))
  .hfunys <- do.call ('substitute', list (body (hfun),
                                          list (z = quote (y / sigma))))
  CDF <- do.call ('substitute', list (.Gfun, list (z = .Hfunys)))
  
  .PDFg <- do.call ('substitute', list (.gfun, list (z = .Hfunys)))
  
  LogPDF <- bquote (log ((1/sigma) * .(.PDFg) * .(.hfunys)))
  
  list (CDF = CDF, LogPDF = LogPDF, G = .Gfun)
}

MakeEGPD <- function (G,
                      Gname = "",
                      y.valid = quote (TRUE),
                      nu.valid = quote (nu > 0),
                      tau.valid = quote (tau > 0),
                      mean = quote (ifelse (nu > 1, mu, NaN)),
                      variance = quote (ifelse (nu > 2, (sigma ^ 2 * nu) / (nu - 2), Inf))) {
  
  .concat <- function (...) paste (..., sep = '')
  Mname   <- .concat ("EGPD", Gname)
  MnameX  <- .concat ("EGPD (", Gname, ")")
  
  Gfrmls  <- formals (G)
  has_nu  <- "nu" %in% names (Gfrmls)
  has_tau <- "tau" %in% names (Gfrmls)
  Frmls   <- c (alist (mu = 0.1, sigma = 1),
                if (has_nu) alist (nu = 2) else alist (),
                if (has_tau) alist (tau = 1) else alist ())
  Mfrmls  <- as.pairlist (c (alist (y = , mu = , sigma = ),
                             if (has_nu) alist (nu = ) else alist (),
                             if (has_tau) alist (tau = ) else alist ()))
  Efrmls  <- as.pairlist (Mfrmls[2:length(Mfrmls)]) #elim `y' from Mfrmls
  Gargs   <- sapply (names (Gfrmls), as.symbol, USE.NAMES = FALSE)
  Margs   <- sapply (names (Mfrmls), as.symbol, USE.NAMES = FALSE)
  Gparams <- Gargs[2:length(Gargs)]
  Mparams <- Margs[2:length(Margs)]
  Gparams_i = c (if (has_nu) quote(nu[i]) else list(),
                 if (has_tau) quote(tau[i]) else list())

  Function     <- function (frmls, body) eval (call ("function", frmls, body))
  ModelFun     <- function (body) Function (Mfrmls, body)
  DefGlobalFun <- function (name, def) eval (bquote (.(as.name (name)) <<- def))
  DefModelFun  <- function (fun, pre_args, post_args, body)
    DefGlobalFun (fun, Function (as.pairlist (c(pre_args, Frmls, post_args)),
                                 body))
  
  EGPDModel <- EGPD_G (G)
  G         <- Function (Gfrmls, bquote ({ .(EGPDModel$G) }))
  CDF       <- ModelFun (bquote ({ .(EGPDModel$CDF) }))
  LogPDF    <- ModelFun (bquote ({ .(EGPDModel$LogPDF) }))
  mean      <- Function (Efrmls, mean)
  variance  <- Function (Efrmls, variance)

  dmu       <- Deriv (LogPDF, "mu", combine="cbind")
  dsigma    <- Deriv (LogPDF, "sigma", combine="cbind")
  if (has_nu) {
    dnu   <- Deriv (LogPDF, "nu", combine="cbind")
  }
  if (has_tau) {
    dtau  <- Deriv (LogPDF, "tau", combine="cbind")
  }
  
  check_pos <- function (b, v)
    bquote (if (.(b) && (any (.(v) <= 0)))
      stop (paste (.(as.character (v)), "must be strictly positive\n")))
  
  check_nu <- check_pos (has_nu, quote (nu))
  check_tau <- check_pos (has_tau, quote (tau))
  
  ##---------------------------------------------------------------------------------------
  dFun <- .concat ("d", Mname)
  DefModelFun (dFun, alist (x = ), alist (log = FALSE), bquote ({
    #if (any(mu < 0))  stop(paste("xi must be positive", "\n", ""))
    if (any(sigma <= 0))  stop(paste("sigma must be strictly positive", "\n", "")) 
    .(check_nu)
    .(check_tau)
    lik <- exp (.(as.call (c(quote (LogPDF), quote (x), Mparams))))
    if (log == FALSE) lik else log (lik)
  }))
  
  ##--------------------------------------------------------------------------------------- 
  pFun <- .concat ("p", Mname)
  DefModelFun (pFun, alist (q = ), alist (), bquote ({
    #if (any(mu < 0))  stop(paste("xi must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be strictly positive", "\n", "")) 
    .(check_nu)
    .(check_tau)
    .(as.call (c(quote (CDF), quote (q), Mparams)))
  }))
  
  ##---------------------------------------------------------------------------------------
  qFun <- .concat ("q", Mname)
  DefModelFun (qFun, alist (p = ), alist (), bquote ({
    #if (any(mu < 0))  stop(paste("xi must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be strictly positive", "\n", "")) 
    .(check_nu)
    .(check_tau)
    lp <- length (p)
    expand <- function (v) rep (v, length = lp)
    if (has_nu) { nu <- expand (nu) }
    if (has_tau) { tau <- expand (tau) }
    loop.uniroot <- function (i) {
      uniroot(function (x) .(as.call (c(quote (G), quote (x), Gparams_i))) - p[i],
              lower = 0, upper = 1, tol = 1e-15)$root
    }
    Ginv <- vapply (seq (along = p), loop.uniroot, numeric (1))
    mu <- expand (mu)
    ifelse (mu != 0,
            (expand (sigma) / mu) * ((1 - Ginv) ^ (- mu) - 1),
            expand (-sigma) * log (1 - Ginv))
  }))

  ##---------------------------------------------------------------------------------------
  ## Random generation function for the model distribution:
  rFun <- .concat ("r", Mname)
  DefModelFun (rFun, alist (n = ), alist (), bquote ({
    #if (any(mu < 0))  stop(paste("xi must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be strictly positive", "\n", "")) 
    .(check_nu)
    .(check_tau)
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))
    .(as.call (c (as.name (qFun), quote (runif (ceiling (n))), Mparams)))
  }))
  
  ##----------------------------------------------------------------------------------------
  ## Set formal arguments /a posteriori/ to avoid a deep/nested `bquote'.
  Ffrmls <- c (alist (   mu.link = "identity", mu.init = 0.15,
                      sigma.link = "log",   sigma.init = 2),
               if (has_nu)  alist ( nu.link = "log",  nu.init = 3) else alist (),
               if (has_tau) alist (tau.link = "log", tau.init = 1) else alist ())

  F <- function () {

    ## Some helpers for definitions below
    .checklink <-
      function (link, s, pos = c("inverse", "log", "identity", "own")) {
        checklink (link, MnameX, s, pos)
      }
    
    d2 <- function (d) ModelFun (bquote ({
      .d <- .(body (d))
      .x <- - .d * .d
      ifelse(.x < -1e-15,  .x, -1e-15)
    }))
    
    dd2 <- function (d, b) ModelFun (bquote ({
      .x <- - .(body (d)) * .(body (b))
      ifelse(.x < -1e-15,  .x, -1e-15)
    }))
    
    ## Define parameter-related lists & functions
    
    mstats <- .checklink ("mu.link",    substitute (mu.link))
    dstats <- .checklink ("sigma.link", substitute (sigma.link))
    parameters <- list (mu = TRUE, sigma = TRUE)
    rqparams <- alist (mu = mu, sigma = sigma)
    
    if (has_nu) {
      parameters["nu"] <- TRUE
      rqparams["nu"] <- expression (nu)
      vstats <- .checklink ("nu.link", substitute (nu.link))
      nu_fields <-
        list (nu.link = as.character (substitute (nu.link)),
              nu.linkfun = vstats$linkfun,
              nu.linkinv = vstats$linkinv,
              nu.dr = vstats$mu.eta,
              nu.valid = eval (bquote (
                function (nu) { all (.(nu.valid)) })),
              nu.initial = eval (bquote (
                expression (nu <- rep (.(nu.init), length (y))))),
              dldv = dnu,
              d2ldv2 = d2 (dnu),
              d2ldmdv = dd2 (dmu, dnu),
              d2ldddv = dd2 (dsigma, dnu))
    } else {
      nu_fields <- list ()
    }

    if (has_tau) {
      parameters["tau"] <- TRUE
      rqparams["tau"] <- expression (tau)
      tstats <- .checklink ("tau.link", substitute (tau.link))
      tau_fields <-
        list (tau.link = as.character (substitute (tau.link)),
              tau.linkfun = tstats$linkfun,
              tau.linkinv = tstats$linkinv,
              tau.dr = tstats$mu.eta,
              tau.valid = eval (bquote (
                function (tau) { all (.(tau.valid)) })),
              tau.initial = eval (bquote (
                expression (tau <- rep (.(tau.init), length (y))))),
              dldt = dtau,
              d2ldt2 = d2 (dtau),
              d2ldmdt = dd2 (dmu, dtau),
              d2ldddt = dd2 (dsigma, dtau))
      tau_fields <- if (has_nu) c (tau_fields, list (d2ldvdt = dd2 (dnu, dtau)))
      else tau_fields
    } else {
      tau_fields <- list ()
    }
    
    ## ---
    
    G.dev.incr <-
      ModelFun (bquote ({
        -2 * .(as.call (c (as.symbol (dFun), Margs, log = TRUE)))
      }))
    
    rqres_call <-
      as.call (c (quote (rqres),
                  pfun = pFun, type = "Continuous", y = quote (y), rqparams))

    fields <-
      list (family = c(Mname, MnameX),
            parameters = parameters,
            nopar = length (parameters),
            type = "Continuous",
            
            mu.link = as.character (substitute (mu.link)),
            mu.linkfun = mstats$linkfun,
            mu.linkinv = mstats$linkinv,
            mu.dr = mstats$mu.eta,
            mu.valid = function (mu) TRUE,
            mu.initial = eval (bquote (
              expression (mu <- rep (.(mu.init), length (y))))),
            
            sigma.link = as.character (substitute (sigma.link)),
            sigma.linkfun = dstats$linkfun,
            sigma.linkinv = dstats$linkinv,
            sigma.dr = dstats$mu.eta,
            sigma.valid = function (sigma) all (sigma > 0),
            sigma.initial = eval (bquote (
              expression (sigma <- rep (.(sigma.init), length (y))))),
            
            dldm = dmu,
            d2ldm2 = d2 (dmu),
            dldd = dsigma,
            d2ldd2 = d2 (dsigma),
            d2ldmdd = dd2 (dmu, dsigma),
            
            G.dev.incr = G.dev.incr,
            rqres = rqres_call,
            
            y.valid = eval (bquote (function (y) { all(.(y.valid)) })),
            mean = mean,
            variance = variance)
    
    S <- structure (c(fields, nu_fields, tau_fields),
                    class = c("gamlss.family", "family"))
    ## NB: (?) whenever returned function is called...
    eval (bquote (.(as.symbol (Mname)) <<- S))
    S
  }
  formals (F) <- as.pairlist (Ffrmls)
  ## ... or once when generic family is declared:
  eval (bquote (.(as.symbol (Mname)) <<- F))
  F
}

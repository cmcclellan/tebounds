
#' tebounds
#'
#' @description This function is a port of the Stata command tebounds (Millimet, McCarthy, and Roy, 2015).  tebounds provides estimates of the bounds of the ATE assuming binary treatment (Y) and a binary outcome (D).  The tebounds command
#' estimates the bounds of the ATE under a range of alternative assumptions regarding the selection process and the extent of measurement error, including assumptions of monotone treatment response and monotone treatment selection.
#' The command also allows for the inclusion of an MIV.  The monotone treatment response, monotone treatment selection, and MIV approaches follow the estimation procedure developed in Manski and Pepper (2000), Kreider and Pepper (2007), Kreider et al. (2012), and Gundersen, Kreider, and Pepper (2012).


#'
#' @param data A data frame or survey design object containing a binary outcome (Y) and a binary treatment (D).
#' The data can optionally include a monotone instrumental variable (Z) and a vector of weights.
#' @param formula A formula indicating the outcome, treatment, and instrument (if any). Instruments are specified using the "|"  at the end of the formula (e.g. Y~D|Z).
#' @param weights Data frame column indicating weight variable, provided without quotation marks.  If survey design with weights is supplied, then this argument is not needed.
#' @param subset  Conditions to subset data.
#' @param erates  Numeric vector, set to (0, 0.5, 0.1, 0.25) by default.  Vector specifying assumed rates of measurement error.
#' @param np Integer, set to 500 by default, specifying number of intervals for grid search in arbitrary measurement-error model
#' @param bootstrap Boolean, set to FALSE by default. Set to TRUE to estimate confidence intervals (CIs) via bootstrap
#' @param boot.reps Integer, set to 100 by default, specifying number of bootstrap replications for confidence interval estimation.
#' @param CI.level Integer from 0 to 100, set to 95 by default.  Level at which to calculate confidence intervals.
#' @param im Boolean, set to FALSE by default.  Set to TRUE to calculate CIs following Imbens and Manski (2004)
#' @param k Integer, set to 100 by default. Number of bootstrap replications for MIV bias correction.
#' @param ncells Integer, set to 5 by default. Number of cells used in the MIV estimator.
#' @importFrom boot boot
#' @import ggplot2
#' @importFrom survey as.svrepdesign
#' @return A list of the class 'tebounds'.  Summary and plot methods are provided to match output of the Stata routine.
#' @export
#'
#'
#' @examples
#'
#'
tebounds <- function(data, formula, weights, subset,
                     erates = c(0, 0.05, 0.1, .25),
                     np = 500,
                     bootstrap = F,
                     boot.reps = 100,
                     CI.level = 95,
                     im = F,
                     k = 100,
                     ncells = 5
                     ) {




  cl <- match.call()
  mf <- match.call(expand.dots = F)
  m <- match(c( "data","formula", "weights", "subset"),
             names(mf), 0L)
  mf <- mf[c(1, m)]

  formula <- Formula::as.Formula(formula)
  mf$formula <- formula

  mf[[1]] <- as.name("model.frame")

  if(any(class(data) %in% c('survey.design', 'survey.design2'))){
    mf$weights <- data$call$weights[[2]]
    svy.design <- data
    data <- svy.design$variables
    mf[[2]] <- data
  }
  mf <- eval(mf, parent.frame())


  Y <- model.response(mf, 'numeric')
  Dx <- terms(formula, rhs = 1)
  D <- model.matrix(Dx, mf)[,2]
  N <- length(Y)

  W <- model.weights(mf)
  if (is.null(W)){

    W <- rep(1,N)

  }


  if(length(formula)[2] < 2L) {
    mtZ <- NULL
    Z <- rep(0,N)
    miv <- F
  } else {
    mtZ <- delete.response(terms(formula, data = data, rhs = 2))
    Z <- model.matrix(mtZ, mf)[,2]
    n.value.z <- length(unique(Z))

    if (n.value.z > ncells) {
      if(exists('svy.design')){
        Z <- cut(Z, c(-Inf,svyquantile(mtZ, svy.design, seq(1/ncells, 1 - 1/ncells, length = ncells - 1), type=2), Inf),include.lowest=TRUE, dig.lab = 5)
        Z[is.na(Z)] <-  levels(Z)[1]
      } else {
        Z <- cut(Z, c(-Inf,quantile(Z, seq(1/ncells, 1 - 1/ncells, length = ncells - 1), type=2), Inf),include.lowest=TRUE, dig.lab = 5 )

      }
    }

    Z <- as.factor(Z)
    p <- model.matrix(~Z-1)
    p_z <- apply(p, 2, function(x) weighted.mean(x, W))
    miv <- T

  }




  d2 <- cbind(Y,D,Z,W)

  if(bootstrap){
    if(exists('svy.design')){
        booted.bounds <- list()
        bsw.dsgn <- survey::as.svrepdesign(svy.design, type = 'subbootstrap', compress=F, replicates = boot.reps)
        bsw <- bsw.dsgn$repweights[complete.cases(data[all.vars(formula)]),]*W

        booted.bounds$t0 <- bound_calc(d2, erates, np)
        t <- matrix(NA, nrow = ncol(bsw), ncol = length(booted.bounds$t0))
        for (b in 1:ncol(bsw)){
          d2[,'W']  <- bsw[,b]
          booted.bound <- bound_calc(d2, erates, np)
          t[b,] <- as.vector(booted.bound)

          d2[,'W']  <- W
        }
        booted.bounds$t <- t
      } else{

        boot_bounds <- function(data, i, ...){
          data2 <- data[i,]
          return(bound_calc(data2, erates, np))
        }

        booted.bounds <- boot::boot(d2, boot_bounds, boot.reps, erates = erates, np=np)
      }

    bounds <- booted.bounds$t0
    colnames(booted.bounds$t) <- paste0(rep(colnames(bounds), each = length(rownames(bounds))),rownames(bounds))
    stdev <- apply(booted.bounds$t,2,sd, na.rm=T)
    nn <- colSums(!is.na(booted.bounds$t))
    if(!im){
      low <-  ((100-CI.level)/2)/100
      high <-  1-low
      l_ci <- apply(booted.bounds$t, 2, quantile, probs = low)
      l_ci <- matrix(l_ci, nrow = 24, byrow = F)
      l_ci <- pmax(l_ci[seq(2,24,2),], -1)

      u_ci <- apply(booted.bounds$t, 2, quantile, probs = high)
      u_ci <- matrix(u_ci, nrow = 24, byrow = F)
      u_ci <- pmin(u_ci[seq(1,23,2),], 1)

    } else {
      delta <- diff(bounds)[seq(1,23,2),]*-1
      stdev <- matrix(stdev, nrow = 24, byrow=F)
      stdev.u <- stdev[seq(1,23,2),]
      stdev.l <- stdev[seq(2,24,2),]
      s <- 1/pmax(stdev.l, stdev.u)
      c <- matrix(1, nrow=nrow(delta), ncol=ncol(delta))

      for ( i in 1:nrow(delta)){
        for (j in 1:ncol(delta)){
          c[i,j] <- qnorm(1-((100-CI.level)/200))

          diff <- 99999
          while (diff>0.00001){
            c[i,j] <- c[i,j]-0.00001
            d <- pnorm(c[i,j]+s[i,j]*delta[i,j])-pnorm(-c[i,j])-(CI.level/100)
            diff <- min(abs(d), diff)
            if (c[i,j]<1){
              warning(paste('Error with I-M CIs: CbarN down to 1 in ', rownames(delta)[i]))
              break
            }
          }
        }
      }

      l_ci <- pmax(bounds[seq(2,24,2),]-c*stdev.l, -1)
      u_ci <- pmin(bounds[seq(1,23,2),]+c*stdev.u,  1)

    }

    l_ci.bounds <- l_ci
    u_ci.bounds <- u_ci

  } else {
    bounds <- bound_calc(d2, erates , np)

  }

  if (miv) {
    if(bootstrap){
      if(exists('svy.design')){
        booted.miv.bounds <- list()
        booted.miv.bounds$t0 <- miv_bounds_calc(d2, erates, k, np)
        t <- matrix(NA, nrow = ncol(bsw), ncol = length(booted.miv.bounds$t0))
        for (b in 1:ncol(bsw)){
          d2[,'W']  <- bsw[,b]
          booted.bound <- miv_bounds_calc(d2, erates, k, np, is.this.a.bootstrap=T)#Flatten matrix to vector to correspond with boot output
          t[b,]<- as.vector(booted.bound)

          d2[,'W']  <- W
        }
        booted.miv.bounds$t <- t
      } else{

        boot_bounds <- function(data, i, ...){
          data2 <- data[i,]
          return(miv_bounds_calc(d2, erates, k, np))
        }
        booted.miv.bounds <- boot::boot(d2, boot_bounds, boot.reps, erates = erates, np=np, k=k, is.this.a.bootstrap=T)
      }
      miv.bounds <- booted.miv.bounds$t0[1:16,]
      bias.miv.bounds <- booted.miv.bounds$t0[17:32,]
      colselect <- rep(1:32, length(erates))
      miv.booted.bounds <- booted.miv.bounds$t[,colselect %in% 1:16]
      #################UPDATE here to match miv_tebounds CI
      stdev <- apply(miv.booted.bounds,2,sd, na.rm=T)
      nn <- colSums(!is.na(miv.booted.bounds))

      if(!im){
        low <-  ((100-CI.level)/2)/100
        high <-  1-low
        l_ci <- apply(miv.booted.bounds, 2, quantile, probs = low, na.rm=T)
        l_ci <- matrix(l_ci, nrow = 16, byrow = F)
        l_ci <- pmax(l_ci[seq(2,16,2),], -1)

        u_ci <- apply(miv.booted.bounds, 2, quantile, probs = low, na.rm=T)
        u_ci <- matrix(u_ci, nrow = 16, byrow = T)
        u_ci <- pmin(u_ci[seq(1,15,2),], 1)

      } else {
        delta <- diff(miv.bounds)[seq(1,15,2),]*-1
        stdev <- matrix(stdev, nrow = 16, byrow=F)
        stdev.u <- stdev[seq(1,15,2),]
        stdev.l <- stdev[seq(2,16,2),]
        s <- 1/pmax(stdev.l, stdev.u)
        c <- matrix(1, nrow=nrow(delta), ncol=ncol(delta))

        for ( i in 1:nrow(delta)){
          for (j in 1:ncol(delta)){
            c[i,j] <- qnorm(1-((100-CI.level)/200))

            diff <- 99999
            while (diff>0.00001){
              c[i,j] <- c[i,j]-0.00001
              d <- pnorm(c[i,j]+s[i,j]*delta[i,j])-pnorm(-c[i,j])-(CI.level/100)
              diff <- min(abs(d), diff)
              if (c[i,j]<1){
                warning(paste('Error with I-M CIs: CbarN down to 1 in ', rownames(delta)[i]))
                break
              }
            }
          }
        }

        l_ci.1 <- pmax(miv.bounds[c(2,4,10,12),]-c[c(1,2,5,6),]*stdev.l[c(1,2,5,6),], -1)
        u_ci.1 <- pmin(miv.bounds[c(1,3, 9,11),]+c[c(1,2,5,6),]*stdev.u[c(1,2,5,6),],  1)

        l_ci.2 <- pmax(miv.bounds[c(6,8,14,16),]-c[c(3,4,7,8),]*stdev.l[c(3,4,7,8),], 0)
        u_ci.2 <- pmin(miv.bounds[c(5,7,13,15),]+c[c(3,4,7,8),]*stdev.u[c(3,4,7,8),], 1)

        l_ci <- rbind(l_ci.1, l_ci.2)
        u_ci <- rbind(u_ci.1, u_ci.2)
      }
      l_ci.miv.bounds <- pmin(l_ci, miv.bounds[c(2,4,6,8,10,12,14,16),])
      u_ci.miv.bounds <- pmax(u_ci, miv.bounds[c(1,3,5,7, 9,11,13,15),])

      #Adjustments to ensure MIV bounds are tighter than MTS/MTR bounds
      miv.bounds[seq(1,15,2)] <- pmin(miv.bounds[seq(1,15,2)], bounds[seq( 9,23,2)])
      miv.bounds[seq(2,16,2)] <- pmax(miv.bounds[seq(2,16,2)], bounds[seq(10,24,2)])


    } else {

      miv.bounds <- miv_bounds_calc(d2, erates, k, np)
      miv.bounds[seq(1,15,2)] <- pmin(miv.bounds[seq(1,15,2)], bounds[seq( 9,23,2)])
      miv.bounds[seq(2,16,2)] <- pmax(miv.bounds[seq(2,16,2)], bounds[seq(10,24,2)])
      if(k>0){

        bias.miv.bounds <- miv.bounds[17:32,]
        miv.bounds <- miv.bounds[1:16,]
      }
    }



  }
  rownames(bounds) <- c('ArbError_ESM_ub','ArbError_ESM_lb',
                        'NFP_ESM_ub', 'NFP_ESM_lb',
                        'ArbError_NMA_ub','ArbError_NMA_lb',
                        'NFP_NMA_ub', 'NFP_NMA_lb',
                        'ArbError_MTS.NS_ub','ArbError_MTS.NS_lb',
                        'NFP_MTS.NS_ub',     'NFP_MTS.NS_lb',
                        'ArbError_MTS.MTR.NS_ub','ArbError_MTS.MTR.NS_lb',
                        'NFP_MTS.MTR.NS_ub',     'NFP_MTS.MTR.NS_lb',

                        'ArbError_MTS.PS_ub','ArbError_MTS.PS_lb',
                        'NFP_MTS.PS_ub',     'NFP_MTS.PS_lb',
                        'ArbError_MTS.MTR.PS_ub','ArbError_MTS.MTR.PS_lb',
                        'NFP_MTS.MTR.PS_ub',     'NFP_MTS.MTR.PS_lb')

  if(bootstrap){
    rownames(l_ci.bounds) <- rownames(bounds)[seq(2,24,2)]
    rownames(u_ci.bounds) <- rownames(bounds)[seq(1,24,2)]
  }

  if (miv) {
    rownames(miv.bounds) <- c('ArbError_MIV.MTS.NS_ub','ArbError_MIV.MTS.NS_lb',
                              'NFP_MIV.MTS.NS_ub',     'NFP_MIV.MTS.NS_lb',
                              'ArbError_MIV.MTS.MTR.NS_ub','ArbError_MIV.MTS.MTR.NS_lb',
                              'NFP_MIV.MTS.MTR.NS_ub',     'NFP_MIV.MTS.MTR.NS_lb',

                              'ArbError_MIV.MTS.PS_ub','ArbError_MIV.MTS.PS_lb',
                              'NFP_MIV.MTS.PS_ub',     'NFP_MIV.MTS.PS_lb',
                              'ArbError_MIV.MTS.MTR.PS_ub','ArbError_MIV.MTS.MTR.PS_lb',
                              'NFP_MIV.MTS.MTR.PS_ub',     'NFP_MIV.MTS.MTR.PS_lb')



    bounds <- rbind(bounds,miv.bounds)

    if(bootstrap){
      rownames(l_ci.miv.bounds) <- rownames(miv.bounds)[seq(2,16,2)]
      rownames(u_ci.miv.bounds) <- rownames(miv.bounds)[seq(1,15,2)]
      l_ci.bounds <- rbind(l_ci.bounds, l_ci.miv.bounds)
      u_ci.bounds <- rbind(u_ci.bounds, l_ci.miv.bounds)
    }

  }
  bound.returns <- list()
  bound.returns$bounds <- bounds
  if(bootstrap){
    bound.returns$lower_ci <- l_ci.bounds
    bound.returns$upper_ci <- u_ci.bounds

  }
  if(miv){
    if(k>0){
      rownames(bias.miv.bounds) <- rownames(miv.bounds)
      bound.returns$miv.bias <- bias.miv.bounds

    }
  }
  bound.returns$call <- cl

  if(bootstrap){
    bound.returns$bsw <- bsw
  }

  class(bound.returns) <- 'tebounds'

  bound.returns
}

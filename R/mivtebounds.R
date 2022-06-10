
miv_bounds_calc <- function(d2, erates, k, np, is.this.a.bootstrap=F ){
    Y <- d2[,'Y']
    D <- d2[,'D']
    W <- d2[,'W']
    Z <- d2[,'Z']

    p_means <- get('p_z', envir = parent.frame())
    names(p_means) <- paste0('z',c(1:length(p_means)))
    data <- get('data', envir = parent.frame())
    formula <- get('formula', envir = parent.frame())

    p_iv_df <- data.frame(D,d1=(Y==1&D==1),
                       d2=(Y==0&D==0),
                       d3=(Y==1&D==0),
                       d4=(Y==0&D==1),
                       Z,
                       W)

    miv.fs.boot <- function(data, i, erates, np, svy.boot=F, is.this.a.bootstrap){
        data2 <- data[i,]

        p_iv <- sapply(split(p_iv_df, data2$Z), function(d) apply(d,2,function(x) weighted.mean(x, w = d$W)))
        p_iv <- p_iv[c('D','d1','d2','d3','d4'),]
        rownames(p_iv) <- c('mDiv','p11iv','p00iv','p10iv','p01iv')
        if(any(p_iv==0)){

          if(all(i == c(1:nrow(data2)))& svy.boot==F  & is.this.a.bootstrap==F) { stop("Error: Some MIV cells do not contain all combinations of {Y,D} in original data. Consider setting fewer cells using the 'ncells' option. ")}
          else {
            warning("Warning: Some MIV cells do not contain all combinations of {Y,D} in a particular bootstrap pseudosample.
                    The program will continue but the number of repetitions will be
                  fewer than the number specified.
                  Consider setting fewer cells using the 'ncells' option.")
            return(NA)
          }
        }

        #Note: There is a difference in results from Stata in the grid search.
        #The Stata module (v1.8) uses a loop to search through the grid, but due to precision
        #issues, it will sometimes not examine the last grid value.  This is because it
        #increments up by 1/np, which can be rounded due to precision and if the rounded value
        #doesn't quite get to the final value, then the loop will skip it.  In this package,
        #instead of looping for the search, the full grid is expanded, then is searched, ensuring
        #all values are included.
        miv.grid.search <- function(x){
          meanD <- x[1]
          p11 <- x[2]
          p00 <- x[3]
          p10 <- x[4]
          p01 <- x[5]


          q <-  length(erates)
          nzr <- (erates!=0)

          t1pos <- pmin(erates, p11)
          t1neg <- pmin(erates, p10)
          t0pos <- pmin(erates, p01)
          t0neg <- pmin(erates, p00)

          #Arbitrary Errors Model: Invoking only A1 (j=1 & i=1)

          ub11_1 <-  rep(-10, q)
          lb11_1 <-  rep(10, q)
          ub11_0 <-  rep(-10, q)
          lb11_0 <-  rep(10, q)
          lb11_1[erates==0] <- ub11_1[erates==0] <- p11/meanD
          lb11_0[erates==0] <- ub11_0[erates==0] <- p10/(1-meanD)

          #Exogenous selection model

          inc_b <- t1pos/(np-1)
          inc_a <- t1neg/(np-1)
          inc_b_temp <- t0neg/(np-1)
          inc_a_temp <- t0pos/(np-1)
          b_list <- mapply(function(x,y)seq(0,x,y), t1pos, inc_b)
          qtempbmax <- mapply(function(x, y, z) pmin(x-y, z), erates, b_list, t0neg)
          tb_list <- mapply(function(a,b)mapply(function(x,y)seq(0,x,y),a, b, SIMPLIFY = F), qtempbmax, inc_b_temp, SIMPLIFY = F)

          r1<- mapply(function(b,tb){
              r1 <- ((p11-b)/(meanD-b+tb))
              r1 <- ifelse(r1>1,NA,r1)
            },
            unlist(b_list), unlist(tb_list, recursive = F)
          )


          r2 <- mapply(function(b,tb){
             r2 <- ((p10+b)/(1-meanD+b-tb))
             r2 <- ifelse(r2>1,NA,r2)
            },
            unlist(b_list), unlist(tb_list, recursive = F)
            )

          ratios <- mapply(function(b,tb){
              r1 <- ((p11-b)/(meanD-b+tb))
              r2 <- ((p10+b)/(1-meanD+b-tb))
              r1 <- ifelse(r1>1,NA,r1)
              r2 <- ifelse(r2>1,NA,r2)
              r <- r1-r2
              r <- ifelse(r< (-1),-1,r)
            },
            unlist(b_list), unlist(tb_list, recursive = F)
          )
          min.index <- lapply(ratios, which.min)
          ratios <- mapply(function(r,i) r[i], ratios, min.index)
          r1 <- mapply(function(r,i) r[i], r1, min.index)
          r2 <- mapply(function(r,i) r[i], r2, min.index)


          if (sum(nzr)!=q)  {
            ratios <- ratios[2:length(ratios)]
            r1 <- r1[2:length(r1)]
            r2 <- r2[2:length(r2)]
          }


          ratios <- split(ratios, cut(seq_along(ratios),sum(nzr),labels=F ))
          r1 <- split(r1, cut(seq_along(r1),sum(nzr),labels=F ))
          r2 <- split(r2, cut(seq_along(r2),sum(nzr),labels=F ))

          min.index <- lapply(ratios, which.min)
          r1 <- mapply(function(r,i) r[i], r1, min.index)
          r2 <- mapply(function(r,i) r[i], r2, min.index)

          lb11_1[nzr] <- r1
          ub11_0[nzr] <- r2

          a_list <- mapply(function(x,y)seq(0,x,y), t1neg, inc_a)
          qtempamax <- mapply(function(x, y, z) pmin(x-y, z), erates, a_list, t0pos)
          ta_list <- mapply(function(a,b)mapply(function(x,y)seq(0,x,y),a, b, SIMPLIFY = F), qtempamax, inc_a_temp, SIMPLIFY = F)

          r1 <- mapply(function(a,ta){
              r1 <- ((p11+a)/(meanD+a-ta))
              r1 <- ifelse(r1>1, NA,r1)
            },
            unlist(a_list), unlist(ta_list, recursive = F)
          )

          r2 <- mapply(function(a,ta){
                r2 <- ((p10-a)/(1-meanD-a+ta))
                r2 <- ifelse(r2>1, NA,r2)
            },
            unlist(a_list), unlist(ta_list, recursive = F)
          )

          ratios <- mapply(function(a,ta){
              r1 <- ((p11+a)/(meanD+a-ta))
              r2 <- ((p10-a)/(1-meanD-a+ta))
              r1 <- ifelse(r1>1, NA,r1)
              r2 <- ifelse(r2>1, NA,r2)
              r <- r1-r2
              r <- ifelse(r>1, 1,r)
            },
            unlist(a_list), unlist(ta_list, recursive = F)
          )

          max.index <- lapply(ratios, which.max)
          ratios <- mapply(function(r,i) r[i], ratios, max.index)
          r1 <- mapply(function(r,i) r[i], r1, max.index)
          r2 <- mapply(function(r,i) r[i], r2, max.index)


          if (sum(nzr)!=q)  {
            ratios <- ratios[2:length(ratios)]
            r1 <- r1[2:length(r1)]
            r2 <- r2[2:length(r2)]
          }

          ratios <- split(ratios, cut(seq_along(ratios),sum(nzr),labels=F ))
          r1 <- split(r1, cut(seq_along(r1),sum(nzr),labels=F ))
          r2 <- split(r2, cut(seq_along(r2),sum(nzr),labels=F ))

          max.index <- lapply(ratios, which.max)
          r1 <- mapply(function(r,i) r[i], r1, max.index)
          r2 <- mapply(function(r,i) r[i], r2, max.index)

          ub11_1[nzr] <- r1
          lb11_0[nzr] <- r2


          #No False Positives Model: Invoking both A1 and A2 => Setting `t1pos' = 0 and `t0pos'= 0 (j=1 & i=2)#####


          ub12_1 <-  ub11_1
          lb12_1 <-  lb11_1
          ub12_0 <-  ub11_0
          lb12_0 <-  lb11_0



          inc_h <- t0neg/(np-1)
          inc_a <- t1neg/(np-1)
          h <-  mapply(function(x,y)seq(0,x,y), t0neg[nzr], inc_h[nzr])
          r1 <- p11/(meanD+h)
          r2 <- p10/(1-meanD-h)
          r1[r1>1] <- NA
          r2[r2>1] <- NA
          ratios <- r1-r2
          ratios[ratios>1] <- NA
          lb12_temp <- apply(ratios,2,which.min)
          lb12_1[nzr] <- r1[cbind(lb12_temp, seq_along(lb12_temp) )]
          ub12_0[nzr] <- r2[cbind(lb12_temp, seq_along(lb12_temp) )]


          #*! Arbitrary Errors Model: Invoking only A1 (j=2 & i=1)######
          ub21_1 <- p11+(1-meanD) + t0pos
          lb21_1 <- p11           - t1pos
          ub21_0 <- p10+meanD     + t0neg
          lb21_0 <- p10           - t1neg

          ub21_1[!nzr] <- p11+(1-meanD)
          lb21_1[!nzr] <- p11
          ub21_0[!nzr] <- p10+meanD
          lb21_0[!nzr] <- p10

          #No False Positives Model: Invoking both A1 and A2 => Setting `t1pos' = 0 and `t0pos'= 0 (j=2 & i=2)#####
          ub22_1 <- p11+(1-meanD) +(t0neg-t0neg)
          lb22_1 <- p11           +(t0neg-t0neg) #add zero vector to make vectors
          ub22_0 <- p10+meanD     + t0neg
          lb22_0 <- p10           - t1neg

          ub22_1[!nzr] <- ub21_1[!nzr]
          lb22_1[!nzr] <- lb21_1[!nzr]
          ub22_0[!nzr] <- ub21_0[!nzr]
          lb22_0[!nzr] <- lb21_0[!nzr]



          #Monotone Treatment Selection (MTSn) Model (j=3)####

          ub31_1 <- pmax(pmin(ub21_1,1),0)
          lb31_1 <- pmin(pmax(lb11_1,0),1)
          ub31_0 <- pmax(pmin(ub11_0,1),0)
          lb31_0 <- pmin(pmax(lb21_0,0),1)

          ub32_1 <- pmax(pmin(ub22_1,1),0)
          lb32_1 <- pmin(pmax(lb12_1,0),1)
          ub32_0 <- pmax(pmin(ub12_0,1),0)
          lb32_0 <- pmin(pmax(lb22_0,0),1)

          # Monotone Treatment Selection (MTSp) Model (j=3p)#####

          ub31p_1 <- pmax(pmin(ub11_1,1),0)
          lb31p_1 <- pmin(pmax(lb21_1,0),1)
          ub31p_0 <- pmax(pmin(ub21_0,1),0)
          lb31p_0 <- pmin(pmax(lb11_0,0),1)

          ub32p_1 <- pmax(pmin(ub12_1,1),0)
          lb32p_1 <- pmin(pmax(lb22_1,0),1)
          ub32p_0 <- pmax(pmin(ub22_0,1),0)
          lb32p_0 <- pmin(pmax(lb12_0,0),1)

          #MTSn and MTR Model (j=5)#####

          ub41_1 <- pmax(pmin(ub21_1,1),0)
          lb41_1 <- pmin(pmax(lb11_1,0),1)
          ub41_0 <- pmax(pmin(ub11_0,1),0)
          lb41_0 <- pmin(pmax(lb21_0,0),1)

          ub42_1 <- pmax(pmin(ub22_1,1),0)
          lb42_1 <- pmin(pmax(lb12_1,0),1)
          ub42_0 <- pmax(pmin(ub12_0,1),0)
          lb42_0 <- pmin(pmax(lb22_0,0),1)
          #MTSp and MTR Model (j=5)#####
          ub41p_1 <- pmax(pmin(ub11_1,1),0)
          lb41p_1 <- pmin(pmax(lb21_1,0),1)
          ub41p_0 <- pmax(pmin(ub21_0,1),0)
          lb41p_0 <- pmin(pmax(lb11_0,0),1)

          ub42p_1 <- pmax(pmin(ub12_1,1),0)
          lb42p_1 <- pmin(pmax(lb22_1,0),1)
          ub42p_0 <- pmax(pmin(ub22_0,1),0)
          lb42p_0 <- pmin(pmax(lb12_0,0),1)



          bnds <- rbind(ub11_1,lb11_1,ub11_0,lb11_0,ub12_1,
                        lb12_1,ub12_0,lb12_0,ub21_1,lb21_1,
                        ub21_0,lb21_0,ub22_1,lb22_1,ub22_0,
                        lb22_0,ub31_1,lb31_1,ub31_0,lb31_0,
                        ub32_1,lb32_1,ub32_0,lb32_0,ub31p_1,
                        lb31p_1,ub31p_0,lb31p_0,ub32p_1,
                        lb32p_1,ub32p_0,lb32p_0,ub41_1,
                        lb41_1, ub41_0, lb41_0, ub42_1,
                        lb42_1, ub42_0, lb42_0, ub41p_1,
                        lb41p_1,ub41p_0,lb41p_0,ub42p_1,
                        lb42p_1,ub42p_0,lb42p_0)

          colnames(bnds) <- paste0('ER',erates)

          return(bnds)
        }
        q <-  length(erates)


        miv.bounds <- apply(p_iv, 2,miv.grid.search)

        bnd.names <- c(
          "ub11_1", "lb11_1","ub11_0", "lb11_0","ub12_1", "lb12_1","ub12_0",
          "lb12_0","ub21_1", "lb21_1","ub21_0", "lb21_0","ub22_1", "lb22_1",
          "ub22_0", "lb22_0","ub31_1", "lb31_1","ub31_0", "lb31_0","ub32_1",
          "lb32_1","ub32_0", "lb32_0","ub31p_1","lb31p_1","ub31p_0","lb31p_0",
          "ub32p_1","lb32p_1","ub32p_0","lb32p_0","ub41_1", "lb41_1", "ub41_0",
          "lb41_0", "ub42_1","lb42_1", "ub42_0","lb42_0", "ub41p_1","lb41p_1","ub41p_0",
          "lb41p_0","ub42p_1","lb42p_1","ub42p_0","lb42p_0")
        miv.bound.names <-     paste(rep(paste0('ER',erates), each = length(bnd.names)), bnd.names, sep = ".")
        rownames(miv.bounds) <- miv.bound.names
        mode(miv.bounds) <- 'numeric'
        miv.bounds.rollmax <- miv.bounds
        miv.bounds.rollmin <- miv.bounds
        for(i in 2:ncol(miv.bounds.rollmax)){
          miv.bounds.rollmax[,i] <- pmax(miv.bounds.rollmax[,i], miv.bounds.rollmax[,i-1])

        }
        for(i in (ncol(miv.bounds.rollmin)-1):1){
          miv.bounds.rollmin[,i] <- pmin(miv.bounds.rollmax[,i], miv.bounds.rollmax[,i+1])

        }
        t.max <- apply(sweep(miv.bounds.rollmax,2, p_means, "*"),1,sum)
        t.min <- apply(sweep(miv.bounds.rollmin,2, p_means, "*"),1,sum)
        t.max <- matrix(t.max, ncol  = q, byrow=F)
        t.min <- matrix(t.min, ncol = q, byrow=F)
        rownames(t.max) <- rownames(t.min) <- bnd.names
        colnames(t.max) <- colnames(t.min) <- paste0('ER',erates)


        t.mat <- t.max[c('lb31_1',  'lb32_1',  'lb41_1','lb42_1',
                         'lb31_0',  'lb32_0',  'lb41_0','lb42_0',
                         'lb31p_1', 'lb32p_1', 'lb41p_1','lb42p_1',
                         'lb31p_0', 'lb32p_0', 'lb41p_0','lb42p_0'),]
        rownames(t.mat) <- c('T11','T21','T31','T41',
                             'T10','T20','T30','T40',
                             'T11p','T21p','T31p','T41p',
                             'T10p','T20p','T30p','T40p')

        u.mat <- t.min[c('ub31_1', 'ub32_1', 'ub41_1','ub42_1',
                         'ub31_0', 'ub32_0', 'ub41_0','ub42_0',
                         'ub31p_1', 'ub32p_1', 'ub41p_1','ub42p_1',
                         'ub31p_0', 'ub32p_0', 'ub41p_0','ub42p_0'),]
        rownames(u.mat) <- c('U11','U21','U31','U41',
                             'U10','U20','U30','U40',
                             'U11p','U21p','U31p','U41p',
                             'U10p','U20p','U30p','U40p')
        first.stage.results <- rbind(t.mat[1:8,], u.mat[1:8,],
                                     t.mat[9:16,],u.mat[9:16,])

        return(first.stage.results)
      }

    #Bias Correction - finish BC, then create a boot function for miv_bounds in main tebounds function. (take out of bounds function)

    if (k>0){
      if(exists('svy.design', parent.frame())){

          bias.bounds <- list()
          svy.design <- get('svy.design', parent.frame())
          bsw.dsgn <- survey::as.svrepdesign(svy.design, type = 'subbootstrap', compress=F, replicates = k)
          bsw <- bsw.dsgn$repweights[complete.cases(data[all.vars(formula)]),]*W
          i <- c(1:nrow(p_iv_df))
          bias.bounds$t0 <- miv.fs.boot(p_iv_df, i=i, erates = erates, np=np, is.this.a.bootstrap = is.this.a.bootstrap)

          t <- matrix(NA, nrow = ncol(bsw), ncol = length(bias.bounds$t0))
          for (b in 1:ncol(bsw)){

            p_iv_df[,'W']<- bsw[,b]
            bias.bound <- miv.fs.boot(p_iv_df, i=i, erates = erates, np=np, svy.boot=T, is.this.a.bootstrap)#Flatten matrix to vector to correspond with boot output
            t[b,] <- as.vector(bias.bound)

            p_iv_df[,'W'] <- W
          }
          bias.bounds$t <- t
          if(all(is.na(t))){
            if(is.this.a.bootstrap){
              return(NA)
            }
            else{
              stop("Error: MIV Bias correction could not be calculated. ")
            }
          }

      } else{

        bias.bounds <- boot::boot(p_iv_df, miv.fs.boot, k, erates = erates, np=np, is.this.a.bootstrap)

      }
      bc.bounds <- matrix(apply(bias.bounds$t, 2, mean, na.rm=T), nrow = nrow(bias.bounds$t0))
      miv.bounds <- 2*bias.bounds$t0 - bc.bounds

      bbias <- -(bias.bounds$t0-bc.bounds)

      bias.ate <- rbind('bias_ub41_ate' = bbias['U11',]-bbias['T10',],
                        'bias_ub42_ate' = bbias['U21',]-bbias['T20',],
                        'bias_ub43_ate' = bbias['U31',]-bbias['T30',],
                        'bias_ub44_ate' = bbias['U41',]-bbias['T40',],

                        'bias_lb41_ate' = bbias['T11',]-bbias['U10',],
                        'bias_lb42_ate' = bbias['T21',]-bbias['U20',],
                        'bias_lb43_ate' = bbias['T31',]-bbias['U30',],
                        'bias_lb44_ate' = bbias['T41',]-bbias['U40',],

                        'bias_ub41p_ate' = bbias['U11p',]-bbias['T10p',],
                        'bias_ub42p_ate' = bbias['U21p',]-bbias['T20p',],
                        'bias_ub43p_ate' = bbias['U31p',]-bbias['T30p',],
                        'bias_ub44p_ate' = bbias['U41p',]-bbias['T40p',],

                        'bias_lb41p_ate' = bbias['T11p',]-bbias['U10p',],
                        'bias_lb42p_ate' = bbias['T21p',]-bbias['U20p',],
                        'bias_lb43p_ate' = bbias['T31p',]-bbias['U30p',],
                        'bias_lb44p_ate' = bbias['T41p',]-bbias['U40p',])
    } else {
      i <- c(1:nrow(p_iv_df))
       miv.bounds<- miv.fs.boot(p_iv_df, i=i, erates = erates, np=np, is.this.a.bootstrap)

    }
    suppressWarnings(if(all( is.na(miv.bounds))){
      if(is.this.a.bootstrap){
        return(NA)
      }
      else{
        stop("Error: Some MIV cells do not contain all combinations of {Y,D} in original data. Consider setting fewer cells using the 'ncells' option. ")
      }
    })
    rownames(miv.bounds) <- c('lb41_1','lb42_1','lb43_1','lb44_1',
                           'lb41_0','lb42_0','lb43_0','lb44_0',
                           'ub41_1','ub42_1','ub43_1','ub44_1',
                           'ub41_0','ub42_0','ub43_0','ub44_0',

                           'lb41p_1','lb42p_1','lb43p_1','lb44p_1',
                           'lb41p_0','lb42p_0','lb43p_0','lb44p_0',
                           'ub41p_1','ub42p_1','ub43p_1','ub44p_1',
                           'ub41p_0','ub42p_0','ub43p_0','ub44p_0')

    miv.bounds.ate <- rbind(

                        'ub41_ate' = pmin(miv.bounds['ub41_1',]-miv.bounds['lb41_0',], 1),
                        'lb41_ate' = pmax(miv.bounds['lb41_1',]-miv.bounds['ub41_0',], -1),
                        'ub42_ate' = pmin(miv.bounds['ub42_1',]-miv.bounds['lb42_0',], 1),
                        'lb42_ate' = pmax(miv.bounds['lb42_1',]-miv.bounds['ub42_0',], -1),
                        'ub43_ate' = pmin(miv.bounds['ub43_1',]-miv.bounds['lb43_0',], 1),
                        'lb43_ate' = pmax(miv.bounds['lb43_1',]-miv.bounds['ub43_0',], 0),
                        'ub44_ate' = pmin(miv.bounds['ub44_1',]-miv.bounds['lb44_0',], 1),
                        'lb44_ate' = pmax(miv.bounds['lb44_1',]-miv.bounds['ub44_0',], 0),


                        'ub41p_ate' = pmin(miv.bounds['ub41p_1',]-miv.bounds['lb41p_0',], 1),
                        'lb41p_ate' = pmax(miv.bounds['lb41p_1',]-miv.bounds['ub41p_0',], -1),
                        'ub42p_ate' = pmin(miv.bounds['ub42p_1',]-miv.bounds['lb42p_0',], 1),
                        'lb42p_ate' = pmax(miv.bounds['lb42p_1',]-miv.bounds['ub42p_0',], -1),
                        'ub43p_ate' = pmin(miv.bounds['ub43p_1',]-miv.bounds['lb43p_0',], 1),
                        'lb43p_ate' = pmax(miv.bounds['lb43p_1',]-miv.bounds['ub43p_0',], 0),
                        'ub44p_ate' = pmin(miv.bounds['ub44p_1',]-miv.bounds['lb44p_0',], 1),
                        'lb44p_ate' = pmax(miv.bounds['lb44p_1',]-miv.bounds['ub44p_0',], 0)

                        )
    if (k>0){
      miv.returns <- rbind(miv.bounds.ate, bias.ate)
    } else{
      miv.returns <- rbind(miv.bounds.ate)
    }
    return(miv.returns)

}

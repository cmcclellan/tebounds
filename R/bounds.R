

bound_calc <- function(dat, erates, np,...){


   Y <- dat[,'Y']
   D <- dat[,'D']
   W <- dat[,'W']
   Z <- dat[,'Z']


   #Sample proportion treated
   meanD <- weighted.mean(D, W)
   #getting the joint probabilities of the various outcomes and the treatment variable
   p11 <- weighted.mean((Y==1&D==1), W)
   p00 <- weighted.mean((Y==0&D==0), W)
   p10 <- weighted.mean((Y==1&D==0), W)
   p01 <- weighted.mean((Y==0&D==1), W)


   #Assumption A1: arbitrary classification errors
   if(!is.numeric(erates)){stop('erates must be numeric')}
   erates[erates>1] <- erates[erates>1]/100
   if(any(erates<0|erates>1)){stop('erates out of range')}

   t1pos <- pmin(erates, p11)
   t1neg <- pmin(erates, p10)
   t0pos <- pmin(erates, p01)
   t0neg <- pmin(erates, p00)


   ####### Getting bounds for ATE under different selection model#########

   #Exogenous selection model
   #Arbitrary Errors Model: Invoking only A1 (j=1 & i=1)
   q <-  length(erates)
   nzr <- (erates!=0)

   ub11_ate <-  rep(-10, q)
   lb11_ate <-  rep(10, q)
   lb11_ate[erates==0] <- ub11_ate[erates==0] <- p11/meanD - p10/(1-meanD)

   # for(i in which(nzr)) {
   #   inc_b <- t1pos[i]/(np-1)
   #   inc_a <- t1neg[i]/(np-1)
   #   inc_b_temp <- t0neg[i]/(np-1)
   #   inc_a_temp <- t0pos[i]/(np-1)
   #   for (b  in seq(0, t1pos[i], inc_b)){
   #     Q_temp_b_max <- min(erates[i]-b, t0neg[i])
   #     for (tb in seq(0, Q_temp_b_max, inc_b_temp)){
   #       ratio1b <-( p11-b)/(meanD-b+tb)
   #       ratio2b <-( p10+b)/(1-meanD+b-tb)
   #       if (ratio1b<=1&ratio2b<=1){
   #         lb11_ate[i] <- min(lb11_ate[i], ratio1b-ratio2b)
   #
   #
   #       }
   #
   #     }
   #   }
   #
   #   for (a  in seq(0, t1neg[i], inc_a)){
   #     Q_temp_a_max <- min(erates[i]-a, t0pos[i])
   #     for (ta in seq(0, Q_temp_a_max, inc_a_temp)){
   #       ratio1a <-( p11+a)/(meanD+a-ta)
   #       ratio2a <-( p10-a)/(1-meanD-a+ta)
   #       if (ratio1a<=1&ratio2a<=1){
   #         ub11_ate[i] <- max(ub11_ate[i], ratio1a-ratio2a)
   #       }
   #     }
   #   }
   # lb11_ate[i] <- max(lb11_ate[i],-1)
   # ub11_ate[i] <- min(ub11_ate[i],1)
   # }

   inc_b <- t1pos/(np-1)
   inc_a <- t1neg/(np-1)
   inc_b_temp <- t0neg/(np-1)
   inc_a_temp <- t0pos/(np-1)

   b_list <- mapply(function(x,y)seq(0,x,y), t1pos, inc_b)
   qtempbmax <- mapply(function(x, y, z) pmin(x-y, z), erates, b_list, t0neg)
   tb_list <- mapply(function(a,b)mapply(function(x,y)seq(0,x,y),a, b, SIMPLIFY = F), qtempbmax, inc_b_temp, SIMPLIFY = F)

   ratios <- mapply(function(b,tb){r1 <- ((p11-b)/(meanD-b+tb))
   r2 <- ((p10+b)/(1-meanD+b-tb))
   r1[r1>1] <- NA
   r2[r2>1] <- NA
   r <- r1-r2
   },
   unlist(b_list), unlist(tb_list, recursive = F))
   ratios <- unlist(lapply(ratios, min, na.rm = T))

   if (sum(nzr)!=q)  {ratios <- ratios[2:length(ratios)]}

   ratios <- split(ratios, cut(seq_along(ratios),sum(nzr),labels=F ))

   lb11_ate[nzr] <-  unlist(pmax(lapply(ratios, min, na.rm = T), -1))


   a_list <- mapply(function(x,y)seq(0,x,y), t1neg, inc_a)
   qtempamax <- mapply(function(x, y, z) pmin(x-y, z), erates, a_list, t0pos)
   ta_list <- mapply(function(a,b)mapply(function(x,y)seq(0,x,y),a, b, SIMPLIFY = F), qtempamax, inc_a_temp, SIMPLIFY = F)

   ratios <- mapply(function(a,ta){r1 <- ((p11+a)/(meanD+a-ta))
     r2 <- ((p10-a)/(1-meanD-a+ta))
     r1[r1>1] <- NA
     r2[r2>1] <- NA
     r <- r1-r2
     },
     unlist(a_list), unlist(ta_list, recursive = F)
    )
   ratios <- suppressWarnings(unlist(lapply(ratios, max, na.rm=T)))

   if (sum(nzr)!=q)  {ratios <- ratios[2:length(ratios)]}

   ratios <- split(ratios, cut(seq_along(ratios),sum(nzr),labels=F ))

   ub11_ate[nzr] <- unlist(pmin(lapply(ratios, max, na.rm=T), 1))



   #for (a  in seq(0, t1neg[i], inc_a)){
   #     Q_temp_a_max <- min(erates[i]-a, t0pos[i])
   #     for (ta in seq(0, Q_temp_a_max, inc_a_temp)){
   #       ratio1a <-( p11+a)/(meanD+a-ta)
   #       ratio2a <-( p10-a)/(1-meanD-a+ta)
   #       if (ratio1a<=1&ratio2a<=1){
   #         ub11_ate[i] <- max(ub11_ate[i], ratio1a-ratio2a)
   #       }
   #     }
   #   }


   #No False Positives Model: Invoking both A1 and A2 => Setting `t1pos' = 0 and `t0pos'= 0 (j=1 & i=2)#####


   ub12_ate <-  ub11_ate
   lb12_ate <-  lb11_ate


   ub12_ate[nzr] <- -10
   lb12_ate[nzr] <- 10
   inc_h <- t0neg/(np-1)
   inc_a <- t1neg/(np-1)
   h <-  mapply(function(x,y)seq(0,x,y), t0neg[nzr], inc_h[nzr])
   r1 <- p11/(meanD+h)
   r2 <- p10/(1-meanD-h)
   r1[r1>1] <- NA
   r2[r2>1] <- NA
   ratios <- r1-r2
   lb12_ate[nzr] <- apply(ratios,2,min, na.rm=T)


   a <-  mapply(function(x,y)seq(0,x,y), t1neg[nzr], inc_a[nzr])
   r1 <- (p11+a)/(meanD+a)
   r2 <- (p10-a)/(1-meanD-a)
   r1[r1>1] <- NA
   r2[r2>1] <- NA
   ratios <- r1-r2
   ub12_ate[nzr] <- apply(ratios,2,max, na.rm=T)

   lb12_ate[nzr] <- pmax(lb12_ate[nzr], -1)
   ub12_ate[nzr] <- pmin(ub12_ate[nzr], 1)



   # ub12_ate <-  ub11_ate
   # lb12_ate <-  lb11_ate
   #
   # for(i in which(nzr)) {
   #   ub12_ate[i] <- -10
   #   lb12_ate[i] <- 10
   #   inc_h <- t0neg[i]/(np-1)
   #   inc_a <- t1neg[i]/(np-1)
   #   for(h in seq(0,t0neg[i],by=inc_h)){
   #     ratio1h <- p11/(meanD+h)
   #     ratio2h <- p10/(1-meanD-h)
   #     if (ratio1h<=1 & ratio2h<=1){
   #       lb12_ate[i] <- min(lb12_ate[i], ratio1h-ratio2h)
   #     }
   #   }
   #   for (a in seq(0, t1neg[i], inc_a)){
   #     ratio1a <- (p11+a)/(meanD+a)
   #     ratio2a <- (p10-a)/(1-meanD-a)
   #     if (ratio1a<=1 &ratio2a<=1){
   #       ub12_ate[i] <- max(ub12_ate[i], ratio1a-ratio2a)
   #     }
   #   }
   #   lb12_ate[i] <- max(lb12_ate[i], -1)
   #   ub12_ate[i] <- min(ub12_ate[i], 1)
   #
   # }
   #
   #*! Worst Case Selection Model (j=2)

   #*! Arbitrary Errors Model: Invoking only A1 (j=2 & i=1)######

   ub21_ate <- p11+(1-meanD)-p10
   lb21_ate <- p11-p10-meanD

   ub21_ate[nzr] <- p11 +(1-meanD) - p10 + pmin(erates[nzr], t0pos[nzr]+t1neg[nzr])
   lb21_ate[nzr] <- p11 - p10 - meanD - pmin(erates[nzr], t1pos[nzr]+t0neg[nzr])

   ub21_ate[nzr] <- pmin(ub21_ate[nzr],  1)
   lb21_ate[nzr] <- pmax(lb21_ate[nzr], -1)

   #No False Positives Model: Invoking both A1 and A2 => Setting `t1pos' = 0 and `t0pos'= 0 (j=2 & i=2)#####
   ub22_ate <- ub21_ate
   lb22_ate <- lb21_ate

   ub22_ate[nzr] <- p11 +(1-meanD) - p10 + t1neg[nzr]
   lb22_ate[nzr] <- p11 - p10 - meanD - t0neg[nzr]

   ub22_ate[nzr] <- pmin(ub22_ate[nzr],  1)
   lb22_ate[nzr] <- pmax(lb22_ate[nzr], -1)

   #Monotone Treatment Selection (MTSn) Model (j=3)####

   ub31_ate <- ub21_ate
   lb31_ate <- lb11_ate
   ub32_ate <- ub22_ate
   lb32_ate <- lb12_ate

   # Monotone Treatment Selection (MTSp) Model (j=3p)#####

   ub31p_ate <- ub11_ate
   lb31p_ate <- lb21_ate
   ub32p_ate <- ub12_ate
   lb32p_ate <- lb22_ate

   #MTSn and MTR Model (j=5)#####

   ub51_ate <- ub31_ate
   lb51_ate <- pmax(lb31_ate, 0)
   ub52_ate <- ub32_ate
   lb52_ate <- pmax(lb32_ate, 0)

   #MTSp and MTR Model (j=5)#####
   ub51p_ate <- ub31p_ate
   lb51p_ate <- pmax(lb31p_ate, 0)
   ub52p_ate <- ub32p_ate
   lb52p_ate <- pmax(lb32p_ate, 0)



   bnds <- rbind(ub11_ate , lb11_ate , ub12_ate , lb12_ate ,
                 ub21_ate , lb21_ate , ub22_ate , lb22_ate ,
                 ub31_ate , lb31_ate , ub32_ate , lb32_ate ,
                 ub51_ate , lb51_ate , ub52_ate , lb52_ate ,
                 ub31p_ate, lb31p_ate, ub32p_ate, lb32p_ate,ub51p_ate,
                 lb51p_ate, ub52p_ate, lb52p_ate)

   colnames(bnds) <- paste0('ER',erates)


   return(bnds)

}






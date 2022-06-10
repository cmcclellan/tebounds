#' Summarize tebound object
#'
#' @param x An object of the class "tebounds".
#' @param ndigits Number of digits to round to.
#'
#' @return A summary of treatment effect bounds
#' @exportS3Method summary
#' @exportS3Method print
#'
#'
#' @keywords internal
#'
summary.tebounds <- function(x, ndigits = 3){



  b.x <- t(x$bounds)
  b.x <- as.data.frame(round(b.x, ndigits))
  b.x$index <- c(1:nrow(b.x))
  b.x[,'Error Rates'] <- gsub('ER', "",rownames(b.x))
  b.x[,' '] <- 'p.e.'

  if ('lower_ci' %in% names(x)){
    #ci.x <- cbind(t(x$lower_ci[seq(2,nrow(x$lower_ci),2),]),t(x$upper_ci[seq(1,nrow(x$upper_ci),2),]))
    ci.x <- cbind(t(x$lower_ci),t(x$upper_ci))
    ci.x <- as.data.frame(round(ci.x, ndigits))

    ci.x$index <- c(1:nrow(ci.x))+.1
    ci.x[,'Error Rates'] <- ""
    ci.x[,' '] <- 'CI'

    b.x <- rbind(b.x, ci.x)
  }

  if ('miv.bias' %in% names(x)){
    bias.x <- (t(x$miv.bias))
    bias.x <- as.data.frame(round(bias.x, ndigits))
    rownames(bias.x) <- NULL
    bias.x$index <- c(1:nrow(bias.x))+.2
    bias.x[,'Error Rates'] <- ""
    bias.x[setdiff(names(b.x), names(bias.x))] <- NA
    bias.x[,' '] <- 'Bias'
    b.x <- rbind(b.x, bias.x)

  }

  b.x <- b.x[order(b.x$index),]


  clean.cell <- function(a){
    arberr_lb <- paste0('ArbError_',a,'_lb')
    arberr_ub <- paste0('ArbError_',a,'_ub')
    nfp_lb <- paste0('NFP_',a,'_lb')
    nfp_ub <- paste0('NFP_',a,'_ub')
    b.x[,"Arbitrary Errors"] <- paste0("[", formatC(b.x[,arberr_lb], format = 'f', flag = '0', digits = ndigits),
                                       ', ', formatC(b.x[,arberr_ub], format = 'f', flag = '0', digits = ndigits),']')
    b.x[,"No False Positives"] <- paste0("[", formatC(b.x[,nfp_lb], format = 'f', flag = '0', digits = ndigits),
                                         ', ', formatC(b.x[,nfp_ub], format = 'f', flag = '0', digits = ndigits),']')

    p.out <- b.x[!is.na(b.x[arberr_lb]),c('Error Rates','Arbitrary Errors','No False Positives', ' ')]

    p.out[grepl('\\[0', p.out[,2]),2] <- gsub('\\[0', '\\[ 0',p.out[grepl('\\[0', p.out[,2]),2])
    p.out[grepl('\\[0', p.out[,3]),3] <- gsub('\\[0', '\\[ 0',p.out[grepl('\\[0', p.out[,3]),3])

    p.out[nchar(p.out[,2])<16,2] <- gsub(',',', ', p.out[nchar(p.out[,2])<16,2])
    p.out[nchar(p.out[,3])<16,3] <- gsub(',',', ', p.out[nchar(p.out[,3])<16,3])


    if (grepl('MIV', a)) {
      p.out[p.out[,4]=='Bias',2] <- gsub("\\[|,|\\]", " ",p.out[p.out[,4]=='Bias',2])
      p.out[p.out[,4]=='Bias',3] <- gsub("\\[|,|\\]", " ",p.out[p.out[,4]=='Bias',3])
    }
    p.out
  }
  #####Exogenous Selection Model####
  p.out <- clean.cell('ESM')

  comment(p.out) <- 'Exogenous Selection Model'

  ESM <- p.out
  #####No Monotonicity Assumptions (Worst Case Selection)####
  p.out <- clean.cell('NMA')

  comment(p.out) <- ('No Monotonicity Assumptions (Worst Case Selection)')

  NMA <- p.out
#####MTS Assumption: Negative Selection####
  p.out <- clean.cell('MTS.NS')

  comment(p.out) <- ('MTS Assumption: Negative Selection')

  MTS.NS <- p.out
#####MTS Assumption: Positive Selection####
  p.out <- clean.cell('MTS.PS')

  comment(p.out) <- ('MTS Assumption: Positive Selection')

  MTS.PS <- p.out
#####MTS and MTR Assumptions: Negative Selection####
  p.out <- clean.cell('MTS.MTR.NS')

  comment(p.out) <- ('MTS and MTR Assumptions: Negative Selection')

  MTS.MTR.NS <- p.out
#####MTS and MTR Assumptions: Positive Selection####
  p.out <- clean.cell('MTS.MTR.PS')

  comment(p.out) <- ('MTS and MTR Assumptions: Positive Selection')
  MTS.MTR.PS <- p.out

  if (any(grepl('MIV', names(b.x)))){

    #####MIV, MTS and MTR Assumptions: Negative Selection####
    p.out <- clean.cell('MIV.MTS.NS')

    comment(p.out) <- ('MIV and MTS Assumptions: Negative Selection')

    MIV.MTS.NS <- p.out
    #####MTS and MTR Assumptions: Positive Selection####
    p.out <- clean.cell('MIV.MTS.PS')

    comment(p.out) <- ('MIV and MTS Assumptions: Positive Selection')

    MIV.MTS.PS <- p.out
#####MIV, MTS and MTR Assumptions: Negative Selection####
    p.out <- clean.cell('MIV.MTS.MTR.NS')

    comment(p.out) <- ('MIV, MTS, and MTR Assumptions: Negative Selection')

    MIV.MTS.MTR.NS <- p.out
#####MTS and MTR Assumptions: Positive Selection####
    p.out <- clean.cell('MIV.MTS.MTR.PS')

    comment(p.out) <- ('MIV, MTS, and MTR Assumptions: Positive Selection')

    MIV.MTS.MTR.PS <- p.out


    rlist <- list('ESM' = ESM, 'NMA'=NMA, 'MTS.PS'= MTS.PS, 'MTS.NS'=MTS.NS,
                  'MTS.MTR.PS' = MTS.MTR.PS, 'MTS.MTR.NS'=MTS.MTR.NS,
                  'MIV.MTS.PS' = MIV.MTS.PS,
                  'MIV.MTS.NS' = MIV.MTS.NS,
                  'MIV.MTS.MTR.PS' = MIV.MTS.MTR.PS,
                  'MIV.MTS.MTR.NS' = MIV.MTS.MTR.NS,
                  call = x$call)

  } else {
    rlist <- list('ESM' = ESM, 'NMA'=NMA, 'MTS.PS'= MTS.PS, 'MTS.NS'=MTS.NS,
                  'MTS.MTR.PS' = MTS.MTR.PS, 'MTS.MTR.NS'=MTS.MTR.NS, call = x$call)
  }

  class(rlist) <- 'summary.tebounds'
  rlist

}

print.summary.tebounds <- function(x){

  cat("\nCall:\n", # S has ' ' instead of '\n'
      paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep = "")
  plist <- names(x)[-grep('call',names(x))]
  for(n in plist){
    t <- x[[n]]
    cat('\n',comment(t),'\n')
    cat(paste(rep("-",64), collapse= ""), '\n')
    print(t, row.names=F)
    cat(paste(rep("-",64), collapse= ""), '\n')
  }


}

combine_bounds <- function(x, select.bounds = c("ESM","NMA", "MTS.PS", "MTS.NS", "MTS.MTR.PS",
                                                "MTS.MTR.NS","MIV.MTS.PS", "MIV.MTS.NS", "MIV.MTS.MTR.PS",
                                                "MIV.MTS.MTR.NS"), Bias = T, CI = T){

  sum <- summary(x)
  sum <- sum[select.bounds]
  panels <- lapply(sum, function(x){
    if (!Bias){
      x <- x[!(grepl('Bias',x[,4])),]
    }
    if (!CI){
      x <- x[!(grepl('CI',x[,4])),]
    }
    rbind(c(comment(x), rep("",ncol(x)-1)), x)})
  panels <- do.call('rbind', panels)
  panels
}


#' Plot Treatment Effect Bounds
#'
#' @param x An object of the class "tebounds"
#'
#' @return A list of plots of treatmetn effect bounds
#' @exportS3Method plot
#'
#'
plot.tebounds <- function(x){
  suppressPackageStartupMessages({
    requireNamespace('ggplot2')
  })

  p.data <- as.data.frame(t(x$bounds))
  p.data[,'ER'] <- as.numeric(gsub('ER', "",rownames(p.data)))


  ####ESM and NMA####
  ESM_NMA <- ggplot(p.data, aes(x = .data$ER)) +
    geom_ribbon(aes(ymin = .data$ArbError_ESM_lb,
                    ymax = .data$ArbError_ESM_ub, color = 'a',linetype = 'a'
                    ),
                alpha = 0) +
    geom_ribbon(aes(ymin = .data$NFP_ESM_lb,
                    ymax = .data$NFP_ESM_ub, linetype = 'b' , color='b'),
                alpha = 0) +
    geom_ribbon(aes(ymin = .data$ArbError_NMA_lb,
                    ymax = .data$ArbError_NMA_ub, linetype = 'c' , color = 'c'),
                alpha = 0) +
    geom_ribbon(aes(ymin = .data$NFP_NMA_lb,
                    ymax = .data$NFP_NMA_ub, linetype = 'd' , color='d'),
                alpha = 0) +
    scale_x_continuous(breaks=p.data$ER)+
  scale_linetype_manual("",labels = c('Exogenous: Arbitrary Errors',
                                        'Exogenous: No False Positives',
                                        'No Selection: Arbitrary Errors',
                                        'No Selection: No False Positives'),
                          values = c("a" = "solid", "b" = "longdash", 'c'='dotdash','d'='dashed'))+
    scale_color_manual(name = "", labels = c('Exogenous: Arbitrary Errors',
                                             'Exogenous: No False Positives',
                                             'No Selection: Arbitrary Errors',
                                             'No Selection: No False Positives'),

                       values = c("a" = "blue4", "b" = "red", 'c'='green4','d'='grey50'))+


    theme_minimal()+
    ggtitle('Exogenous and No Selection Assumptions')+
    xlab('Maximum Allowed Degree of Misclassification')+
    ylab('ATE')+
    theme(legend.position = 'bottom',
          panel.grid.minor = element_blank())+
    guides(linetype=guide_legend(nrow=2,byrow=TRUE))







  ####MTS.NS####
  AE_MTS.NS <- ggplot(p.data, aes(x = .data$ER)) +
    geom_ribbon(aes(ymin = .data$ArbError_MTS.NS_lb,
                    ymax = .data$ArbError_MTS.NS_ub, color = 'a',linetype = 'a'
    ),
    alpha = 0)
  if("ArbError_MIV.MTS.NS_lb" %in% colnames(p.data)){
    AE_MTS.NS <- AE_MTS.NS +
      geom_ribbon(aes(ymin = .data$ArbError_MIV.MTS.NS_lb,
                      ymax = .data$ArbError_MIV.MTS.NS_ub, color = 'b',linetype = 'b'),
      alpha = 0)+
      scale_linetype_manual("",labels = c('MTS Alone',
                                          'Joint MTS & MIV'),
                            values = c("a" = "dotted", "b" = "longdash"))+
      scale_color_manual(name = "", labels = c('MTS Alone',
                                               'Joint MTS & MIV'),

                         values = c("a" = "blue4", "b" = "red"))
  } else {
    AE_MTS.NS <- AE_MTS.NS +
      scale_linetype_manual("",labels = c('MTS Alone'),
                            values = c("a" = "solid"))+
      scale_color_manual(name = "", labels = c('MTS Alone'),

                         values = c("a" = "blue4"))
  }

  AE_MTS.NS <- AE_MTS.NS +
    scale_x_continuous(breaks=p.data$ER)+
    theme_minimal()+
    ggtitle(label = 'MTSn & MTSn-IV Assumptions', subtitle = 'Arbitrary Errors')+

    xlab('Maximum Allowed Degree of Misclassification')+
    ylab('ATE')+
    theme(legend.position = 'bottom',
          panel.grid.minor = element_blank())+
    guides(linetype=guide_legend(nrow=1,byrow=TRUE))




  NFP_MTS.NS <- ggplot(p.data, aes(x = .data$ER)) +
    geom_ribbon(aes(ymin = .data$NFP_MTS.NS_lb,
                    ymax = .data$NFP_MTS.NS_ub, color = 'a',linetype = 'a'
    ),
    alpha = 0)
  if("ArbError_MIV.MTS.NS_lb" %in% colnames(p.data)){
    NFP_MTS.NS <- NFP_MTS.NS +
      geom_ribbon(aes(ymin = .data$NFP_MIV.MTS.NS_lb,
                      ymax = .data$NFP_MIV.MTS.NS_ub, color = 'b',linetype = 'b'),
                  alpha = 0)+
      scale_linetype_manual("",labels = c('MTS Alone',
                                          'Joint MTS & MIV'),
                            values = c("a" = "dotted", "b" = "longdash"))+
      scale_color_manual(name = "", labels = c('MTS Alone',
                                               'Joint MTS & MIV'),

                         values = c("a" = "blue4", "b" = "red"))
  } else {
    NFP_MTS.NS <- NFP_MTS.NS +
      scale_linetype_manual("",labels = c('MTS Alone'),
                            values = c("a" = "solid"))+
      scale_color_manual(name = "", labels = c('MTS Alone'),

                         values = c("a" = "blue4"))
  }

  NFP_MTS.NS <- NFP_MTS.NS +
    scale_x_continuous(breaks=p.data$ER)+
    theme_minimal()+
    ggtitle(label = 'MTSn & MTSn-IV Assumptions', subtitle = 'No False Positives')+

    xlab('Maximum Allowed Degree of Misclassification')+
    ylab('ATE')+
    theme(legend.position = 'bottom',
          panel.grid.minor = element_blank())+
    guides(linetype=guide_legend(nrow=1,byrow=TRUE))


  ####MTS.MTR.NS####
  AE_MTS.MTR.NS <- ggplot(p.data, aes(x = .data$ER)) +
    geom_ribbon(aes(ymin = .data$ArbError_MTS.MTR.NS_lb,
                    ymax = .data$ArbError_MTS.MTR.NS_ub, color = 'a',linetype = 'a'
    ),
    alpha = 0)
  if("ArbError_MIV.MTS.MTR.NS_lb" %in% colnames(p.data)){
    AE_MTS.MTR.NS <- AE_MTS.MTR.NS +
      geom_ribbon(aes(ymin = .data$ArbError_MIV.MTS.MTR.NS_lb,
                      ymax = .data$ArbError_MIV.MTS.MTR.NS_ub, color = 'b',linetype = 'b'),
                  alpha = 0)+
      scale_linetype_manual("",labels = c('Joint MTS & MTR',
                                          'Joint MTS, MTR, & MIV'),
                            values = c("a" = "dotted", "b" = "longdash"))+
      scale_color_manual(name = "", labels = c('Joint MTS & MTR',
                                               'Joint MTS, MTR, & MIV'),

                         values = c("a" = "blue4", "b" = "red"))
  } else {
    AE_MTS.MTR.NS <- AE_MTS.MTR.NS +
      scale_linetype_manual("",labels = c('Joint MTS & MTR'),
                            values = c("a" = "solid"))+
      scale_color_manual(name = "", labels = c('Joint MTS & MTR'),

                         values = c("a" = "blue4"))
  }

  AE_MTS.MTR.NS <- AE_MTS.MTR.NS +
    scale_x_continuous(breaks=p.data$ER)+
    theme_minimal()+
    ggtitle(label = 'MTSn-MTR & MTSn-MTR-IV Assumptions', subtitle = 'Arbitrary Errors')+

    xlab('Maximum Allowed Degree of Misclassification')+
    ylab('ATE')+
    theme(legend.position = 'bottom',
          panel.grid.minor = element_blank())+
    guides(linetype=guide_legend(nrow=1,byrow=TRUE))




  NFP_MTS.MTR.NS <- ggplot(p.data, aes(x = .data$ER)) +
    geom_ribbon(aes(ymin = .data$NFP_MTS.MTR.NS_lb,
                    ymax = .data$NFP_MTS.MTR.NS_ub, color = 'a',linetype = 'a'
    ),
    alpha = 0)
  if("ArbError_MIV.MTS.MTR.NS_lb" %in% colnames(p.data)){
    NFP_MTS.MTR.NS <- NFP_MTS.MTR.NS +
      geom_ribbon(aes(ymin = .data$NFP_MIV.MTS.MTR.NS_lb,
                      ymax = .data$NFP_MIV.MTS.MTR.NS_ub, color = 'b',linetype = 'b'),
                  alpha = 0)+
      scale_linetype_manual("",labels = c('Joint MTS & MTR',
                                          'Joint MTS, MTR, & MIV'),
                            values = c("a" = "dotted", "b" = "longdash"))+
      scale_color_manual(name = "", labels = c('Joint MTS & MTR',
                                               'Joint MTS, MTR, & MIV'),

                         values = c("a" = "blue4", "b" = "red"))
  } else {
    NFP_MTS.MTR.NS <- NFP_MTS.MTR.NS +
      scale_linetype_manual("",labels = c('Joint MTS & MTR'),
                            values = c("a" = "solid"))+
      scale_color_manual(name = "", labels = c('Joint MTS & MTR'),

                         values = c("a" = "blue4"))
  }

  NFP_MTS.MTR.NS <- NFP_MTS.MTR.NS +
    scale_x_continuous(breaks=p.data$ER)+
    theme_minimal()+
    ggtitle(label = 'MTSn, MTR, & MTSn-IV Assumptions', subtitle = 'No False Positives')+

    xlab('Maximum Allowed Degree of Misclassification')+
    ylab('ATE')+
    theme(legend.position = 'bottom',
          panel.grid.minor = element_blank())+
    guides(linetype=guide_legend(nrow=1,byrow=TRUE))





  ####MTS.PS####
  AE_MTS.PS <- ggplot(p.data, aes(x = .data$ER)) +
    geom_ribbon(aes(ymin = .data$ArbError_MTS.PS_lb,
                    ymax = .data$ArbError_MTS.PS_ub, color = 'a',linetype = 'a'
    ),
    alpha = 0)
  if("ArbError_MIV.MTS.PS_lb" %in% colnames(p.data)){
    AE_MTS.PS <- AE_MTS.PS +
      geom_ribbon(aes(ymin = .data$ArbError_MIV.MTS.PS_lb,
                      ymax = .data$ArbError_MIV.MTS.PS_ub, color = 'b',linetype = 'b'),
                  alpha = 0)+
      scale_linetype_manual("",labels = c('MTS Alone',
                                          'Joint MTS & MIV'),
                            values = c("a" = "dotted", "b" = "longdash"))+
      scale_color_manual(name = "", labels = c('MTS Alone',
                                               'Joint MTS & MIV'),

                         values = c("a" = "blue4", "b" = "red"))
  } else {
    AE_MTS.PS <- AE_MTS.PS +
      scale_linetype_manual("",labels = c('MTS Alone'),
                            values = c("a" = "solid"))+
      scale_color_manual(name = "", labels = c('MTS Alone'),

                         values = c("a" = "blue4"))
  }

  AE_MTS.PS <- AE_MTS.PS +
    scale_x_continuous(breaks=p.data$ER)+
    theme_minimal()+
    ggtitle(label = 'MTSp & MTSp-IV Assumptions', subtitle = 'Arbitrary Errors')+

    xlab('Maximum Allowed Degree of Misclassification')+
    ylab('ATE')+
    theme(legend.position = 'bottom',
          panel.grid.minor = element_blank())+
    guides(linetype=guide_legend(nrow=1,byrow=TRUE))




  NFP_MTS.PS <- ggplot(p.data, aes(x = .data$ER)) +
    geom_ribbon(aes(ymin = .data$NFP_MTS.PS_lb,
                    ymax = .data$NFP_MTS.PS_ub, color = 'a',linetype = 'a'
    ),
    alpha = 0)
  if("ArbError_MIV.MTS.PS_lb" %in% colnames(p.data)){
    NFP_MTS.PS <- NFP_MTS.PS +
      geom_ribbon(aes(ymin = .data$NFP_MIV.MTS.PS_lb,
                      ymax = .data$NFP_MIV.MTS.PS_ub, color = 'b',linetype = 'b'),
                  alpha = 0)+
      scale_linetype_manual("",labels = c('MTS Alone',
                                          'Joint MTS & MIV'),
                            values = c("a" = "dotted", "b" = "longdash"))+
      scale_color_manual(name = "", labels = c('MTS Alone',
                                               'Joint MTS & MIV'),

                         values = c("a" = "blue4", "b" = "red"))
  } else {
    NFP_MTS.PS <- NFP_MTS.PS +
      scale_linetype_manual("",labels = c('MTS Alone'),
                            values = c("a" = "solid"))+
      scale_color_manual(name = "", labels = c('MTS Alone'),

                         values = c("a" = "blue4"))
  }

  NFP_MTS.PS <- NFP_MTS.PS +
    scale_x_continuous(breaks=p.data$ER)+
    theme_minimal()+
    ggtitle(label = 'MTSp & MTSp-IV Assumptions', subtitle = 'No False Positives')+

    xlab('Maximum Allowed Degree of Misclassification')+
    ylab('ATE')+
    theme(legend.position = 'bottom',
          panel.grid.minor = element_blank())+
    guides(linetype=guide_legend(nrow=1,byrow=TRUE))


  ####MTS.MTR.PS####
  AE_MTS.MTR.PS <- ggplot(p.data, aes(x = .data$ER)) +
    geom_ribbon(aes(ymin = .data$ArbError_MTS.MTR.PS_lb,
                    ymax = .data$ArbError_MTS.MTR.PS_ub, color = 'a',linetype = 'a'
    ),
    alpha = 0)
  if("ArbError_MIV.MTS.MTR.PS_lb" %in% colnames(p.data)){
    AE_MTS.MTR.PS <- AE_MTS.MTR.PS +
      geom_ribbon(aes(ymin = .data$ArbError_MIV.MTS.MTR.PS_lb,
                      ymax = .data$ArbError_MIV.MTS.MTR.PS_ub, color = 'b',linetype = 'b'),
                  alpha = 0)+
      scale_linetype_manual("",labels = c('Joint MTS & MTR',
                                          'Joint MTS, MTR, & MIV'),
                            values = c("a" = "dotted", "b" = "longdash"))+
      scale_color_manual(name = "", labels = c('Joint MTS & MTR',
                                               'Joint MTS, MTR, & MIV'),

                         values = c("a" = "blue4", "b" = "red"))
  } else {
    AE_MTS.MTR.PS <- AE_MTS.MTR.PS +
      scale_linetype_manual("",labels = c('Joint MTS & MTR'),
                            values = c("a" = "solid"))+
      scale_color_manual(name = "", labels = c('Joint MTS & MTR'),

                         values = c("a" = "blue4"))
  }

  AE_MTS.MTR.PS <- AE_MTS.MTR.PS +
    scale_x_continuous(breaks=p.data$ER)+
    theme_minimal()+
    ggtitle(label = 'MTSp-MTR & MTSp-MTR-IV Assumptions', subtitle = 'Arbitrary Errors')+

    xlab('Maximum Allowed Degree of Misclassification')+
    ylab('ATE')+
    theme(legend.position = 'bottom',
          panel.grid.minor = element_blank())+
    guides(linetype=guide_legend(nrow=1,byrow=TRUE))




  NFP_MTS.MTR.PS <- ggplot(p.data, aes(x = .data$ER)) +
    geom_ribbon(aes(ymin = .data$NFP_MTS.MTR.PS_lb,
                    ymax = .data$NFP_MTS.MTR.PS_ub, color = 'a',linetype = 'a'
    ),
    alpha = 0)
  if("ArbError_MIV.MTS.MTR.PS_lb" %in% colnames(p.data)){
    NFP_MTS.MTR.PS <- NFP_MTS.MTR.PS +
      geom_ribbon(aes(ymin = .data$NFP_MIV.MTS.MTR.PS_lb,
                      ymax = .data$NFP_MIV.MTS.MTR.PS_ub, color = 'b',linetype = 'b'),
                  alpha = 0)+
      scale_linetype_manual("",labels = c('Joint MTS & MTR',
                                          'Joint MTS, MTR, & MIV'),
                            values = c("a" = "dotted", "b" = "longdash"))+
      scale_color_manual(name = "", labels = c('Joint MTS & MTR',
                                               'Joint MTS, MTR, & MIV'),

                         values = c("a" = "blue4", "b" = "red"))
  } else {
    NFP_MTS.MTR.PS <- NFP_MTS.MTR.PS +
      scale_linetype_manual("",labels = c('Joint MTS & MTR'),
                            values = c("a" = "solid"))+
      scale_color_manual(name = "", labels = c('Joint MTS & MTR'),

                         values = c("a" = "blue4"))
  }

  NFP_MTS.MTR.PS <- NFP_MTS.MTR.PS +
    scale_x_continuous(breaks=p.data$ER)+
    theme_minimal()+
    ggtitle(label = 'MTSp, MTR, & MTSp-IV Assumptions', subtitle = 'No False Positives')+

    xlab('Maximum Allowed Degree of Misclassification')+
    ylab('ATE')+
    theme(legend.position = 'bottom',
          panel.grid.minor = element_blank())+
    guides(linetype=guide_legend(nrow=1,byrow=TRUE))





  out <- list("Exogenous and No Selection Assumptions" = ESM_NMA,
              "Monotone sel.- Arb. Error - NS" = AE_MTS.NS,
              "Monotone sel.- NFP - NS" = NFP_MTS.NS,
              "Monotone sel.- Arb. Error - PS" = AE_MTS.PS,
              "Monotone sel.- NFP - PS" = NFP_MTS.PS,
              "Monotone sel. & treat.- Arb. Error - NS" = AE_MTS.MTR.NS,
              "Monotone sel. & treat.- NFP - NS" = NFP_MTS.MTR.NS,
              "Monotone sel. & treat.- Arb. Error - PS" = AE_MTS.MTR.PS,
              "Monotone sel. & treat.- NFP - PS" = NFP_MTS.MTR.PS)

  return(out)
}
#
#

#
# 'ArbError_MTS.NS_ub','ArbError_MTS.NSM_lb',
# 'NFP_MTS.NS_ub',     'NFP_MTS.NS_lb',
# 'ArbError_MTS.PS_ub','ArbError_MTS.PSM_lb',
# 'NFP_MTS.PS_ub',     'NFP_MTS.PS_lb',
#
# 'ArbError_MTS.MTR.NS_ub','ArbError_MTS.MTR.NSM_lb',
# 'NFP_MTS.MTR.NS_ub',     'NFP_MTS.MTR.NS_lb',
# 'ArbError_MTS.MTR.PS_ub','ArbError_MTS.MTR.PSM_lb',
# 'NFP_MTS.MTR.PS_ub',     'NFP_MTS.MTR.PS_lb'
#
# 'ArbError_MIV.MTS.NS_ub','ArbError_MIV.MTS.NSM_lb',
# 'NFP_MIV.MTS.NS_ub',     'NFP_MIV.MTS.NS_lb',
# 'ArbError_MIV.MTS.PS_ub','ArbError_MIV.MTS.PSM_lb',
# 'NFP_MIV.MTS.PS_ub',     'NFP_MIV.MTS.PS_lb',
#
# 'ArbError_MIV.MTS.MTR.NS_ub','ArbError_MIV.MTS.MTR.NSM_lb',
# 'NFP_MIV.MTS.MTR.NS_ub',     'NFP_MIV.MTS.MTR.NS_lb',
# 'ArbError_MIV.MTS.MTR.PS_ub','ArbError_MIV.MTS.MTR.PSM_lb',
# 'NFP_MIV.MTS.MTR.PS_ub',     'NFP_MIV.MTS.MTR.PS_lb'

##########################################################################################
#
#  TTITLE:  Probability of observing the hypothesized treatment effect at the end
#           of study conditional on the interim look results
#  AUTHOR:  Shamil Sadikhov
#
##########################################################################################


require(ggplot2)
require(gridExtra)
require(mvtnorm)
require(mnormt)


cprob <- function(delta                 # true effect of drug on ADAS-Cog (mu.trt-mu.pl)
                  , d.IA                # observed mean adas effect at IA
                  , d.FA                # observed mean adas effect at FA
                  , n.trt.IA            # number of subjects in trt arm at IA
                  , n.trt.FA            # number of subjects in trt arm at FA
                  , n.pl.IA             # number of subjects in placebo arm at IA
                  , n.pl.FA             # number of subjects in trt arm at FA
                  , stddev.trt          # sd for change from baseline in ADAS-Cog in active treatment arm
                  , stddev.pl           # sd for change from baseline in ADAS-Cog for indiv in placebo arm
                  ) {          
  
  t1 <- (d.FA-delta)
  t2 <- (d.IA-delta)
  t3 <- (stddev.trt^2/n.trt.FA  + stddev.pl^2/n.pl.FA)
  t4 <- (stddev.trt^2/n.trt.IA  + stddev.pl^2/n.pl.IA)

    # probability that the FA effect conditional on the IA effect will exceed the FA effect
  tstat <- ifelse(t1 >= t2, abs(t1 - (t3/t4)*t2)/sqrt(t3-t3^2/t4),
                           t1 - (t3/t4)*t2/sqrt(t3-t3^2/t4))
  pnorm(tstat, lower.tail = FALSE)
}


##########################################################################################
#
#                              ADAS-COG calculations
#
##########################################################################################

cprob(delta = -2.23
      , d.IA = -2.23
      , d.FA = -2.23
      , n.trt.IA = 40
      , n.pl.IA = 40
      , n.trt.FA = 118
      , n.pl.FA = 118
      , stddev.trt = 8.03
      , stddev.pl = 8.03
      )


cprob(delta = -2.23
      , d.IA = -1.36
      , d.FA = -2.23
      , n.trt.IA = 40
      , n.pl.IA = 40
      , n.trt.FA = 118
      , n.pl.FA = 118
      , stddev.trt = 8.03
      , stddev.pl = 8.03
      )

cprob(delta = -2.23
      , d.IA = 0
      , d.FA = -2.23
      , n.trt.IA = 40
      , n.pl.IA = 40
      , n.trt.FA = 118
      , n.pl.FA = 118
      , stddev.trt = 8.03
      , stddev.pl = 8.03
      )

cprob(delta = 2.23
      , d.IA = 1.36
      , d.FA = 2.23
      , n.trt.IA = 40
      , n.pl.IA = 40
      , n.trt.FA = 118
      , n.pl.FA = 118
      , stddev.trt = 8.03
      , stddev.pl = 8.03
      )

cprob(delta = 2.4
      , d.IA = 2
      , d.FA = 2.4
      , n.trt.IA = 118
      , n.pl.IA = 118
      , n.trt.FA = 400
      , n.pl.FA = 400
      , stddev.trt = 8
      , stddev.pl = 8
      )

# plot the conditional probability as a function of the observed effect at IA given the information fraction, and hypothesized effect and sd at the FA.

cprob.cog.ia.eff <- function(x) { cprob(delta = -2.23
                                  , d.IA = x
                                  , d.FA = -2.23
                                  , n.trt.IA = 40
                                  , n.pl.IA = 40
                                  , n.trt.FA = 118
                                  , n.pl.FA = 118
                                  , stddev.trt = 8.03
                                  , stddev.pl = 8.03
                                  )}

curve(cprob.cog.ia.eff, from=-3, to=1, xlab = 'Observed treatment effect on ADAS-COG at interim look', ylab = 'Conditional probability')

# or

cprob.cog.plot1 <- ggplot(data.frame(x=c(-2.23, 1)), aes(x)) + stat_function(fun = cprob.cog.ia.eff) + 
                      scale_x_continuous('Observed treatment effect on ADAS-COG at interim look') +
                      scale_y_continuous('Conditional probability')  


cprob.cog.plot1 <- arrangeGrob(cprob.cog.plot1, sub = textGrob("* Conditional proability of observing at least the hypothesized treatment effect given the interim look data", x = 0.1, hjust = 0.1, vjust= 0.4, gp = gpar(fontface = "italic", fontsize = 8)))

cprob.cog.plot1

# qplot(c(-2.23, 1), fun= cprob.ia.eff, stat="function", geom="line")

##########################################################################################
#
#                               ADCS-ADL Calculations
#
##########################################################################################
cprob.adl.ia.eff <- function(x) { cprob(delta = 2.26
                                  , d.IA = x
                                  , d.FA = 2.26
                                  , n.trt.IA = 40
                                  , n.pl.IA = 40
                                  , n.trt.FA = 118
                                  , n.pl.FA = 118
                                  , stddev.trt = 12
                                  , stddev.pl = 12
                                  )}

cprob.adl.ia.eff <- function(x) { cprob(delta = 2.26
                                  , d.IA = x
                                  , d.FA = 2.26
                                  , n.trt.IA = 40
                                  , n.pl.IA = 40
                                  , n.trt.FA = 118
                                  , n.pl.FA = 118
                                  , stddev.trt = 12
                                  , stddev.pl = 12
                                  )}

cprob.adl.plot1 <- ggplot(data.frame(x=c(0, 6)), aes(x)) + stat_function(fun = cprob.adl.ia.eff) + 
                      scale_x_continuous('Observed treatment effect on ADAS-ADL at interim look') +
                      scale_y_continuous('Conditional probability')  


cprob.adl.plot1 <- arrangeGrob(cprob.adl.plot1, sub = textGrob("* Conditional proability of observing at least the hypothesized treatment effect given the interim look data", x = 0.1, hjust = 0.1, vjust= 0.4, gp = gpar(fontface = "italic", fontsize = 8)))




##########################################################################################
#
#                               Joint probability using Bivariate Normal
#
##########################################################################################

jcprob <- function(delta                 # vector of true effects of drug (mu.trt-mu.pl)
                   , d.IA                # vector of observed mean diffs - effect at IA
                   , d.FA                # vector of observed mean diffs - effect at FA 
                   , n.trt.IA            # number of subjects in trt arm at IA
                   , n.trt.FA            # number of subjects in trt arm at FA
                   , n.pl.IA             # number of subjects in placebo arm at IA
                   , n.pl.FA             # number of subjects in trt arm at FA
                   , stddev1.trt         # sd for change from baseline in ADAS-Cog in active treatment arm
                   , stddev1.pl          # sd for change from baseline in ADAS-Cog indiv in placebo arm
                   , stddev2.trt         # sd for change from baseline in ADCS-ADL in active treatment arm
                   , stddev2.pl          # sd for change from baseline in ADCS-ADL indiv in placebo arm                   
                   , rho                 # correlation between endpoints
                   ) {
    
    var1.FA <- (stddev1.trt^2/n.trt.FA  + stddev1.pl^2/n.pl.FA)
    var1.IA <- (stddev1.trt^2/n.trt.IA  + stddev1.pl^2/n.pl.IA)
    var2.FA <- (stddev2.trt^2/n.trt.FA  + stddev2.pl^2/n.pl.FA)
    var2.IA <- (stddev2.trt^2/n.trt.IA  + stddev2.pl^2/n.pl.IA)

    dbar.IA <- d.IA - delta 
    dbar.FA <- d.FA - delta

    # FA correlation matrix
    S.11 <- matrix(c(var1.FA, rho,
                     rho, var2.FA),
                          nrow = 2, byrow=TRUE)
    # IA correlation matrix
    S.22 <- matrix(c(var1.IA, rho,
                     rho, var2.IA),
                          nrow = 2, byrow=TRUE)
    # Covariances
    S.12 <-  matrix(c(var1.FA, rho*sqrt(var1.FA*var2.FA*var1.FA/var2.IA),
                      rho*sqrt(var1.FA*var2.FA*var1.FA/var2.IA), var2.FA),
                           nrow = 2, byrow=TRUE)
    S.21 <- t(S.12)

                                        # conditional bivariate distribution of FA|IA
    E.mean  <- as.vector(dbar.FA + S.11%*%solve(S.22)%*%(dbar.IA))
    E.sigma <- S.11 - S.12 %*% solve(S.22) %*% S.21 


                                        # multivariate normal probability
    res <- pmvnorm(lower = rep(0, 2), upper = rep(Inf, 2), mean = E.mean, sigma = E.sigma)                       
                                        # pmvt(lower = rep(-Inf, 2), upper = rep(Inf, 2), delta = E.mean, sigma = E.sigma, df=0)
                                        # res <- mnormt::pmnorm(x = c(0,0), mean = E.mean, varcov = E.sigma)

    return(res)    

}

jcprob(delta = c(2.26, 2.26)
       , d.IA = c(2.26, 2.26)
       , d.FA = c(2.26, 2.26)
       , n.trt.IA = 50
       , n.trt.FA = 200
       , n.pl.IA = 50
       , n.pl.FA = 200
       , stddev1.trt = 8         
       , stddev1.pl = 8         
       , stddev2.trt = 12       
       , stddev2.pl = 12                          
       , rho = 0.0  )

jcprob(delta = c(2.26, 2.26)
       , d.IA = c(2.26, 2.26)
       , d.FA = c(2.26, 2.26)
       , n.trt.IA = 50
       , n.trt.FA = 200
       , n.pl.IA = 50
       , n.pl.FA = 200
       , stddev1.trt = 8         
       , stddev1.pl = 8         
       , stddev2.trt = 12       
       , stddev2.pl = 12                          
       , rho = 0.5  )


jcprob(delta = c(2.26, 2.26)
       , d.IA = c(1.26, 1.26)
       , d.FA = c(1.26, 1.26)
       , n.trt.IA = 50
       , n.trt.FA = 200
       , n.pl.IA = 50
       , n.pl.FA = 200
       , stddev1.trt = 8         
       , stddev1.pl = 8         
       , stddev2.trt = 12       
       , stddev2.pl = 12                          
       , rho = 0.0  )

jcprob(delta = c(2.26, 2.26)
       , d.IA = c(1.26, 1.26)
       , d.FA = c(2.26, 2.26)
       , n.trt.IA = 50
       , n.trt.FA = 200
       , n.pl.IA = 50
       , n.pl.FA = 200
       , stddev1.trt = 8         
       , stddev1.pl = 8         
       , stddev2.trt = 12       
       , stddev2.pl = 12                          
       , rho = 0.0  )

jcprob(delta = c(2.26, 2.26)
       , d.IA = c(1.26, 1.26)
       , d.FA = c(2.26, 2.26)
       , n.trt.IA = 50
       , n.trt.FA = 200
       , n.pl.IA = 50
       , n.pl.FA = 200
       , stddev1.trt = 8         
       , stddev1.pl = 8         
       , stddev2.trt = 12       
       , stddev2.pl = 12                          
       , rho = 0.3  )


### Plot ellipse for bivariate normal
require(MASS)
library(car)

   k<-1
   Sigma <- E.sigma
   rho <- Sigma[1,2]/sqrt(Sigma[1,1]*Sigma[2,2])
   rho
   eta1<-replicate(300,mvrnorm(k, mu= E.mean, Sigma))


   dataEllipse(eta1[1,],eta1[2,], levels=c(0.05, 0.95))

##########################################################################################
#
#  TTITLE:  Probability of observing the hypothesized treatment effect at the end
#           of study conditional on the interim look results
#  AUTHOR:  Shamil Sadikhov
#
##########################################################################################


require(ggplot2)
require(gridExtra)


cprob <- function(delta                 # true effect of drug on ADAS-Cog (mu.trt-mu.pl)
                  , d.IA           # observed mean adas effect at IA
                  , d.FA           # observed mean adas effect at FA
                  , n.trt.IA             # number of subjects in trt arm at IA
                  , n.trt.FA             # number of subjects in trt arm at FA
                  , n.pl.IA              # number of subjects in placebo arm at IA
                  , n.pl.FA              # number of subjects in trt arm at FA
                  , stddev.trt          # number of subjects in placebo arm at IA
                  , stddev.pl           # sd for change from baseline in ADAS-Cog for indiv in placebo arm
                  ) {          
  
  t1 <- (d.FA-delta)
  t2 <- (d.IA-delta)
  t3 <- (stddev.trt^2/n.trt.FA  + stddev.pl^2/n.pl.FA)
  t4 <- (stddev.trt^2/n.trt.IA  + stddev.pl^2/n.pl.IA)

    # probability that the FA effect conditional on the IA effect will exceed the FA effect
  pnorm(abs(t1 - (t3/t4)*t2)/sqrt(t3-t3^2/t4), lower.tail = FALSE)
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

curve(cprob.cog.ia.eff, from=-2.23, to=1, xlab = 'Observed treatment effect on ADAS-COG at interim look', ylab = 'Conditional probability')

# or

cprob.cog.plot1 <- ggplot(data.frame(x=c(-2.23, 1)), aes(x)) + stat_function(fun = cprob.cog.ia.eff) + 
                      scale_x_continuous('Observed treatment effect on ADAS-COG at interim look') +
                      scale_y_continuous('Conditional probability')  


cprob.cog.plot1 <- arrangeGrob(cprob.cog.plot1, sub = textGrob("* Conditional proability of observing at least the hypothesized treatment effect given the interim look data", x = 0.1, hjust = 0.1, vjust= 0.4, gp = gpar(fontface = "italic", fontsize = 8)))

cprob.cog.plot1

# qplot(c(-2.23, 1), fun= cprob.ia.eff, stat="function", geom="line")

##########################################################################################
#
#                               ADAS-ADL Calculations
#
##########################################################################################
cprob.adl.ia.eff <- function(x) { cprob(delta = -3.6
                                  , d.IA = x
                                  , d.FA = -3.6
                                  , n.trt.IA = 40
                                  , n.pl.IA = 40
                                  , n.trt.FA = 118
                                  , n.pl.FA = 118
                                  , stddev.trt = 12
                                  , stddev.pl = 12
                                  )}

cprob.adl.plot1 <- ggplot(data.frame(x=c(-3.6, 1)), aes(x)) + stat_function(fun = cprob.adl.ia.eff) + 
                      scale_x_continuous('Observed treatment effect on ADAS-ADL at interim look') +
                      scale_y_continuous('Conditional probability')  


cprob.adl.plot1 <- arrangeGrob(cprob.adl.plot1, sub = textGrob("* Conditional proability of observing at least the hypothesized treatment effect given the interim look data", x = 0.1, hjust = 0.1, vjust= 0.4, gp = gpar(fontface = "italic", fontsize = 8)))

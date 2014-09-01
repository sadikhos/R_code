##########################################################################################
#
#  Conditional power and predictive power calculations for MAOB Mayflower Road interim analysis
#  Author: Shamil Sadikhov
#  Date:
#  Reference: Proschan, Lan, Wittes. Statsitical Monitoring of Clinical Trials
# 
##########################################################################################

# rm(list=objects())

# setwd()
source('U:\\My Documents\\Alzheimers\\MAO_B\\BP28248\\IA\\cpower1.R')


######################################################################
#
#                       Predictive Power
#
######################################################################
# Note: Predictive power averages conditional power over the posterior distribution of the drift parameter given the data (at interim analysis).

## Note: to give 10% weight to the prior variance information relative to the data, solve (1/sigma_0^2) / (1+ 1/sigma.0^2) = 0.10


## # Assuming the prior on true drift parameter theta is under alternative (z.alpha.half + z.beta)

# prior var of 20 gives approximately 5% weight to the prior on variance of theta
ppowerf(type1 = 0.1
        ,type2 = 0.2
        ,delta.t=2.23
        ,sigma.t=8.03
        ,N=118
        ,n.a=40
        ,n.b=40
        ,theta.0= qnorm(0.1, lower.tail = FALSE) + qnorm(1-0.2, lower.tail=TRUE)
        ,sigma.0= 20
        ,two.sided=FALSE)

# At the beginning of the trial
ppowerf(type1 = 0.1
        ,type2 = 0.2
        ,delta.t=2.23
        ,sigma.t=8.03
        ,N=118
        ,n.a=0.0000000000000001
        ,n.b=0.0000000000000001
        ,theta.0= qnorm(0.1, lower.tail = FALSE) + qnorm(1-0.2, lower.tail=TRUE)
        ,sigma.0=20
        ,two.sided=FALSE)

# predictive power = conditional power when prior variance on theta = 0
ppowerf(type1 = 0.1
        ,type2 = 0.2
        ,delta.t=2.23
        ,sigma.t=8.03
        ,N=118
        ,n.a=0.0000000000000001
        ,n.b=0.0000000000000001
        ,theta.0= qnorm(0.1, lower.tail = FALSE) + qnorm(1-0.2, lower.tail=TRUE)
        ,sigma.0=0
        ,two.sided=FALSE)

ppowerf(type1 = 0.1
        ,type2 = 0.2
        ,delta.t=2.23
        ,sigma.t=8.03
        ,N=118
        ,n.a=40
        ,n.b=40
        ,theta.0= qnorm(0.1, lower.tail = FALSE) + qnorm(1-0.2, lower.tail=TRUE)
        ,sigma.0=0
        ,two.sided=FALSE)

ppowerf(type1 = 0.1
        ,type2 = 0.2
        ,delta.t=2.23
        ,sigma.t=8.03
        ,N=118
        ,n.a=40
        ,n.b=40
        ,theta.0= qnorm(0.1, lower.tail = FALSE) + qnorm(1-0.2, lower.tail=TRUE)
        ,sigma.0=0
        ,two.sided=FALSE) 

######################################################################
#
#                       Conditional Power
#
######################################################################

# Assumption: under null hypothesis
cpower(type1 = 0.1
       ,type2 =  0.2     # type 2 error rate beta
       ,delta = 2.23
       ,delta.t = 2.23    # observed z score at interim
       ,sigma.t = 8.03    # observed variance of the z-score at interim
       ,N = 118           # per arm sample size
       ,n.a = 40          # n at interim in arm 'a'
       ,n.b = 40          # n at interim in arm 'b'
       ,two.sided = FALSE)

cpower(type1 = 0.1
       ,type2 =  0.2     # type 2 error rate beta
       ,delta = 2.23
       ,delta.t = 2.23    # observed z score at interim
       ,sigma.t = 8.03    # observed variance of the z-score at interim
       ,N = 118           # per arm sample size
       ,n.a = 1          # n at interim in arm 'a'
       ,n.b = 1          # n at interim in arm 'b'
       ,two.sided = FALSE)

# Assmuption: Original treatment size, but observed variability at interim
cpower(type1 = 0.1
       ,type2 =  0.2      # type 2 error rate beta
       ,delta = 2.23
       ,delta.t = 2.23    # observed z score at interim
       ,sigma.t = 7.24    # observed variance of the z-score at interim
       ,N = 118           # per arm sample size
       ,n.a = 40          # n at interim in arm 'a'
       ,n.b = 40          # n at interim in arm 'b'
       ,two.sided = FALSE)

cpower(type1 = 0.2
       ,type2 =  0.2      # type 2 error rate beta
       ,delta = 2.23
       ,delta.t = 1.34    # observed z score at interim
       ,sigma.t = 8.03    # observed variance of the z-score at interim
       ,N = 118           # per arm sample size
       ,n.a = 40          # n at interim in arm 'a'
       ,n.b = 40          # n at interim in arm 'b'
       )


cpower(type1 = 0.1
       ,type2 =  0.1      # type 2 error rate beta
       ,delta = 2.23
       ,delta.t = 2.23    # observed z score at interim
       ,sigma.t = 8.03    # observed variance of the z-score at interim
       ,N = 118           # per arm sample size
       ,n.a = 118          # n at interim in arm 'a'
       ,n.b = 118          # n at interim in arm 'b'
       )

cpower(type1 = 0.1
       ,type2 =  0.1      # type 2 error rate beta
       ,delta = 2.23
       ,delta.t = 2.23    # observed z score at interim
       ,sigma.t = 8.03    # observed variance of the z-score at interim
       ,N = 118           # per arm sample size
       ,n.a = 1          # n at interim in arm 'a'
       ,n.b = 1          # n at interim in arm 'b'
       )




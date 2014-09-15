##########################################################################################
#
#      Validation of IA outputs using ADaM data sets from Sept 2014
#    Author: Shamil Sadikhov
#    Date:   2 Sept 2014
#
##########################################################################################

require(Hmisc)
require(nlme)
require(lsmeans)
require(pbkrtest)
require(multcomp)
require(multcompView)

options( width = 90, digits = 5 )
Sys.setlocale("LC_TIME", "English") # with incorrect locale as.Date does not recognize months, months appear in German
setwd("u:/My Documents/Alzheimers/MAO_B/BP28248/IA/IA Code and Data/Validation")

################################################################################
#                          Import data
################################################################################

maobdata <- sasxport.get("H:\\cdt4983b\\i28248b\\libraries\\csvdata_12SEP2014", method='csv')
# names(maobdata)



##########################################################################################
#
#       Validation of the descriptive statistics and models
#
##########################################################################################

################################################################################
#                          Data Management
################################################################################
### ADAS-COG
adqscog <- upData(maobdata$adqscog)
# names(adqscog)
# View(adqscog)

adascogtot <- droplevels(subset(unique(adqscog), paramcd=='ADCSCORE' & avisit %in% c('Baseline Day -1 Predose', 'Week 12', 'Week 24', 'Week 36', 'Week 52 EOT')))
                                        #, select = c(usubjid, visit, visitnum, agegr, agegr2, sex, race, randdt, arm, hippcat, backgcat, region, edlevel, astdy, paramcd, aval, chg))
# View(adascogtot)

# order the data by subject and time

adascogtot <- adascogtot[order(adascogtot$usubjid, adascogtot$astdy ), ]

# Visit date vs visit
with(adascogtot, summary.formula(astdy ~ avisit, overall = FALSE))


# get subjects with multiple records at each visit
sum(with(adascogtot, tapply(visitnum, list(usubjid), function(x) length(x) > 5)))
sum(with(adascogtot, tapply(avisit, list(usubjid), function(x) length(x) > 5)))
sum(with(adascogtot, tapply(visitnum, list(usubjid), function(x) length(x) > 4)))
sum(with(adascogtot, tapply(avisit, list(usubjid), function(x) length(x) > 4)))
which(with(adascogtot, tapply(visitnum, list(usubjid), function(x) length(x) > 5)))

# subject with more than one record at each visit
which(with(subset(adascogtot, anl01fl=='Y'), tapply(visitnum, list(usubjid, avisit), function(x) length(x) > 1)))


View(adascogtot[adascogtot$usubjid=='BP28248-250562-250562006', ])
View(adascogtot[adascogtot$usubjid=='BP28248-250311-250311009', ])  # major outlier

# Drop multiple records at each visit, if they exist (unclean data)
adascogtot2 <- do.call("rbind", by(data = adascogtot, INDICES = list(adascogtot$usubjid, adascogtot$avisit), FUN = function(x) x[1, ]))
                       
# Grouped Data Object
adascogtot.grp <- groupedData(aval ~ astdy | usubjid, data = adascogtot2)
# gsummary(adascogtot.grp, invariantsOnly = TRUE)

### BEHAVE-AD
adqsbhv <- upData(maobdata$adqsbhv)

# total score calculation
adqsbhv2 <- droplevels(subset(adqsbhv, avisit %in% c('Baseline Day -1 Predose', 'Week 12', 'Week 24', 'Week 52 EOT')))
adqsbhv.ts <- do.call("rbind", by(adqsbhv2, list(adqsbhv2$usubjid, adqsbhv2$avisit),
                                  function(x){
                                      res <- 0
                                      for (i in 1:25){
                                         val <- x[x$paramn==2*i-1, ]$aval 
                                         freq <- x[x$paramn==2*i, ]$aval
                                         if (length(val) > 0) {
                                             if (length(freq) > 0)
                                                 {res <- res + val*freq}
                                             else {res <- res + val}
                                         }
                                     }
                                      res + x[x$paramn==51, ]$aval
                                  }))





##########################################################################################
#
#                               Summary statistics
# 
##########################################################################################
# get means for ADAS-Cog 11  by visit and arm
with(adascogtot, tapply(aval, list(arm, visit), function(x) mean(x, na.rm=TRUE)))
with(adascogtot, tapply(aval, list(arm, visit), function(x) sum(!is.na(x))))
with(adascogtot, tapply(aval, list(arm, visit), function(x) sd(x, na.rm = TRUE)))
with(adascogtot, tapply(aval, list(visit), function(x) sd(x, na.rm = TRUE)))

# Statistics for change from baseline in the ITT population:
sum.adas.cog <- Hmisc::summary.formula(chg ~ avisit + stratify(arm), data = subset(adascogtot, ittfl=='Y'))
print(sum.adas.cog)

# plot(adascogtot.grp[1:100,], outer ~ arm, aspect = 3)
# plot(adascogtot.grp)

xyplot( aval ~ astdy | usubjid, data = subset(adascogtot, ittfl=='Y'), type = 'b', groups = arm )

##########################################################################################
#
#                                    Models
# 
##########################################################################################

############## Marginal Model with unstructured covariance for each subject##############

# create an interger time variable for the corrSymm() general covariance structure
adascogtot3 <- droplevels(subset(adascogtot2, visit != 'Baseline Day -1 Predose' &  anl01fl == 'Y' & ittfl=='Y'))
adascogtot3$timeint <- as.numeric(as.factor(adascogtot3$avisit))

# set reference level for treatmetn arm to placebo
adascogtot3$arm <- relevel(adascogtot3$arm, ref = 'PLACEBO')

# order the data by subject and time/visit
adascogtot3 <- adascogtot3[order(adascogtot3$usubjid, adascogtot3$visit ), ]
mod1 <- gls(chg ~  avisit* arm, correlation = corSymm(form = ~ timeint | usubjid),
            weights = varIdent(form = ~ 1 | timeint),
            data = subset(adascogtot3, anl01fl == 'Y' & ittfl=='Y'), na.action = na.omit)

# corSymm(form = ~ 1 | SUBJECT), gls() will assume that the first measurement for a particular subject corresponds to time = 1. If that's not the case (and that can of course happen if rows with missing data have been deleted), then things get misalligned. Then you need the time covariate, so that things can be matched up properly.
# This should work fine with the data properly sorted
# mod1 <- gls(chg ~  avisit * arm, correlation = corSymm(form = ~ 1 | usubjid), weights = varIdent(form = ~ 1 | timeint), data = adascogtot3, na.action = na.omit)

summary(mod1)

# approx intervals
nlme:::intervals.gls(mod1)

# Diagnostic plots
nlme:::plot.gls(mod1)
nlme:::qqnorm.gls(mod1)
#nlme:::residuals.gls(mod1)
 
# Variance Components

# get correlation matrices for all subjects
corMatrix(obj$modelStruct$corStruct)

var.weights <- varWeights(mod1$modelStruct$varStruct)

# LS Means
lsmeans(mod1, specs = pairwise ~ avisit:arm|avisit)
lsmeans(mod1, specs = trt.vs.ctrl ~ avisit:arm|avisit)


############## Random intercept and categorical visit - AR(1) correlations within subjects model #################
# assuming equally spaced visits -- almost 
mod2a.ar1 <- lme(chg ~ avisit*arm, random = ~ 1|usubjid, correlation = corAR1(),  method = "REML", data = adascogtot3, na.action = na.omit)
summary(mod2a.ar1)


############## Random coefficients model with continuous time (visit date) #################
                                        # random intercept
mod2a <- lme(aval ~ astdy + arm + astdy:arm, random = ~ 1|usubjid, method = "REML", data=adascogtot.grp)
mod2a.ml <- update(mod2a, method = "ML")
                                        # random intercept and slope
mod2b <- lme(aval ~ astdy + arm + astdy:arm, random = ~ (1 + astdy)|usubjid, method = "REML", data=adascogtot.grp)
mod2b.ml <- update(mod2b, method = "ML")
                                        # random intercept and slope
mod2b.quad <- lme(aval ~ astdy + I(astdy^2) + arm + astdy:arm, random = ~ (1 + astdy)|usubjid, method = "REML", data=adascogtot.grp)
mod2b.quad.ml <- update(mod2b.quad, method = "ML")
                                        # random intercept plus continuous AR(1) correlation within subject
mod2a.car1 <- lme(aval ~ astdy + arm + astdy:arm, random = ~ 1|usubjid, correlation = corCAR1(form = ~ astdy), method = "REML", data=adascogtot.grp)
summary(mod2a.car1)
mod2a.car1.ml <- update(mod2a.car1, method = "ML")

#Comparison of models
anova(mod2a.ml, mod2b.ml, mod2b.quad.ml, mod2a.car1.ml)

# Diagnostic plots
plot(augPred(mod2b, level = 0:1), layout = c(4,3,1))
plot(augPred(mod2b.quad, level = 0:1), layout = c(4,3,1))

# conditional residuals vs. conditional predicted values
plot(mod2b, resid(., type = 'p') ~ fitted(.),
     layout = c(1,1), aspect=2, abline = 0, id = 0.01)

fixef(mod2b.quad)
# 1mg vs. placebo trt effect:
#-0.34818 + (-0.000888889)*52*7

##########################################################################################
#
#              Deriving ADAS-Cog total score from individual items
#
##########################################################################################

# View(subset(adqscog, visit=='Baseline Day -1 Predose'))

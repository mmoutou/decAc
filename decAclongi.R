### ... Dropbox/FIL_aux/R_scripts/decAclongi.R
### Various functions and script bits for correlating decision acuity 
### estimates, as derived outside here, with psychiatric measures and also 
### londitudinal analyses of same

# orient and load basics ---------------------------------------------------
try(source('C:/Users/mmpsy/Dropbox/FIL_aux/R_scripts/dirs4r.R'))
try(source('~/Dropbox/FIL_aux/R_scripts/dirs4r.R'))
try(source('/media/sf_mmpsy/Dropbox/FIL_aux/R_scripts/dirs4r.R'))

dirs <- dirs4r()
setwd(dirs$uchjan17on)
### was: symfacdec <- importCSV(paste(dirs$uchjan17on,"symDispFactDecAcIQ.csv",sep=''))
symfacdec <- importCSV(paste(dirs$uchjan17on,"symfacdec.csv",sep=''))
### -------------------------------------------------------------------------

### Simple lme - from https://rpsychologist.com/r-guide-longitudinal-lme-lmer  which is the most  
###                   amazingly awesome list at increasing levels of hierarchical complexity EVER
###              and https://stats.stackexchange.com/questions/40647/lme-error-iteration-limit-reached
###                  hence the ctrl <- lmeControl(opt='optim') option

library(nlme)    # needed for 'lmeControl' etc.
library(dplyr)   # needed for 'arrange' and %>%
### was as follows, changed so that rather than using extrapolated sl5 etc. scores,
### use interpolated decAc4b scores further down below
###
# v2do <- c('nspnID', 'sex',                      # REM here M=1, F=2
#          c("decAc4b.bsl","decAc4b.fu2"),
#          c("age_iua1","age_iua2"),
#          c("sl5_general_1","sl5genext2"),
#          c("sl5_sf1_1","sl5sf1ext2"),                # sl5_sf1 is something like 'self-confidence'
#          c("sl5_sf2_1","sl5sf2ext2"),                # sl5_sf2 is something like 'antisocial'
#          c("sl5_sf3_1","sl5sf3ext2"),                # sl5_sf3 is something like 'worry'
#          c("sl5_sf4_1","sl5sf4ext2"),                # sl5_sf4 is something like 'aberrant thinking'
#          c("sl5_sf5_1","sl5sf5ext2"),                # sl5_sf5 is something like 'mood'
#          c("IQmatrix.raw.bsl","IQmatrix.raw.fu2"),
#          c("IQvocab.raw.bsl","IQvocab.raw.fu2"),
#          c("dispGeneral_Factor_1", "dispSpecific_factor1_1","dispSpecific_factor2_1","dispSpecific_factor3_1","dispSpecific_factor4_1")
# )
 
v2do <- c('nspnID', 'sex',                      # REM here M=1, F=2
          c("decAc4b.bsl","decAchqp2est"),
          c("age_iua1","age_hqp2_returned"),
          c("sl5_general_1","sl5_general_2"),
          c("sl5_sf1_1","sl5_sf1_2"),                # sl5_sf1 is something like 'self-confidence'
          c("sl5_sf2_1","sl5_sf2_2"),                # sl5_sf2 is something like 'antisocial'
          c("sl5_sf3_1","sl5_sf3_2"),                # sl5_sf3 is something like 'worry'
          c("sl5_sf4_1","sl5_sf4_2"),                # sl5_sf4 is something like 'aberrant thinking'
          c("sl5_sf5_1","sl5_sf5_2"),                # sl5_sf5 is something like 'mood'
          c("IQmatrix.raw.bsl","IQmatr.raw.hqp2est"),
          c("IQvocab.raw.bsl","IQvoc.raw.hqp2est"),
          c("dispGeneral_Factor_1", "dispSpecific_factor1_1","dispSpecific_factor2_1","dispSpecific_factor3_1","dispSpecific_factor4_1")
)


        
yl <- symfacdec[,v2do];
yl$dt <- yl$age_hqp2_returned - yl$age_iua1
# reshape selected for long format
yl <- reshape(yl,
              direction='long',
              varying= list(c("decAc4b.bsl","decAchqp2est"),
                            c("age_iua1","age_hqp2_returned"),
                            c("sl5_general_1","sl5_general_2"),
                            c("sl5_sf1_1","sl5_sf1_2"),
                            c("sl5_sf2_1","sl5_sf2_2"),
                            c("sl5_sf3_1","sl5_sf3_2"),
                            c("sl5_sf4_1","sl5_sf4_2"),
                            c("sl5_sf5_1","sl5_sf5_2"),
                            c("IQmatrix.raw.bsl","IQmatr.raw.hqp2est"),
                            c("IQvocab.raw.bsl","IQvoc.raw.hqp2est")
              ),
              idvar = "nspnID",
              timevar = 'assessment'
)
### was:
# yl <- reshape(yl,
#               direction='long',
#               varying= list(c("decAc4b.bsl","decAc4b.fu2"),
#                             c("age_iua1","age_iua2"),
#                             c("sl5_general_1","sl5genext2"),
#                             c("sl5_sf1_1","sl5sf1ext2"),
#                             c("sl5_sf2_1","sl5sf2ext2"),
#                             c("sl5_sf3_1","sl5sf3ext2"),
#                             c("sl5_sf4_1","sl5sf4ext2"),
#                             c("sl5_sf5_1","sl5sf5ext2"),
#                             c("IQmatrix.raw.bsl","IQmatrix.raw.fu2"),
#                             c("IQvocab.raw.bsl","IQvocab.raw.fu2")
#               ),
#               idvar = "nspnID",
#               timevar = 'iua'
# )


yl <- {yl %>% arrange(-desc(nspnID))}
# label nicely
colnames(yl) <- c('nspnID','sexM1F2',
                  'dispo_general', "dispo_sf1", "dispo_sf2", "dispo_sf3", "dispo_sf4",
                  'dt','assessment','decAc4b', 'age', 
                  'sl5_general',  'sl5_sf1', 'sl5_sf2', 'sl5_sf3', 'sl5_sf4',  'sl5_sf5',
                  'IQmatrix.raw', 'IQvocab.raw')              

### BIC etc for this big model when symptoms were extrapolated was: 
# Linear mixed-effects model fit by REML
# AIC      BIC    logLik
# 5047.977 5135.685 -2506.989
ctrl <- lmeControl(opt='optim'); l1 <- lme(decAc4b ~ age + sl5_general + dispo_general + 
                                             sl5_sf1 + sl5_sf2 + sl5_sf3 + sl5_sf4 + sl5_sf5 + 
                                             dispo_sf1 + dispo_sf2 + dispo_sf3 + dispo_sf4, 
                                           random =  ~age | nspnID, data=yl, na.action='na.omit', control = ctrl); 
summary(l1)
### BIC etc for this big model when symptoms were extrapolated was much worse: 
###        AIC      BIC    logLik
###     5596.012 5683.917 -2781.006

### Just as good or better if no random slope or specific dispos:
###   AIC      BIC    logLik
###   5025.773 5082.559 -2501.887
ctrl <- lmeControl(opt='optim'); l2 <- lme(decAc4b ~ age + sl5_general + dispo_general + 
                                             sl5_sf1 + sl5_sf2 + sl5_sf3 + sl5_sf4 + sl5_sf5+IQmatrix.raw+IQvocab.raw, 
                                           random =  ~1 | nspnID, data=yl, na.action='na.omit', 
                                           control = ctrl); 
summary(l2)

### pruned down, no dispos:
ctrl <- lmeControl(opt='optim'); l3 <- lme(decAc4b ~ age + sl5_general + 
                                             sl5_sf1 + sl5_sf2 + sl5_sf3 + sl5_sf4 + sl5_sf5, 
                                           random =  ~1 | nspnID, data=yl, na.action='na.omit', 
                                           control = ctrl); 
summary(l3)



# -------------------------------------------------------------------------
## Work based on first draft of Moutoussis-Garzon / Guitart draft
#  See scipsyplan.txt 18/02/2019 - 20/02/2019
y <- symfacdec  # restore in case section above has changed it.
nodep <- !vecTRUE(y[,'cohort'] =='Depression')

# sf2,3,4 are ANTISOC, WORRY, ABERRANT THINKING :
lr = summary(lm(y[nodep,"decAc4b.bsl"] ~ (y[nodep,"sl5_general_1"] + 
                                            y[nodep,"sl5_sf1_1"] + 
                                            y[nodep,"sl5_sf2_1"] + 
                                            y[nodep,"sl5_sf3_1"] + 
                                            y[nodep,"sl5_sf4_1"] + 
                                            y[nodep,"sl5_sf5_1"] + 
                                            y[nodep,"dispGeneral_Factor_1"] + 
                                            y[nodep,"dispSpecific_factor1_1"] + 
                                            y[nodep,"dispSpecific_factor2_1"]+ 
                                            y[nodep,"dispSpecific_factor1_1"]+ 
                                            y[nodep,"dispSpecific_factor4_1"] ))) ; 
lr

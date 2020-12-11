# Factor analyses over UCHANGE parameters, Sept 18 on.

# #########  Provide basic paths depending on machine - hunt out dirs4r and source it ##############
remove(list=ls())
bigreload = 0;   # Whether to reload the big RData file uch_Jan17on.RData with 
                 # AllD18 but also work so far etc.
repartitiondata = 0  # redo the data partion - GENERALLY KEEP AT ZERO UNLESS V. THOROUGHLY BACKED UP!!
                     # ... and ... NO FISHING !
redotraintest = 0;   # If we are to re-prepare the train/test data.
repreprocbsl =0;     # If we want to redo the baseline preprocessing.
repreprocfu1 =0;     # If we want to redo stability preprocessing.
repreprocfu2 =0;     # If we want to redo main follow-up preprocessing.
redefinemissing = 0; # Redo various structures with missing pts. Benjamin communicated.
save2disk = 0;       # Only save if things worked ... 

print('IGNORE \'cannot change working directory\' errors below (switch statement) unless further crash later:',quote='F')
switch(Sys.info()[['sysname']],
       Linux  = {try(setwd('/media/sf_D_DRIVE/mmoutou/Dropbox/FIL_aux/R_scripts/'));   # Kikidi VM
         try(setwd('~/Dropbox/FIL_aux/R_scripts/'));                           # Usual linux w Dbox
         try(setwd('/share/scratch/mmoutou/R/R4hal/scripts/'));                # Hal jury rig
         try(setwd('~/R/R4hal/scripts/'));                                  # if Hal home functions.
         try(setwd("/media/sf_mmpsy/Dropbox/FIL_aux/R_scripts"))},          # G7 ubuntu VM  
       Windows= {try(setwd('D:/mmoutou/Dropbox/FIL_aux/R_scripts/'));       #  Cedric Windoze
         try(setwd('C:/Users/mmpsy/Dropbox/FIL_aux/R_scripts/'))}           #  Dell G7 ubuntu VM
) ;
source('dirs4r.R');
wherearewe <- dirs4r();  # provides paths for main R scripts, a home dir, and optionally for LD work.

#       #########   Load key data       ##############
if (bigreload){
   # load(paste(wherearewe$uchjan17on,'semPavD.RData',sep=''))
   load(paste(wherearewe$uchjan17on,'uch_Jan17on.RData',sep=''))  # No need to be stingy !
}

#  File below created as follows: save(AllD,AllDhd,IDs,file='AllD.RData')
#  So contains basically big data array AllD, but also handy lists
#  with catagorzations of participants and variables into waves, tasks etc. :
# load(paste(wherearewe$uch17on,'AllD18.RData',sep=''))

#                         #########  Preprocessing parameters  ##############

# .  Select manageable number, e.g. 30-40
# .  Initial transforms, e.g. power
# .  Orthogonalize most highly correl. params, incl. fishes task Cost parameters into 
#    their sum and difference
# .  z-scoring

if (repreprocbsl){ # if we are to do the baseline preprocessing, 
    
  # 1. Go-NoGo  baseline preprocessing ###########################################################
  print(hd18$gng) # print(AllDhd$gng)
  # Next var., bslcog4fa, to be augmented as we process each task.
  bslcog4fa <- c("nspnID",
                 "G2W.lt.RT.bsl", "G2AP.lt.RT.bsl",  "invG2W.bsl", "invG2AP.bsl",         
                 "beta.b2a.bsl",  "lrnR_Appet.b2a.bsl", "lrnR_Aver.b2a.bsl",  
                 "PavBias.b2a.bsl", "irNoiseXi.b2a.bsl", "GoBias.b2a.bsl")
  bsltrgng <- c("nspnID",
                "meanLnRT.bsl", "difLnRT.bsl",           
                "ln.beta.bsl",  "lrnRap.to4th.bsl", "lrnRav.to4th.bsl",  "ln.PavBias.bsl", "ln.irNoise.bsl", "actionBias.bsl")
  bsltrcog <- bsltrgng;  # we will redo this incrementally below. 
  
  # Reaction time orthogonalization :
  AllD18$meanLnRT.bsl <- (log(AllD18$G2W.lt.RT.bsl) + log(AllD18$G2AP.lt.RT.bsl))/2
  AllD18$difLnRT.bsl  <- (log(AllD18$G2AP.lt.RT.bsl)- log( AllD18$G2W.lt.RT.bsl))
  
  # transf. of the others:
  AllD18$ln.beta.bsl <- log(AllD18$beta.b2a.bsl)
  AllD18$lrnRap.to4th.bsl <- AllD18$lrnR_Appet.b2a.bsl^(1/4)
  AllD18$lrnRav.to4th.bsl <- AllD18$lrnR_Aver.b2a.bsl^(1/4)
  AllD18$ln.PavBias.bsl <- log(AllD18$PavBias.b2a.bsl)
  AllD18$ln.irNoise.bsl <- log(AllD18$irNoiseXi.b2a.bsl)
  AllD18$actionBias.bsl <- AllD18$GoBias.b2a.bsl
  
  # Z score all  the new ones, but keep a record of their mean and SD first:
  gng4faM  = colMeans(na.omit(AllD18[,c(bsltrgng[2:length(bsltrgng)] )])          )
  gng4faSD = apply(   na.omit(AllD18[,c(bsltrgng[2:length(bsltrgng)] )]),2,sqrtVar)
  AllD18[,c('nspnID', bsltrgng[2:length(bsltrgng)])] <- zScore(AllD18[,c('nspnID',bsltrgng[2:length(bsltrgng)] )],IDcol=1)
  
  # 2. Roulette/MVS  baseline preprocessing #####################################################
  
  bslcog4fa <- c(bslcog4fa, hd18$bgmvs[1:4])  # Augment bslcog4fa
  bsltrbgmvs <- c("basGambl.bsl" , "evSens.bsl",   "varSens.bsl",  "skewSens.bsl") 
  bsltrcog <- c(bsltrcog,bsltrbgmvs)
  
  mvs4faM  = colMeans(na.omit(AllD18[,hd18$bgmvs[1:4]])          )          
  mvs4faSD = apply(   na.omit(AllD18[,hd18$bgmvs[1:4]]),2,sqrtVar)
  AllD18[,c("basGambl.bsl","evSens.bsl","varSens.bsl","skewSens.bsl")] <- zScore(AllD18[,hd18$bgmvs[1:4]]) ;
  
  # 3. fishes  baseline preprocessing ###########################################################
  
  # Many more baseline measures are derived, these are the very basics:
  bslcog4fa <- c(bslcog4fa,c("vaniLnTmax.bsl",  "vaniLnCSmax.bsl", "decrLnTmax.bsl",  "decrLnCSmax.bsl"));
  # Headings for z-scored columns:
  bsltrfish <-c("vaniLnTz.bsl",  "vaniLnCSz.bsl", "decrLnTz.bsl",  "decrLnCSz.bsl"); 
  bsltrcog <- c(bsltrcog,bsltrfish)
  
  fish4faM  = colMeans(na.omit(AllD18[,c("vaniLnTmax.bsl",  "vaniLnCSmax.bsl", "decrLnTmax.bsl",  "decrLnCSmax.bsl")])          )          
  fish4faSD = apply(   na.omit(AllD18[,c("vaniLnTmax.bsl",  "vaniLnCSmax.bsl", "decrLnTmax.bsl",  "decrLnCSmax.bsl")]),2,sqrtVar)
  AllD18[,bsltrfish] <- zScore(AllD18[,c("vaniLnTmax.bsl",  "vaniLnCSmax.bsl", "decrLnTmax.bsl",  "decrLnCSmax.bsl")])
  
  # AllDhd$bslcog4fa <- bslcog4fa;
  
  # 4. Bach approach-avoidance preprocessing ####################################################
  # aka robber, predator
  # Factor scores for baseline from Dominik, new version of Dec 19, which are NOT z-scored.
  # Deliberately exclude the 'performance' factor too!
  bslcog4fa <- c(bslcog4fa, c("f1threats.bsl", "f2losss.bsl"));    # ,     "f3predperf.bsl" ));
  bsltrpred <-  c("f1threatsz.bsl", "f2losssz.bsl");                 # ,     "f3predperf.bsl" );
  bsltrcog <- c(bsltrcog,bsltrpred)

  # Mean and SD as usual, then normalize:
  pred4faM  = colMeans(na.omit(AllD18[,c("f1threats.bsl",  "f2losss.bsl")]) )
  pred4faSD = apply(   na.omit(AllD18[,c("f1threats.bsl",  "f2losss.bsl")]),2,sqrtVar)
  AllD18[,bsltrpred] <- zScore(AllD18[,c("f1threats.bsl",  "f2losss.bsl")])
  
  # AllDhd$bslcog4fa <- bslcog4fa;
  
  # 5. Delegated intertemporal discounting preprocessing ###########################################
  bslcog4fa <- c(bslcog4fa, c("Km_mLnK.bsl", "Ku_sLnK.bsl", "Xi.bsl",        "sqrt.S11.bsl",    "Tso.bsl"))
  bsltrdid  <-              c('kmz.bsl',    'sqrtkuz.bsl', 'xididto7thz.bsl', 'sqrtsrefz.bsl','sqrttsoz.bsl');
  bsltrcog <- c(bsltrcog,bsltrdid)
  
  # Do the transformations step by step for clarity:
  # Copy over:
  AllD18[,bsltrdid] <- AllD18[,c("Km_mLnK.bsl", "Ku_sLnK.bsl", "Xi.bsl", "sqrt.S11.bsl", "Tso.bsl")]
  # square-root transform the appropriate ones:
  AllD18[,c('sqrtkuz.bsl','sqrtsrefz.bsl','sqrttsoz.bsl')] <- AllD18[,c('sqrtkuz.bsl','sqrtsrefz.bsl','sqrttsoz.bsl')]^(1/2)
  # lapse rate to the one-seventh:
  AllD18[,'xididto7thz.bsl'] <- AllD18[,'xididto7thz.bsl']^(1/7)
  # Record means and SDs, then z-Score all. NOTE IN THE OTHER VARS WE WOULDN'T DO IT LIKE THIS.
  did4faM  = colMeans(na.omit(AllD18[,bsltrdid])          )
  did4faSD = apply(   na.omit(AllD18[,bsltrdid]),2,sqrtVar)
  AllD18[,bsltrdid] <- zScore(AllD18[,bsltrdid])
  
  # 6. Investor-Trustee preprocessing ##############################################################
  #    Again really basic compared to possbilities ...
  bslcog4fa <- c(bslcog4fa, c("invOpen.bsl", "respMag.bsl", "respAng.bsl"))
  bsltrtrust <- c("initinvz.bsl", "respmagz.bsl", "respangz.bsl")
  bsltrcog <- c(bsltrcog,bsltrtrust)

  # Record means and SDs, then z-Score all:
  trust4faM  = colMeans(na.omit(AllD18[,c("invOpen.bsl", "respMag.bsl", "respAng.bsl")])          )
  trust4faSD = apply(   na.omit(AllD18[,c("invOpen.bsl", "respMag.bsl", "respAng.bsl")]),2,sqrtVar)
  AllD18[,bsltrtrust] <- zScore(AllD18[,c("invOpen.bsl", "respMag.bsl", "respAng.bsl")])
  
  
  # 7. Two-step task preprocessing #################################################################
  #    As per 4 Oct on , using mfx estimates of 5 param model - but his may evolve.
  bslcog4fa <- c(bslcog4fa, c("TSTlrnrate5mxf.bsl", "TSTbeta5mxf.bsl", 
                              "TSTelig5mxf.bsl", "TSTw5mxf.bsl", "TSTprsv5mxf.bsl" ))
  bsltrtst <- c("sqrttstlrnr5.bsl", "logtstbeta5.bsl", "tstelig5.bsl", "sqrttstw5.bsl", "tstprsv5.bsl" )
  bsltrcog <- c(bsltrcog,bsltrtst)
  
  # copy over
  AllD18[,bsltrtst] <- AllD18[,c("TSTlrnrate5mxf.bsl", "TSTbeta5mxf.bsl", 
                             "TSTelig5mxf.bsl", "TSTw5mxf.bsl", "TSTprsv5mxf.bsl" )]
  # sqrt transformed:
  AllD18[,c("sqrttstlrnr5.bsl", "sqrttstw5.bsl" )] <- AllD18[,c("sqrttstlrnr5.bsl", "sqrttstw5.bsl" )]^(1/2)
  # log transf:
  AllD18[,"logtstbeta5.bsl"] <- log(AllD18[,"logtstbeta5.bsl"])
  
  # Record means and SDs and zScore
  tst4faM  = colMeans(na.omit(AllD18[,bsltrtst])          )
  tst4faSD = apply(   na.omit(AllD18[,bsltrtst]),2,sqrtVar)
  AllD18[,bsltrtst] <- zScore( AllD18[,bsltrtst])
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Record var names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  hd18$bslcog4fa <- bslcog4fa; hd18$bsltrcog <- bsltrcog

} # end if (re)do the baseline preprocessing, if (repreprocbsl) ...


## ###########################  MAIN FOLLOW-UP PREPROC ####################################### ##
if (repreprocfu2){
  
  # 8. Go-NoGo main fu (fu2) preprocessing #######################################################
  print(hd18$gng) # print(AllDhd$gng)
  # Next var., fu2cog4fa, to be augmented as we process each task.
  fu2cog4fa <- c("nspnID",
                 "G2W.lt.RT.fu2", "G2AP.lt.RT.fu2",  "invG2W.fu2", "invG2AP.fu2",         
                 "beta.fu2",  "lrnR_Appet.fu2", "lrnR_Aver.fu2",  "PavBias.fu2", "irNoiseXi.fu2", "GoBias.fu2")
  fu2trgng <- c("nspnID",
                "meanLnRT.fu2", "difLnRT.fu2",           
                "ln.beta.fu2",  "lrnRap.to4th.fu2", "lrnRav.to4th.fu2",  "ln.PavBias.fu2", "ln.irNoise.fu2", "actionBias.fu2")
  fu2trcog <- fu2trgng;  # we will redo this incrementally below. 
    
  # Reaction time orthogonalization :
  AllD18$meanLnRT.fu2 <- (log(AllD18$G2W.lt.RT.fu2) + log(AllD18$G2AP.lt.RT.fu2))/2
  AllD18$difLnRT.fu2  <- (log(AllD18$G2AP.lt.RT.fu2)- log( AllD18$G2W.lt.RT.fu2))
  
  # transf. of the others:
  AllD18$ln.beta.fu2 <- log(AllD18$beta.b2a.fu2)
  AllD18$lrnRap.to4th.fu2 <- AllD18$lrnR_Appet.b2a.fu2^(1/4)
  AllD18$lrnRav.to4th.fu2 <- AllD18$lrnR_Aver.b2a.fu2^(1/4)
  AllD18$ln.PavBias.fu2 <- log(AllD18$PavBias.b2a.fu2)
  AllD18$ln.irNoise.fu2 <- log(AllD18$irNoiseXi.b2a.fu2)
  AllD18$actionBias.fu2 <- AllD18$GoBias.b2a.fu2
  
  # Instead of z-scoring w.r.t. the fu2 set, normalize w.r.t the baseline score:
  y <- AllD18[,fu2trgng[2:length(fu2trgng)]]; 
  for (k in 1:length(gng4faM)){  y[,k] <- (y[,k] - gng4faM[k] ) / gng4faSD[k] }
  AllD18[,fu2trgng[2:length(fu2trgng)]] <- y; 
  
  # 9. fishes  fu2 preprocessing ###########################################################
  
  # Many more baseline measures are derived, these are the very basics:
  fu2cog4fa <- c(fu2cog4fa,c("vaniLnTmax.fu2",  "vaniLnCSmax.fu2", "decrLnTmax.fu2",  "decrLnCSmax.fu2"));
  # Headings for z-scored columns:
  fu2trfish <-c("vaniLnTz.fu2",  "vaniLnCSz.fu2", "decrLnTz.fu2",  "decrLnCSz.fu2"); 
  fu2trcog <- c(fu2trcog,fu2trfish)
  
  # Normalize w.r.t. baseline measures:
  y <- AllD18[,c("vaniLnTmax.fu2",  "vaniLnCSmax.fu2", "decrLnTmax.fu2",  "decrLnCSmax.fu2")]; 
  for (k in 1:length(fish4faM)){  y[,k] <- (y[,k] - fish4faM[k] ) / fish4faSD[k] }
  AllD18[,fu2trfish] <- y; 
  
   # 10. Bach approach-avoidance fu2 preprocessing ###################################################
  # *** MOVED BELOW *** 
  
  # 11. Delegated intertemporal discounting fu2 preprocessing ##########################################
  fu2cog4fa <- c(fu2cog4fa, c("Km_mLnK.fu2", "Ku_sLnK.fu2", "Xi.fu2",        "sqrt.S11.fu2",    "Tso.fu2"))
  fu2trdid  <-              c('kmz.fu2',    'sqrtkuz.fu2', 'xididto7thz.fu2', 'sqrtsrefz.fu2','sqrttsoz.fu2');
  fu2trcog <- c(fu2trcog,fu2trdid)
  
  # Do the transformations step by step for clarity:
  # Copy over:
  AllD18[,fu2trdid] <- AllD18[,c("Km_mLnK.fu2", "Ku_sLnK.fu2", "Xi.fu2", "sqrt.S11.fu2", "Tso.fu2")]
  # square-root transform the appropriate ones:
  AllD18[,c('sqrtkuz.fu2','sqrtsrefz.fu2','sqrttsoz.fu2')] <- AllD18[,c('sqrtkuz.fu2','sqrtsrefz.fu2','sqrttsoz.fu2')]^(1/2)
  # lapse rate to the one-seventh:
  AllD18[,'xididto7thz.fu2'] <- AllD18[,'xididto7thz.fu2']^(1/7)
  
  
  # Normalize w.r.t. baseline measures:
  y <- AllD18[,fu2trdid]; 
  for (k in 1:length(did4faM)){  y[,k] <- (y[,k] - did4faM[k] ) / did4faSD[k] }
  AllD18[,fu2trdid] <- y; 
  
  # 12. Investor-Trustee fu2 preprocessing ##############################################################
  #    Again really basic compared to possbilities ...
  fu2cog4fa <- c(fu2cog4fa, c("invOpen.fu2", "respMag.fu2", "respAng.fu2"))
  fu2trtrust <- c("initinvz.fu2", "respmagz.fu2", "respangz.fu2")
  fu2trcog <- c(fu2trcog,fu2trtrust)
  
  # Normalize w.r.t. baseline measures:
  y <- AllD18[,c("invOpen.fu2", "respMag.fu2", "respAng.fu2")]; 
  for (k in 1:length(trust4faM)){  y[,k] <- (y[,k] - trust4faM[k] ) / trust4faSD[k] }
  AllD18[,fu2trtrust] <- y; 
  
  # hd18$fu2cog4fa <- fu2cog4fa; hd18$fu2trcog <- fu2trcog
  
  # 13. Two-step task preprocessing #################################################################
  #    As per 4 Oct, using mfx estimates of 5 param model - but his may evolve.
  fu2cog4fa <- c(fu2cog4fa, c("TSTlrnrate5mxf.fu2", "TSTbeta5mxf.fu2", 
                              "TSTelig5mxf.fu2", "TSTw5mxf.fu2", "TSTprsv5mxf.fu2" ))
  fu2trtst <- c("sqrttstlrnr5.fu2", "logtstbeta5.fu2", "tstelig5.fu2", "sqrttstw5.fu2", "tstprsv5.fu2" )
  fu2trcog <- c(fu2trcog,fu2trtst)
  
  # copy over
  AllD18[,fu2trtst] <- AllD18[,c("TSTlrnrate5mxf.fu2", "TSTbeta5mxf.fu2", 
                                 "TSTelig5mxf.fu2", "TSTw5mxf.fu2", "TSTprsv5mxf.fu2" )]
  # sqrt transformed:
  AllD18[,c("sqrttstlrnr5.fu2", "sqrttstw5.fu2" )] <- AllD18[,c("sqrttstlrnr5.fu2", "sqrttstw5.fu2" )]^(1/2)
  # log transf:
  AllD18[,"logtstbeta5.fu2"] <- log(AllD18[,"logtstbeta5.fu2"])
  
  
  # Normalize w.r.t. baseline measures:
  y <- AllD18[,fu2trtst]; 
  for (k in 1:length(tst4faM)){  y[,k] <- (y[,k] - tst4faM[k] ) / tst4faSD[k] }
  AllD18[,fu2trtst] <- y; 
  
  
  # 10. Bach approach-avoidance fu2 preprocessing ###################################################
  # aka robber, predator

  # Factor scores for baseline from Dominik, new version of Dec 19, which are NOT z-scored.
  # Deliberately exclude the 'performance' factor too!
  fu2cog4fa <- c(fu2cog4fa, c("f1threats.fu2", "f2losss.fu2"));    # ,     "f3predperf.fu2" ));
  fu2trpred <-  c("f1threatsz.fu2", "f2losssz.fu2");                 # ,     "f3predperf.fu2" );
  fu2trcog <- c(fu2trcog,fu2trpred)
  
  # normalize w.r.t. baseline: 
  y <- AllD18[,c("f1threats.fu2", "f2losss.fu2")]; 
  for (k in 1:length(pred4faM)){  y[,k] <- (y[,k] - pred4faM[k] ) / pred4faSD[k] }
  AllD18[,fu2trpred] <- y; 
  
  # NOT TO DO: 14. Happiness fu2 preprocessing #####################################################################
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Record var names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  hd18$fu2cog4fa <- fu2cog4fa; hd18$fu2trcog <- fu2trcog
  
}  # end if (repreprocfu2)

if (redefinemissing){
  # ~~~~~~~~~~~~~~~~~~~~~~ Objects to keep track of missing files etc ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Missing as of Jan 18 - Benjamin wrote:
  missingJan18 <- list();
  missingJan18$decAcNB = 'For these subjects Benjamin had imaging but not decision acuity data'
  missingJan18$decAc = c(13771, 18150, 26740, 27433, 27722, 28647, 33373, 34025, 34421, 34843, 
                         35501, 37903, 40758, 40980, 41053, 41095, 41806, 41814, 41897, 43943, 
                         46383, 46730, 46839, 46888, 46912, 46920, 46987, 47159, 47563, 47571, 
                         47720, 47738, 47829, 437, 47860, 47910, 47928, 48132, 48199, 48231, 
                         48348, 48371, 48439, 48454, 48504, 48512, 48520, 59238);
  missingJan18$IQsubscNB = 'For these subjects had imaging but not raw matrix or vocab IQ scores'
  missingJan18$IQsubsc = c(13771, 27722, 33373, 34843, 46383, 47571, 47720, 48199, 48371, 59238);
  missingJan18$IQfullNB = 'For these subjects had imaging but not IQ scores (wasi_zz_iq_full2_iq)'
  missingJan18$IQfull = c(13771, 27722, 28142, 33373, 33845, 34843, 41095, 46383, 46920, 47571, 47720, 48199, 48348,
                          48371, 48512, 59238);
  missingJan18$no.task.dataNB = 'no task data wha tsoever as per Compendium 2 Mar 18'
  missingJan18$no.task.data = c(13771,33373,46383,46987);   
  missingJan18$foundNB = 'were not in Compendium but have been recovered (by Gita Prabhu team)'
  missingJan18$found = c(41095); 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}

# DECFA18 <- function(){   # Dummy function header (and closing, below )

#library(Rcmdr)
library(psych)
library(missMDA) # This loads also required package FactoMineR
library(nFactors)
library(lavaan)
library(sem)
library(corrplot)
library(caret)
library(car)
library(ppcor)

## Interpolate missing and optionally create new data partiton - 
#  First, here, with 2K, non-Depression, iua baseline data
if (redotraintest){
  #  including interpolation of missing values:
  bsl2k4fa <- AllD18[AllD18$nspnID %in% IDs18$iuabsl[IDs18$iuabsl[,2]=='2K_Cohort',1], hd18$bsltrcog]; 
  totvarn <- dim(bsl2k4fa)[2]-1 
  # find number of principal components, up to a max. to be checked later, 
  # by cross-validation:
  NPC1 <- estim_ncpPCA(bsl2k4fa[,2:(totvarn+1)],ncp.min=0,ncp.max=15,method.cv="Kfold",nbsim=50)  # try ncp.max = 15 ??
  # Use this to impute missing values:
  impbsl2k <- imputePCA(bsl2k4fa[,2:(totvarn+1)],ncp=NPC1$ncp)
  zParD1 <- impbsl2k$completeObs;
  
  prcPar1  <- prcomp(zParD1,center=TRUE,scale=TRUE); summary(prcPar1)
  UPar1 <- prcPar1$rotation
  # Display 5 components nicely:
  for (nv in 5:1){
      load1 <- prcPar1$rotation; 
      sorted.load1 <- load1[order(load1[,nv]),nv];  
      dotchart(abs(sorted.load1),main=paste('|loadings| for cognitive PC ',nv),xlab='|loadings|')
  }
  zParD1 <- data.frame(bsl2k4fa$nspnID,zParD1)  # put the IDs back in the dataframe
  colnames(zParD1)[1] <- 'nspnID'
  
  if (repartitiondata){
     trainData <- createDataPartition(zParD1$nspnID,  p=.50, list=F)
  }
  D.Train <- zParD1[trainData,2:dim(bsl2k4fa)[2]]
  D.Test <- zParD1[-trainData,2:dim(bsl2k4fa)[2]]

} # end  if redotraintest

## Determine the number of factors to extract 

library(nFactors)
ev <- eigen(cor(D.Train)) # get eigenvalues
ap <- parallel(subject=nrow(D.Train),var=ncol(D.Train),
  rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS, main='Scree Test solutions - non-Depr. train set'); 

## Run the factor analysis on no. of factors suggested by the scree plot
factN <- nS[[1]]$noc; # see nS[[1]]
nocfit<-fa(r=cor(D.Train), nfactors=factN,fm="minres"); nocfit$loadings
# Minimal EFA with just one factor:
fit1<-fa(r=cor(D.Train), nfactors=1,fm="minres"); fit1$loadings

## The FA on whole sample as suggested by 
#  the SEM analysis below:
fit3 <- fa(r= zParD1[,2:dim(zParD1)[2]], 
             nfactors=3,
             fm="minres",  scores="Bartlett",  rotate="varimax"); # fit3$loadings
fit4 <- fa(r= zParD1[,2:dim(zParD1)[2]], 
           nfactors=4,
           fm="minres",  scores="Bartlett",  rotate="varimax"); # fit4$loadings

# Then the fully exploratory, 3 and 4 component for whole sample. 
# At least the3 last two are typically very similar.
for (nv in factN:1) {
  load1 <- nocfit$loadings[,];
  for (nv in factN:1) {
    sorted.load1 <- load1[order(load1[,nv]),nv];  
    dotchart(abs(sorted.load1),main=paste('TRAIN DATA |loadings|, ECFA component ',nv),xlab='|loadings|')
    abline(v=0.2)
  }
  if (nv <=3){
     fload1 <- fit3$loadings[,];
     sorted.fload1 <- fload1[order(fload1[,nv]),nv];  
     dotchart(abs(sorted.fload1),main=paste('|loadings| for 3-factor component ',nv),xlab='|loadings|')
  }
  if (nv <= 4) {
     fload1 <- fit4$loadings[,];
     sorted.fload1 <- fload1[order(fload1[,nv]),nv];  
     dotchart(abs(sorted.fload1),main=paste('|loadings| for 4-factor component ',nv),xlab='|loadings|')
     abline(v=0.2)
  }
}


# Tidy up the scores for the 3 and 4 factor fits:
fa3decAc <- as.data.frame(fit3$scores); 
fa3decAc$nspnID <- zParD1$nspnID
colnames(fa3decAc)[1:3] <- c('decAc.bsl', 'f2FA3.bsl', 'f3FA3.bsl');
fa4decAc <- as.data.frame(fit4$scores); 
fa4decAc$nspnID <- zParD1$nspnID
colnames(fa4decAc)[1:4] <- c('decAc.bsl', 'f2FA4.bsl', 'f3FA4.bsl', 'f3FA4.bsl');

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  CFA  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Will put together headings so as to make models for fitting with
#  different combinations of variables, using adjustable thresholds
#  (hyperparameter) from the ECFA

## Old example done by hand:
#.model <- 
#   c('Factor.1: eligTraceTST, ethnicity, gender, Inv_Temp.1, ln.bet_AppetZ, ln.bet_AverZ, ln.G2AP.lt.RTZ, ln.G2W.lt.RTZ, ln.gngXiZ, ln.KuZ, ln.lrnRateZ, ln.PavBiasZ, MVS_prefSkZ, respMagZ, sepCSvani_jtc09Z, TinvTST, XiTST')
#  .model <- cfa(file=textConnection(.model), reference.indicators=FALSE)
#   .Data <- D.Test[, c('eligTraceTST', 'ethnicity', 'gender', 'Inv_Temp.1', 'ln.bet_AppetZ', 'ln.bet_AverZ', 'ln.G2AP.lt.RTZ', 
#   'ln.G2W.lt.RTZ', 'ln.gngXiZ', 'ln.KuZ', 'ln.lrnRateZ', 'ln.PavBiasZ', 'MVS_prefSkZ', 'respMagZ', 'sepCSvani_jtc09Z', 'TinvTST', 
#   'XiTST')] 
#   .Data <- as.data.frame(scale(.Data))
#   summary(sem(.model, data=.Data), robust=FALSE, fit.indices=c("AIC","BIC","RMSEA","NFI","NNFI","CFI"))

thr <- 0.25;   # can try 0.2
strF <- list();  hdF <- list(); hdIncr <- list();
for (nv in 1:factN){
  hdF[[nv]] <- rownames(load1)[abs(load1[,nv])>= thr];
  if (nv==1){
    hdIncr[[nv]] <- hdF[[nv]];
  } else {
    hdIncr[[nv]] <- union(hdF[[nv]],hdIncr[[nv-1]]);
  }
  strF[[nv]] <- paste(hdF[[nv]],collapse=', ');
  strF[[nv]] <- paste('Factor.',nv,': ',strF[[nv]],sep='');
}

FAsem <- list();

## One factor 
model1 <- c(strF[[1]])
model1 <- cfa(file=textConnection(model1), reference.indicators=FALSE)
.Data <- D.Test[, hdIncr[[1]]] 
.Data <- as.data.frame(scale(.Data))
FAsem[[1]] <- sem(model1, data=.Data);
summary(FAsem[[1]], robust=FALSE, fit.indices=c("AIC","BIC","RMSEA","NFI","NNFI","CFI"))


## Two factors 
model2 <- c(strF[[1]],strF[[2]])
   model2 <- cfa(file=textConnection(model2), reference.indicators=FALSE)
   .Data <- D.Test[, hdIncr[[2]]] 
   .Data <- as.data.frame(scale(.Data))
   FAsem[[2]] <- sem(model2, data=.Data);
   summary(FAsem[[2]], robust=FALSE, fit.indices=c("AIC","BIC","RMSEA","NFI","NNFI","CFI"))

# Three factors 
model3 <-   c(strF[[1]],strF[[2]], strF[[3]]) # , strF[[4]])
model3 <- cfa(file=textConnection(model3), reference.indicators=FALSE)
.Data <- D.Test[, hdIncr[[3]] ] #   union(hdF[[1]] ... hdF[[3]]))] 
.Data <- as.data.frame(scale(.Data))
FAsem[[3]] <- sem(model3, data=.Data);
summary(FAsem[[3]], robust=FALSE, fit.indices=c("AIC","BIC","RMSEA","NFI","NNFI","CFI"))

# Four factors 
model4 <-   c(strF[[1]],strF[[2]], strF[[3]], strF[[4]]) # , strF[[4]])
model4 <- cfa(file=textConnection(model4), reference.indicators=FALSE)
.Data <- D.Test[, hdIncr[[4]] ] #   union(hdF[[1]] ... hdF[[3]]))] 
.Data <- as.data.frame(scale(.Data))
FAsem[[4]] <- sem(model4, data=.Data);
summary(FAsem[[4]], robust=FALSE, fit.indices=c("AIC","BIC","RMSEA","NFI","NNFI","CFI"))

# Five factors 
model5 <-   c(strF[[1]],strF[[2]], strF[[3]], strF[[4]], strF[[5]] ) # , strF[[4]])
model5 <- cfa(file=textConnection(model5), reference.indicators=FALSE)
.Data <- D.Test[, hdIncr[[5]] ] #   union(hdF[[1]] ... hdF[[5]]))] 
.Data <- as.data.frame(scale(.Data))
FAsem[[5]] <- sem(model5, data=.Data);
summary(FAsem[[5]], robust=FALSE, fit.indices=c("AIC","BIC","RMSEA","NFI","NNFI","CFI"))

# Eight factors 
model8 <-   c(strF[[1]],strF[[2]], strF[[3]], strF[[4]], strF[[5]], strF[[6]], strF[[7]], strF[[8]] ) # , strF[[4]])
model8 <- cfa(file=textConnection(model8), reference.indicators=FALSE)
.Data <- D.Test[, hdIncr[[8]] ] #   union(hdF[[1]] ... hdF[[8]]))] 
.Data <- as.data.frame(scale(.Data))
FAsem[[8]] <- sem(model8, data=.Data);
summary(FAsem[[8]], robust=FALSE, fit.indices=c("AIC","BIC","RMSEA","NFI","NNFI","CFI"))


# ~~~~~~~~~~~~~~~~~ Exploration not excluding the 'Depression' pts ~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~ first, imputation for bsl: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bsliua4fa <- AllD18[AllD18$nspnID %in% IDs18$iuabsl[,1], hd18$bsltrcog]; 
totvarn <- dim(bsliua4fa)[2]-1 
# find number of principal components, up to a max. to be checked later, 
# by cross-validation:
NPCallbsl <- estim_ncpPCA(bsliua4fa[,2:(totvarn+1)],ncp.min=0,ncp.max=15,method.cv="Kfold",nbsim=50)  # try ncp.max = 10 ??
# Use this to impute missing values:
impbsliua <- imputePCA(bsliua4fa[,2:(totvarn+1)],ncp=NPCallbsl$ncp)
zParDall <- impbsliua$completeObs;
## ~~~~~~~~~~~~~~~~~~~~~~~~ then, imputation for fu2: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fu2iua4fa <- AllD18[AllD18$nspnID %in% IDs18$iuafu2[,1], hd18$fu2trcog]; 
fu2varn <- dim(fu2iua4fa)[2]-1 
# find number of principal components, up to a max. to be checked later, 
# by cross-validation:
NPCallfu2 <- estim_ncpPCA(fu2iua4fa[,2:(fu2varn+1)],ncp.min=0,ncp.max=15,method.cv="Kfold",nbsim=50)  # try ncp.max = 10 ??
# Use this to impute missing values:
impfu2iua <- imputePCA(fu2iua4fa[,2:(fu2varn+1)],ncp=NPCallfu2$ncp)
z6allfu2 <- impfu2iua$completeObs;
z6allfu2 <- data.frame(fu2iua4fa$nspnID,z6allfu2); colnames(z6allfu2)[1] <- 'nspnID'
##          ~~~~~~~~~~~~~~~~~~~~~~~~  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assemble in same order as baseline and add dummies, as will be heading order for z6allbsl further down,
# in    ### Scores excluding Roulette task ... :
test <- z6allfu2[,1:9]
test[,10:13] <- 0; colnames(test)[10:13] <- c("basGambl.fu2", "evSens.fu2"  , "varSens.fu2" , "skewSens.fu2")  # dummies
test[,14:32] <- z6allfu2[,c("vaniLnTz.fu2",     "vaniLnCSz.fu2",    "decrLnTz.fu2" ,    "decrLnCSz.fu2" , 
                            "f1threatsz.fu2" ,  "f2losssz.fu2" ,  
                            "kmz.fu2" ,         "sqrtkuz.fu2"  ,    "xididto7thz.fu2" , "sqrtsrefz.fu2"  , "sqrttsoz.fu2"  ,   
                            "initinvz.fu2" ,    "respmagz.fu2"   ,  "respangz.fu2"  ,   
                            "sqrttstlrnr5.fu2",   "logtstbeta5.fu2",  "tstelig5.fu2"  ,   "sqrttstw5.fu2"   , "tstprsv5.fu2" )]
z6all <- merge(z6allbsl,test,by.x='nspnID',by.y='nspnID',all=TRUE) # mostly for reference ... 
if (save2disk){exportCSV(paste(wherearewe$uchjan17on,'z6allbslfu2.csv',sep=''),z6all)}
z6allfu2 <- test; remove(test)

# Baseline 
prcPar1  <- prcomp(zParDall,center=TRUE,scale=TRUE); summary(prcPar1)
UPar1 <- prcPar1$rotation
# Display 5 components nicely:
for (nv in 5:1){
  load1a <- prcPar1$rotation; 
  sorted.load1a <- load1a[order(load1a[,nv]),nv];  
  dotchart(abs(sorted.load1a),main=paste('|loadings| for cog. PC incl. Dep',nv),xlab='|loadings|')
}
zParDall <- data.frame(bsliua4fa$nspnID,zParDall)  # put the IDs back in the dataframe
colnames(zParDall)[1] <- 'nspnID'

# ####  (re) create train-test partition(s) if need be  ####
if (redotraintest){
  if (repartitiondata){
      trainDall <- createDataPartition(zParDall$nspnID,  p=.50, list=F)
  }
  Dall.Train <- zParDall[trainDall,2:dim(bsliua4fa)[2]]
  Dall.Test <- zParDall[-trainDall,2:dim(bsliua4fa)[2]]
}

## Determine the number of factors to extract if Depr. included,
#  From training set :
ev <- eigen(cor(Dall.Train[,2:dim(Dall.Train)[2]])) # get eigenvalues
ap <- parallel(subject=nrow(Dall.Train[,2:dim(Dall.Train)[2]]),
               var=ncol(Dall.Train[,2:dim(Dall.Train)[2]]),
               rep=100,cent=.05)
nSall<- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nSall,main='Training half of full cohort')
# Compare with same for whole:
ev <- eigen(cor( zParDall[,2:dim(bsliua4fa)[2]] )) # get eigenvalues
ap <- parallel(subject=nrow( zParDall[,2:dim(bsliua4fa)[2]]  ),
               var=ncol( zParDall[,2:dim(bsliua4fa)[2]] ),
               rep=100,cent=.05)
nSwhole<- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nSwhole,main='Whole dataset Incl. Depr. cohort')


## Run the factor analysis on no. of factors suggested by the scree plot
factN <- nSall[[1]]$noc; # see nS[[1]]
nocfitall<-fa(r=cor(Dall.Train[,2:dim(Dall.Train)[2]]), 
              nfactors=factN,fm="minres"); nocfitall$loadings
# EFAs on the training set:
fit1atr<-fa(r=cor(Dall.Train[,2:dim(Dall.Train)[2]]),
            nfactors=1,fm="minres"); fit1atr$loadings

## The FA on whole sample for 2-6 factors. The idea is to
# demonstrate robustness to hyperparams eventually.
fit2a <- fa(r= zParDall[,2:dim(zParDall)[2]], 
            nfactors=2,
            fm="minres",  scores="Bartlett",  rotate="varimax"); 
fit3a <- fa(r= zParDall[,2:dim(zParDall)[2]], 
           nfactors=3,
           fm="minres",  scores="Bartlett",  rotate="varimax"); 
fit4a <- fa(r= zParDall[,2:dim(zParDall)[2]], 
           nfactors=4,
           fm="minres",  scores="Bartlett",  rotate="varimax"); 
fit5a <- fa(r= zParDall[,2:dim(zParDall)[2]], 
            nfactors=5,
            fm="minres",  scores="Bartlett",  rotate="varimax"); 
fit6a <- fa(r= zParDall[,2:dim(zParDall)[2]], 
            nfactors=6,
            fm="minres",  scores="Bartlett",  rotate="varimax"); 
fit7a <- fa(r= zParDall[,2:dim(zParDall)[2]], 
            nfactors=7,
            fm="minres",  scores="Bartlett",  rotate="varimax"); 

# See how much diffence it would make to the decAc estimate if we
# didn't use the 4-factor but another solution:
y = data.frame(fit4a$scores[,1], fit2a$scores[,1]); co=pcor(na.omit(y)); 
plot(y,xlab='decAc 4-factor',ylab='decAc 2-factor', pch=21, col='black',
     main=paste('r=',round(co$est[1,2],4))); abline(0,1,col='cyan3',lwd=3)
y = data.frame(fit4a$scores[,1], fit3a$scores[,1]); co=pcor(na.omit(y)); 
plot(y,xlab='decAc 4-factor',ylab='decAc 3-factor',pch=21, col='black',
     main=paste('r=',round(co$est[1,2],4))); abline(0,1,col='cyan4',lwd=3)
y = data.frame(fit4a$scores[,1], fit5a$scores[,1]); co=pcor(na.omit(y)); 
plot(y,xlab='decAc 4-factor',ylab='decAc 5-factor',pch=21, col='black',
     main=paste('r=',round(co$est[1,2],4))); abline(0,1,col='blue2',lwd=3)
y = data.frame(fit4a$scores[,1], -fit6a$scores[,1]); co=pcor(na.omit(y)); 
plot(y,xlab='decAc 4-factor',ylab='decAc 6-factor',pch=21, col='black',
     main=paste('r=',round(co$est[1,2],4))); abline(0,1,col='blue4',lwd=3)
y = data.frame(fit4a$scores[,1], -fit7a$scores[,1]); co=pcor(na.omit(y)); 
plot(y,xlab='decAc 4-factor',ylab='decAc 7-factor',pch=21, col='black',
     main=paste('r=',round(co$est[1,2],4))); abline(0,1,col='blue4',lwd=3)


# Plot the fully exploratory components for whole sample, and the training one. 
# These should be similar. The idea is to put one analysis in the main paper, and
# provide evidence that the structure changes little, and the scores very little,
# if alternatives are used. 
for (nv in 6:1) {
  if (nv <= factN){
    load1a <- nocfitall$loadings[,];
    sorted.load1a <- load1a[order(load1a[,nv]),nv];  
    dotchart(abs(sorted.load1a),main=paste('TRAIN DATA |loadings|, ECFA component ',nv),xlab='|loadings|')
    abline(v=0.25, col='gray')
  }
  if (nv <=3){
    fload1a <- fit3a$loadings[,];
    sorted.fload1a <- fload1a[order(fload1a[,nv]),nv];  
    dotchart(abs(sorted.fload1a),main=paste('|loadings| for 3-factor component ',nv),xlab='|loadings|')
    abline(v=0.25, col='gray')
  }
  if (nv <= 4) {
    fload1a <- fit4a$loadings[,];
    sorted.fload1a <- fload1a[order(fload1a[,nv]),nv];  
    dotchart(abs(sorted.fload1a),main=paste('|loadings| for 4-factor component ',nv),xlab='|loadings|')
    abline(v=0.25, col='gray')
  }
  if (nv <= 5) {
    fload1a <- fit5a$loadings[,];
    sorted.fload1a <- fload1a[order(fload1a[,nv]),nv];  
    dotchart(abs(sorted.fload1a),main=paste('|loadings| for 5-factor component ',nv),xlab='|loadings|')
    abline(v=0.25, col='gray')
  }
}


# Tidy up the scores for the 3 and 4 factor fits:
fa3decAcAll <- as.data.frame(fit3a$scores); 
fa3decAcAll$nspnID <- zParDall$nspnID
colnames(fa3decAcAll)[1:3] <- c('decAcAll.bsl', 'f2FA3All.bsl', 'f3FA3All.bsl');
fa4decAcAll <- as.data.frame(fit4a$scores); 
fa4decAcAll$nspnID <- zParDall$nspnID
colnames(fa4decAcAll)[1:4] <- c('decAc4All.bsl', 'f2FA4All.bsl', 'f3FA4All.bsl', 'f4FA4All.bsl');


# #####################  CFA #####################
load1a <- nocfitall$loadings[,];
thr <- 0.25;   # can try 0.2
strF <- list();  hdF <- list(); hdIncr <- list();
for (nv in 1:factN){
  hdF[[nv]] <- rownames(load1a)[abs(load1a[,nv])>= thr];
  if (nv==1){
    hdIncr[[nv]] <- hdF[[nv]];
  } else {
    hdIncr[[nv]] <- union(hdF[[nv]],hdIncr[[nv-1]]);
  }
  strF[[nv]] <- paste(hdF[[nv]],collapse=', ');
  strF[[nv]] <- paste('Factor.',nv,': ',strF[[nv]],sep='');
}

FAsemAll <- list();

## One factor 
allcfa1 <- c(strF[[1]])
allcfa1 <- cfa(file=textConnection(allcfa1), reference.indicators=FALSE)
.Data <- Dall.Test[, hdIncr[[1]]] 
.Data <- as.data.frame(scale(.Data))
FAsemAll[[1]] <- sem(allcfa1, data=.Data);
summary(FAsemAll[[1]], robust=FALSE, fit.indices=c("AIC","BIC","RMSEA","NFI","NNFI","CFI"))

## Two factors
allcfa2 <- c(strF[[1]],strF[[2]])
allcfa2 <- cfa(file=textConnection(allcfa2), reference.indicators=FALSE)
.Data <- Dall.Test[, hdIncr[[2]]] 
.Data <- as.data.frame(scale(.Data))
FAsemAll[[2]] <- sem(allcfa2, data=.Data);
summary(FAsemAll[[2]], robust=FALSE, fit.indices=c("AIC","BIC","RMSEA","NFI","NNFI","CFI"))

## Three factors
allcfa3 <- c(strF[[1]],strF[[2]],strF[[3]])
allcfa3 <- cfa(file=textConnection(allcfa3), reference.indicators=FALSE)
.Data <- Dall.Test[, hdIncr[[3]]] 
.Data <- as.data.frame(scale(.Data))
FAsemAll[[3]] <- sem(allcfa3, data=.Data);
summary(FAsemAll[[3]], robust=FALSE, fit.indices=c("AIC","BIC","RMSEA","NFI","NNFI","CFI"))

## Four factors
allcfa4 <- c(strF[[1]],strF[[2]],strF[[3]],strF[[4]] )
allcfa4 <- cfa(file=textConnection(allcfa4), reference.indicators=FALSE)
.Data <- Dall.Test[, hdIncr[[4]]] 
.Data <- as.data.frame(scale(.Data))
FAsemAll[[4]] <- sem(allcfa4, data=.Data);
summary(FAsemAll[[4]], robust=FALSE, fit.indices=c("AIC","BIC","RMSEA","NFI","NNFI","CFI"))

## Five factors
allcfa5 <- c(strF[[1]],strF[[2]],strF[[3]],strF[[4]],
             strF[[5]])
allcfa5 <- cfa(file=textConnection(allcfa5), reference.indicators=FALSE)
.Data <- Dall.Test[, hdIncr[[5]]] 
.Data <- as.data.frame(scale(.Data))
FAsemAll[[5]] <- sem(allcfa5, data=.Data);
summary(FAsemAll[[5]], robust=FALSE, fit.indices=c("AIC","BIC","RMSEA","NFI","NNFI","CFI"))

## Six factors
allcfa6 <- c(strF[[1]],strF[[2]],strF[[3]],strF[[4]],
             strF[[5]], strF[[6]])
allcfa6 <- cfa(file=textConnection(allcfa6), reference.indicators=FALSE)
.Data <- Dall.Test[, hdIncr[[6]]] 
.Data <- as.data.frame(scale(.Data))
FAsemAll[[6]] <- sem(allcfa6, data=.Data);
summary(FAsemAll[[6]], robust=FALSE, fit.indices=c("AIC","BIC","RMSEA","NFI","NNFI","CFI"))

## Seven factors
allcfa7 <- c(strF[[1]],strF[[2]],strF[[3]],strF[[4]],
             strF[[5]],strF[[6]],strF[[7]])
allcfa7 <- cfa(file=textConnection(allcfa7), reference.indicators=FALSE)
.Data <- Dall.Test[, hdIncr[[7]]] 
.Data <- as.data.frame(scale(.Data))
FAsemAll[[7]] <- sem(allcfa7, data=.Data);
summary(FAsemAll[[7]], robust=FALSE, fit.indices=c("AIC","BIC","RMSEA","NFI","NNFI","CFI"))

# ##################### Scores excluding Roulette task, all pts but 6 tasks only: ######################
z6allbsl <- zParDall
z6allbsl[,bsltrbgmvs] <- 0
# analogous to fa4decAcAll :
fa4decAcAl6 <- list();
fa4decAcAl6[[1]] <- as.matrix(z6allbsl[,2:dim(z6allbsl)[2]]) %*% fit4a$loadings
fa4decAcAl6[[1]] <- data.frame(z6allbsl$nspnID, fa4decAcAl6[[1]])
colnames(fa4decAcAl6[[1]])[1] <- 'nspnID'
for (v in 5:2){
   co=pcor(data.frame(fa4decAcAll[,v-1],fa4decAcAl6[[1]][,v])); 
   plot(fa4decAcAll[,v-1],fa4decAcAl6[[1]][,v], 
        main=paste('Component',v-1, 'with vs. without MVS; r=', round(co$est[2,1],3)));  
   abline(0,1)
}

# Now using the same weights for follow-up:
fa4decAcAl6[[2]] <- as.matrix(z6allfu2[,2:dim(z6allfu2)[2]]) %*% fit4a$loadings
fa4decAcAl6[[2]] <- data.frame(z6allfu2$nspnID, fa4decAcAl6[[2]]);
colnames(fa4decAcAl6[[2]])[1] <- 'nspnID'
fa4decAcAl6[[3]] <- merge(fa4decAcAl6[[1]],fa4decAcAl6[[2]],by.x='nspnID',by.y='nspnID',all=TRUE)
colnames(fa4decAcAl6[[3]]) <- c('nspnID','decAc4b.bsl','f2b.bsl','f3b.bsl','f4b.bsl','decAc4b.fu2','f2b.fu2','f3b.fu2','f4b.fu2')
# Now check convention to make sure decAc 'the right way around':
# This should be r < 0 :
co <- pcor(na.omit(data.frame(fa4decAcAl6[[3]][,'decAc4b.bsl'],zParDall[,"vaniLnTz.bsl"])));
if (co$est[2,1] > 0){ fa4decAcAl6[[3]][,2:5] <- -fa4decAcAl6[[3]][,2:5]}
co <- pcor(na.omit(data.frame(fa4decAcAl6[[3]][,'decAc4b.fu2'],z6all[,"vaniLnTz.fu2"])));
if (co$est[2,1] > 0){ fa4decAcAl6[[3]][,6:9] <- -fa4decAcAl6[[3]][,6:9]}
co <- pcor(na.omit(data.frame(fa4decAcAll[,'decAc4All.bsl'],z6all[,"vaniLnTz.fu2"])));
if (co$est[2,1] > 0){ 
  fa4decAcAll[,c('decAc4All.bsl', 'f2FA4All.bsl', 'f3FA4All.bsl', 
                 'f4FA4All.bsl')] <- -fa4decAcAll[,c('decAc4All.bsl', 'f2FA4All.bsl', 'f3FA4All.bsl', 
                                                     'f4FA4All.bsl')] 
}
if (save2disk){ exportCSV(paste(wherearewe$uchjan17on,'decAc4All.bsl.csv',sep=''),fa4decAcAll) }

# Simple display of reliabilities etc:
for (k in 5:2){
  co=pcor(na.omit(fa4decAcAl6[[3]][,c(k, k+4)])); 
  print(co)
  plot( fa4decAcAl6[[3]][,k], fa4decAcAl6[[3]][,k+4],
       main=paste('Component',k-1, '(no MVS), bsl vs. fu2; r=', round(co$est[2,1],3)));  
  abline(0,1)
  
}

if (save2disk){ exportCSV(paste(wherearewe$uchjan17on,'decAc4b.csv',sep=''),fa4decAcAl6[[3]]) }

#                     #####################  ######################


# Invariance by sex and age
# First div up the data:
d4invAll <- merge(AllD18[,hd18$cms[c(1,4)]], pseudoage[pseudoage[,2]=='iua_baseline_arm_1',c(1,6)], by='nspnID')
d4invAll <- merge(d4invAll, zParDall, by='nspnID')
dwom <- d4invAll[d4invAll[,2]=='Female',4:dim(d4invAll)[2]]
dman <- d4invAll[d4invAll[,2]=='Male',4:dim(d4invAll)[2]]
medpsage <- median(na.omit(d4invAll[,"cog_pseudo_age"]))
dyoung <- d4invAll[vecTRUE(d4invAll[,3] < medpsage),4:dim(d4invAll)[2]]
dold <- d4invAll[vecTRUE(d4invAll[,3] >= medpsage),4:dim(d4invAll)[2]]

dtrain = Dall.Train[,2:dim(Dall.Train)[2]];
dtest = Dall.Test[ ,2:dim(Dall.Test )[2]]

faTrain<-fa(r=as.matrix(dtrain), nfactors=4,fm="minres",  rotate="varimax", scores="components"); 
faTest <-fa(r=as.matrix(dtest), nfactors=4,fm="minres",  rotate="varimax", scores="components"); 

faWom <- fa(r= as.matrix(dwom), nfactors=4, fm="minres",  rotate="varimax", scores="components"); 
faMan <- fa(r= as.matrix(dman), nfactors=4, fm="minres",  rotate="varimax", scores="components"); 

faYoung <- fa(r= as.matrix(dyoung), nfactors=4, fm="minres",  rotate="varimax", scores="components"); 
faOld <- fa(r= as.matrix(dold), nfactors=4, fm="minres",  rotate="varimax", scores="components"); 

# As a baseline, see how much the train and test sets, which should not
# differ, do differ:
test = (dds2w %*% faTest$loadings)
co = pcor(data.frame(test[,1],faTrain$scores[,1]))
plot(test[,1],faTrain$scores[,1],main=paste('r=',round(co$est[1,2],4)),
     xlab='Train based on Test', ylab='own Train'); abline(0,1)
# Now same for sex:
test = (as.matrix(dman) %*% faWom$loadings)
co = pcor(data.frame(test[,1],-faMan$scores[,1]))
plot(test[,1],-faMan$scores[,1],main=paste('r=',round(co$est[1,2],4)),
     xlab='Male based on Female', ylab='Male own'); abline(0,1)
test = (as.matrix(dwom) %*% faMan$loadings)
co = pcor(data.frame(test[,1],-faWom$scores[,1]))
plot(test[,1],-faWom$scores[,1],main=paste('r=',round(co$est[1,2],4)),
     xlab='W based on M', ylab='own'); abline(0,1)
# This to make nicer:
barplot(t(as.matrix(data.frame(faWom$loadings[,1],-faMan$loadings[,1]))),
        beside = T)
# Invariance for age:
test = (as.matrix(dyoung) %*% faOld$loadings)
co = pcor(data.frame(test[,1],faYoung$scores[,1]))
plot(test[,1],faYoung$scores[,1],main=paste('r=',round(co$est[1,2],4)),
     xlab='Young based on Old', ylab='Young own'); abline(0,1)

######## PLOTS ##########
# Labels to use:
paramlabs = c("meanLnRT.bsl",       "mean.log.RT.GoNoGo",
              "difLnRT.bsl",        "diff.log.RT.GoNoGo",
              "ln.beta.bsl",        "Beta.tr.GoNoGo", 
              "lrnRap.to4th.bsl",   "Appe.lrnRate.tr.GoNoGo",
              "lrnRav.to4th.bsl",   "Aver.lrnRate.tr.GoNoGo",
              "ln.PavBias.bsl",     "Pavl.Bias.tr.GoNoGo",
              "ln.irNoise.bsl",     "LapseRate.tr.GoNoGo",
              "actionBias.bsl",     "Action.Bias.tr.GoNoGo",
              "basGambl.bsl",       "Basic.Gambling.EconRisk",
              "evSens.bsl",         "EV.sens.EconRisk",    
              "skewSens.bsl",       "Skewness.sens.EconRisk",  
              "varSens.bsl",        "Risk.sens.EconRisk",  
              "vaniLnTz.bsl",       "T.uncost.tr.InfoGath",
              "vaniLnCSz.bsl",      "SubjCost.uncost.tr.InfoGath",
              "decrLnTz.bsl",       "T.costed.tr.InfoGath",  
              "decrLnCSz.bsl",      "SubjCost.costed.tr.InfoGath", 
              "f1threats.bsl",      "ApprAvoid.threat.sens",
              "f2losss.bsl",        "ApprAvoid.loss.sens",
              # "predF1.bsl",         "ApproachAvoid.F1",
              # "predF2.bsl",         "ApproachAvoid.F2",  
              # "predF3.bsl",         "ApproachAvoid.F3",   
              "kmz.bsl",            "Intertemp.Disc.DID",
              "sqrtkuz.bsl",        "TasteUncert.tr.DID",
              "xididto7thz.bsl",    "LapseRate.tr.DID",
              "sqrtsrefz.bsl",      "EpiTrust.tr.DID",
              "sqrttsoz.bsl",       "T.choose4other.tr.DID",
              "initinvz.bsl",       "Init.invest.InvTrust",
              "respmagz.bsl",       "Reactiveness.InvTrust",
              "respangz.bsl",       "Coop.Responding.InvTrust",
              "sqrttstlrnr5.bsl",   "lrnRate.tr.TwoStep",
              "logtstbeta5.bsl",    "Beta.tr.TwoStep",
              "tstelig5.bsl",       "Eligibility.TwoStep",
              "sqrttstw5.bsl",      "Modelbasedness.tr.TwoStep",
              "tstprsv5.bsl",       "Persev.TwoStep"
)
paramlabs = t(matrix(paramlabs,2,31)); 


# Plot weights
opar <- par();   # At the end restore by par(opar)
fload1a <- fit4a$loadings[,];  
rownames(fload1a) <-  paramlabs[,2]
sorted.fload1a <- fload1a[order(fload1a[,1]),1];  
dotchart((sorted.fload1a),main=paste('Decison Acuity loadings (4-factor sol.)'),xlab='loadings')
abline(v=0.25, col='gray',lwd=2); abline(v=0); abline(v=-0.25, col='gray',lwd=2)
par(opar)

# plot uniquenesses
opar <- par();   # At the end restore by par(opar)
y = fit4a$uniquenesses; names(y) = paramlabs[,2]
par(las=1) # horizontal labels
par(mar=c(5.1,16,4.1,2.1))  # D, L, U, R margins in cm
barplot(-sort(-y), hor=TRUE, main='4-factor Uniquenesses',xlim=c(0,1))
grid(nx = NULL, ny = NULL,col="gray20",lty = "dotted", lwd = par("lwd"), equilogs=TRUE)
par(opar); 


# Save: Basic training data & analysis, Test data & analyses,
#       
save(Dall.Train, nSall, nocfitall,
     Dall.Test, hdIncr, 
     allcfa1,allcfa2,allcfa3,allcfa4, FAsemAll,
     fa4decAcAll,fa4decAcAl6,
     file=paste(wherearewe$uchjan17on,'decAc4_train-test_working.RData',sep='')
)

# } # end of dummy function DECFA18

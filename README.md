# decAc
Decision Acuity work for Moutoussis - Garzon et al.

Data files:

AllD18.RData has all decision-making measures that were transformed and factor-analyzed.

symfacdeciq.csv includes symptom and disposition factor scores, decision acuity scores, demographic variables and IQ scores needed for linear mixed effects longitudinal analyses.

decAc4_train-test_singleShot.RData contains the single-shot discovery vs. test split of the cognitive data, used to select the factor-analytic model that was was used to derive Decision Acuity scores. The structure of Decision Acuity was then fixed for all subsequent analyses, both LME analyses as above but also, crucially, imaging analyses.

Script files:

CFA-decAc.R Derivation of decision acuity factors, including Exploratory and Confirmatory analyses. Note that it has the ability to re-partition the data into test and train tests, which should not routinely be done.

decAclong.R Longitudinal analyses of decision acuity.

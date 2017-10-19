## step 1: import data and understanding the meanings
Days <- read.csv(file="data.csv",header=TRUE,sep=",")
names(Days)
head(Days)
n <- nrow(Days)
# remove "registerd" and "cnt" responses 
Days$registered <- NULL;Days$cnt <- NULL
# remove "instant", "dteday" predictors
Days$instant <- NULL;Days$dteday <- NULL
# check the class of the predictor variables
apply(Days, 2, class)
# understand the meaning of all the variables ...(takes a while). Now it's done.
# some of the predictors should be converted to categorical ones
Days$season <- as.factor(Days$season)
Days$yr <- as.factor(Days$yr) # since there are only 2 yrs, no need to convert this
Days$mnth <- as.factor(Days$mnth)
Days$holiday <- as.factor(Days$holiday) # similar situation with 'yr' variable
Days$weekday <- as.factor(Days$weekday)
Days$workingday <- as.factor(Days$workingday) # similar situation with 'yr' variable
Days$weathersit <- as.factor(Days$weathersit)

# # log transformation of Y
# Days$casual <- log(Days$casual)

# step 1.0: should first convert categorical predictors to dummy predictors
Days_dum <- model.matrix(casual ~ ., data = Days)
Days_dum <- as.data.frame(Days_dum)
#names(Days_dum)
#head(Days_dum)
Days_dum$`(Intercept)` <- NULL
Days_dum$casual <- Days$casual

## step 2: simple check of variable relations
pairs(Days) # check relations visually by scatter plot matrix
# We cannot do it for the dummy variable table because it's too large

cor(Days_dum[,c(2,22,23,26,27)]) # check correlation coefficients

# Conclusion: 

## step 3: choose the best linear model using automatic model selection methods

# step 3.1: using stepwise method WITH DIFFERENT CRITERIA (F-test & AIC)
lm0 <- lm(casual ~ 1, data= Days_dum)
lmfull <- lm(casual ~., data = Days_dum)
upper <- formula(lmfull)

lm_forward_F <- step(lm0, scope = upper, direction = "both",test="F") # use F-test crit
lm_forward_A <- step(lm0, scope = upper, direction = "both") # use AIC crit
lm_forward_F$call;lm_forward_A$call # compare the two results (they are the same)

lm_backward_F <- step(lmfull, direction = "both",test="F")
lm_backward_A <- step(lmfull, direction = "both")
lm_backward_F$call;lm_backward_A$call # compare the two results (they are the same)

lm_forward_A$call; lm_backward_A$call # compare the two results (they are NOT same)

# Conclusion: with different criteria, best model is the same.
#             with different direction, best model is not the same. 

# step 3.2: using bestsubset method (since # of predictors is small)
library(leaps)
# Previously we know temp and atemp are highly correlated, since both forward and backward stepwise ...
# ... methods included temp but not atemp, we could manually exclude atemp for subsets to stablize the run.
best_sub_regs <- regsubsets(casual ~., data = Days_dum, method = "exhaustive",nvmax = 22, nbest =2)
best_sub_sum <- summary(best_sub_regs)
temp_pr <- data.frame(num_pre = row.names(best_sub_sum$which),
           Cp = best_sub_sum$cp, 
           adjr2 = best_sub_sum$adjr2,
           aic = best_sub_sum$bic+(2-log(n))*(as.numeric(row.names(best_sub_sum$which))+1),
           bic = best_sub_sum$bic,
           best_sub_sum$which)#,
           #row.names = row.names(best_sub_sum$which))
temp_pr[25:44,] # several best models based on different criteria are identified !!!
pool_len <- nrow(temp_pr)

Cp_temp <- temp_pr$Cp
Cp_whichmin1 <- which.min(Cp_temp)
Cp_temp[Cp_whichmin1] <- 1e8
Cp_whichmin2 <- which.min(Cp_temp)

adjr2_temp <- temp_pr$adjr2
adjr2_whichmax1 <- which.max(adjr2_temp)
adjr2_temp[adjr2_whichmax1] <- -1e8
adjr2_whichmax2 <- which.max(adjr2_temp)

aic_temp <- temp_pr$aic
aic_whichmin1 <- which.min(aic_temp)
aic_temp[aic_whichmin1] <- 1e8
aic_whichmin2 <- which.min(aic_temp)

bic_temp <- temp_pr$bic
bic_whichmin1 <- which.min(bic_temp)
bic_temp[bic_whichmin1] <- 1e8
bic_whichmin2 <- which.min(bic_temp)
all_best_ind <- c(Cp_whichmin1,Cp_whichmin2,adjr2_whichmax1,adjr2_whichmax2,
                           aic_whichmin1,aic_whichmin2,bic_whichmin1,bic_whichmin2)
best_model_ind <- unique(all_best_ind)

# save all the best models in best subset
best_subset_models <- list()
all_pred_names <- names(temp_pr[,7:35])
library(pracma)
for (i in 1:length(best_model_ind)) {
  best_ind <- best_model_ind[i]
  best_logical_ind <- as.logical(temp_pr[best_ind,7:35])
  model_names <- all_pred_names[best_logical_ind]
  formu_string <- paste('casual~',paste(model_names,collapse = '+'),sep = '')
  tmp = strcat("model_", num2str(i, 0))
  best_subset_models[[tmp]] <- as.formula(formu_string)
}
# attach the previous two best models from stepwise method
tmp = strcat("model_", num2str(length(best_model_ind)+1, 0))
best_subset_models[[tmp]] <- as.formula(lm_forward_A$call)
tmp = strcat("model_", num2str(length(best_model_ind)+2, 0))
best_subset_models[[tmp]] <- as.formula(lm_backward_A$call)
## step 4: Manually determine the best linear model (with PRESSp and 10-CV)

## step 4.1: using PRESSp
temp_len <- length(best_model_ind)+2
PRESSp <- matrix(0,1,temp_len)
for (i in 1:temp_len) {
  tmp = strcat("model_", num2str(i, 0))
  lmfit <- lm(best_subset_models[[tmp]], data = Days_dum)
  res <- resid(lmfit)
  pr <- res/(1-lm.influence(lmfit)$hat)
  PRESSp[i] <- sum(pr^2)
}
cat(PRESSp)
# result is [93990820 94156987 94100101 94980465 95825395 94620656 94236365]. The smaller the better.

## step 4.2: using 10-fold CV
#### function for creating list of K index sets for K-fold CV ############
CVInd <- function(n,K)	{
  m <- floor(n/K) 	# approximate size of each part
  r <- n-m*K
  I <- sample(n,n) 	# random reordering of the indices
  Ind <- list() 	# will be list of indices for all K parts
  length(Ind) <- K
  for (k in 1:K) {
    if (k <= r) kpart <- ((m+1)*(k-1)+1):((m+1)*k)
    else kpart <- ((m+1)*r+m*(k-r-1)+1):((m+1)*r+m*(k-r))
    Ind[[k]] <- I[kpart]	# indices for kth part of data
  }
  Ind
}
##########################################################################
SEED <- 123
set.seed(SEED)
Nrep <- 200	# number of replicates of CV
K <- 10		# K-fold CV on each replicate
temp_len <- length(best_model_ind)+2
CV_10fold_SSE <- matrix(0,Nrep,temp_len)
y <- Days_dum$casual
for (i in 1:temp_len) {
  tmp = strcat("model_", num2str(i, 0))
  for (j in 1:Nrep) {
    Ind <- CVInd(n,K)
    yhat <- y
    for (k in 1:K) {
      lmfit <- lm(best_subset_models[[tmp]], data = Days_dum[-Ind[[k]],])
      yhat[Ind[[k]]] <- as.numeric(predict(lmfit, Days_dum[Ind[[k]],]))
    }	# end of k loop
    CV_10fold_SSE[j,i]=sum((y-yhat)^2)
  }
  
}
CV_10fold <- apply(CV_10fold_SSE,2,mean)
cat(CV_10fold)
# result is [94287472 94420869 94374413 95232309 96036132 94987072 94560306]. The smaller the better.

# Conclusion: The 1st of the 7 best models is the best, although the 2nd and 3rd and 7th are very close.

## step 5: Check influential observations (hopefully none. If any, then we need to go back to step 3)
# fit the best linear model
lm_best <- lm(best_subset_models$model_1, data = Days_dum)
p <- length(lm_best$coefficients)
inf_table <- influence.measures(lm_best)
inf_mat <- inf_table$infmat
dfbetas_data <- inf_mat[,2:20];dffits_data <- inf_mat[,21];cooks_data <- inf_mat[,23];
dfbeta_max <- max(abs(dfbetas_data)); cat(dfbeta_max, '\t', 1) # 2*sqrt(1/n) is too conservative
dffits_max <- max(abs(dffits_data)); cat(dffits_max, '\t', 2*sqrt(p/n))
pf(max(cooks_data),p,n-p)

# conclusion: different influential detections yield different conclusions. So ...
# ... we should test residuals and perhaps adopt new models before we go back to this step.

## step 6: Check model assumption (normal error assumption) 
# (If residuals are non-normal, then we need to consider higher order terms)

res <- lm_best$residuals
plot(1:n,res, pch = 19, xlab = 'index', ylab = 'residuals')
qqnorm(res) # pretty normal, just with a bit heavy tails
qqline(res) # pretty normal, just with a bit heavy tails

par(mfrow=c(2,2))
plot(Days_dum$temp, res, pch = 18, xlab = 'temp', ylab = 'residuals') # no trend (?)
plot(Days_dum$hum, res, pch = 18, xlab = 'humidity', ylab = 'residuals') # no trend
plot(Days_dum$windspeed, res, pch = 18, xlab = 'windspeed', ylab = 'residuals') # no trend
plot(Days_dum$weekday6*Days_dum$temp, res, pch = 18, xlab = 'weekday6*temp', ylab = 'residuals') # quadratic trend

## step 7: add a nonlinear term to the best model
formu_string <- as.character(best_subset_models$model_1)
formu_string <- paste('casual~',formu_string[3],'+I(weekday6*temp)+I(weekday6*temp^2)',sep = '')
nonl_model <- as.formula(formu_string)
nonl_fit <- lm(nonl_model, data = Days_dum)
# compare the nonlinear model with the best linear model via 10-fold CV
SEED <- 123
set.seed(SEED)
Nrep <- 200	# number of replicates of CV
K <- 10		# K-fold CV on each replicate
CV_10fold_SSE <- matrix(0,Nrep,2)
y <- Days_dum$casual
for (j in 1:Nrep) {
  Ind <- CVInd(n,K)
  yhat1 <- y
  yhat2 <- y
  for (k in 1:K) {
    lmfit1 <- lm(best_subset_models$model_1, data = Days_dum[-Ind[[k]],])
    lmfit2 <- lm(nonl_model, data = Days_dum[-Ind[[k]],])
    yhat1[Ind[[k]]] <- as.numeric(predict(lmfit1, Days_dum[Ind[[k]],]))
    yhat2[Ind[[k]]] <- as.numeric(predict(lmfit2, Days_dum[Ind[[k]],]))
  }	# end of k loop
  CV_10fold_SSE[j,]=c(sum((y-yhat1)^2),sum((y-yhat2)^2))
}
CV_10fold <- apply(CV_10fold_SSE,2,mean)
cat(CV_10fold)
# the result is [94287472 84155325] which means the nonlinear model is appreciably better
# but there are so many nonlinear terms to find and we cannot afford doing it manually.
# all 2nd-order terms = 29*28/2=406 (!); 3rd-order terms = 29*28*27/6 = 3654 (!!)
# we need some automated procedure or try other (advanced) types of models

# also try the updated qq plot
par(mfrow = c(1,1))
res <- nonl_fit$residuals
plot(1:n,res, pch = 19, xlab = 'index', ylab = 'residuals')
qqnorm(res) # pretty normal, just with a bit heavy tails
qqline(res) # pretty normal, just with a bit heavy tails

# Step 8: interpret the model parameters
summary(nonl_fit) # adjr2 is 0.7705
formu_string <- as.character(nonl_model)
formu_string <- paste(formu_string[c(2,1,3)], collapse = "")
nonl_model2 <- paste(formu_string, '-mnth8','-weekday1',sep = '')
nonl_fit2 <- lm(nonl_model2, data = Days_dum)
summary(nonl_fit2) # adjr2 is 0.7693 (smaller)

# Step 9: predict a future event
# see the doc for explanations

# step 10: use the powerful model gbm in comparison with the best model above
library(gbm)
gbm_model <- gbm(casual ~., data = Days_dum, distribution = 'gaussian', 
                 n.trees = 200, interaction.depth = 5, shrinkage = 0.05,
                 train.fraction = 0.8)

# compare the nonlinear model with the best linear model via 10-fold CV
SEED <- 123
set.seed(SEED)
Nrep <- 200	# number of replicates of CV
K <- 10		# K-fold CV on each replicate
CV_10fold_SSE <- matrix(0,Nrep,2)
y <- Days_dum$casual
for (j in 1:Nrep) {
  Ind <- CVInd(n,K)
  yhat1 <- y
  yhat2 <- y
  for (k in 1:K) {
    lmfit1 <- gbm(casual ~., data = Days_dum[-Ind[[k]],], distribution = 'gaussian', 
                  n.trees = 200, interaction.depth = 5, shrinkage = 0.05,
                  train.fraction = 0.8)
    lmfit2 <- lm(nonl_model, data = Days_dum[-Ind[[k]],])
    yhat1[Ind[[k]]] <- as.numeric(predict(lmfit1, Days_dum[Ind[[k]],]))
    yhat2[Ind[[k]]] <- as.numeric(predict(lmfit2, Days_dum[Ind[[k]],]))
  }	# end of k loop
  CV_10fold_SSE[j,]=c(sum((y-yhat1)^2),sum((y-yhat2)^2))
}
CV_10fold <- apply(CV_10fold_SSE,2,mean)
cat(CV_10fold) # result is 57245383 84189808

# compare the predictions between the nonl model and gbm model
gbm_yhat <- predict(gbm_model, Days_dum)
hist(Days_dum$casual-gbm_yhat, pch = 18, 
     main = "Histogram of gbm model residuals",
     xlab = 'Residuals') # histogram of residuals
plot(Days_dum$casual, gbm_yhat) # yhat vs y

nonl_yhat <- nonl_fit$fitted.values
hist(Days_dum$casual-nonl_yhat, pch = 18, 
     main = "Histogram of polynomial model residuals",
     xlab = 'Residuals') # histogram of residuals
plot(Days_dum$casual, nonl_yhat) # yhat vs y
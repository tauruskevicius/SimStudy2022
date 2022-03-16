library(mvtnorm)
library(MatchIt)
library(mice)
#set.seed(1234)

#0 is MCAR
#1 is MAR1
missingness_index <- 0

# Create a dataframe to store estimated treatment effects
resultsData <- data.frame(matrix(ncol = 8, nrow = 0,
                          dimnames = list(NULL, c('Complete Data', 'CCA', 'W_MI5', 'W_MI20', 'W_MI50',
                                                  'A_MI5', 'A_MI20', 'A_MI50'))))
for (j in 1:10){
   n <- 1000
   
   # Generate covariates
   Sigma <- matrix(0.5, ncol = 5, nrow = 5)
   diag(Sigma) <- 5
   mu_X <- rep(10,5)
   
   X_raw <- X <- rmvnorm(n, mu_X, Sigma)
   colnames(X_raw) <- c('X1', 'X2', 'X3', 'X4', 'X5')
   colnames(X) <- c('X1', 'X2', 'X3', 'X4', 'X5')
   
   # Generate treatments
   beta_trt <- rep(-.25,5)
   beta0 <- 10
   xb <- beta0 + X %*% beta_trt
   p <- exp(xb) / (1 + exp(xb))
   trt <- rbinom(n, 1, p)
   boxplot(X[,4]~trt)
   sum(trt)
   
   # Generate responses 
   beta <- rep(1,1,5)
   beta0 <- 0
   mu1 <- beta0 + 5 + X %*% beta
   mu0 <- beta0 + X %*% beta
   sigma <- 1
   
   # Treatment
   y1 <- rnorm(n, mu1, sigma)
   # Control
   y0 <- rnorm(n, mu0, sigma)
   
   # Remove the unobserved treatments
   y <- y0
   y[trt == 1] <- y1[trt == 1]
   
   #Check baseline bias
   #True effect
   delta <- mean(y1) - mean(y0)
   resultsData[j, 'Complete.Data'] = delta     
   
   #Biased effect.  
   mean(y1[trt == 1]) - mean(y0[trt == 0])
   boxplot(y0[trt == 1],y0[trt == 0])
   
   #Add missingness in X1
   #MCAR
   if (missingness_index == 0){
     ind <- sample(1:nrow(X), nrow(X)/2,replace = FALSE)
     X[ind,1] <- NA
   }
   
   #MAR1
    if (missingness_index == 1){
      logi <- -20 + 2*X[,2]
      p <- exp(logi) / (1+ exp(logi))
      R <- rbinom(nrow(X),1,p)
      X[R == 1,1] <- NA
    }
   #MAR2
   #MNAR
   
   # Create a dataframe with missingness
   dat <- data.frame(y, trt, X)
   
   #CCA
   # Remove rows with NAN
   datCCA <- na.omit(dat)
   # Do matching
   m.out <- matchit(trt ~ X1 + X2 + X3 + X4 + X5 , data = datCCA, method = 'nearest',
                    ratio = 1) 
   dataMatched <- match.data(m.out)
   # Fit a model and store estimate treatment effect
   lm1 <- lm(y ~ trt + X1 + X2 + X3 + X4 + X5 , data = dataMatched)
   resultsData[j, 'CCA'] <- summary(lm1)$coef[2,1]                               
   
   # Multiple Imputation - 5 datasets 
   m = 5
   withinResults <- c()
   propensityScores <- data.frame(matrix(ncol=0, nrow = n))
   mids <- mice(dat[,3:7], m = m)
   
   for(k in 1:m) {
      # Access imputed dataset
      dataImputed <- complete(mids, k)
      # Append y and trt
      dataImputed$trt = dat$trt
      dataImputed$y = dat$y
      # Do matching
      m.out <- matchit(trt ~ X1 + X2 + X3 + X4 + X5 , data = dataImputed, method = 'nearest',
                       ratio = 1)
      dataMatched <- match.data(m.out)
      # Append propensity scores of each obs in a dataset
      propensityScores <- cbind(propensityScores, data.frame(m.out$distance))
      # Fit a model and store estimate treatment effect
      lm0 <- lm(y ~ trt + X1 + X2 + X3 + X4 + X5, data = dataMatched)
      withinResults[k] <- summary(lm0)$coef[2,1]
   }
   resultsData[j, 'W_MI5'] <- mean(withinResults)  
   
   # Assign avg propensity scores to m.out and match
   avg_prop <- rowMeans(propensityScores)
   m.out$distance <- avg_prop
   dataMatched <- match.data(m.out)
   # Fit a model and store estimate treatment effect
   lm1 <- lm(y ~ trt + X1 + X2 + X3 + X4 + X5, data = dataMatched)
   resultsData[j, 'A_MI5'] <- summary(lm1)$coefficients[2,1]
   
   # Multiple Imputation - 20 datasets
   m = 20
   withinResults <- c()
   propensityScores <- data.frame(matrix(ncol=0, nrow = n))
   mids <- mice(dat[,3:7], m = m)
   
   for(k in 1:m) {
      # Impute
      dataImputed <- complete(mids, k)
      # Append y and trt
      dataImputed$trt = dat$trt
      dataImputed$y = dat$y
      # Do matching
      m.out <- matchit(trt ~ X1 + X2 + X3 + X4 + X5 , data = dataImputed, method = 'nearest',
                       ratio = 1)
      dataMatched <- match.data(m.out)
      # Append propensity scores of each obs in a dataset
      propensityScores <- cbind(propensityScores, data.frame(m.out$distance))
      # Fit a model and store estimate treatment effect
      lm0 <- lm(y ~ trt + X1 + X2 + X3 + X4 + X5, data = dataMatched)
      withinResults[k] <- summary(lm0)$coef[2,1]
   }
   resultsData[j, 'W_MI20'] <- mean(withinResults)
   
   # Assign avg propensity scores to m.out and match
   avg_prop <- rowMeans(propensityScores)
   m.out$distance <- avg_prop
   dataMatched <- match.data(m.out)
   # Fit a model and store estimate treatment effect
   lm1 <- lm(y ~ trt + X1 + X2 + X3 + X4 + X5, data = dataMatched)
   resultsData[j, 'A_MI20'] <- summary(lm1)$coefficients[2,1]
   
   # Multiple Imputation - 50 datasets
   m = 50
   withinResults <- c()
   propensityScores <- data.frame(matrix(ncol=0, nrow = n))
   mids <- mice(dat[,3:7], m = m)
   
   for(k in 1:m) {
      # Impute
      dataImputed <- complete(mids, k)
      # Append y and trt
      dataImputed$trt = dat$trt
      dataImputed$y = dat$y
      # Do matching
      m.out <- matchit(trt ~ X1 + X2 + X3 + X4 + X5 , data = dataImputed, method = 'nearest',
                       ratio = 1)
      dataMatched <- match.data(m.out)
      # Append propensity scores of each obs in a dataset
      propensityScores <- cbind(propensityScores, data.frame(m.out$distance))
      # Fit a model and store estimate treatment effect
      lm0 <- lm(y ~ trt + X1 + X2 + X3 + X4 + X5, data = dataMatched)
      withinResults[k] <- summary(lm0)$coef[2,1]
   }
   resultsData[j, 'W_MI50'] <- mean(withinResults)
   
   # Assign avg propensity scores to m.out and match
   avg_prop <- rowMeans(propensityScores)
   m.out$distance <- avg_prop
   dataMatched <- match.data(m.out)
   # Fit a model and store estimate treatment effect
   lm1 <- lm(y ~ trt + X1 + X2 + X3 + X4 + X5, data = dataMatched)
   resultsData[j, 'A_MI50'] <- summary(lm1)$coefficients[2,1]
}

results <- data.frame(colMeans(resultsData))
t(results)

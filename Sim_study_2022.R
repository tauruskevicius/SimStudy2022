library(mvtnorm)
library(MatchIt)
library(Zelig)
set.seed(1234)

results0 <- results1 <- c()
#for (j in 1:500){print(j)

n <- 1000

Sigma <- matrix(0.5, ncol = 5, nrow = 5)
diag(Sigma) <- 5
mu_X <- rep(10,5)

# Covariates
X <- rmvnorm(n, mu_X, Sigma)

# Generate the treatments
beta_trt <- rep(-.25,5)
beta0 <- 10
xb <- beta0 + X %*% beta_trt
p <- exp(xb) / (1 + exp(xb))
trt <- rbinom(n, 1, p)
boxplot(X[,4]~trt)
sum(trt)



# Generate the responses 
beta <- rep(1,1,5)
beta0 <- 0
mu1 <- beta0 + 5 + X %*% beta
mu0 <- beta0 + X %*% beta
sigma <- 1
# Treatment
y1 <- rnorm(n, mu1, sigma)
# Control
y0 <- rnorm(n, mu0, sigma)

mean(y1 - y0)

# Remove the unobserved treatments
y <- y0
y[trt == 1] <- y1[trt == 1]


#Check baseline bias
#True effect
delta <- mean(y1) - mean(y0)
#Biased effect.  
mean(y1[trt == 1]) - mean(y0[trt == 0])
boxplot(y0[trt == 1],y0[trt == 0])


dat <- data.frame(y, trt, X)

lm0 <- lm(y ~trt + X1 + X2 , data = dat)
results0[j] <- summary(lm0)$coef[2,1]


# Propensity score matching
m.out <- matchit(trt ~ X1 + X2 , data = dat, method = 'nearest',
                 ratio = 1) # Matching methods: exact, subclass, nearest, optimal, genetic, cem
                            # Other options: NNM wiht or without replacement, NNM with or without caliper
#summary(m.out, standardize = TRUE)
#plot(m.out, type='jitter')
#plot(m.out, type='hist')

# does not seem to be a good match, it is suggested that SMD are below either 0.25 or 0.1.
# Things to try: bigger sample, matching NNM with caliper (smaller matched dataset)
dataMatched <- match.data(m.out)

lm1 <- lm(y ~ trt + X1 + X2 , data = dataMatched)
#treatment est
results1[j] <- summary(lm1)$coef[2,1]




#}


boxplot(results0,results1)
abline(h = 1, col = "red")
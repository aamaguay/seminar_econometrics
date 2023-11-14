library(foreign)
library(AER)
library(caret)
library(MASS)
library(parallel)
library(dplyr)

set.seed(123)


# establish true coefficients
beta0 <- 1
beta1 <- 3
beta2 <- 5
beta3 <- 2
beta4 <- 4
betas.vector <- c(beta0, beta1, beta2, beta3, beta4)
list.true.coefficients <- list()
for (i in 1:length(betas.vector)) list.true.coefficients[[paste("beta", i-1, sep = "")]] <- betas.vector[i]

# E(.) for x's, z's, and \epsilon
mu_eps <- 0
mu_z1 <- 0
mu_z2 <- 0
mu_x1 <- 0
mu_x2 <- 0
mu_x3 <- 0
mu_x4 <- 0
mu_vector <- c(mu_x1, mu_z1, mu_z2, mu_x2, mu_x3, mu_x4, mu_eps)

# covariances between features
cov_x1_z1 <- 0.7
cov_x1_z2 <- 0.1
cov_x1_x2 <- 0.2
cov_x1_x3 <- 0.3
cov_x1_x4 <- 0.2
cov_x1_eps <- 0.6

cov_x2_z1 <- 0.1
cov_x2_z2 <- 0.1
cov_x2_x3 <- 0.2
cov_x2_x4 <- 0.3
cov_x2_eps <- 0

cov_x3_z1 <- 0.1
cov_x3_z2 <- 0.1
cov_x3_x4 <- 0.3
cov_x3_eps <- 0

cov_x4_z1 <- 0.1
cov_x4_z2 <- 0.1
cov_x4_eps <- 0

cov_z1_z2 <- 0.2
cov_z1_eps <- 0
cov_z2_eps <- 0

# var-covariance matrix
var.cor <- matrix(c(1 , cov_x1_z1, cov_x1_z2, cov_x1_x2, cov_x1_x3, cov_x1_x4, cov_x1_eps,
            cov_x1_z1 ,1, cov_z1_z2 ,cov_x2_z1, cov_x3_z1, cov_x4_z1, cov_z1_eps,
            cov_x1_z2, cov_z1_z2, 1, cov_x2_z2, cov_x3_z2, cov_x4_z2, cov_z2_eps,
            cov_x1_x2, cov_x2_z1, cov_x2_z2, 1, cov_x2_x3, cov_x2_x4, cov_x2_eps,
            cov_x1_x3, cov_x3_z1, cov_x3_z2, cov_x2_x3, 1, cov_x3_x4, cov_x3_eps,
            cov_x1_x4, cov_x4_z1, cov_x4_z2, cov_x2_x4, cov_x3_x4, 1, cov_x4_eps,
            cov_x1_eps, cov_z1_eps, cov_z2_eps, cov_x2_eps, cov_x3_eps, cov_x4_eps, 1),
            nrow=7)
rownames(var.cor)<-colnames(var.cor)<-c("x1","z1","z2","x2", "x3", "x4", "eps")
print(var.cor)
# mvrnorm(n = n_obs, mu = mu_vector, Sigma = var.cor)
# eigen(var.cor)$values


# Monte Carlo simulations from IV regr.
estimate_iv_simulation_1 <- function(n_obs, mu_vector, varCovMatrix, ls.true.coefficients,
                                     used.coefficients, instrum.var, endog.var, is.ivmodel){
  
  # prepare data set for the Monte Carlo runs, regressors and errors
  sample = mvrnorm(n = n_obs, mu = mu_vector, Sigma = varCovMatrix)
  full.ds = as.data.frame(sample)
  # compute "y" using the simulated "X's regressors" and errors
  filter.true.coefficients <- list()
  for (i in used.coefficients) filter.true.coefficients[[i]] <- ls.true.coefficients[[i]]
  filter.true.coefficients <- as.vector(unlist(filter.true.coefficients))
  x_vector_columns <- sub("beta", "x", used.coefficients)
  x_vector_columns <- x_vector_columns[x_vector_columns != "x0"]
  
  mx.regressors <- cbind(1, full.ds %>% select(all_of(c(x_vector_columns,"eps") ))) %>% as.matrix()
  y <- mx.regressors %*% matrix(c(filter.true.coefficients,1), 
                                nrow = length(filter.true.coefficients)+1 )
  full.ds$y <- y
  
  # run IV regression
  if (is.ivmodel == TRUE){
    formula.iv.reg <- formula( paste( c( paste("y ~", paste(x_vector_columns, collapse = " + ") ),
                                       paste(c(instrum.var, x_vector_columns[x_vector_columns != endog.var]), 
                                             collapse = " + ") ), collapse = " | " )
                               )
    iv.model <- ivreg(formula = formula.iv.reg, data = full.ds)
    betas.estimated <- iv.model$coefficient
  } else{
    
    formula.ols <- paste("y ~", paste(x_vector_columns, collapse = " + ") )
    ols.model <- lm(formula = formula.ols , data = full.ds)
    betas.estimated <- ols.model$coefficient
  }

  
  # n_obs = 30
  # varCovMatrix = var.cor
  # ls.true.coefficients = list.true.coefficients
  # used.coefficients = c("beta0", "beta1", "beta2")
  # library(broom)
  # library(knitr)
  # first.ols <- lm( x1 ~ z1 + x2 + x3 + x4, data=as.data.frame(full.ds))
  # kable(tidy(first.ols), digits=4, align='c',caption=
  #         "First stage in the 2SLS model for the equation")
  # 
  # firstHat <- fitted(first.ols)
  # scd.2sls <- lm(y~ firstHat + x2+ x3+x4, data=as.data.frame(full.ds)  )
  # kable(tidy(scd.2sls), digits=4, align='c',caption=
  #         "Second stage in the 2SLS model for the equation")
  # 
  # model.iv <- ivreg( y ~ x1 + x2 + x3 + x4|
  #                      x2 + x3 + x4 + z1, data=as.data.frame(full.ds))
  # 
  # kable(tidy(model.iv), digits=4, align='c',caption=
  #         "full model for the IV equation")
  # summary(model.iv)
  # summary(model.iv, vcov = sandwich, df = Inf, diagnostics = TRUE)
  
  # beta1_IVz1 <- ivreg(y ~ x1 + x2 | z1, data = full.ds)$coefficient[2] # IV with instrument being valid
  # beta0_IVz1 <- ivreg(y ~ x1 | z1, data = full.ds)$coefficient[1]
  # beta1_x2 <- ivreg(y ~ z | x2,data=myData)$coefficient[2] # IV with instrument being not valid
  # beta1_x1x2 <- ivreg(y ~ z | x1+x2,data=myData)$coefficient[2] # IV with one of two instruments being
  return(as.vector(betas.estimated))
}


# Setting up the parameters
# numb. of Monte Carlo runs
N <- 500
# numb. of data points
n <- 1000
simulation.resuts.mc.1 <- mclapply(1:N, 
                                   FUN = function(i) estimate_iv_simulation_1(n, 
                                                                              mu_vector, var.cor, 
                                                                              list.true.coefficients,
                                                                              c("beta0","beta1","beta2","beta3","beta4"), 
                                                                              c("z1"), 
                                                                              c("x1")), TRUE)

simulation.resuts.mc.1 <- t(do.call(rbind, simulation.resuts.mc.1))
dim(simulation.resuts.mc.1)

hist(simulation.resuts.mc.1[5,])




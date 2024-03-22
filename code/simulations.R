library(foreign)
library(AER)
library(caret)
library(MASS)
library(parallel)
library(dplyr)
library(corrplot)
library(ggplot2)
library(cowplot)
library(matrixStats)

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

corrplot(var.cor, method = "color", type="upper", number.digits=2,
         tl.offset = 0.5, addCoef.col = "black")


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
  return(as.vector(betas.estimated))
}

# ******************************************************************************
# Setting up the parameters
# numb. of Monte Carlo runs
N <- 100
# numb. of observations
n <- 100

# case 1, ols without instrument
simulation.resuts.mc.olsCase1 <- mclapply(1:N, 
                                   FUN = function(i) estimate_iv_simulation_1(n, 
                                                                              mu_vector, var.cor, 
                                                                              list.true.coefficients,
                                                                              c("beta0","beta1","beta2"), 
                                                                              c(""), 
                                                                              c(""), FALSE ) )
simulation.resuts.mc.olsCase1 <- t(do.call(rbind, simulation.resuts.mc.olsCase1))
mean.estimation.olsCase1 <- as.matrix(rowMeans(simulation.resuts.mc.olsCase1))

# case 1, IV with instrument (z1)
simulation.resuts.mc.ivCase1 <- mclapply(1:N, 
                                         FUN = function(i) estimate_iv_simulation_1(n, 
                                                                              mu_vector, var.cor, 
                                                                              list.true.coefficients,
                                                                              c("beta0","beta1","beta2"), 
                                                                              c("z1"), 
                                                                              c("x1"), TRUE) )
simulation.resuts.mc.ivCase1 <- t(do.call(rbind, simulation.resuts.mc.ivCase1))
mean.estimation.ivCase1 <- as.matrix(rowMeans(simulation.resuts.mc.ivCase1))

# comparative analysis without IV and Z1
results.case1.z1 <- cbind(mean.estimation.olsCase1,mean.estimation.ivCase1,
                          rowSds(simulation.resuts.mc.olsCase1),
                          rowSds(simulation.resuts.mc.ivCase1) )
rownames(results.case1.z1) <- c("beta0", "beta1", "beta2")
dimnames(results.case1.z1)[[2]] <- c("ols_N_100_n_100", "iv_z1_N_100_n_100", "ols_std", "iv_z1_std")

df.ols.z1 <- data.frame(method = c(rep("OLS",300),rep("IV(z1)",300)),
                        value = c(rbind(matrix(t(simulation.resuts.mc.olsCase1),ncol=1),
                                        matrix(t(simulation.resuts.mc.ivCase1),ncol=1) ) ),
                        parameter = rep(c(rep("beta0",100),rep("beta1",100),rep("beta2",100)),2)
)

#plot ols vs iv (z1)
plot_grid(nrow = 1, ncol = 3,
ggplot(df.ols.z1 %>% filter(parameter == "beta0"), aes(x=value, fill=method)) +
  geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) + 
  geom_vline(aes(xintercept = list.true.coefficients$beta0), colour="black", linetype = "dashed") + 
  ggtitle("Distribution of beta0")+theme(plot.title = element_text(hjust = 0.5))
,
ggplot(df.ols.z1 %>% filter(parameter == "beta1"), aes(x=value, fill=method)) +
  geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
  geom_vline(aes(xintercept = list.true.coefficients$beta1), colour="black", linetype = "dashed") +
  ggtitle("Distribution of beta1")+theme(plot.title = element_text(hjust = 0.5))
,
ggplot(df.ols.z1 %>% filter(parameter == "beta2"), aes(x=value, fill=method)) +
  geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
  geom_vline(aes(xintercept = list.true.coefficients$beta2), colour="black", linetype = "dashed") +
  ggtitle("Distribution of beta2")+theme(plot.title = element_text(hjust = 0.5))
)

# case 1, IV with instrument (z2)
simulation.resuts.mc.ivCase1.z2 <- mclapply(1:N, 
                                         FUN = function(i) estimate_iv_simulation_1(n, 
                                                                                    mu_vector, var.cor, 
                                                                                    list.true.coefficients,
                                                                                    c("beta0","beta1","beta2"), 
                                                                                    c("z2"), 
                                                                                    c("x1"), TRUE) )
simulation.resuts.mc.ivCase1.z2 <- t(do.call(rbind, simulation.resuts.mc.ivCase1.z2))
mean.estimation.ivCase1.z2 <- as.matrix(rowMeans(simulation.resuts.mc.ivCase1.z2))

# comparative analysis without IV and Z2
results.case1.z2 <- cbind(mean.estimation.olsCase1, mean.estimation.ivCase1.z2,
                          rowSds(simulation.resuts.mc.olsCase1),
                          rowSds(simulation.resuts.mc.ivCase1.z2) )
rownames(results.case1.z2) <- c("beta0", "beta1","beta2")
dimnames(results.case1.z2)[[2]] <- c("ols_N_100_n_100", "iv_z2_N_100_n_100", "ols_std", "iv_z2_std")


df.ols.z2 <- data.frame(method = c(rep("OLS",300),rep("IV(z2)",300)),
                        value = c(rbind(matrix(t(simulation.resuts.mc.olsCase1),ncol=1),
                                        matrix(t(simulation.resuts.mc.ivCase1.z2),ncol=1) ) ),
                        parameter = rep(c(rep("beta0",100),rep("beta1",100),rep("beta2",100)),2)
)


#plot ols vs iv (z2)
plot_grid(nrow = 1, ncol = 3,
  ggplot(df.ols.z2 %>% filter(parameter == "beta0"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) + 
    geom_vline(aes(xintercept = list.true.coefficients$beta0), colour="black", linetype = "dashed") + 
    ggtitle("Distribution of beta0")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z2 %>% filter(parameter == "beta1"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta1), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta1")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z2 %>% filter(parameter == "beta2"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta2), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta2")+theme(plot.title = element_text(hjust = 0.5))
)

# concat results, Case #1
results.case1 <- cbind(mean.estimation.olsCase1, mean.estimation.ivCase1, mean.estimation.ivCase1.z2)
rownames(results.case1) <- c("beta0", "beta1","beta2")
dimnames(results.case1)[[2]] <- c("ols_N_100_n_100", "iv_z1_N_100_n_100", "iv_z2_N_100_n_100")

# ******************************************************************************
# numb. of Monte Carlo runs
N <- 1000
# numb. of observations
n <- 1000

# case 1, ols without instrument and increase N and n
simulation.resuts.mc.olsCase1.increaseNn <- mclapply(1:N, 
                                          FUN = function(i) estimate_iv_simulation_1(n, 
                                                                                     mu_vector, var.cor, 
                                                                                     list.true.coefficients,
                                                                                     c("beta0","beta1","beta2"), 
                                                                                     c(""), 
                                                                                     c(""), FALSE ) )
simulation.resuts.mc.olsCase1.increaseNn <- t(do.call(rbind, simulation.resuts.mc.olsCase1.increaseNn))
mean.estimation.olsCase1.increaseNn <- as.matrix(rowMeans(simulation.resuts.mc.olsCase1.increaseNn))

# case 1, IV with instrument (z1) and increase N and n
simulation.resuts.mc.ivCase1.increaseNn <- mclapply(1:N, 
                                         FUN = function(i) estimate_iv_simulation_1(n, 
                                                                                    mu_vector, var.cor, 
                                                                                    list.true.coefficients,
                                                                                    c("beta0","beta1","beta2"), 
                                                                                    c("z1"), 
                                                                                    c("x1"), TRUE) )
simulation.resuts.mc.ivCase1.increaseNn <- t(do.call(rbind, simulation.resuts.mc.ivCase1.increaseNn))
mean.estimation.ivCase1.increaseNn <- as.matrix(rowMeans(simulation.resuts.mc.ivCase1.increaseNn))

# comparative analysis without IV and Z1
results.case1.increaseNn.z1 <- cbind(mean.estimation.olsCase1.increaseNn,mean.estimation.ivCase1.increaseNn,
                          rowSds(simulation.resuts.mc.olsCase1.increaseNn),
                          rowSds(simulation.resuts.mc.ivCase1.increaseNn) )
rownames(results.case1.increaseNn.z1) <- c("beta0", "beta1", "beta2")
dimnames(results.case1.increaseNn.z1)[[2]] <- c("ols_N_1000_n_1000", "iv_z1_N_1000_n_1000", "ols_std", "iv_z1_std")

df.ols.z1.increaseNn <- data.frame(method = c(rep("OLS",3000),rep("IV(z1)",3000)),
                        value = c(rbind(matrix(t(simulation.resuts.mc.olsCase1.increaseNn),ncol=1),
                                        matrix(t(simulation.resuts.mc.ivCase1.increaseNn),ncol=1) ) ),
                        parameter = rep(c(rep("beta0",1000),rep("beta1",1000),rep("beta2",1000)),2)
)

#plot ols vs iv (z1), .increaseNn
plot_grid(nrow = 1, ncol = 3,
  ggplot(df.ols.z1.increaseNn %>% filter(parameter == "beta0"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) + 
    geom_vline(aes(xintercept = list.true.coefficients$beta0), colour="black", linetype = "dashed") + 
    ggtitle("Distribution of beta0")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z1.increaseNn %>% filter(parameter == "beta1"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta1), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta1")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z1.increaseNn %>% filter(parameter == "beta2"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta2), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta2")+theme(plot.title = element_text(hjust = 0.5))
)



# case 1, IV with instrument (z2) and increase N and n
simulation.resuts.mc.ivCase1.z2.increaseNn <- mclapply(1:N, 
                                                    FUN = function(i) estimate_iv_simulation_1(n, 
                                                                                               mu_vector, var.cor, 
                                                                                               list.true.coefficients,
                                                                                               c("beta0","beta1","beta2"), 
                                                                                               c("z2"), 
                                                                                               c("x1"), TRUE) )
simulation.resuts.mc.ivCase1.z2.increaseNn <- t(do.call(rbind, simulation.resuts.mc.ivCase1.z2.increaseNn))
mean.estimation.ivCase1.z2.increaseNn <- as.matrix(rowMeans(simulation.resuts.mc.ivCase1.z2.increaseNn))

# comparative analysis without IV and Z2
results.case1.z2.increaseNn <- cbind(mean.estimation.olsCase1.increaseNn,mean.estimation.ivCase1.z2.increaseNn,
                          rowSds(simulation.resuts.mc.olsCase1.increaseNn),
                          rowSds(simulation.resuts.mc.ivCase1.z2.increaseNn) )
rownames(results.case1.z2.increaseNn) <- c("beta0", "beta1","beta2")
dimnames(results.case1.z2.increaseNn)[[2]] <- c("ols_N_1000_n_1000", "iv_z2_N_1000_n_1000", "ols_std", "iv_z2_std")

df.ols.z2.increaseNn <- data.frame(method = c(rep("OLS",3000),rep("IV(z2)",3000)),
                        value = c(rbind(matrix(t(simulation.resuts.mc.olsCase1.increaseNn),ncol=1),
                                        matrix(t(simulation.resuts.mc.ivCase1.z2.increaseNn),ncol=1) ) ),
                        parameter = rep(c(rep("beta0",1000),rep("beta1",1000),rep("beta2",1000)),2)
)


#plot ols vs iv (z2)
plot_grid(ncol = 3, nrow = 1,
  ggplot(df.ols.z2.increaseNn %>% filter(parameter == "beta0"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) + 
    geom_vline(aes(xintercept = list.true.coefficients$beta0), colour="black", linetype = "dashed") + 
    ggtitle("Distribution of beta0")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z2.increaseNn %>% filter(parameter == "beta1"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta1), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta1")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z2.increaseNn %>% filter(parameter == "beta2"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta2), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta2")+theme(plot.title = element_text(hjust = 0.5))
)


# concat results, Case #1
results.case1.increaseNn <- cbind(mean.estimation.olsCase1.increaseNn, mean.estimation.ivCase1.increaseNn, mean.estimation.ivCase1.z2.increaseNn)
rownames(results.case1.increaseNn) <- c("beta0", "beta1")
dimnames(results.case1.increaseNn)[[2]] <- c("ols_N_1000_n_1000", "iv_z2_N_1000_n_1000", "iv_z2_N_1000_n_1000")


# SECOND CASE
# ******************************************************************************
# Setting up the parameters
# numb. of Monte Carlo runs
N <- 100
# numb. of observations
n <- 100

# case 2, ols without instrument, many regressors and small n and N
simulation.resuts.mc.olsCase2 <- mclapply(1:N, 
                                   FUN = function(i) estimate_iv_simulation_1(n, 
                                                                              mu_vector, var.cor, 
                                                                              list.true.coefficients,
                                                                              c("beta0","beta1","beta2","beta3","beta4"), 
                                                                              c(""), 
                                                                              c(""), FALSE) )

simulation.resuts.mc.olsCase2 <- t(do.call(rbind, simulation.resuts.mc.olsCase2))
mean.estimation.olsCase2 <- as.matrix(rowMeans(simulation.resuts.mc.olsCase2))

# case 2, IV with instrument (z1), many regressors and small n and N
simulation.resuts.mc.ivCase2 <- mclapply(1:N, 
                                          FUN = function(i) estimate_iv_simulation_1(n, 
                                                                                     mu_vector, var.cor, 
                                                                                     list.true.coefficients,
                                                                                     c("beta0","beta1","beta2","beta3","beta4"), 
                                                                                     c("z1"), 
                                                                                     c("x1"), TRUE) )
simulation.resuts.mc.ivCase2 <- t(do.call(rbind, simulation.resuts.mc.ivCase2))
mean.estimation.ivCase2 <- as.matrix(rowMeans(simulation.resuts.mc.ivCase2))



# comparative analysis without IV and Z1
results.case2.z1 <- cbind(mean.estimation.olsCase2, mean.estimation.ivCase2,
                          rowSds(simulation.resuts.mc.olsCase2),
                          rowSds(simulation.resuts.mc.ivCase2) )
rownames(results.case2.z1) <- c("beta0", "beta1", "beta2", "beta3", "beta4")
dimnames(results.case2.z1)[[2]] <- c("ols_N_100_n_100", "iv_z1_N_100_n_100", "ols_std", "iv_z1_std")

df.ols.z1.case2 <- data.frame(method = c(rep("OLS",100*5),rep("IV(z1)",100*5)),
                        value = c(rbind(matrix(t(simulation.resuts.mc.olsCase2),ncol=1),
                                        matrix(t(simulation.resuts.mc.ivCase2),ncol=1) ) ),
                        parameter = rep(c(rep("beta0",100), rep("beta1",100), rep("beta2",100),
                                          rep("beta3",100), rep("beta4",100)
                                          ),2)
)

#plot ols vs iv (z1)
plot_grid(
  ggplot(df.ols.z1.case2 %>% filter(parameter == "beta0"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) + 
    geom_vline(aes(xintercept = list.true.coefficients$beta0), colour="black", linetype = "dashed") + 
    ggtitle("Distribution of beta0")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z1.case2 %>% filter(parameter == "beta1"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta1), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta1")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z1.case2 %>% filter(parameter == "beta2"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta2), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta2")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z1.case2 %>% filter(parameter == "beta3"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta3), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta3")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z1.case2 %>% filter(parameter == "beta4"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta4), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta4")+theme(plot.title = element_text(hjust = 0.5))
)



# case 2, IV with instrument (z2), many regressors and small n and N
simulation.resuts.mc.ivCase2.z2 <- mclapply(1:N, 
                                         FUN = function(i) estimate_iv_simulation_1(n, 
                                                                                    mu_vector, var.cor, 
                                                                                    list.true.coefficients,
                                                                                    c("beta0","beta1","beta2","beta3","beta4"), 
                                                                                    c("z2"), 
                                                                                    c("x1"), TRUE) )
simulation.resuts.mc.ivCase2.z2 <- t(do.call(rbind, simulation.resuts.mc.ivCase2.z2))
mean.estimation.ivCase2.z2 <- as.matrix(rowMeans(simulation.resuts.mc.ivCase2.z2))




# comparative analysis without IV and Z2
results.case2.z2 <- cbind(mean.estimation.olsCase2, mean.estimation.ivCase2.z2,
                          rowSds(simulation.resuts.mc.olsCase2),
                          rowSds(simulation.resuts.mc.ivCase2.z2) )
rownames(results.case2.z2) <- c("beta0", "beta1", "beta2", "beta3", "beta4")
dimnames(results.case2.z2)[[2]] <- c("ols_N_100_n_100", "iv_z2_N_100_n_100", "ols_std", "iv_z2_std")


df.ols.z2.case2 <- data.frame(method = c(rep("OLS",100*5),rep("IV(z1)",100*5)),
                              value = c(rbind(matrix(t(simulation.resuts.mc.olsCase2),ncol=1),
                                              matrix(t(simulation.resuts.mc.ivCase2.z2),ncol=1) ) ),
                              parameter = rep(c(rep("beta0",100), rep("beta1",100), rep("beta2",100),
                                                rep("beta3",100), rep("beta4",100)
                              ),2)
)


#plot ols vs iv (z2)
plot_grid(
  ggplot(df.ols.z2.case2 %>% filter(parameter == "beta0"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) + 
    geom_vline(aes(xintercept = list.true.coefficients$beta0), colour="black", linetype = "dashed") + 
    ggtitle("Distribution of beta0")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z2.case2 %>% filter(parameter == "beta1"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta1), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta1")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z2.case2 %>% filter(parameter == "beta2"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta2), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta2")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z2.case2 %>% filter(parameter == "beta3"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta3), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta3")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z2.case2 %>% filter(parameter == "beta4"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta4), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta4")+theme(plot.title = element_text(hjust = 0.5))
)



# concat results, Case #2
results.case2 <- cbind(mean.estimation.olsCase2, mean.estimation.ivCase2, mean.estimation.ivCase2.z2 )
rownames(results.case2) <- c("beta0", "beta1", "beta2", "beta3", "beta4")
dimnames(results.case2)[[2]] <- c("ols_fullvars_N_100_n_100", "iv_z1_fullvars_N_100_n_100", "iv_z2_fullvars_N_100_n_100")

# ******************************************************************************
# Setting up the parameters
# numb. of Monte Carlo runs
N <- 1000
# numb. of observations
n <- 1000

# case 2, ols without instrument, many regressors and large n and N
simulation.resuts.mc.olsCase2.increaseNn <- mclapply(1:N, 
                                          FUN = function(i) estimate_iv_simulation_1(n, 
                                                                                     mu_vector, var.cor, 
                                                                                     list.true.coefficients,
                                                                                     c("beta0","beta1","beta2","beta3","beta4"), 
                                                                                     c(""), 
                                                                                     c(""), FALSE) )

simulation.resuts.mc.olsCase2.increaseNn <- t(do.call(rbind, simulation.resuts.mc.olsCase2.increaseNn))
mean.estimation.olsCase2.increaseNn <- as.matrix(rowMeans(simulation.resuts.mc.olsCase2.increaseNn))

# case 2, IV with instrument (z1), many regressors and large n and N
simulation.resuts.mc.ivCase2.increaseNn <- mclapply(1:N, 
                                                    FUN = function(i) estimate_iv_simulation_1(n, 
                                                                                               mu_vector, var.cor, 
                                                                                               list.true.coefficients,
                                                                                               c("beta0","beta1","beta2","beta3","beta4"), 
                                                                                               c("z1"), 
                                                                                               c("x1"), TRUE) )
simulation.resuts.mc.ivCase2.increaseNn <- t(do.call(rbind, simulation.resuts.mc.ivCase2.increaseNn))
mean.estimation.ivCase2.increaseNn <- as.matrix(rowMeans(simulation.resuts.mc.ivCase2.increaseNn))

# comparative analysis without IV and Z1
results.case2.increaseNn.z1 <- cbind(mean.estimation.olsCase2.increaseNn,mean.estimation.ivCase2.increaseNn,
                                     rowSds(simulation.resuts.mc.olsCase2.increaseNn),
                                     rowSds(simulation.resuts.mc.ivCase2.increaseNn) )
rownames(results.case2.increaseNn.z1) <- c("beta0", "beta1", "beta2", "beta3", "beta4")
dimnames(results.case2.increaseNn.z1)[[2]] <- c("ols_N_1000_n_1000", "iv_z1_N_1000_n_1000", "ols_std", "iv_z1_std")

df.ols.z1.increaseNn.case2 <- data.frame(method = c(rep("OLS",1000*5),rep("IV(z1)",1000*5)),
                                   value = c(rbind(matrix(t(simulation.resuts.mc.olsCase2.increaseNn),ncol=1),
                                                   matrix(t(simulation.resuts.mc.ivCase2.increaseNn),ncol=1) ) ),
                                   parameter = rep(c(rep("beta0",1000), rep("beta1",1000), rep("beta2",1000),
                                                     rep("beta3",1000), rep("beta4",1000)
                                   ),2)
)

#plot ols vs iv (z1), .increaseNn
plot_grid(
  ggplot(df.ols.z1.case2 %>% filter(parameter == "beta0"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) + 
    geom_vline(aes(xintercept = list.true.coefficients$beta0), colour="black", linetype = "dashed") + 
    ggtitle("Distribution of beta0")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z1.case2 %>% filter(parameter == "beta1"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta1), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta1")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z1.case2 %>% filter(parameter == "beta2"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta2), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta2")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z1.case2 %>% filter(parameter == "beta3"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta3), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta3")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z1.case2 %>% filter(parameter == "beta4"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta4), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta4")+theme(plot.title = element_text(hjust = 0.5))
)






# case 2, IV with instrument (z2), many regressors and large n and N
simulation.resuts.mc.ivCase2.z2.increaseNn <- mclapply(1:N, 
                                            FUN = function(i) estimate_iv_simulation_1(n, 
                                                                                       mu_vector, var.cor, 
                                                                                       list.true.coefficients,
                                                                                       c("beta0","beta1","beta2","beta3","beta4"), 
                                                                                       c("z2"), 
                                                                                       c("x1"), TRUE) )
simulation.resuts.mc.ivCase2.z2.increaseNn <- t(do.call(rbind, simulation.resuts.mc.ivCase2.z2.increaseNn))
mean.estimation.ivCase2.z2.increaseNn <- as.matrix(rowMeans(simulation.resuts.mc.ivCase2.z2.increaseNn))


# comparative analysis without IV and Z1
results.case2.increaseNn.z1 <- cbind(mean.estimation.olsCase2.increaseNn,mean.estimation.ivCase2.increaseNn,
                                     rowSds(simulation.resuts.mc.olsCase2.increaseNn),
                                     rowSds(simulation.resuts.mc.ivCase2.increaseNn) )
rownames(results.case2.increaseNn.z1) <- c("beta0", "beta1", "beta2", "beta3", "beta4")
dimnames(results.case2.increaseNn.z1)[[2]] <- c("ols_N_1000_n_1000", "iv_z1_N_1000_n_1000", "ols_std", "iv_z1_std")

df.ols.z1.increaseNn.case2 <- data.frame(method = c(rep("OLS",1000*5),rep("IV(z1)",1000*5)),
                                         value = c(rbind(matrix(t(simulation.resuts.mc.olsCase2.increaseNn),ncol=1),
                                                         matrix(t(simulation.resuts.mc.ivCase2.increaseNn),ncol=1) ) ),
                                         parameter = rep(c(rep("beta0",1000), rep("beta1",1000), rep("beta2",1000),
                                                           rep("beta3",1000), rep("beta4",1000)
                                         ),2)
)

#plot ols vs iv (z1), .increaseNn
plot_grid(
  ggplot(df.ols.z1.case2 %>% filter(parameter == "beta0"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) + 
    geom_vline(aes(xintercept = list.true.coefficients$beta0), colour="black", linetype = "dashed") + 
    ggtitle("Distribution of beta0")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z1.case2 %>% filter(parameter == "beta1"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta1), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta1")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z1.case2 %>% filter(parameter == "beta2"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta2), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta2")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z1.case2 %>% filter(parameter == "beta3"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta3), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta3")+theme(plot.title = element_text(hjust = 0.5))
  ,
  ggplot(df.ols.z1.case2 %>% filter(parameter == "beta4"), aes(x=value, fill=method)) +
    geom_histogram( color='#e9ecef', alpha=0.5, position='identity' ) +
    geom_vline(aes(xintercept = list.true.coefficients$beta4), colour="black", linetype = "dashed") +
    ggtitle("Distribution of beta4")+theme(plot.title = element_text(hjust = 0.5))
)




# concat results, Case #2 large n and N
results.case2.increaseNn <- cbind(mean.estimation.olsCase2.increaseNn, mean.estimation.ivCase2.increaseNn, 
                                  mean.estimation.ivCase2.z2.increaseNn )
rownames(results.case2.increaseNn) <- c("beta0", "beta1", "beta2", "beta3", "beta4")
dimnames(results.case2.increaseNn)[[2]] <- c("ols_fullvars_N_1000_n_1000", "iv_z1_fullvars_N_1000_n_1000", 
                                             "iv_z2_fullvars_N_1000_n_1000")




set.seed(123)
n <- 1000
exp_samples <- rexp(n, rate = 3)

estimate_ecdf <- function(data, x) {
  # Count observations less than or equal to each x
  ecdf_values <- sapply(x, function(val) sum(data <= val) / length(data))
  return(ecdf_values)
}

x_values <- seq(-1, 4, by = 0.01)
ecdf_values_30 <- estimate_ecdf(exp_samples[1:30], x_values)
ecdf_values_60 <- estimate_ecdf(exp_samples[1:60], x_values)
ecdf_values_100 <- estimate_ecdf(exp_samples[1:100], x_values)
ecdf_values_1000 <- estimate_ecdf(exp_samples, x_values)


plot(x_values, pexp(x_values, rate = 3), type = 'l', col = 'blue', lty = 2,
     ylab = 'Probability', xlab = 'x', main = 'Estimated ECDF vs Actual CDF')
lines(x_values, ecdf_values_30, col = 'red')
lines(x_values, ecdf_values_60, col = 'green')
lines(x_values, ecdf_values_100, col = 'purple')
lines(x_values, ecdf_values_1000, col = 'orange')
legend("right", legend = c("Actual CDF", "n=30", "n=60", "n=100", "n=1000"), col = c("blue", "red", "green", "purple", "orange"), lty = c(2, 1, 1, 1, 1))



sup_diff_30 <- max(abs(ecdf_values_30 - pexp(x_values, rate = 3)))
sup_diff_60 <- max(abs(ecdf_values_60 - pexp(x_values, rate = 3)))
sup_diff_100 <- max(abs(ecdf_values_100 - pexp(x_values, rate = 3)))
sup_diff_1000 <- max(abs(ecdf_values_1000 - pexp(x_values, rate = 3)))

print(cat("n=30:", sup_diff_30, "n=60:", sup_diff_60, "n=100:", sup_diff_100, "n=1000:", sup_diff_1000, sep = "  "))

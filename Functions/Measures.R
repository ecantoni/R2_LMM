###################
# Load the packages
library(Hmisc) # for rcorr.cens()
library(lme4)  # for lmer()
library(MASS)  # for mvrnorm() and ginv()
library(rstan) # for stan()
library(cAIC4) # for cAIC()

###########
# Functions

########################
# Snijders & Bosker 1994

# Within-groups variance
WGV <- function(variable, groupvariable)
{
  numtotal <- 0
  denomtotal <- 0
  for (j in 1:J)
  {  nj <- length(variable[groupvariable==j])
  varj <- var(variable[groupvariable==j])
  if(length(variable[groupvariable==j])==1){varj=0}
  numj <- (nj - 1)*varj
  denomj <- nj
  numtotal <- numtotal + numj
  denomtotal <- denomtotal + denomj
  }
  # calculate the within-groups variance
  Vw <- numtotal / (denomtotal - J)
  return(Vw)
}

# Between-groups variance
BGV <- function(variable, groupvariable)
{
  grandmean <- mean(variable)
  
  numtotal <- 0
  for (j in 1:J)
  {
    nj <- length(variable[groupvariable==j])
    meanj <- mean(variable[groupvariable==j])
    numj <- nj * ((meanj - grandmean)^2)
    numtotal <- numtotal + numj
  }
  # calculate the between-groups variance
  Vb <- numtotal / (J - 1)
  Vb <- Vb[[1]]
  return(Vb)
}

# Measures of Snijders & Bosker 1994
SB <- function(model, corr){
  
  # The measures are proposed for balanced data, as the data are unbalanced,
  # I choose a representative value which is the median or 30% trimmed mean
  # of the n_i
  nSB <- 5
  
  # Within- and between-groups variance
  Vw <- WGV(variable=x, groupvariable=county)
  Vb <- BGV(variable=x, groupvariable=county)
  
  # The null model corresponds to m2
  sigma2.y.2 <- summary(m2)$sigma^2
  sigma2.a.2 <- VarCorr(m2)$county[1]
  
  sigma2.y <- summary(model)$sigma^2
  sigma2.a <- VarCorr(model)$county[1]
  sigma.a <- sqrt(sigma2.a)
  
  if(corr==0){
    sigma2.b <- VarCorr(model)$county.1[1]
    sigma.b <- sqrt(sigma2.b)
    rho <- 0 
  }
  
  if(corr!=0){
    sigma2.b <- VarCorr(model)$county[4]
    sigma.b <- sqrt(sigma2.b)
    rho <- attr(VarCorr(model)$county, "correlation")[1,2]
  }
  
  num.1 <- sigma2.a + 2*mean(x)*rho*sigma.a*sigma.b + sigma2.b*((mean(x)^2)+Vb+Vw)
           + sigma2.y
  denom.1 <- sigma2.y.2 + sigma2.a.2
  num.2 <- sigma2.a + 2*mean(x)*rho*sigma.a*sigma.b
           + sigma2.b*((mean(x)^2)+Vb+(Vw/nSB)) + sigma2.y/nSB
  denom.2 <- (sigma2.y.2/nSB) + sigma2.a.2
  
  R2.1 <- 1 - num.1/denom.1
  R2.2 <- 1 - num.2/denom.2
  
  return(c(R2.1, R2.2))
}

####################
# Vonesh et al. 1996

ccc <- function(y, y.hat, adj)
{
  N <- length(y)
  ind <- rep(1, N)
  y.bar <- mean(y)
  y.hat.bar <- mean(y.hat)
  
  rc <- 1 - (t(y-y.hat)%*%(y-y.hat)) /
        ((t(y-y.bar*ind)%*%(y-y.bar*ind)) + 
         (t(y.hat-y.hat.bar*ind)%*%(y.hat-y.hat.bar*ind)) +
          N*((y.bar-y.hat.bar)^2))
  
  k <- N/(N-adj)
  
  rcadj <- 1-k*(1-rc)
  
  return(c(rc, rcadj))
}

# adj = number of fixed effects = p

##########################
# Vonesh & Chinchilli 1996

R2VC <- function(y, y.hat, model0, adj)
{
  N <- length(y)
  
  sigma2.y.0 <- (summary(model0)$sigma)^2
  y.hat.0 <- fitted(model0)
  V <- sigma2.y.0*diag(N)
  
  R2 <- 1 - (t(y-y.hat)%*%solve(V)%*%(y-y.hat))/
       (t(y-y.hat.0)%*%solve(V)%*%(y-y.hat.0))
  
  k <- N/(N-adj)
  
  R2adj <- 1-k*(1-R2)
  
  return(c(R2, R2adj))
}

# adj = number of fixed effects = p

############
# Zheng 2000

Drand <- function(y, y.hat, model0)
{
  D0 <- deviance(model0)
  
  Drand <- 1 - sum((y-y.hat)^2)/D0
  Drand
}

Prand <- function(y, y.hat, model, model0, lm, corr)
{
  D0 <- deviance(model0)
  
  # function lm() is used
  if(lm==1){
    dev <- sum((y-y.hat)^2)
    phi <- summary(model)$sigma^2
    
    Prand <- 1 - (dev/(2*phi)) / (D0/(2*phi))
  }
  
  # function lmer() is used
  if(lm==0){
    
    # random intercept
    if(ncol(ranef(model)$county)==1)
    {
      dev <- sum((y-y.hat)^2)
      phi <- summary(model)$sigma^2
      D <- VarCorr(model)$county[1]
      re <- matrix(ranef(model)$county[,1], ncol=1)
      
      Prand <- 1 - (dev/(2*phi) + (t(re)%*%ginv(D%x%diag(J))%*%re)/2) / (D0/(2*phi))
    }
    
    # random intercept and slope
    if(ncol(ranef(model)$county)==2)
    {  
      # uncorrelated random effects
      if(corr==0){
        D <- matrix(c(VarCorr(model)$county[1], 0, 0, VarCorr(model)$county.1[1]),
                    ncol=2)
      }
      
      # correlated random effects
      if(corr!=0){
        D <- matrix(c(VarCorr(model)$county[1], VarCorr(model)$county[2],
                      VarCorr(model)$county[3], VarCorr(model)$county[4]), ncol=2)
      }
      
      dev <- sum((y-y.hat)^2)      
      phi <- summary(model)$sigma^2
      re <- matrix(c(ranef(model)$county[,1], ranef(model)$county[,2]), ncol=2)
      re.vec <- as.vector(re)
      
      Prand <- 1 - (dev/(2*phi) + (t(re.vec)%*%ginv(D%x%diag(J))%*%re.vec)/2) 
               / (D0/(2*phi))
    }
  }
  
  Prand 
}

#########
# Xu 2003

Xu <- function(y, y.hat, model, model0)
{
  N <- length(y)
  
  sigma2.y.0 <- (summary(model0)$sigma)^2
  
  sigma2.y <- (summary(model)$sigma)^2
  
  r2 <- 1 - (sigma2.y/sigma2.y.0)
  
  RSS <- sum((y-y.hat)^2)
  RSS0 <- sum((y-fitted(model0))^2)
  R2 <- 1 - (RSS/RSS0)
  
  rho2 <- 1 - (sigma2.y/sigma2.y.0)*(exp((RSS/(N*sigma2.y))
               - (RSS0/(N*sigma2.y.0))))
  
  return(c(r2, R2, rho2))
}

#################
# Liu et al. 2008

Liu <- function(y, model, lm, adj, corr, problem){
  
  N <- length(y)
  
  ind <- rep(1, N)
  
  # function lm() is used
  if(lm==1){
    X <- matrix(c(rep(1, N), x), ncol=2, nrow=N)
    
    fe <- as.matrix(coef(model))
    
    # As there is no random effect, all measures are equal:
    R2F <- 1 - t(y-X%*%fe)%*%(y-X%*%fe)/
      t(y-mean(y)*ind)%*%(y-mean(y)*ind)
    R2T <- R2F
    R2TF <- R2F
    
    kF <- N/(N-adj)
    R2Fa <- 1 - kF*(1-R2F)
    
    kTF <- N/(N-qr(X)$rank)
    R2TFa <- 1 - kTF*(1-R2TF)
  }
  
  # function lmer() is used
  if(lm==0){
    
    X <- as.matrix(getME(model,"X"))
    
    # uncorrelated random effects
    if(corr==0){
      Z <- as.matrix(getME(model,"Z"))
    } 
    
    # correlated random effects
    if(corr==1){
      Z <- cbind(t(as.matrix(getME(model,"Ztlist")$`county.(Intercept)`)),
                 t(as.matrix(getME(model,"Ztlist")$county.x)))
    }  
    
    XZ <- cbind(X, Z)
    
    fe <- fixef(model)
    
    # random intercept
    if(ncol(ranef(model)$county)==1){
      re <- matrix(ranef(model)$county[1:J,], ncol=1)
    }
    
    # random intercept and slope
    if(ncol(ranef(model)$county)==2){
      re <- matrix(c(ranef(model)$county[1:J,1], ranef(model)$county[1:J,2]),
            ncol=1)
    }
    
    # the matrix (t(XZ)%*%XZ) is not invertible for the model 3
    # eta cannot be computed
    if(problem==1){
      
      R2F <- 1 - t(y-X%*%fe)%*%(y-X%*%fe)/
        t(y-mean(y)*ind)%*%(y-mean(y)*ind)
      
      R2T <- 1 - t(y-X%*%fe-(Z%*%re))%*%(y-X%*%fe-(Z%*%re))/
        t(y-mean(y)*ind)%*%(y-mean(y)*ind)
      
      R2TF <- NA
      
      kF <- N/(N-adj)
      R2Fa <- 1 - kF*(1-R2F)
      
      kTF <- N/(N-qr(XZ)$rank)
      R2TFa <- 1 - kTF*(1-R2TF) 
    }
    
    if(problem==0){
      eta <- ginv(t(XZ)%*%XZ)%*%t(XZ)%*%y
      
      R2F <- 1 - t(y-X%*%fe)%*%(y-X%*%fe)/
        t(y-mean(y)*ind)%*%(y-mean(y)*ind)
      
      R2T <- 1 - t(y-X%*%fe-(Z%*%re))%*%(y-X%*%fe-(Z%*%re))/
        t(y-mean(y)*ind)%*%(y-mean(y)*ind)
      
      R2TF <- 1 - t(y-XZ%*%eta)%*%(y-XZ%*%eta)/
        t(y-mean(y)*ind)%*%(y-mean(y)*ind)
      
      kF <- N/(N-adj)
      R2Fa <- 1 - kF*(1-R2F)
      
      kTF <- N/(N-qr(XZ)$rank)
      R2TFa <- 1 - kTF*(1-R2TF)
    }
  }
  
  return(c(R2F, R2Fa, R2T, R2TF, R2TFa))
}

# adj = (number of fixed effects) + (number of random effects) = p + q

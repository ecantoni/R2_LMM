#########################################################
#                                                       #
# R code for the analysis of the home radon levels data #
#                                                       #
#########################################################

#######################
# Set up the radon data

srrs2 <- read.table("srrs2.dat", header=T, sep=",")

# To have only the Minnesota state
mn <- srrs2$state=="MN"
radon <- srrs2$activity[mn]

log.radon <- log (ifelse (radon==0, .1, radon))
floor <- srrs2$floor[mn]                     # 0 for basement, 1 for first floor
basement <- c(rep(1,length(floor))) - floor  # 1 for basement, 0 for first floor

n <- length(radon)
y <- log.radon
x <- basement

# Get county index variable
county.name <- as.vector(srrs2$county[mn])

# Remove duplicate elements
uniq <- unique(county.name)
J <- length(uniq)
county <- rep (NA, J)
for (j in 1:J){
  county[county.name==uniq[j]] <- j
}

# Get the county-level predictor
srrs2.fips <- srrs2$stfips*1000 + srrs2$cntyfips
cty <- read.table ("cty.dat", header=T, sep=",")
usa.fips <- 1000*cty[,"stfips"] + cty[,"ctfips"]
usa.rows <- match (unique(srrs2.fips[mn]), usa.fips)
uranium <- cty[usa.rows,"Uppm"]
u <- log (uranium)  # length 85
u.full <- u[county] # length 919

#############################
# Matrix to stock the results

res <- matrix(rep(NA, 43*8), nrow=43, ncol=8)
colnames(res) <- c("Measure", "m0", "m1", "m2", "m3", "m4", "m5", "m6")

############################
# Estimate the models 0 to 6

# Frequentist estimation

m0 <- lm(y ~ 1)

m1 <- lm(y ~ x)

m2ML <- lmer(y ~ 1 + (1|county), REML=FALSE)
m2 <- lmer(y ~ 1 + (1|county))

m3ML <- lmer(y ~ x + (1 + x|county), REML=FALSE)
m3 <- lmer(y ~ x + (1 + x|county))

m4ML <- lmer(y ~ x + u.full + (1 + x|county), REML=FALSE)
m4 <- lmer(y ~ x + u.full + (1 + x|county))

m5ML <- lmer(y ~ x + u.full + x:u.full + (1|county) + (0 + x|county), REML=FALSE)
m5 <- lmer(y ~ x + u.full + x:u.full + (1|county) + (0 + x|county))

m6ML <- lmer(y ~ x + u.full + x:u.full + (1 + x|county), REML=FALSE)
m6 <- lmer(y ~ x + u.full + x:u.full + (1 + x|county))

# Bayesian estimation

radon.data0 <- list(n=n, y=y)
M0stan <- stan(file="M0stan.stan", data=radon.data0, iter=40000, chains=3, thin=30)
e0 <- extract(M0stan)

radon.data1 <- list(n=n, y=y, x=x)
M1stan <- stan(file="M1stan.stan", data=radon.data1, iter=40000, chains=3, thin=30)
e1 <- extract(M1stan)

radon.data2 <- list(n=n, J=J, y=y, county=county)
M2stan <- stan(file="M2stan.stan", data=radon.data2, iter=40000, chains=3, thin=30)
e2 <- extract(M2stan)

radon.data3 <- list(n=n, J=J, y=y, x=x, county=county)
M3stan <- stan(file="M3stan.stan", data=radon.data3, iter=40000, chains=3, thin=30)
e3 <- extract(M3stan)

radon.data4 <- list(n=n, J=J, y=y, x=x, u=u, county=county)
M4stan <- stan(file="M4stan.stan", data=radon.data4, iter=40000, chains=3, thin=30)
e4 <- extract(M4stan)

radon.data5 <- list(n=n, J=J, y=y, x=x, u=u, county=county)
M5stan <- stan(file="M5stan.stan", data=radon.data5, iter=40000, chains=3, thin=30)
e5 <- extract(M5stan)

radon.data6 <- list(n=n, J=J, y=y, x=x, u=u, county=county)
M6stan <- stan(file="M6stan.stan", data=radon.data6, iter=40000, chains=3,
               thin=30)
e6 <- extract(M6stan)

######################
# Compute the measures

######################
# Gelman & Pardoe 2006

# rsquared.y0 and lambda.y0 are always equal to 0

rsquared.y1 <- 1 - mean (apply (e1$e_y, 1, var)) / var (y)
lambda.y1 <- 1 - var (apply (e1$e_y, 2, mean)) / mean (apply (e1$e_y, 1, var))

rsquared.y2 <- 1 - mean (apply (e2$e_y, 1, var)) / var (y)
lambda.y2 <- 1 - var (apply (e2$e_y, 2, mean)) / mean (apply (e2$e_y, 1, var))
rsquared.a2 <- 1 - mean (apply (e2$e_a, 1, var)) / mean (apply (e2$a, 1, var))
lambda.a2 <- 1 - var (apply (e2$e_a, 2, mean)) / mean (apply (e2$e_a, 1, var))

rsquared.y3 <- 1 - mean (apply (e3$e_y, 1, var)) / var (y)
lambda.y3 <- 1 - var (apply (e3$e_y, 2, mean)) / mean (apply (e3$e_y, 1, var))
rsquared.a3 <- 1 - mean (apply (e3$e_a, 1, var)) / mean (apply (e3$a, 1, var))
lambda.a3 <- 1 - var (apply (e3$e_a, 2, mean)) / mean (apply (e3$e_a, 1, var))
rsquared.b3 <- 1 - mean (apply (e3$e_b, 1, var)) / mean (apply (e3$b, 1, var))
lambda.b3 <- 1 - var (apply (e3$e_b, 2, mean)) / mean (apply (e3$e_b, 1, var))

rsquared.y4 <- 1 - mean (apply (e4$e_y, 1, var)) / var (y)
lambda.y4 <- 1 - var (apply (e4$e_y, 2, mean)) / mean (apply (e4$e_y, 1, var))
rsquared.a4 <- 1 - mean (apply (e4$e_a, 1, var)) / mean (apply (e4$a, 1, var))
lambda.a4 <- 1 - var (apply (e4$e_a, 2, mean)) / mean (apply (e4$e_a, 1, var))
rsquared.b4 <- 1 - mean (apply (e4$e_b, 1, var)) / mean (apply (e4$b, 1, var))
lambda.b4 <- 1 - var (apply (e4$e_b, 2, mean)) / mean (apply (e4$e_b, 1, var))

rsquared.y5 <- 1 - mean (apply (e5$e_y, 1, var)) / var (y)
lambda.y5 <- 1 - var (apply (e5$e_y, 2, mean)) / mean (apply (e5$e_y, 1, var))
rsquared.a5 <- 1 - mean (apply (e5$e_a, 1, var)) / mean (apply (e5$a, 1, var))
lambda.a5 <- 1 - var (apply (e5$e_a, 2, mean)) / mean (apply (e5$e_a, 1, var))
rsquared.b5 <- 1 - mean (apply (e5$e_b, 1, var)) / mean (apply (e5$b, 1, var))
lambda.b5 <- 1 - var (apply (e5$e_b, 2, mean)) / mean (apply (e5$e_b, 1, var))

rsquared.y6 <- 1 - mean (apply (e6$e_y, 1, var)) / var (y)
lambda.y6 <- 1 - var (apply (e6$e_y, 2, mean)) / mean (apply (e6$e_y, 1, var))
rsquared.a6 <- 1 - mean (apply (e6$e_a, 1, var)) / mean (apply (e6$a, 1, var))
lambda.a6 <- 1 - var (apply (e6$e_a, 2, mean)) / mean (apply (e6$e_a, 1, var))
rsquared.b6 <- 1 - mean (apply (e6$e_b, 1, var)) / mean (apply (e6$b, 1, var))
lambda.b6 <- 1 - var (apply (e6$e_b, 2, mean)) / mean (apply (e6$e_b, 1, var))

########################
# Snijders & Bosker 1994

SB3.1 <- SB(model=m3, corr=1)[1]
SB3.2 <- SB(model=m3, corr=1)[2]
SB4.1 <- SB(model=m4, corr=1)[1]
SB4.2 <- SB(model=m4, corr=1)[2]
SB5.1 <- SB(model=m5, corr=0)[1]
SB5.2 <- SB(model=m5, corr=0)[2]
SB6.1 <- SB(model=m6, corr=1)[1]
SB6.2 <- SB(model=m6, corr=1)[2]

####################
# Vonesh et al. 1996

# adj = number of fixed effects = p

# Marginal
ccc.marg1 <-     ccc(y=y, y.hat=coef(m1)[1] + coef(m1)[2]*x, adj=2)[1]
ccc.marg1.adj <- ccc(y=y, y.hat=coef(m1)[1] + coef(m1)[2]*x, adj=2)[2]
ccc.marg2 <-     ccc(y=y, y.hat=fixef(m2), adj=1)[1]
ccc.marg2.adj <- ccc(y=y, y.hat=fixef(m2), adj=1)[2]
ccc.marg3 <-     ccc(y=y, y.hat=fixef(m3)[1] + fixef(m3)[2]*x, adj=2)[1]
ccc.marg3.adj <- ccc(y=y, y.hat=fixef(m3)[1] + fixef(m3)[2]*x, adj=2)[2]
ccc.marg4 <-     ccc(y=y, y.hat=fixef(m4)[1] + fixef(m4)[2]*x
                     + fixef(m4)[3]*u.full, adj=3)[1]
ccc.marg4.adj <- ccc(y=y, y.hat=fixef(m4)[1] + fixef(m4)[2]*x
                     + fixef(m4)[3]*u.full, adj=3)[2]
ccc.marg5 <-     ccc(y=y, y.hat=fixef(m5)[1] + fixef(m5)[2]*x
                     + fixef(m5)[3]*u.full + fixef(m5)[4]*x*u.full, adj=4)[1]
ccc.marg5.adj <- ccc(y=y, y.hat=fixef(m5)[1] + fixef(m5)[2]*x
                     + fixef(m5)[3]*u.full + fixef(m5)[4]*x*u.full, adj=4)[2]
ccc.marg6 <-     ccc(y=y, y.hat=fixef(m6)[1] + fixef(m6)[2]*x
                     + fixef(m6)[3]*u.full + fixef(m6)[4]*x*u.full, adj=4)[1]
ccc.marg6.adj <- ccc(y=y, y.hat=fixef(m6)[1] + fixef(m6)[2]*x
                     + fixef(m6)[3]*u.full + fixef(m6)[4]*x*u.full, adj=4)[2]

# Conditional
ccc.cond1 <-     ccc(y=y, y.hat=fitted(m1), adj=2)[1]
ccc.cond1.adj <- ccc(y=y, y.hat=fitted(m1), adj=2)[2]
ccc.cond2 <-     ccc(y=y, y.hat=fitted(m2), adj=1)[1]
ccc.cond2.adj <- ccc(y=y, y.hat=fitted(m2), adj=1)[2]
ccc.cond3 <-     ccc(y=y, y.hat=fitted(m3), adj=2)[1]
ccc.cond3.adj <- ccc(y=y, y.hat=fitted(m3), adj=2)[2]
ccc.cond4 <-     ccc(y=y, y.hat=fitted(m4), adj=3)[1]
ccc.cond4.adj <- ccc(y=y, y.hat=fitted(m4), adj=3)[2]
ccc.cond5 <-     ccc(y=y, y.hat=fitted(m5), adj=4)[1]
ccc.cond5.adj <- ccc(y=y, y.hat=fitted(m5), adj=4)[2]
ccc.cond6 <-     ccc(y=y, y.hat=fitted(m6), adj=4)[1]
ccc.cond6.adj <- ccc(y=y, y.hat=fitted(m6), adj=4)[2]

##########################
# Vonesh & Chinchilli 1996

# adj = number of fixed effects = p

# Marginal
R2.marg1.1 <-     R2VC(y=y, y.hat=coef(m1)[1] + coef(m1)[2]*x, model0=m0, adj=2)[1]
R2.marg1.1.adj <- R2VC(y=y, y.hat=coef(m1)[1] + coef(m1)[2]*x, model0=m0, adj=2)[2]
R2.marg2.1 <-     R2VC(y=y, y.hat=fixef(m2), model0=m0, adj=1)[1]
R2.marg2.1.adj <- R2VC(y=y, y.hat=fixef(m2), model0=m0, adj=1)[2]
R2.marg3.1 <-     R2VC(y=y, y.hat=fixef(m3)[1] + fixef(m3)[2]*x, model0=m0,
                       adj=2)[1]
R2.marg3.1.adj <- R2VC(y=y, y.hat=fixef(m3)[1] + fixef(m3)[2]*x, model0=m0,
                       adj=2)[2]
R2.marg4.1 <-     R2VC(y=y, y.hat=fixef(m4)[1] + fixef(m4)[2]*x
                       + fixef(m4)[3]*u.full, model0=m0,adj=3)[1]
R2.marg4.1.adj <- R2VC(y=y, y.hat=fixef(m4)[1] + fixef(m4)[2]*x
                       + fixef(m4)[3]*u.full, model0=m0, adj=3)[2]
R2.marg5.1 <-     R2VC(y=y, y.hat=fixef(m5)[1] + fixef(m5)[2]*x
                       + fixef(m5)[3]*u.full + fixef(m5)[4]*x*u.full, model0=m0,
                       adj=4)[1]
R2.marg5.1.adj <- R2VC(y=y, y.hat=fixef(m5)[1] + fixef(m5)[2]*x
                       + fixef(m5)[3]*u.full + fixef(m5)[4]*x*u.full, model0=m0,
                       adj=4)[2]
R2.marg6.1 <-     R2VC(y=y, y.hat=fixef(m6)[1] + fixef(m6)[2]*x
                       + fixef(m6)[3]*u.full + fixef(m6)[4]*x*u.full, model0=m0,
                       adj=4)[1]
R2.marg6.1.adj <- R2VC(y=y, y.hat=fixef(m6)[1] + fixef(m6)[2]*x
                       + fixef(m6)[3]*u.full + fixef(m6)[4]*x*u.full, model0=m0,
                       adj=4)[2]

# Conditional
R2.cond1.1 <-     R2VC(y=y, y.hat=fitted(m1), model0=m0, adj=2)[1]
R2.cond1.1.adj <- R2VC(y=y, y.hat=fitted(m1), model0=m0, adj=2)[2]

R2.cond2.1 <-     R2VC(y=y, y.hat=fitted(m2), model0=m0, adj=1)[1]
R2.cond2.1.adj <- R2VC(y=y, y.hat=fitted(m2), model0=m0, adj=1)[2]

R2.cond3.1 <-     R2VC(y=y, y.hat=fitted(m3), model0=m0, adj=2)[1]
R2.cond3.1.adj <- R2VC(y=y, y.hat=fitted(m3), model0=m0, adj=2)[2]
R2.cond3.2 <-     R2VC(y=y, y.hat=fitted(m3), model0=m2, adj=2)[1]
R2.cond3.2.adj <- R2VC(y=y, y.hat=fitted(m3), model0=m2, adj=2)[2]

R2.cond4.1 <-     R2VC(y=y, y.hat=fitted(m4), model0=m0, adj=3)[1]
R2.cond4.1.adj <- R2VC(y=y, y.hat=fitted(m4), model0=m0, adj=3)[2]
R2.cond4.2 <-     R2VC(y=y, y.hat=fitted(m4), model0=m2, adj=3)[1]
R2.cond4.2.adj <- R2VC(y=y, y.hat=fitted(m4), model0=m2, adj=3)[2]

R2.cond5.1 <-     R2VC(y=y, y.hat=fitted(m5), model0=m0, adj=4)[1]
R2.cond5.1.adj <- R2VC(y=y, y.hat=fitted(m5), model0=m0, adj=4)[2]
R2.cond5.2 <-     R2VC(y=y, y.hat=fitted(m5), model0=m2, adj=4)[1]
R2.cond5.2.adj <- R2VC(y=y, y.hat=fitted(m5), model0=m2, adj=4)[2]

R2.cond6.1 <-     R2VC(y=y, y.hat=fitted(m6), model0=m0, adj=4)[1]
R2.cond6.1.adj <- R2VC(y=y, y.hat=fitted(m6), model0=m0, adj=4)[2]
R2.cond6.2 <-     R2VC(y=y, y.hat=fitted(m6), model0=m2, adj=4)[1]
R2.cond6.2.adj <- R2VC(y=y, y.hat=fitted(m6), model0=m2, adj=4)[2]

############
# Zheng 2000

# marginal Drand
Drand1.marg <- Drand(y=y, y.hat=coef(m1)[1] + coef(m1)[2]*x, model0=m0)
Drand2.marg <- Drand(y=y, y.hat=fixef(m2), model0=m0)
Drand3.marg <- Drand(y=y, y.hat=fixef(m3)[1] + fixef(m3)[2]*x, model0=m0)
Drand4.marg <- Drand(y=y, y.hat=fixef(m4)[1] + fixef(m4)[2]*x
                     + fixef(m4)[3]*u.full, model0=m0)
Drand5.marg <- Drand(y=y, y.hat=fixef(m5)[1] + fixef(m5)[2]*x
                     + fixef(m5)[3]*u.full + fixef(m5)[4]*x*u.full, model0=m0)
Drand6.marg <- Drand(y=y, y.hat=fixef(m6)[1] + fixef(m6)[2]*x
                     + fixef(m6)[3]*u.full + fixef(m6)[4]*x*u.full, model0=m0)

# Drand
Drand1 <- Drand(y=y, y.hat=fitted(m1), model0=m0)
Drand2 <- Drand(y=y, y.hat=fitted(m2), model0=m0)
Drand3 <- Drand(y=y, y.hat=fitted(m3), model0=m0)
Drand4 <- Drand(y=y, y.hat=fitted(m4), model0=m0)
Drand5 <- Drand(y=y, y.hat=fitted(m5), model0=m0)
Drand6 <- Drand(y=y, y.hat=fitted(m6), model0=m0)

# marginal Prand
Prand1.marg <- Drand1.marg
Prand2.marg <- Drand2.marg
Prand3.marg <- Drand3.marg
Prand4.marg <- Drand4.marg
Prand5.marg <- Drand5.marg
Prand6.marg <- Drand6.marg

# Prand
Prand1 <- Prand(y=y, y.hat=fitted(m1), model=m1, model0=m0, lm=1, corr=0)
Prand2 <- Prand(y=y, y.hat=fitted(m2), model=m2, model0=m0, lm=0, corr=0)
Prand3 <- Prand(y=y, y.hat=fitted(m3), model=m3, model0=m0, lm=0, corr=1)
Prand4 <- Prand(y=y, y.hat=fitted(m4), model=m4, model0=m0, lm=0, corr=1)
Prand5 <- Prand(y=y, y.hat=fitted(m5), model=m5, model0=m0, lm=0, corr=0)
Prand6 <- Prand(y=y, y.hat=fitted(m6), model=m6, model0=m0, lm=0, corr=1)

# marginal Concordance index c
c1.marg <- rcorr.cens((coef(m1)[1] + coef(m1)[2]*x), y)[1]
c2.marg <- rcorr.cens((fixef(m2)*rep(1,length(y))), y)[1]
c3.marg <- rcorr.cens((fixef(m3)[1] + fixef(m3)[2]*x), y)[1]
c4.marg <- rcorr.cens((fixef(m4)[1] + fixef(m4)[2]*x + fixef(m4)[3]*u.full), y)[1]
c5.marg <- rcorr.cens((fixef(m5)[1] + fixef(m5)[2]*x + fixef(m5)[3]*u.full
                      + fixef(m5)[4]*x*u.full), y)[1]
c6.marg <- rcorr.cens((fixef(m6)[1] + fixef(m6)[2]*x + fixef(m6)[3]*u.full
                      + fixef(m6)[4]*x*u.full), y)[1]

# Concordance index c
c1 <- rcorr.cens(fitted(m1), y)[1]
c2 <- rcorr.cens(fitted(m2), y)[1]
c3 <- rcorr.cens(fitted(m3), y)[1]
c4 <- rcorr.cens(fitted(m4), y)[1]
c5 <- rcorr.cens(fitted(m5), y)[1]
c6 <- rcorr.cens(fitted(m6), y)[1]

#########
# Xu 2003

# Marginal
Xu1.1.R2.marg <-   Xu(y=y, y.hat=coef(m1)[1] + coef(m1)[2]*x, model=m1,
                      model0=m0)[2]
Xu2.1.R2.marg <-   Xu(y=y, y.hat=fixef(m2), model=m2, model0=m0)[2]
Xu3.1.R2.marg <-   Xu(y=y, y.hat=fixef(m3)[1] + fixef(m3)[2]*x, model=m3,
                      model0=m0)[2]
Xu4.1.R2.marg <-   Xu(y=y, y.hat=fixef(m4)[1] + fixef(m4)[2]*x
                      + fixef(m4)[3]*u.full, model=m4, model0=m0)[2]
Xu5.1.R2.marg <-   Xu(y=y, y.hat=fixef(m5)[1] + fixef(m5)[2]*x
                      + fixef(m5)[3]*u.full + fixef(m5)[4]*x*u.full, model=m5,
                      model0=m0)[2]
Xu6.1.R2.marg <-   Xu(y=y, y.hat=fixef(m6)[1] + fixef(m6)[2]*x
                      + fixef(m6)[3]*u.full + fixef(m6)[4]*x*u.full, model=m6,
                      model0=m0)[2]

# Conditional
Xu1.1.r2 <-   Xu(y=y, y.hat=fitted(m1), model=m1, model0=m0)[1]
Xu1.1.R2 <-   Xu(y=y, y.hat=fitted(m1), model=m1, model0=m0)[2]
Xu1.1.rho2 <- Xu(y=y, y.hat=fitted(m1), model=m1, model0=m0)[3]

Xu2.1.r2 <-   Xu(y=y, y.hat=fitted(m2), model=m2, model0=m0)[1]
Xu2.1.R2 <-   Xu(y=y, y.hat=fitted(m2), model=m2, model0=m0)[2]
Xu2.1.rho2 <- Xu(y=y, y.hat=fitted(m2), model=m2, model0=m0)[3]

Xu3.1.r2 <-   Xu(y=y, y.hat=fitted(m3), model=m3, model0=m0)[1]
Xu3.1.R2 <-   Xu(y=y, y.hat=fitted(m3), model=m3, model0=m0)[2]
Xu3.1.rho2 <- Xu(y=y, y.hat=fitted(m3), model=m3, model0=m0)[3]
Xu3.2.r2 <-   Xu(y=y, y.hat=fitted(m3), model=m3, model0=m2)[1]
Xu3.2.R2 <-   Xu(y=y, y.hat=fitted(m3), model=m3, model0=m2)[2]
Xu3.2.rho2 <- Xu(y=y, y.hat=fitted(m3), model=m3, model0=m2)[3]

Xu4.1.r2 <-   Xu(y=y, y.hat=fitted(m4), model=m4, model0=m0)[1]
Xu4.1.R2 <-   Xu(y=y, y.hat=fitted(m4), model=m4, model0=m0)[2]
Xu4.1.rho2 <- Xu(y=y, y.hat=fitted(m4), model=m4, model0=m0)[3]
Xu4.2.r2 <-   Xu(y=y, y.hat=fitted(m4), model=m4, model0=m2)[1]
Xu4.2.R2 <-   Xu(y=y, y.hat=fitted(m4), model=m4, model0=m2)[2]
Xu4.2.rho2 <- Xu(y=y, y.hat=fitted(m4), model=m4, model0=m2)[3]

Xu5.1.r2 <-   Xu(y=y, y.hat=fitted(m5), model=m5, model0=m0)[1]
Xu5.1.R2 <-   Xu(y=y, y.hat=fitted(m5), model=m5, model0=m0)[2]
Xu5.1.rho2 <- Xu(y=y, y.hat=fitted(m5), model=m5, model0=m0)[3]
Xu5.2.r2 <-   Xu(y=y, y.hat=fitted(m5), model=m5, model0=m2)[1]
Xu5.2.R2 <-   Xu(y=y, y.hat=fitted(m5), model=m5, model0=m2)[2]
Xu5.2.rho2 <- Xu(y=y, y.hat=fitted(m5), model=m5, model0=m2)[3]

Xu6.1.r2 <-   Xu(y=y, y.hat=fitted(m6), model=m6, model0=m0)[1]
Xu6.1.R2 <-   Xu(y=y, y.hat=fitted(m6), model=m6, model0=m0)[2]
Xu6.1.rho2 <- Xu(y=y, y.hat=fitted(m6), model=m6, model0=m0)[3]
Xu6.2.r2 <-   Xu(y=y, y.hat=fitted(m6), model=m6, model0=m2)[1]
Xu6.2.R2 <-   Xu(y=y, y.hat=fitted(m6), model=m6, model0=m2)[2]
Xu6.2.rho2 <- Xu(y=y, y.hat=fitted(m6), model=m6, model0=m2)[3]


#################
# Liu et al. 2008

# adj = (number of fixed effects) + (number of random effects) = p + q

Liu1.R2F <-   Liu(y=y, model=m1, lm=1, adj=2+0, corr=0, problem=0)[1]
Liu1.R2Fa <-  Liu(y=y, model=m1, lm=1, adj=2+0, corr=0, problem=0)[2]
Liu1.R2T <-   Liu(y=y, model=m1, lm=1, adj=2+0, corr=0, problem=0)[3]
Liu1.R2TF <-  Liu(y=y, model=m1, lm=1, adj=2+0, corr=0, problem=0)[4]
Liu1.R2TFa <- Liu(y=y, model=m1, lm=1, adj=2+0, corr=0, problem=0)[5]

Liu2.R2F <-   Liu(y=y, model=m2, lm=0, adj=1+1, corr=0, problem=0)[1]
Liu2.R2Fa <-  Liu(y=y, model=m2, lm=0, adj=1+1, corr=0, problem=0)[2]
Liu2.R2T <-   Liu(y=y, model=m2, lm=0, adj=1+1, corr=0, problem=0)[3]
Liu2.R2TF <-  Liu(y=y, model=m2, lm=0, adj=1+1, corr=0, problem=0)[4]
Liu2.R2TFa <- Liu(y=y, model=m2, lm=0, adj=1+1, corr=0, problem=0)[5]

Liu3.R2F <-   Liu(y=y, model=m3, lm=0, adj=2+2, corr=1, problem=1)[1]
Liu3.R2Fa <-  Liu(y=y, model=m3, lm=0, adj=2+2, corr=1, problem=1)[2]
Liu3.R2T <-   Liu(y=y, model=m3, lm=0, adj=2+2, corr=1, problem=1)[3]
Liu3.R2TF <-  Liu(y=y, model=m3, lm=0, adj=2+2, corr=1, problem=1)[4]
Liu3.R2TFa <- Liu(y=y, model=m3, lm=0, adj=2+2, corr=1, problem=1)[5]

Liu4.R2F <-   Liu(y=y, model=m4, lm=0, adj=3+2, corr=1, problem=0)[1]
Liu4.R2Fa <-  Liu(y=y, model=m4, lm=0, adj=3+2, corr=1, problem=0)[2]
Liu4.R2T <-   Liu(y=y, model=m4, lm=0, adj=3+2, corr=1, problem=0)[3]
Liu4.R2TF <-  Liu(y=y, model=m4, lm=0, adj=3+2, corr=1, problem=0)[4]
Liu4.R2TFa <- Liu(y=y, model=m4, lm=0, adj=3+2, corr=1, problem=0)[5]

Liu5.R2F <-   Liu(y=y, model=m5, lm=0, adj=4+2, corr=0, problem=0)[1]
Liu5.R2Fa <-  Liu(y=y, model=m5, lm=0, adj=4+2, corr=0, problem=0)[2]
Liu5.R2T <-   Liu(y=y, model=m5, lm=0, adj=4+2, corr=0, problem=0)[3]
Liu5.R2TF <-  Liu(y=y, model=m5, lm=0, adj=4+2, corr=0, problem=0)[4]
Liu5.R2TFa <- Liu(y=y, model=m5, lm=0, adj=4+2, corr=0, problem=0)[5]

Liu6.R2F <-   Liu(y=y, model=m6, lm=0, adj=4+2, corr=1, problem=0)[1]
Liu6.R2Fa <-  Liu(y=y, model=m6, lm=0, adj=4+2, corr=1, problem=0)[2]
Liu6.R2T <-   Liu(y=y, model=m6, lm=0, adj=4+2, corr=1, problem=0)[3]
Liu6.R2TF <-  Liu(y=y, model=m6, lm=0, adj=4+2, corr=1, problem=0)[4]
Liu6.R2TFa <- Liu(y=y, model=m6, lm=0, adj=4+2, corr=1, problem=0)[5]

#####
# DIC

Dbar0 <- mean(e0$dev)
sigma_ybar0 <- mean(e0$sigma_y)
abar0 <- mean(e0$a)
y_hatbar0 <- abar0
Dhat0 <- -2*sum(dnorm(y , mean=y_hatbar0 , sd=sigma_ybar0, log=TRUE))
pD0 <- Dbar0 - Dhat0
DIC0 <- Dbar0 + pD0

Dbar1 <- mean(e1$dev)
sigma_ybar1 <- mean(e1$sigma_y)
abar1 <- mean(e1$a)
bbar1 <- mean(e1$b)
y_hatbar1 <- abar1 + bbar1*x
Dhat1 <- -2*sum(dnorm(y , mean=y_hatbar1, sd=sigma_ybar1 , log=TRUE))
pD1 <- Dbar1 - Dhat1
DIC1 <- Dbar1 + pD1

Dbar2 <- mean(e2$dev)
sigma_ybar2 <- mean(e2$sigma_y)
abar2 <- apply(e2$a, 2, mean)
y_hatbar2 <- abar2[county]
Dhat2 <- -2*sum(dnorm(y , mean=y_hatbar2, sd=sigma_ybar2, log=TRUE))
pD2 <- Dbar2 - Dhat2
DIC2 <- Dbar2 + pD2

Dbar3 <- mean(e3$dev)
sigma_ybar3 <- mean(e3$sigma_y)
abar3 <- apply(e3$a, 2, mean)
bbar3 <- apply(e3$b, 2, mean)
y_hatbar3 <- abar3[county] + bbar3[county]*x
Dhat3 <- -2*sum( dnorm(y , mean=y_hatbar3, sd=sigma_ybar3, log=TRUE))
pD3 <- Dbar3 - Dhat3
DIC3 <- Dbar3 + pD3

Dbar4 <- mean(e4$dev)
sigma_ybar4 <- mean(e4$sigma_y)
abar4 <- apply(e4$a, 2, mean)
bbar4 <- apply(e4$b, 2, mean)
y_hatbar4 <- abar4[county] + bbar4[county]*x
Dhat4 <- -2*sum( dnorm(y , mean=y_hatbar4, sd=sigma_ybar4, log=TRUE))
pD4 <- Dbar4 - Dhat4
DIC4 <- Dbar4 + pD4

Dbar5 <- mean(e5$dev)
sigma_ybar5 <- mean(e5$sigma_y)
abar5 <- apply(e5$a, 2, mean)
bbar5 <- apply(e5$b, 2, mean)
y_hatbar5 <- abar5[county] + bbar5[county]*x
Dhat5 <- -2*sum( dnorm(y , mean=y_hatbar5, sd=sigma_ybar5, log=TRUE))
pD5 <- Dbar5 - Dhat5
DIC5 <- Dbar5 + pD5

Dbar6 <- mean(e6$dev)
sigma_ybar6 <- mean(e6$sigma_y)
abar6 <- apply(e6$a, 2, mean)
bbar6 <- apply(e6$b, 2, mean)
y_hatbar6 <- abar6[county] + bbar6[county]*x
Dhat6 <- -2*sum( dnorm(y , mean=y_hatbar6, sd=sigma_ybar6, log=TRUE))
pD6 <- Dbar6 - Dhat6
DIC6 <- Dbar6 + pD6

######
# -2LL

LL0 <- -2*logLik(m0)[1]
LL1 <- -2*logLik(m1)[1]

LL2ML <- -2*logLik(m2ML)[1]
LL3ML <- -2*logLik(m3ML)[1]
LL4ML <- -2*logLik(m4ML)[1]
LL5ML <- -2*logLik(m5ML)[1]
LL6ML <- -2*logLik(m6ML)[1]

LL2REML <- as.numeric(-2*logLik(m2)[1])
LL3REML <- as.numeric(-2*logLik(m3)[1])
LL4REML <- as.numeric(-2*logLik(m4)[1])
LL5REML <- as.numeric(-2*logLik(m5)[1])
LL6REML <- as.numeric(-2*logLik(m6)[1])

#####
# AIC

# extractAIC() gives negative values with lm()
AIC0 <- AIC(m0)
AIC1 <- AIC(m1)
# for REML estimation, AIC() gives the AIC value based on -2LL REML and
# extractAIC() gives the AIC value based on -2LL ML -> prefer extractAIC()
AIC2 <- extractAIC(m2ML)[2]
AIC3 <- extractAIC(m3ML)[2]
AIC4 <- extractAIC(m4ML)[2]
AIC5 <- extractAIC(m5ML)[2]
AIC6 <- extractAIC(m6ML)[2]

#####
# BIC

BIC0 <- AIC(m0, k=log(n))
BIC1 <- AIC(m1, k=log(n))
BIC2 <- extractAIC(m2ML, k=log(n))[2]
BIC3 <- extractAIC(m3ML, k=log(n))[2]
BIC4 <- extractAIC(m4ML, k=log(n))[2]
BIC5 <- extractAIC(m5ML, k=log(n))[2]
BIC6 <- extractAIC(m6ML, k=log(n))[2]

######
# cAIC

# ML
cAIC2ML <- cAIC(m2ML)$caic
cAIC3ML <- cAIC(m3ML)$caic
cAIC4ML <- cAIC(m4ML)$caic
cAIC5ML <- cAIC(m5ML)$caic
cAIC6ML <- cAIC(m6ML)$caic

# REML
cAIC2REML <- cAIC(m2)$caic
cAIC3REML <- cAIC(m3)$caic
cAIC4REML <- cAIC(m4)$caic
cAIC5REML <- cAIC(m5)$caic
cAIC6REML <- cAIC(m6)$caic

##########################
# Stock the results in res

res[1,] <- c("rsquared.y", NA, rsquared.y1, rsquared.y2, rsquared.y3,
             rsquared.y4, rsquared.y5, rsquared.y6)
res[2,] <- c("lambda.y", NA, lambda.y1, lambda.y2, lambda.y3,
             lambda.y4, lambda.y5, lambda.y6)
res[3,] <- c("rsquared.a", NA, NA, rsquared.a2, rsquared.a3,
             rsquared.a4, rsquared.a5, rsquared.a6)
res[4,] <- c("lambda.a", NA, NA, lambda.a2, lambda.a3,
             lambda.a4, lambda.a5, lambda.a6)
res[5,] <- c("rsquared.b", NA, NA, NA, rsquared.b3,
             rsquared.b4, rsquared.b5, rsquared.b6)
res[6,] <- c("lambda.b", NA, NA, NA, lambda.b3,
             lambda.b4, lambda.b5, lambda.b6)

res[7,] <- c("SB.1", NA, NA, NA, SB3.1, SB4.1, SB5.1, SB6.1)
res[8,] <- c("SB.2", NA, NA, NA, SB3.2, SB4.2, SB5.2, SB6.2)

res[9,] <- c("ccc.marg", NA, ccc.marg1, ccc.marg2, ccc.marg3,
             ccc.marg4, ccc.marg5, ccc.marg6)
res[10,] <- c("ccc.marg.adj", NA, ccc.marg1.adj, ccc.marg2.adj, ccc.marg3.adj,
              ccc.marg4.adj, ccc.marg5.adj, ccc.marg6.adj)
res[11,] <- c("ccc.cond", NA, ccc.cond1, ccc.cond2, ccc.cond3,
              ccc.cond4, ccc.cond5, ccc.cond6)
res[12,] <- c("ccc.cond.adj", NA, ccc.cond1.adj, ccc.cond2.adj, ccc.cond3.adj,
              ccc.cond4.adj, ccc.cond5.adj, ccc.cond6.adj)

res[13,] <- c("R2.marg", NA, R2.marg1.1, R2.marg2.1, R2.marg3.1,
              R2.marg4.1, R2.marg5.1, R2.marg6.1)
res[14,] <- c("R2.marg.adj", NA, R2.marg1.1.adj, R2.marg2.1.adj, R2.marg3.1.adj,
              R2.marg4.1.adj, R2.marg5.1.adj, R2.marg6.1.adj)
res[15,] <- c("R2.cond.1", NA, R2.cond1.1, R2.cond2.1, R2.cond3.1,
              R2.cond4.1, R2.cond5.1, R2.cond6.1)
res[16,] <- c("R2.cond.1.adj", NA, R2.cond1.1.adj, R2.cond2.1.adj, R2.cond3.1.adj,
              R2.cond4.1.adj, R2.cond5.1.adj, R2.cond6.1.adj)
res[17,] <- c("R2.cond.2", NA, NA, NA, R2.cond3.2, R2.cond4.2,
              R2.cond5.2, R2.cond6.2)
res[18,] <- c("R2.cond.2.adj", NA, NA, NA, R2.cond3.2.adj, R2.cond4.2.adj,
              R2.cond5.2.adj, R2.cond6.2.adj)

res[19,] <- c("Drand.marg", NA, Drand1.marg, Drand2.marg, Drand3.marg,
              Drand4.marg, Drand5.marg, Drand6.marg)
res[20,] <- c("Prand.marg", NA, Prand1.marg, Prand2.marg, Prand3.marg,
              Prand4.marg, Prand5.marg, Prand6.marg)
res[21,] <- c("c.marg", NA, c1.marg, c2.marg, c3.marg, c4.marg, c5.marg, c6.marg)
res[22,] <- c("Drand", NA, Drand1, Drand2, Drand3, Drand4, Drand5, Drand6)
res[23,] <- c("Prand", NA, Prand1, Prand2, Prand3, Prand4, Prand5, Prand6)
res[24,] <- c("c", NA, c1, c2, c3, c4, c5, c6)

res[25,] <- c("Xu.1.R2.marg", NA, Xu1.1.R2.marg, Xu2.1.R2.marg, Xu3.1.R2.marg,
              Xu4.1.R2.marg, Xu5.1.R2.marg, Xu6.1.R2.marg)
res[26,] <- c("Xu.1.r2", NA, Xu1.1.r2, Xu2.1.r2, Xu3.1.r2,
              Xu4.1.r2, Xu5.1.r2, Xu6.1.r2)
res[27,] <- c("Xu.1.R2", NA, Xu1.1.R2, Xu2.1.R2, Xu3.1.R2,
              Xu4.1.R2, Xu5.1.R2, Xu6.1.R2)
res[28,] <- c("Xu.1.rho2", NA, Xu1.1.rho2, Xu2.1.rho2, Xu3.1.rho2,
              Xu4.1.rho2, Xu5.1.rho2, Xu6.1.rho2)
res[29,] <- c("Xu.2.r2", NA, NA, NA, Xu3.2.r2, Xu4.2.r2, Xu5.2.r2, Xu6.2.r2)
res[30,] <- c("Xu.2.R2", NA, NA, NA, Xu3.2.R2, Xu4.2.R2, Xu5.2.R2, Xu6.2.R2)
res[31,] <- c("Xu.2.rho2", NA, NA, NA, Xu3.2.rho2, Xu4.2.rho2,
              Xu5.2.rho2, Xu6.2.rho2)

res[32,] <- c("Liu.R2F", NA, Liu1.R2F, Liu2.R2F, Liu3.R2F,
              Liu4.R2F, Liu5.R2F, Liu6.R2F)
res[33,] <- c("Liu.R2Fa", NA, Liu1.R2Fa, Liu2.R2Fa, Liu3.R2Fa,
              Liu4.R2Fa, Liu5.R2Fa, Liu6.R2Fa)
res[34,] <- c("Liu.R2T", NA, Liu1.R2T, Liu2.R2T, Liu3.R2T,
              Liu4.R2T, Liu5.R2T, Liu6.R2T)
res[35,] <- c("Liu.R2TF", NA, Liu1.R2TF, Liu2.R2TF, Liu3.R2TF,
              Liu4.R2TF, Liu5.R2TF, Liu6.R2TF)
res[36,] <- c("Liu.R2TFa", NA, Liu1.R2TFa, Liu2.R2TFa, Liu3.R2TFa,
              Liu4.R2TFa, Liu5.R2TFa, Liu6.R2TFa)

res[37,] <- c("DIC", DIC0, DIC1, DIC2, DIC3, DIC4, DIC5, DIC6)
res[38,] <- c("LLML", LL0, LL1, LL2ML, LL3ML, LL4ML, LL5ML, LL6ML)
res[39,] <- c("LLREML", NA, NA, LL2REML, LL3REML, LL4REML, LL5REML, LL6REML)
res[40,] <- c("AIC", AIC0, AIC1, AIC2, AIC3, AIC4, AIC5, AIC6)
res[41,] <- c("BIC", BIC0, BIC1, BIC2, BIC3, BIC4, BIC5, BIC6)
res[42,] <- c("cAICML", NA, NA, cAIC2ML, cAIC3ML, cAIC4ML, cAIC5ML, cAIC6ML)
res[43,] <- c("cAICREML", NA, NA, cAIC2REML, cAIC3REML, cAIC4REML,
              cAIC5REML, cAIC6REML)

write.table(res, file="radon_example.csv", sep=";", row.names=FALSE)

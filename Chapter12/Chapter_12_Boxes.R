# R Script for Growth Models and Inference Chapter - July 4, 2017

## A logical to identify if TIFFs of the figures should be made
##   Set to TRUE to make the TIFFs
##   Note any line with makeTIFFFigs is not shown in boxes
makeTIFFFigs <- FALSE

## Figure options
lwd <- 2
figdim <- 5
res <- 600
ptsz <- 12

par(mar=c(3.05,3.05,0.65,0.65),mgp=c(1.9,0.3,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)

#
#
# #####################################################################
# =====================================================================
# ---------------------------------------------------------------------
#
# Box 12.1 (must be run for all other boxes)
#
# ---------------------------------------------------------------------
# =====================================================================

library(FSA)        # for filterD, residPlot, confint, predict (Ogle 2017b)
library(FSAdata)    # for WalleyeErie2, BluegillIL, WalleyeML (Ogle 2017c)
library(MASS)       # for confint, profile (Venables and Ripley 2002)
library(nlstools)   # for nlsBoot (Baty et al. 2015)
library(minpack.lm) # for nlsLM (Elzhov 2015)
library(cvTools)    # for cvFit (Alfons 2012)
library(AICcmodavg) # for aictab (Mazerolle 2016)
library(nlme)       # for nlme and related (Pinheiro et al. 2016)
library(lattice)    # for xyplot (Sarkar 2008)
library(msm)        # for deltamethod (Jackson 2011)
library(R2jags)     # for jags and related (Su and Yajima 2015)
# R2jags loads coda package for traceplot & densplot (Plummer et al. 2006)

options(show.signif.stars=FALSE)

data(WalleyeErie2)
str(WalleyeErie2)
head(WalleyeErie2,n=3)

data(BluegillIL)
str(BluegillIL)
head(BluegillIL,n=3)

data(WalleyeML)
str(WalleyeML)
head(WalleyeML,n=3)

#
#
# #####################################################################
# =====================================================================
# ---------------------------------------------------------------------
#
# Box 12.2 (requires running Box 12.1 script)
#
# ---------------------------------------------------------------------
# =====================================================================

wf14T <- filterD(WalleyeErie2,sex=="female",year==2014,loc==1)

vbT <- function(T,Linf,K=NULL,t0=NULL) {
  if (length(Linf)==3) {
    t0 <- Linf[[3]]
    K <- Linf[[2]]
    Linf <- Linf[[1]]
  }
  Linf*(1-exp(-K*(T-t0)))
}

vbT(3,Linf=650,K=0.3,t0=-2)  # one age & separate parameter arguments
vbT(2:4,c(650,0.3,-2))       # three ages & parameters in one vector

plot(tl~age,data=wf14T)

if (makeTIFFFigs) tiff("Figs/Box12_2_StartingValues.tif",
                       width=figdim,height=figdim,
                       units="in",pointsize=ptsz,family="sans",res=res)
par(mar=c(3.05,3.05,0.65,0.65),mgp=c(1.9,0.3,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
plot(tl~age,data=wf14T)
curve(vbT(x,650,0.3,0),from=0,to=11,add=TRUE,lty=2) # first guess
curve(vbT(x,650,0.3,-2),from=0,to=11,add=TRUE)      # second guess
if (makeTIFFFigs) dev.off()

svV <- list(Linf=650,K=0.3,t0=-2)

fitV <- nls(tl~vbT(age,Linf,K,t0),data=wf14T,start=svV)

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
cfT <- coef(fitV)
sumT <- summary(fitV)
#######################################################################

coef(fitV)

summary(fitV,correlation=TRUE)

deviance(fitV) # minimum RSS
logLik(fitV)   # maximum log-likelihood

#
#
# #####################################################################
# =====================================================================
# ---------------------------------------------------------------------
#
# Box 12.3 (requires running Boxes 12.1 and 12.2 scripts)
#
# ---------------------------------------------------------------------
# =====================================================================

if (makeTIFFFigs) tiff("Figs/Box12_3_A_ResidualPlot.tif",
                       width=2*figdim,height=figdim,
                       units="in",pointsize=ptsz,family="sans",res=res)
par(mar=c(3.05,3.05,0.65,0.65),mgp=c(1.9,0.3,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
residPlot(fitV,resid.type="standardized",loess=FALSE)
if (makeTIFFFigs) dev.off()

if (makeTIFFFigs) tiff("Figs/Box12_3_B_FittedPlot.tif",
                       width=figdim,height=figdim,
                       units="in",pointsize=ptsz,family="sans",res=res)
par(mar=c(3.05,3.05,0.65,0.65),mgp=c(1.9,0.3,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
plot(tl~age,data=wf14T,xlab="Age (years)",ylab="Total Length (mm)")
curve(vbT(x,coef(fitV)),from=0,to=11,lwd=2,add=TRUE)
if (makeTIFFFigs) dev.off()

svV2 <- list(Linf=500,K=0.3,t0=0)
fitV2 <- nls(tl~vbT(age,Linf,K,t0),data=wf14T,start=svV2)
svV3 <- list(Linf=800,K=0.2,t0=1)
fitV3 <- nls(tl~vbT(age,Linf,K,t0),data=wf14T,start=svV3)
fitV4 <- nls(tl~vbT(age,Linf,K,t0),data=wf14T,start=svV,algorithm="port")
fitV5 <- nlsLM(tl~vbT(age,Linf,K,t0),data=wf14T,start=svV)
rbind(coef(fitV),coef(fitV2),coef(fitV3),coef(fitV4),coef(fitV5))

#
#
# #####################################################################
# =====================================================================
# ---------------------------------------------------------------------
#
# Box 12.4 (requires running Boxes 12.1 and 12.2 scripts)
#
# ---------------------------------------------------------------------
# =====================================================================

#######################################################################
# This code sets the random seed for analyses in this box. This is not
#   required for these analyses but will cause the same results to be
#   returned as in the printed chapter if this script is sourced.
set.seed(14354454)
#######################################################################

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
plh.conf <- confint(fitV)
#######################################################################

confint(fitV)

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
plh.conf
#######################################################################

bootT <- nlsBoot(fitV)
head(bootT$coefboot,n=3)

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
boot.conf <- confint(bootT)
#######################################################################

confint(bootT)

vbT(4,coef(fitV))

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
p4T <- predict(bootT,vbT,T=4)
#######################################################################

predict(bootT,vbT,T=4)

predict(bootT,vbT,T=3:7)

#
#
# #####################################################################
# =====================================================================
# ---------------------------------------------------------------------
#
# Box 12.5 (requires running Boxes 12.1 and 12.2 scripts)
#
# ---------------------------------------------------------------------
# =====================================================================

cvFit(fitV,data=wf14T,y=wf14T$tl,K=nrow(wf14T))

predL <- vbT(wf14T$age,coef(fitV))
tmp <- lm(wf14T$tl~predL)
summary(tmp)$r.squared

#
#
# #####################################################################
# =====================================================================
# ---------------------------------------------------------------------
#
# Box 12.6 (requires running Box 12.1 script)
#
# ---------------------------------------------------------------------
# =====================================================================

#######################################################################
# This code sets the random seed for analyses in this box. This is not
#   required for these analyses but will cause the same results to be
#   returned as in the printed chapter if this script is sourced.
set.seed(14354454)
#######################################################################

vbW <- function(Lm,dt,Linf,K=NULL,beta=NULL) {
  if (length(Linf)==3) {
    beta <- Linf[[3]]
    K <- Linf[[2]]
    Linf <- Linf[[1]]
  }
  (Linf+beta*(Lm-mean(Lm))-Lm)*(1-exp(-K*dt))
}

svW <- list(Linf=250,K=0.3,beta=0)
fitW <- nls(deltaLen~vbW(lenMark,deltaTime,Linf,K,beta),
            data=BluegillIL,start=svW)

if (makeTIFFFigs) tiff("Figs/Box12_6_A_ResidualPlot.tif",
                       width=2*figdim,height=figdim,
                       units="in",pointsize=ptsz,family="sans",res=res)
par(mar=c(3.05,3.05,0.65,0.65),mgp=c(1.9,0.3,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
residPlot(fitW,resid.type="standardized",loess=FALSE)
if (makeTIFFFigs) dev.off()

fitW2 <- nls(deltaLen~vbW(lenMark,deltaTime,Linf,K,beta),data=BluegillIL,
             start=list(Linf=400,K=0.5,beta=0.5))
fitW3 <- nls(deltaLen~vbW(lenMark,deltaTime,Linf,K,beta),data=BluegillIL,
             start=list(Linf=200,K=0.2,beta=0.5))
fitW4 <- nls(deltaLen~vbW(lenMark,deltaTime,Linf,K,beta),data=BluegillIL,
             start=svW,algorithm="port")
fitW5 <- nlsLM(deltaLen~vbW(lenMark,deltaTime,Linf,K,beta),data=BluegillIL,
               start=svW)
rbind(coef(fitW),coef(fitW2),coef(fitW3),coef(fitW4),coef(fitW5))

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
bootW <- nlsBoot(fitW)
conf.boot <- confint(bootW)
#######################################################################

coef(fitW)

bootW <- nlsBoot(fitW)
confint(bootW)

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
conf.boot
#######################################################################

if (makeTIFFFigs) tiff("Figs/Box12_6_B_FittedPlot.tif",
                       width=figdim,height=figdim,
                       units="in",pointsize=ptsz,family="sans",res=res)
par(mar=c(3.05,3.05,0.65,0.65),mgp=c(1.9,0.3,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
plot(NA,xlim=c(0,3),ylim=c(0,85),xlab="Change in time (years)",
     ylab="Change in total length (mm)",xaxt="n")
axis(side=1,at=0:3)
tmp.dts <- seq(0,3,0.1)
for (i in seq(150,250,25)) lines(tmp.dts,vbW(i,tmp.dts,coef(fitW)),lwd=2)
text(2.2,82,"Lm=150 mm",cex=0.9)
text(2.3,3,"Lm=250 mm",cex=0.9)
if (makeTIFFFigs) dev.off()

#
#
# #####################################################################
# =====================================================================
# ---------------------------------------------------------------------
#
# Box 12.7 (requires running Boxes 12.1 and 12.2 scripts)
#
# ---------------------------------------------------------------------
# =====================================================================

w14T1 <- filterD(WalleyeErie2,year==2014,loc==1)

vbLKT <- tl~Linf[sex]*(1-exp(-K[sex]*(age-t0[sex])))
vbLK <- tl~Linf[sex]*(1-exp(-K[sex]*(age-t0)))
vbLT <- tl~Linf[sex]*(1-exp(-K*(age-t0[sex])))
vbKT <- tl~Linf*(1-exp(-K[sex]*(age-t0[sex])))
vbL <- tl~Linf[sex]*(1-exp(-K*(age-t0)))
vbK <- tl~Linf*(1-exp(-K[sex]*(age-t0)))
vbT <- tl~Linf*(1-exp(-K*(age-t0[sex])))
vbO <- tl~Linf*(1-exp(-K*(age-t0)))

# Starting values if all parameters differ
svLKT <- list(Linf=c(550,625),K=c(0.3,0.3),t0=c(0,0))
# Starting values if two parameters differ
svLK <- list(Linf=c(550,625),K=c(0.3,0.3),t0=0)
svLT <- list(Linf=c(550,625),K=0.3,t0=c(0,0))
svKT <- list(Linf=590,K=c(0.3,0.3),t0=c(0,0))
# Starting values if one parameter differs
svL <- list(Linf=c(550,625),K=0.3,t0=0)
svK <- list(Linf=590,K=c(0.3,0.3),t0=0)
svT <- list(Linf=590,K=0.3,t0=c(0,0))
# Starting values if no parameters differ
svO <- list(Linf=590,K=0.3,t0=0)

fitLKT <- nls(vbLKT,data=w14T1,start=svLKT)

if (makeTIFFFigs) tiff("Figs/Box12_7_ResidualPlot.tif",
                       width=2*figdim,height=figdim,
                       units="in",pointsize=ptsz,family="sans",res=res)
par(mar=c(3.05,3.05,0.65,0.65),mgp=c(1.9,0.3,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
residPlot(fitLKT,resid.type="standardized",loess=FALSE)
if (makeTIFFFigs) dev.off()

fitLK <- nls(vbLK,data=w14T1,start=svLK)
fitLT <- nls(vbLT,data=w14T1,start=svLT)
fitKT <- nls(vbKT,data=w14T1,start=svKT)
fitL <- nls(vbL,data=w14T1,start=svL)
fitK <- nls(vbK,data=w14T1,start=svK)
fitT <- nls(vbT,data=w14T1,start=svT)
fitO <- nls(vbO,data=w14T1,start=svO)

ms <- list(fitO,fitL,fitK,fitT,fitLK,fitLT,fitKT,fitLKT)
mnames <- c("None differ","Linf differs","K differs","t0 differs",
            "Linf & K differ","Linf & t0 differ","K & t0 differ",
            "All differ")
aictab(ms,mnames)

confint(fitLT)

#
#
# #####################################################################
# =====================================================================
# ---------------------------------------------------------------------
#
# Box 12.8 (requires running Boxes 12.1 and 12.2 scripts)
#
# ---------------------------------------------------------------------
# =====================================================================

gomp1 <- function(T,Linf,tstar=NULL,gstar=NULL) {
  if (length(Linf)==3) {
    gstar <- Linf[[3]]
    tstar <- Linf[[2]]
    Linf <- Linf[[1]]
  }
  Linf*exp(-exp(-gstar*(T-tstar)))
}

rich1 <-  function(T,Linf,k=NULL,tstar=NULL,b=NULL) {
  if (length(Linf)==4) {
    b <- Linf[[4]]
    tstar <- Linf[[3]]
    k <- Linf[[2]]
    Linf <- Linf[[1]]
  }
  Linf*(1-(1/b)*exp(-k*(T-tstar)))^b
}

plot(tl~age,data=wf14T)
curve(gomp1(x,650,0.2,0.5),from=0,to=11,add=TRUE,lwd=2,lty=2)
curve(rich1(x,650,0.6,0.3,3),from=0,to=11,add=TRUE,lwd=2)
svG <- list(Linf=650,tstar=0.2,gstar=0.5)
svR <- list(Linf=650,k=0.6,tstar=0.3,b=3)

fitG <- nls(tl~gomp1(age,Linf,tstar,gstar),data=wf14T,start=svG)
fitR <- nls(tl~rich1(age,Linf,k,tstar,b),data=wf14T,start=svR)

aictab(list(fitV,fitG,fitR),c("von Bertalanffy","Gompertz","Richards"))

cbind(AIC(fitV,fitG,fitR),BIC(fitV,fitG,fitR))

#
#
# #####################################################################
# =====================================================================
# ---------------------------------------------------------------------
#
# Box 12.9 (requires running Boxes 12.1 and 12.2 scripts)
#
# ---------------------------------------------------------------------
# =====================================================================

#######################################################################
# This code sets the random seed for analyses in this box. This is not
#   required for these analyses but will cause the same results to be
#   returned as in the printed chapter if this script is sourced.
set.seed(14354454)
#######################################################################

wf14T.list <- list(tl=wf14T$tl,age=wf14T$age,N=nrow(wf14T))

wf14T.params <- c("Linf","K","t0","sigma")

wf14T.inits <- list(list(logLinf=log(500),logK=log(0.1),t0=-4,sigma= 5),
                    list(logLinf=log(600),logK=log(0.3),t0=-2,sigma=50),
                    list(logLinf=log(700),logK=log(0.5),t0= 0,sigma=95))

VBGF.model <- function() {
  # 1. Sets priors
  logLinf ~ dnorm(0,1/1000)
  logK ~ dnorm(0,1/1000)
  t0 ~ dnorm(0,1/1000)
  sigma ~ dunif(0,100)
  # 2. Exponentiate parameters to be characterized
  Linf <- exp(logLinf)
  K <- exp(logK)
  # 3. Compute model likelihood for each individual
 	for (i in 1:N) {
 		mu[i] <- Linf*(1-exp(-K*(age[i]-t0)))
  	tl[i] ~ dnorm(mu[i],1/(sigma^2))
 	}
}

#######################################################################
# The random number seed is set here to assure that the user can
#   reliably reproduce the results in the chapter. Run this code
#   before each run of jags() below.
set.seed(123456)
#
# The running of jags() below may take several hours to complete.
#   Consider setting n.iter and n.burnin to small values until you
#   are sure that the code runs successfully.
#######################################################################

wf14T.jags <- jags(wf14T.list,inits= wf14T.inits,
                   parameters.to.save=wf14T.params,model.file=VBGF.model,
                   n.chains=3,n.iter=500000,n.burnin=100000,n.thin=4)

wf14T.mcmc <- as.mcmc(wf14T.jags)

if (makeTIFFFigs) tiff("Figs/Box12_9_A_TracePlot.tif",
                       width=2*figdim,height=1.5*figdim,
                       units="in",pointsize=ptsz,family="sans",res=res)
oldpar <- par(mfrow=c(3,2),mar=c(3,3,2,1),mgp=c(1.7,0.4,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
traceplot(wf14T.mcmc,col=c("gray25","gray50","gray75"))
par(oldpar)
if (makeTIFFFigs) dev.off()

gelman.diag(wf14T.mcmc)

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
tmpLinf <- wf14T.jags$BUGSoutput$sims.list$Linf
tmpsigma <- wf14T.jags$BUGSoutput$sims.list$sigma
#######################################################################

wf14T.jags

if (makeTIFFFigs) tiff("Figs/Box12_9_B_DensityPlot.tif",
                       width=2*figdim,height=1.5*figdim,
                       units="in",pointsize=ptsz,family="sans",res=res)
oldpar <- par(mfrow=c(3,2),mar=c(3,3,2,1),mgp=c(1.7,0.4,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
densplot(wf14T.mcmc)
par(oldpar)
if (makeTIFFFigs) dev.off()

#
#
# #####################################################################
# =====================================================================
# ---------------------------------------------------------------------
#
# Box 12.10 (requires running Box 12.1 script)
#
# ---------------------------------------------------------------------
# =====================================================================

#######################################################################
# This code sets the random seed for analyses in this box. This is not
#   required for these analyses but will cause the same results to be
#   returned as in the printed chapter if this script is sourced.
set.seed(14354454)
#######################################################################

wml.01 <- filterD(WalleyeML,Year==2001)

vbTlog <- function(T,logLinf,K=NULL,t0=NULL) {
  if (length(logLinf)==3) {
    t0 <- logLinf[[3]]
    K <- logLinf[[2]]
    logLinf <- logLinf[[1]]
  }
  exp(logLinf)*(1-exp(-K*(T-t0)))
}

#######################################################################
# This code is needed but was not shown in the printed chapter
gomplog <- function(T,logLinf,tstar=NULL,gstar=NULL) {
  if (length(logLinf)==3) {
    gstar <- logLinf[[3]]
    tstar <- logLinf[[2]]
    logLinf <- logLinf[[1]]
  }
  exp(logLinf)*exp(-exp(-gstar*(T-tstar)))
}

lgstclog <- function(T,logLinf,tstar=NULL,gninf=NULL) {
  if (length(logLinf)==3) {
    gninf <- logLinf[[3]]
    tstar <- logLinf[[2]]
    logLinf <- logLinf[[1]]
  }
  exp(logLinf)/(1+exp(-gninf*(T-tstar)))
}
#######################################################################

fitVB.01 <- nlme(BC.Len~vbTlog(BC.Age,logLinf,K,t0),data=wml.01,
                 fixed=list(logLinf~1,K~1,t0~1),
                 random=logLinf+K+t0~1|ID,
                 start=list(fixed=c(logLinf=log(600),K=0.3,t0=0)))

if (makeTIFFFigs) tiff("Figs/Box12_10_A_ResidualPlot.tif",
                       width=2*figdim,height=figdim,
                       units="in",pointsize=ptsz,family="sans",res=res)
par(mar=c(3.05,3.05,0.65,0.65),mgp=c(1.9,0.3,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
residPlot(fitVB.01,resid.type="standardized",loess=FALSE,
          col=col2rgbt("black",1/5))
if (makeTIFFFigs) dev.off()

#######################################################################
# This code is needed but was not shown in the printed chapter
fitG.01 <- nlme(BC.Len~gomplog(BC.Age,logLinf,tstar,gstar),data=wml.01,
                fixed=list(logLinf~1,tstar~1,gstar~1),
                random=logLinf+tstar+gstar~1|ID,
                start=list(fixed=c(logLinf=log(600),tstar=2,gstar=0.5)))

fitL.01 <- nlme(BC.Len~lgstclog(BC.Age,logLinf,tstar,gninf),data=wml.01,
                fixed=list(logLinf~1,tstar~1,gninf~1),
                random=logLinf+tstar+gninf~1|ID,
                start=list(fixed=c(logLinf=log(600),tstar=2,gninf=0.8)))
#######################################################################

fitVB.01p <- nlme(BC.Len~vbTlog(BC.Age,logLinf,K,t0),data=wml.01,
                  fixed=list(logLinf~1,K~1,t0~1),
                  random=logLinf+K+t0~1|ID,
                  start=list(fixed=c(logLinf=log(600),K=0.3,t0=0)),
                  weights=varPower(form=~BC.Age))

#######################################################################
# This code is needed but was not shown in the printed chapter
fitG.01p <- nlme(BC.Len~gomplog(BC.Age,logLinf,tstar,gstar),data=wml.01,
                 fixed=list(logLinf~1,tstar~1,gstar~1),
                 random=logLinf+tstar+gstar~1|ID,
                 start=list(fixed=c(logLinf=log(600),tstar=2,gstar=0.5)),
                 weights=varPower(form=~BC.Age))

fitL.01p <- nlme(BC.Len~lgstclog(BC.Age,logLinf,tstar,gninf),data=wml.01,
                 fixed=list(logLinf~1,tstar~1,gninf~1),
                 random=logLinf+tstar+gninf~1|ID,
                 start=list(fixed=c(logLinf=log(600),tstar=2,gninf=0.8)),
                 weights=varPower(form=~BC.Age))
#######################################################################

aictab(list(fitVB.01,fitG.01,fitL.01,
            fitVB.01p,fitG.01p,fitL.01p),
       c("VBGF","Gompertz","logistic",
         "VBGF hetero.","Gompertz hetero.","logistic hetero."))

if (makeTIFFFigs) tiff("Figs/Box12_10_B_ResidualPlot.tif",
                       width=2*figdim,height=figdim,
                       units="in",pointsize=ptsz,family="sans",res=res)
par(mar=c(3.05,3.05,0.65,0.65),mgp=c(1.9,0.3,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
residPlot(fitG.01p,resid.type="standardized",loess=FALSE,
          col=col2rgbt("black",1/5))
if (makeTIFFFigs) dev.off()

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
fe <- fixef(fitG.01p)
sds <- VarCorr(fitG.01p)
#######################################################################

fixef(fitG.01p)

intervals(fitG.01p,which="fixed")

exp(fixef(fitG.01p)[[1]])

deltamethod(~exp(x1),fixef(fitG.01p),vcov(fitG.01p))

VarCorr(fitG.01p)

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
delta <- summary(fitG.01p)$modelStruct$varStruct[1]
vrnce <- as.numeric(sds["Residual","Variance"])
#######################################################################

ranef(fitG.01p)  # only first 3 shown

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
ranef(fitG.01p)[1:3,]
#######################################################################

coef(fitG.01p)  # only first 3 shown

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
coef(fitG.01p)[1:3,]
#######################################################################

# Base plot
if (makeTIFFFigs) tiff("Figs/Box12_10_C_FittedPlot.tif",
                       width=figdim,height=figdim,
                       units="in",pointsize=ptsz,family="sans",res=res)
par(mar=c(3.05,3.05,0.65,0.65),mgp=c(1.9,0.3,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
plot(NA,xlab="Age (years)",ylab="Total Length (mm)",
     xlim=c(0,20),ylim=c(0,800))
# Add individual best-fit curves
for (i in 1:length(unique(wml.01$ID))) {
  curve(gomplog(x,coef(fitG.01p)[i,]),from=0,to=20,
        col=col2rgbt("black",1/10),add=TRUE)
}
# Add the population-average curve
curve(gomplog(x,fixef(fitG.01p)),from=0,to=20,add=TRUE,
      lty="dashed",lwd=3,col="gray90")
if (makeTIFFFigs) dev.off()

#
#
# #####################################################################
# =====================================================================
# ---------------------------------------------------------------------
#
# Box 12.11 (requires running Boxes 12.1 and 12.10 scripts)
#
# ---------------------------------------------------------------------
# =====================================================================

#######################################################################
# This code sets the random seed for analyses in this box. This is not
#   required for these analyses but will cause the same results to be
#   returned as in the printed chapter if this script is sourced.
set.seed(14354454)
#######################################################################

fitG.01p2 <- nlme(BC.Len~gomplog(BC.Age,logLinf,tstar,gstar),data=wml.01,
                  fixed=list(logLinf~Sex,tstar~Sex,gstar~Sex),
                  random= logLinf+tstar+gstar~1|ID,
                  start=list(fixed=c(log(600),0,2,0,0.5,0)),
                  weights=varPower(form=~BC.Age))

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
fe2 <- fixef(fitG.01p2)
sds2 <- VarCorr(fitG.01p2)
#######################################################################

intervals(fitG.01p2,which="fixed")

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
ses <- deltamethod(list(~exp(x1),~exp(x1+x2),~x3+x4,~x5+x6),
                   fixef(fitG.01p2),vcov(fitG.01p2))
#######################################################################

deltamethod(list(~exp(x1),~exp(x1+x2),~x3+x4,~x5+x6),
            fixef(fitG.01p2),vcov(fitG.01p2))

VarCorr(fitG.01p2)

ranef(fitG.01p2)  # only first 3 shown

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
ranef(fitG.01p2)[1:3,]
#######################################################################

coef(fitG.01p2)  # only first 3 shown

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
coef(fitG.01p2)[1:3,]
#######################################################################

anova(fitG.01p,fitG.01p2)

# Isolate coefficients for females and males
tmp <- coef(fitG.01p2)
coefF <- filterD(tmp,grepl("F",rownames(tmp)))
coefM <- filterD(tmp,grepl("M",rownames(tmp)))
# create a base plot
if (makeTIFFFigs) tiff("Figs/Box12_11_FittedPlot.tif",
                       width=figdim,height=figdim,
                       units="in",pointsize=ptsz,family="sans",res=res)
par(mar=c(3.05,3.05,0.65,0.65),mgp=c(1.9,0.3,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
plot(NA,xlab="Age (years)",ylab="Total Length (mm)",
     xlim=c(0,20),ylim=c(0,800))
# set point color and symbols
clrs <- c("gray","black")
# add trajectories for individual females (in gray)
for (i in 1:nrow(coefF)) {
  curve(gomplog(x,coefF[i,c(1,3,5)]),from=0,to=20,add=TRUE,col=clrs[1])
}
# add trajectories for individual males (in black)
for (i in 1:nrow(coefM)) {
  cfM <- coefM[i,c(1,3,5)]+coefM[i,c(2,4,6)] # construct coefficients
  curve(gomplog(x,cfM),from=0,to=20,col=clrs[2],add=TRUE)
}
# add legend
legend("topleft",c("Female","Male"),pch=NA,col=clrs,bty="n",lty=1)
if (makeTIFFFigs) dev.off()

#
#
# #####################################################################
# =====================================================================
# ---------------------------------------------------------------------
#
# Box 12.12 (requires running Box 12.1 and Box 12.10 scripts)
#
# ---------------------------------------------------------------------
# =====================================================================

#######################################################################
# This code sets the random seed for analyses in this box. This is not
#   required for these analyses but will cause the same results to be
#   returned as in the printed chapter if this script is sourced.
set.seed(14354454)
#######################################################################

wml.f0106 <- filterD(WalleyeML,Sex=="F",Year>=2001,Year<=2006)
wml.f0106$Year <- factor(wml.f0106$Year)

if (makeTIFFFigs) tiff("Figs/Box12_12_A_XYPlot.tif",
                       width=1.5*figdim,height=figdim,
                       units="in",pointsize=ptsz,family="sans",res=res)
par(mar=c(3.05,3.05,0.65,0.65),mgp=c(1.9,0.3,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
xyplot(BC.Len~BC.Age|Year,data=wml.f0106,col=col2rgbt("black",1/5),
       xlab="Age (years)",ylab="Total Length (mm)")
if (makeTIFFFigs) dev.off()

fitVB.f0106 <- nlme(BC.Len~vbTlog(BC.Age,logLinf,K,t0),data=wml.f0106,
                    fixed=list(logLinf~1,K~1,t0~1),
                    random= list(Year=logLinf+K+t0~1,ID=logLinf+K+t0~1),
                    start=list(fixed=c(log(600),0.3,0.0)))

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
fe3 <- fixef(fitVB.f0106)
sds3 <- VarCorr(fitVB.f0106)
re3 <- ranef(fitVB.f0106)
#######################################################################

intervals(fitVB.f0106,which="fixed")
VarCorr(fitVB.f0106)

( fe <- fixef(fitVB.f0106) )
re <- ranef(fitVB.f0106)
re$Year

re$ID  # only first 3 shown

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
re$ID[1:3,]
#######################################################################

( logLinf.byYear <- fe["logLinf"]+re$Year[,"logLinf"] )

( cfs.byID <- coef(fitVB.f0106) ) # only first three shown

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
coef(fitVB.f0106)[1:3,]
#######################################################################

year <- as.numeric(substr(rownames(cfs.byID),1,4))
if (makeTIFFFigs) tiff("Figs/Box12_12_B_MeansPlot.tif",
                       width=figdim,height=figdim,
                       units="in",pointsize=ptsz,family="sans",res=res)
par(mar=c(3.05,3.05,0.65,0.65),mgp=c(1.9,0.3,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
plot(exp(cfs.byID$logLinf)~year,col=col2rgbt("black",1/4),
     xlab="Year",ylab=expression(L[infinity]))
year <- as.numeric(rownames(re$Year))
lines(exp(logLinf.byYear)~year,lwd=2,col="gray70")
if (makeTIFFFigs) dev.off()

#
#
# #####################################################################
# =====================================================================
# ---------------------------------------------------------------------
#
# Box 12.13 (requires running Boxes 12.1 and 12.10 scripts)
#
# ---------------------------------------------------------------------
# =====================================================================

npar <- 3
NFish <- length(unique(wml.01$ID))
wml.01.list <- list(Age=wml.01$BC.Age,Len=wml.01$BC.Len,ID=wml.01$ID,
                    N=nrow(wml.01),NFish=NFish,npar=npar,
                    scale=diag(npar),prec=diag(npar)/1000,mn=rep(0,npar))

wml.01.params <- c("pa.Linf","pa.tstar","pa.gstar",
                   "Linf.1","tstar.1","gstar.1",
                   "sigma","delta")

tmp <- log(c(500,0.1,0.1))
inits1 <- list(pa.log=tmp,
               params.log=matrix(rep(tmp,each=NFish),ncol=npar),
               TAU=diag(npar),sigma=5,delta=-2)

tmp <- log(c(650,2,0.5))
inits2 <- list(pa.log=tmp,
               params.log=matrix(rep(tmp,each=NFish),ncol=npar),
               TAU=diag(npar),sigma=50,delta=0)
tmp <- log(c(800,4,0.9))
inits3 <- list(pa.log=tmp,
               params.log=matrix(rep(tmp,each=NFish),ncol=npar),
               TAU=diag(npar),sigma=95,delta=2)
wml.01.inits <- list(inits1,inits2,inits3)

gomp.model <- function(){
  # 1. Set hyperpriors
  pa.log ~ dmnorm(mn,prec)
  TAU ~ dwish(scale,npar+1)
  # 2. Set priors for deviations of individual fish
  for (j in 1:NFish) {
    params.log[j,1:npar] ~ dmnorm(pa.log,TAU)
    Linf[j] <- exp(params.log[j,1])
    tstar[j] <- exp(params.log[j,2])
    gstar[j] <- exp(params.log[j,3])
  }
  # 3. Set priors for sigma and delta
  sigma ~ dunif(0,100)
  delta ~ dnorm(0,1/1000)
  # 4. Isolate params to characterize posterior
  #    population averages
  pa.Linf <- exp(pa.log[1])
  pa.tstar <- exp(pa.log[2])
  pa.gstar <- exp(pa.log[3])
  #    for the first fish
  Linf.1 <- Linf[1]
  tstar.1 <- tstar[1]
  gstar.1 <- gstar[1]
  # 5. Likelihood for each observation
  for(i in 1:N) {
    mu[i] <- Linf[ID[i]]*exp(-exp(-gstar[ID[i]]*(Age[i]-tstar[ID[i]])))
    tau[i] <- 1/((sigma^2)*Age[i]^(2*delta))
    Len[i] ~ dnorm(mu[i],tau[i])
  }
}

#######################################################################
# The random number seed is set here to assure that the user can
#   reliably reproduce the results in the chapter. Run this code
#   before each run of jags() below.
set.seed(123456)
#
# The running of jags() below may take several hours to complete.
#   Consider setting n.iter and n.burnin to small values until you
#   are sure that the code runs successfully.
#######################################################################

wml.01.jags <- jags(wml.01.list,inits=wml.01.inits,
                    parameters.to.save=wml.01.params,model.file=gomp.model,
                    n.chains=3,n.iter=200000,n.burnin=40000,n.thin=4)
wml.01.mcmc <- as.mcmc(wml.01.jags)

if (makeTIFFFigs) tiff("Figs/Box12_13_TracePlot.tif",
                       width=2*figdim,height=2.5*figdim,
                       units="in",pointsize=ptsz,family="sans",res=res)
oldpar <- par(mfrow=c(5,2),mar=c(3,3,2,1),mgp=c(1.7,0.4,0),tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
traceplot(wml.01.mcmc,col=c("gray25","gray50","gray75"))
par(oldpar)
if (makeTIFFFigs) dev.off()

gelman.diag(wml.01.mcmc)

#######################################################################
# This code was used only to typeset the chapter; it is not needed (can be deleted)
tmpLinf <- wml.01.jags$BUGSoutput$sims.list$pa.Linf
tmpsigma <- wml.01.jags$BUGSoutput$sims.list$sigma
tmpdelta <- wml.01.jags$BUGSoutput$sims.list$delta
tmpLinf1 <- wml.01.jags$BUGSoutput$sims.list$Linf.1
#######################################################################

wml.01.jags


# Script created at 2017-07-04 11:10:08

## crocodile population projection
## 
## accompanies paper: Bradshaw et al. 2006. Incorporating known sources of uncertainty
## to determine precautionary harvests of saltwater crocodiles. Ecol. Applic. 16: 1436-1448
## https://doi.org/10.1890/1051-0761(2006)016[1436:IKSOUT]2.0.CO;2
##
## ----------------------------------------------------
## base age model - density independent; deterministic
## survival values fitted from Webb & Manolis 1993;
## fecundity adjusted so N varies from 4500 to ~35000
## females in 20 years (1978 - 1998)
## ----------------------------------------------------
## CJA Bradshaw 2006

## Removes everything in memory
rm(list = ls())

## Add libraries
library(base)
library(MASS)

## define max age
age.max <- 13

## define survival
## seh <- 0.30 ## survival of egg to hatchling (Messell)
## seh <- 0.295 ## survival of egg to hatchling for C. johnstoni (Smith & Webb 1985) +/- 0.15 error

seh <- 0.314 ## survival of egg to hatchling (Webb et al. 1983)

## s0 <- 0.2 ## 0.20 to 0.80 (DD) (Webb & Manolis 1991; Webb et al. 2000)
s0 <- mean(c(0.20,0.80))
s.vec <- c(s0,0.54,0.66,0.66,0.70,0.85,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95)

## Fit survival values
## Create data frame
age.vec <- seq(0,length(s.vec)-1)
surv.fit.data <- data.frame(age.vec,s.vec)

## model formula (exponential rise to maximum): yhat = y0 + (a*(1-exp((-b)*x)))
fit.s.vec <- nls(surv.fit.data$s.vec ~ y0.coeff + (a.coeff*(1-(exp((-b.coeff)*surv.fit.data$age.vec)))),
		data = surv.fit.data,
		start = list(y0.coeff = 0.1, a.coeff = 1, b.coeff = 1),
		trace = TRUE)
summary(fit.s.vec)
sf.vec <- as.numeric(fitted(fit.s.vec))

## fecundity parameters
primi <- 13 ## age at maturity (females)

## ef <- 0.906 ## egg fertility
ef <- 1

clutch = 53.1 ## Webb et al. 1983 (combined means) with 12.23 SD

r13.hi <- 0.80 ## proportion females breeding (high) Webb et al. 1984
r13.lo <- 0.25 ## proportion females breeding (low) Webb, pers. comm.
r13 <- mean(c(r13.hi,r13.lo))
r12 <- 0.44 ## proportion 12-year olds breeding
r11 <- 0.33 ## proportion 11-year olds breeding
r10 <- 0.20 ## proportion 10-year olds breeding
r9 <- 0.1 ## proportion 9-year olds breeding
r8 <- 0.02 ## proportion 8-year olds breeding
r7 <- 0 ## proportion 7-year olds breeding

r.vec <- c(0,0,0,0,0,0,r7,r8,r9,r10,r11,r12,r13,r13,r13,r13,r13,r13)
r.age.vec <- seq(1,(length(r.vec)))

plot(r.age.vec,r.vec,xlim=range(c(1,19)),ylim=range(c(0,1)))

## Create data frame
r.fit.data <- data.frame(r.age.vec,r.vec)

## model formula 3-parameter logistic: yhat = a/(1+(x/x0)^b)

fit.r.vec <- nls(r.fit.data$r.vec ~ a.coeff / (1 + (r.fit.data$r.age.vec/x0.coeff)^b.coeff),
		data = r.fit.data,
		start = list(a.coeff = 0.5, b.coeff = -11, x0.coeff = 10),
		trace = TRUE)
summary(fit.r.vec)
rf.vec <- as.numeric(fitted(fit.r.vec))

plot(r.age.vec,rf.vec,xlim=range(c(1,19)),ylim=range(c(0,1)),type="l",xlab="age (years)",ylab="proportion females breeding")

## Sex ratios
hsr.avg <- 0.50 ## average hatchling sex ratio (prop female)
sr <- 0.50 ## adult sex ratio

## construct fertility vector
##m.vec <- c(0,0,0,0,0,0,0,seh*s.vec[1]*hsr.avg*clutch*ef*rf.vec[8],seh*s.vec[1]*hsr.avg*clutch*ef*r9,seh*s.vec[1]*hsr.avg*clutch*ef*r10,seh*s.vec[1]*hsr.avg*clutch*ef*r11,seh*s.vec[1]*hsr.avg*clutch*ef*r12,seh*s.vec[1]*hsr.avg*clutch*ef*r13)

m.vec <- rf.vec*s.vec[1]*seh*clutch*hsr.avg

## initial population sizes
Navg <- 9000 ## Webb et al. 2000

##total females in population
f <- Navg*sr

##total males in population
mal <- Navg*(1-sr)

##Initial population size vector
N <- f + mal

## Matrix size
k <- 2*age.max

## Create matrix shell
a <- matrix(data<-0,nrow<-k,ncol<-k)

## Add survival vectors to matrix
diag(a[2:age.max,1:(age.max-1)]) <- sf.vec[2:age.max] ## adds sf.vec to female quadrant
a[age.max,age.max] <- sf.vec[age.max+1] ## put final survival value in female age.max/ag.max cell

diag(a[(age.max+2):k,(age.max+1):(k-1)]) <- sf.vec[2:age.max] ## adds sf.vec to male quadrant
a[k,k] <- sf.vec[age.max+1] ## put final survival value in male age.max/ag.max cell

a[1,1:(age.max)] <- m.vec[1:(age.max)] ## adds female fecundity
a[(age.max+1),1:(age.max)] <- m.vec[1:(age.max)] ## adds male fecundity

## Maximum lambda function
max.lambda <- function(x) Re((eigen(x)$values)[1]) ## where 'x' is a Leslie matrix

## Maximum r function
max.r <- function(x) log(Re((eigen(x)$values)[1])) ## where 'x' is a Leslie matrix

## Stable stage distribution
stable.stage.dist <- function(x) ((x %*% (Re((eigen(x)$vectors)[,1])))/(sum((x %*% (Re((eigen(x)$vectors)[,1]))))))[,1]

## Generation length function

R.vals <- function(X,age.max) ## reproductive value (R0) where X = Leslie matrix; age.max = maximum age of females
{		
		## define the transition matrix ****Change according to Matrix type****
		T <- X[1:age.max,1:age.max]
		T[1,2:(age.max)] <- 0

		## define the fertility matrix ****Change according to Matrix type****
		F <- X[1:age.max,1:age.max]
		F[2:age.max,1:(age.max)] <- 0

		## define the identity matrix
		I <- matrix(data<-0,nrow<-age.max,ncol<-age.max)
		diag(I) <- 1

		## define the fundamental matrix
		library(MASS)
		N.fund <- ginv(I-T)

		## define the reproductive matrix
		R <- F %*% N.fund

		## define R0 (number of female offspring produced per female during lifetime)
		R0 <- Re((eigen(R)$values)[1])

		## Mean age at death (life expectancy)
		LE <- colSums(N.fund)[1]
		LE.var <- ((colSums((2*(N.fund %*% N.fund)) - N.fund)) - (colSums(N.fund)*colSums(N.fund)))[1]
		LE.sd <- sqrt(LE.var)

		G <- (log(R0))/(log(Re((eigen(X)$values)[1])))

		## Mean age of parents of the offspring produced by a cohort of its lifetime (MU1)
		MU1 <- ((F %*% (N.fund %*% N.fund))[1])/R0

		## Mean age of the parents of the offspring produced by a population at the stable age distribution (Abar)
		Abar <- (2*G) - MU1 ## ****May not give appropriate value for age-classified model****

		## output
		print("_______________________________________________________________________________________")
		print("number of female offspring produced per female during its lifetime")
		print(R0)
		print("_______________________________________________________________________________________")
		print("mean generation time (T) = time required for population to increase by a factor of R0")
		print(G)
		print("_______________________________________________________________________________________")
		print("Mean age of parents of the offspring produced by a cohort of its lifetime")
		print(MU1)
		print("_______________________________________________________________________________________")
		print("Mean age of the parents of the offspring produced by a population at the stable age distribution")
		print(Abar)
		print("_______________________________________________________________________________________")
		print("mean age at death (life expectancy")
		print(LE)
		print("_______________________________________________________________________________________")
		print("sd of mean age at death (life expectancy")
		print(LE.sd)

		return(data.frame(R0=R0,G=G,MU1=MU1,Abar=Abar,LE=LE,LE.sd=LE.sd))
}


## max lambda
maxlambda <- max.lambda(a)

## left and right eigenvectors & eigenvalues
w <- eigen(a)
v <- Conj(ginv(w$vectors))

## left and right eigenvalues
wsd <- Re(w$vectors[,1])
vsd <- Re(v[1,])

## Reproductive values
rep.vals <- R.vals(a,age.max)

## r
r <- max.r(a)

## Stable stage distribution
ssd <- stable.stage.dist(a)

## ssd classes
    ## female
    ssd.juvf <- sum(ssd[1:primi-1]); ssd.adf <- sum(ssd[primi:age.max])

    ## male
    ssd.juvm <- sum(ssd[(age.max+1):(age.max+primi-1)]); ssd.adm <- sum(ssd[(age.max+primi):k])

## sensitivities and elasticities
senmat <- (vsd %o% wsd)
emat <- (a * senmat) / maxlambda

## Average age of mothers at stability
G.age.mean <- weighted.mean(age.vec[age.max],ssd[age.max])
    
## pick initial vectors
n <- ssd*N

## age vector
age.vec <- rep(0,age.max)
	for (j in 1:(age.max-1)) {
	    age.vec[j+1] <- j
	}
age.vec <- age.vec+1

##Calculate Quasi-extinction times
thresh <- 50 ## Define quasi-exinction threshold

Q <- (log(thresh/sum(n)))/log(maxlambda)
	if (Q < 0) {
	    Q <- "infinity"
	}
	Q <- Q

## do a simulation
## first specify the initial condition and length of simulation
tlimit <- 20

##set population size year step vector
pop.vec <- rep(1,tlimit+1)
pop.vec[1] <- sum(n)

## set age class vectors
n1.vec <- rep(1,tlimit+1); n1.vec[1] <- n[1]
n2.vec <- rep(1,tlimit+1); n2.vec[1] <- n[2]
n3.vec <- rep(1,tlimit+1); n3.vec[1] <- n[3]
n4.vec <- rep(1,tlimit+1); n4.vec[1] <- n[4]
n5.vec <- rep(1,tlimit+1); n5.vec[1] <- n[5]
n6.vec <- rep(1,tlimit+1); n6.vec[1] <- n[6]
n7.vec <- rep(1,tlimit+1); n7.vec[1] <- n[7]
n8.vec <- rep(1,tlimit+1); n8.vec[1] <- n[8]
n9.vec <- rep(1,tlimit+1); n9.vec[1] <- n[9]
n10.vec <- rep(1,tlimit+1); n10.vec[1] <- n[10]
n11.vec <- rep(1,tlimit+1); n11.vec[1] <- n[11]
n12.vec <- rep(1,tlimit+1); n12.vec[1] <- n[12]

##set year step vector
yr.vec <- rep(1,tlimit+1)
yr.vec[1] <- 0

	for (j in 1:tlimit) {
    		yr.vec[j+1] <- j
	}

## then iterate

	for (ii in 1:tlimit) { 
	   n <- a %*% n 
	   pop.vec[ii+1] <- sum(n)
	   n1.vec[ii+1] <- n[1,1]; n2.vec[ii+1] <- n[2,1]; n3.vec[ii+1] <- n[3,1]; n4.vec[ii+1] <- n[4,1]; n5.vec[ii+1] <- n[5,1]; n6.vec[ii+1] <- n[6,1]
	   n7.vec[ii+1] <- n[7,1]; n8.vec[ii+1] <- n[8,1]; n9.vec[ii+1] <- n[9,1]; n10.vec[ii+1] <- n[10,1]; n11.vec[ii+1] <- n[11,1]; n12.vec[ii+1] <- n[12,1]
	} 

log.pop.vec <- log10(pop.vec)

##total population size after 'tlimit' years
pop.st <- N
pop.end <- sum(n)
tlimit
maxlambda
r
Q
ssd.juv <- ssd.juvf + ssd.juvm
ssd.ad <- ssd.adf + ssd.adm

##Make density independent plots
row <- 2
col <- 2
par(mfrow=c(row,col))
plot(yr.vec,pop.vec,xlab="year",ylab="N",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(pop.vec)+(0.25*(max(pop.vec))))),type='l')

##Make survival,fecundity and age-class proportion plots
plot(age.vec,sf.vec[1:(age.max)],xlab="age (years)",ylab="P(survival)",xlim=range(0.5:(age.max+0.5)),ylim=range(0:1),type='l')
plot(age.vec,m.vec[1:(age.max)],xlab="age (years)",ylab="m",xlim=range(0.5:(age.max+0.5)),ylim=range(0:(max(m.vec))),type='l')

plot(age.vec,ssd[1:age.max],xlab="age (years)",ylab="proportion in class",xlim=range(0.5:age.max+0.5),cex = 0.5, type='h')
par(mfrow=c(1,2))

## Make age-specific plots (n vector by year)
rowa <- 3
cola <- 4
par(mfrow=c(rowa,cola))
plot(yr.vec,n1.vec,xlab="year",ylab="1-year olds",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(n1.vec)+(0.25*(max(n1.vec))))),type='l')
plot(yr.vec,n2.vec,xlab="year",ylab="2-year olds",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(n2.vec)+(0.25*(max(n2.vec))))),type='l')
plot(yr.vec,n3.vec,xlab="year",ylab="3-year olds",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(n3.vec)+(0.25*(max(n3.vec))))),type='l')
plot(yr.vec,n4.vec,xlab="year",ylab="4-year olds",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(n4.vec)+(0.25*(max(n4.vec))))),type='l')
plot(yr.vec,n5.vec,xlab="year",ylab="5-year olds",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(n5.vec)+(0.25*(max(n5.vec))))),type='l')
plot(yr.vec,n6.vec,xlab="year",ylab="6-year olds",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(n6.vec)+(0.25*(max(n6.vec))))),type='l')
plot(yr.vec,n7.vec,xlab="year",ylab="7-year olds",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(n7.vec)+(0.25*(max(n7.vec))))),type='l')
plot(yr.vec,n8.vec,xlab="year",ylab="8-year olds",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(n8.vec)+(0.25*(max(n8.vec))))),type='l')
plot(yr.vec,n9.vec,xlab="year",ylab="9-year olds",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(n9.vec)+(0.25*(max(n9.vec))))),type='l')
plot(yr.vec,n10.vec,xlab="year",ylab="10-year olds",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(n10.vec)+(0.25*(max(n10.vec))))),type='l')
plot(yr.vec,n11.vec,xlab="year",ylab="11-year olds",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(n11.vec)+(0.25*(max(n11.vec))))),type='l')
plot(yr.vec,n12.vec,xlab="year",ylab="12-year olds",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(n12.vec)+(0.25*(max(n12.vec))))),type='l')
par(mfrow=c(1,2))



############################
############################
## Growth analysis
## size vs. age
############################
############################

## Female HL versus age up to 0.80 SVL (Webb et al. 1978)
Smax.smf <- 28.82 ## theoretical maximum HL in cm
Sinit.smf <- 4.6 ## initial size HL in cm
tau.smf <- 1274 ## time constant for small animals
day.vec.smf <- c(365,548,730,913,1095,1278,1460,1643,1825,2008)
gr.yr.vec.smf <- day.vec.smf/365
HL.pred.smf <- ((Smax.smf - Sinit.smf) * (1 - (exp(-day.vec.smf/tau.smf)))) + Sinit.smf
SVL.pred.smf <- (3.534*HL.pred.smf) - 3.03
TL.pred.smf <- ((HL.pred.smf - 0.84)/0.137)/100 ## predicted total length in metres

## Female HL versus age over 0.80 SVL (Webb et al. 1978)
Smax.lgf <- 51 ## theoretical maximum HL in cm
Sinit.lgf <- 10.3 ## initial size HL in cm
tau.lgf <- 4721 ## time constant for small animals
day.vec.lgf <- c(2190,2373,2555,2738,2920,3103,3285,3468,3650,3833,4015,4198,4380,4563,4745,4928,5110,5293,5475,5658,5840,6023,6205,6388,6570,6753,6935,7118,7300,7483,7665,7848,8030,8213,8395,8578,8760,8943,9125)
gr.yr.vec.lgf <- day.vec.lgf/365
HL.pred.lgf <- ((Smax.lgf - Sinit.lgf) * (1 - (exp(-day.vec.lgf/tau.lgf)))) + Sinit.lgf
SVL.pred.lgf <- (3.534*HL.pred.lgf) - 3.03
TL.pred.lgf <- ((HL.pred.lgf - 0.84)/0.137)/100 ## predicted total length in metres

par(mfrow=c(1,1))
plot(gr.yr.vec.lgf,TL.pred.lgf,xlab="age (years)",ylab="total length (m)",xlim=range(c(0,max(gr.yr.vec.lgf))),ylim=range(c(0,(max(TL.pred.lgf)))),col="red",type="p",pch=19)
points(gr.yr.vec.smf,TL.pred.smf,col="black",pch=19)
par(mfrow=c(1,1))

## Male HL versus age up to 0.80 SVL (Webb et al. 1978)
Smax.smm <- 33.17 ## theoretical maximum HL in cm
Sinit.smm <- 4.6 ## initial size HL in cm
tau.smm <- 1429 ## time constant for small animals
day.vec.smm <- c(365,548,730,913,1095,1278,1460,1643,1825,2008,2190)
gr.yr.vec.smm <- day.vec.smm/365
HL.pred.smm <- ((Smax.smm - Sinit.smm) * (1 - (exp(-day.vec.smm/tau.smm)))) + Sinit.smm
SVL.pred.smm <- (3.534*HL.pred.smm) - 3.03
TL.pred.smm <- ((HL.pred.smm - 0.84)/0.137)/100 ## predicted total length in metres

## Male HL versus age over 0.80 SVL (Webb et al. 1978)
Smax.lgm <- 65 ## theoretical maximum HL in cm
Sinit.lgm <- 9 ## initial size HL in cm
tau.lgm <- 5051 ## time constant for small animals
day.vec.lgm <- c(2373,2555,2738,2920,3103,3285,3468,3650,3833,4015,4198,4380,4563,4745,4928,5110,5293,5475,5658,5840,6023,6205,6388,6570,6753,6935,7118,7300,7483,7665,7848,8030,8213,8395,8578,8760,8943,9125)
gr.yr.vec.lgm <- day.vec.lgm/365
HL.pred.lgm <- ((Smax.lgm - Sinit.lgm) * (1 - (exp(-day.vec.lgm/tau.lgm)))) + Sinit.lgm
SVL.pred.lgm <- (3.534*HL.pred.lgm) - 3.03
TL.pred.lgm <- ((HL.pred.lgm - 0.84)/0.137)/100 ## predicted total length in metres

par(mfrow=c(1,1))
plot(gr.yr.vec.lgm,TL.pred.lgm,xlab="age (years)",ylab="total length (m)",xlim=range(c(0,max(gr.yr.vec.lgm))),ylim=range(c(0,(max(TL.pred.lgm)))),col="red",type="p",pch=19)
points(gr.yr.vec.smm,TL.pred.smm,col="black",pch=19)
par(mfrow=c(1,1))


## Fit von Bertalanffy growth function M(t) = Mmax - (Mmax - M0) exp(-kt) to combined dataset

## Growth functions
f.TL.vec <- c(TL.pred.smf,TL.pred.lgf); m.TL.vec <- c(TL.pred.smm,TL.pred.lgm)
f.age.vec <- c(gr.yr.vec.smf,gr.yr.vec.lgf); m.age.vec <- c(gr.yr.vec.smm,gr.yr.vec.lgm)

## Create data frame
gr.data <- data.frame(f.TL.vec,f.age.vec,m.TL.vec,m.age.vec)

## model formula
## Von Bertalanffy growth function: M(t) = Mmax - (Mmax - M0) exp(-kt)
## M0 = mean start mass
## Mmax = mean maximum mass
## k = rate constant per year

f.max.l <- 4 ## (Webb et al. 1978)
m.max.l <- 5 ## (Webb et al. 1978)

## integer age vector
age.vec.int <- seq(1,25,1)

fit.gr.f <- nls(gr.data$f.TL.vec ~ f.max.l - (f.max.l - gr.data$f.TL.vec[1]) * exp(-k.coeff*gr.data$f.age.vec),
		data = gr.data,
		start = list(k.coeff = 0.1),
		trace = TRUE)
sum.fit.gr.f <- summary(fit.gr.f)

fit.gr.m <- nls(gr.data$m.TL.vec ~ m.max.l - (m.max.l - gr.data$m.TL.vec[1]) * exp(-k.coeff*gr.data$m.age.vec),
		data = gr.data,
		start = list(k.coeff = 0.1),
		trace = TRUE)
sum.fit.gr.m <- summary(fit.gr.m)

## Coefficients from fit
coeff.fit.gr.f <- as.numeric(sum.fit.gr.f$parameters)
coeff.fit.gr.m <- as.numeric(sum.fit.gr.m$parameters)

## Predict new q vector with integer age inputs
grf <- f.max.l - (f.max.l - gr.data$f.TL.vec[1]) * exp(-coeff.fit.gr.f[1]*age.vec.int)
grm <- m.max.l - (m.max.l - gr.data$m.TL.vec[1]) * exp(-coeff.fit.gr.m[1]*age.vec.int)

f.mat.lo <- 2
f.mat.up <- 3

m.mat <- 3.30

f.mat.age.lo <- (log((f.max.l - f.mat.lo)/(f.max.l-gr.data$f.TL.vec[1])))/-coeff.fit.gr.f[1]
f.mat.age.up <- (log((f.max.l - f.mat.up)/(f.max.l-gr.data$f.TL.vec[1])))/-coeff.fit.gr.f[1]

f.mat.age.mean <- mean(c(f.mat.age.lo,f.mat.age.up))

m.mat.age <- (log((m.max.l - m.mat)/(m.max.l-gr.data$m.TL.vec[1])))/-coeff.fit.gr.m[1]

row <- 2
col <- 1
par(mfrow=c(row,col))
plot(age.vec.int,grf,xlab="age (years)",ylab="total length (m)",ylim=range(c(0,5)),type="l",col="red")
points(gr.yr.vec.smf,TL.pred.smf,col="black",pch=20,cex=0.4)
points(gr.yr.vec.lgf,TL.pred.lgf,col="black",pch=20,cex=0.4)
abline(h=f.mat.lo,col="black")
abline(h=f.mat.up,col="black")
abline(v=f.mat.age.lo,col="black")
abline(v=f.mat.age.up,col="black")
plot(age.vec.int,grm,col="red",xlab="age (years)",ylab="total length (m)",ylim=range(c(0,5)),type="l")
points(gr.yr.vec.smm,TL.pred.smm,col="black",pch=20,cex=0.4)
points(gr.yr.vec.lgm,TL.pred.lgm,col="black",pch=20,cex=0.4)
abline(h=m.mat,col="black")
abline(v=m.mat.age,col="black")
##title(main="von Bertalanffy Growth Functions",sub="(male vs. female)")
##legend(20,1,c("male","female"),c("red","black"), c(19,19))
par(mfrow=c(1,2))

################################
## Manuscript figure - growth ##
################################

## metres to feet
grf.ft <- grf * (1/0.3048)
grm.ft <- grm * (1/0.3048)

######################################################
## **************************************************#
## Density-dependence in survival & fertility        #
## - stochasticity in nest flooding events (user set)#
## - set harvest rates (user set)                    #
## **************************************************#
######################################################

## Start user-defined section
## ***********************************************************************************************
##K population size vector
Nd.k <- 140000

## Set initial population vector (where vital rates are maximal)
Nd.st <- 9000

## set simulation initial population size (set to Nd.st if starting from recovery)
Nd.init <- 50000 ## Nd.st ## 50000

## identify start year for plotting purposes
st.yr <- 2004 ## 1978

## Set whether to use start (1978) or recent (1996-1998) age distributions for start vector
## 1 = use start proportions (for retrospective analysis)
## 2 = use recent proportions (for projections)
sd.init <- 2

## Use observed (average 1996-1998) recent or predicted proportions
## 0 = observed
## 1 = predicted
sd.source <- 0

## Set time limit for projection
tlimit <- 30

## Set min and max s0 values
s0.min <- 0.8 ## 0.8
s0.max <- 0.2 ## 0.2

## Set min and max s1 values
s1.min <- sf.vec[2]*1.6 ## 0.80
s1.max <- sf.vec[2]*0.4 ## 0.45

## Set min and max s2 values
s2.min <- 0.95 ## sf.vec[3]*1.6 ## 0.85
s2.max <- sf.vec[3]*0.4 ## 0.50

## Set min and max s3 values
s3.min <- 0.95
s3.max <- sf.vec[4]*0.4 ## 0.60

## Set min and max proportion of mature females breeding
pfb.min <- rf.vec[13]*1.6 ## 0.80
pfb.max <- rf.vec[13]*0.4 ## 0.20

## Set density-dependence scenario (1 = single dd relationship over all regions; 2 = separate dd relationships per region)
dd.select <- 2

## Set whether stochastic flooding events occur (vary survival of egg to hatching [seh])
## Set min and max seh
seh.min <- 0.24
seh.max <- 0.36

## 1 = deterministic; 2 = stochastic
stoch.select <- 2

## Set number of iterations to estimate mean & confidence interval for projections based on flooding stochasticity
perm <- 1000 ## number of permutations
if (stoch.select == 1) perm <- 1 else perm <- perm

## Egg harvest
## For retrospective analysis, set proportions of eggs harvested per region
## Average values from 1999-2002
pharv.1 <- 0.0130
pharv.2 <- 0.0014
pharv.3 <- 0.0075

## For projected populations, set proportion of total eggs harvested per region
ptharv.1 <- 0.1549
ptharv.2 <- 0.0149
ptharv.3 <- 0.1251

## Model harvest from aboriginal subsistence (150/yr), problem (200/yr), incidental fisheries catches (500/yr) & landowner harvests (250/yr)
## include historic harvests at the following rates (1 = yes; 0 = no)
harv.hist.incl <- 1

harv.as <- 150
harv.pr <- 200
harv.if <- 500
harv.lo <- 250
harv.hist <- harv.as + harv.pr + harv.if + harv.lo

## Set kill type
## 0 = no kill
## 1 = adult males > 3 m (12+ -year olds)
## 2 = indiscrimant mature adult kills (> 2.5 m)
## 3 = constant number per year (as set in NT Croc MP) - must define kill.num values
kill.type <- 3

##################################################################
## Set constant annual kill rates according to NT Management Plan
##################################################################
kill.num.ad <- 600
kill.num.juv <- 500

## Set kill rate (per cent/year)
kill.pc <- 0.0 ## % of [start only OR time-step; see below] population per year

########################################################
## NB: kill rate is constant; if adaptive (i.e. varies
##     according to current population size, set below
########################################################
## 0 = constant kill rate
## 1 = varies according to time step population size
kill.vary <- 1

########################################################
## NB: kill is divided equally among Regions by default
##     to set individual-region kill rates, set below
########################################################
## 0 = equally divided among regions
## 1 = user-set rates varying among regions
kill.dist <- 0

## Set kill rates (percentages)
kill.pc1 <- 0 ## user-defined kill rate for Region 1
kill.pc2 <- 0 ## user-defined kill rate for Region 2
kill.pc3 <- 0 ## user-defined kill rate for Region 3

## ************************************************************************************************
## End user control


## Define start vectors for each region
################################################################################
## Region 1
## Adelaide/Mary/Wildman/W Alligator/S Alligator/E Alligator/Coopers/Murganella
################################################################################
##st.sd1.tmp <- (c(0.22,0.12,0.19,0.13,0.08,0.045,0.040,0.0292,0.0292,0.0292,0.0292,0.0291,0.0291))/2
st.sd1.tmp <- (c(0.2954,0.1279,0.1435,0.1056,0.0720,0.0320,0.0320,0.0320,0.0320,0.0319,0.0319,0.0319,0.0319))/2
st.sd1 <- c(st.sd1.tmp,st.sd1.tmp) ## Region 1 start vector
rowa <- 3
cola <- 1
par(mfrow=c(rowa,cola))
plot(st.sd1[1:age.max],type="l",xlab="age class",ylab="proportion",ylim=range(c(0,0.25)))
title(main="Region 1")

##########################################################################
## Region 2
## Victoria/Fitzmaurice/Moyle/Daly/Reynolds/Finniss/
## Hart/Rose/Roper/Towns/Nathan/McArthur/Johnson/Wearyan/Robinson/Calvert
##########################################################################
##st.sd2.tmp <- (c(0.11,0.09,0.16,0.14,0.10,0.06,0.06,0.047,0.047,0.047,0.047,0.046,0.046))/2
st.sd2.tmp <- (c(0.1118,0.0944,0.1623,0.1475,0.1017,0.0478,0.0478,0.0478,0.0478,0.0478,0.0478,0.0478,0.0477))/2
st.sd2 <- c(st.sd2.tmp,st.sd2.tmp) ## Region 2 start vector
plot(st.sd2[1:age.max],type="l",xlab="age class",ylab="proportion",ylim=range(c(0,0.25)))
title(main="Region 2")

###############################################################################################################
## Region 3
## King/Goomadeer/Liverpool/Tomkinson/Cadell/Blyth/Glyde/Goyder/Woolen/Kalarwoi/
## Buckingham/Warawurowoi/Kurala/Habgood/Darwarunga/Baralminar/Gobalpa/Goromuru/Cato/Peter John/Burungbirinung
###############################################################################################################
##st.sd3.tmp <- (c(0.48,0.16,0.17,0.09,0.03,0.01,0.01,0.01,0.01,0.01,0.01,0.005,0.005))/2
st.sd3.tmp <- (c(0.4966,0.1334,0.1494,0.0817,0.0383,0.01258,0.01258,0.01258,0.01258,0.01258,0.01258,0.01258,0.01254))/2
st.sd3 <- c(st.sd3.tmp,st.sd3.tmp)
plot(st.sd3[1:age.max],type="l",xlab="age class",ylab="proportion",ylim=range(c(0,0.25)))
title(main="Region 3")
par(mfrow=c(1,2))

## Proportion of total initial population in each region
num.1 <- 160*8; num.2 <- 40*16; num.3 <- 80*21; num <- num.1 + num.2 + num.3
dens.1 <- 1.9; dens.2 <- 1.2; dens.3 <- 3.2
km.1 <- 600; km.2 <- 1400; km.3 <- 500
km.dens.1 <- km.1*dens.1; km.dens.2 <- km.2*dens.2; km.dens.3 <- km.3*dens.3
nprop1 = km.dens.1 / (km.dens.1+km.dens.2+km.dens.3); nprop2 = km.dens.2 / (km.dens.1+km.dens.2+km.dens.3); nprop3 = km.dens.3 / (km.dens.1+km.dens.2+km.dens.3)

## Egg harvest check
if (sd.init == 1) ptharv.1 <- 0
if (sd.init == 1) ptharv.2 <- 0
if (sd.init == 1) ptharv.3 <- 0

## Historic harvest values
harv.hist.age <- harv.hist/k/3 ## uniform over all 3 regions
harv.hist.vec <- rep(harv.hist.age,k)
if (harv.hist.incl == 1) harv.hist <- harv.hist else harv.hist <- 0

## Kill-related parameter calculation
kill <- Nd.init * (kill.pc/100)

if (kill.type == 0) scenario <- "no kill"
if (kill.type == 1) scenario <- "males >3 m killed"
if (kill.type == 2) scenario <- "males & females >2.5 m killed"
if (kill.vary == 0) kill.vary.scenario <- "(constant)"
if (kill.vary == 1) kill.vary.scenario <- "(f(Ni))"
if (kill.dist == 0) {
	kill.1 <- kill/3; kill.2 <- kill/3; kill.3 <- kill/3
}
	kill.1 <- ((kill.pc1/100)*Nd.init*nprop1); kill.2 <- ((kill.pc2/100)*Nd.init*nprop2); kill.3 <- ((kill.pc3/100)*Nd.init*nprop3)
if (kill.dist == 0) kill <- kill else kill <- kill.1 + kill.2 + kill.3
if (kill.type == 0) kill <- 0
if (kill.type == 0) kill.1 <- 0
if (kill.type == 0) kill.2 <- 0
if (kill.type == 0) kill.3 <- 0
if (kill.type == 0) kill.pc <- 0
if (kill.pc == 0) kill.vary.scenario <- "(constant)"

## ************************************************************************************************
## Re-calculate logistic functions for survival & fecundity based on initial population vector

library(stats)
library(base)

## Vector based on total population (Region 1)
pop.k.1 <- Nd.k * nprop1; pop.min.1 <- Nd.st * nprop1; pop.init.1 <- Nd.init * nprop1; pop.step <- 50; xmid.1 <- (pop.k.1 + pop.min.1)/2; num.vec.1 <- seq(pop.min.1,pop.k.1,((pop.k.1-pop.min.1)/pop.step))

## Vector based on total population (Region 2)
pop.k.2 <- Nd.k * nprop2; pop.min.2 <- Nd.st * nprop2; pop.init.2 <- Nd.init * nprop2; pop.step <- 50; xmid.2 <- (pop.k.2 + pop.min.2)/2; num.vec.2 <- seq(pop.min.2,pop.k.2,((pop.k.2-pop.min.2)/pop.step))

## Vector based on total population (Region 3)
pop.k.3 <- Nd.k * nprop3; pop.min.3 <- Nd.st * nprop3; pop.init.3 <- Nd.init * nprop3; pop.step <- 50; xmid.3 <- (pop.k.3 + pop.min.3)/2; num.vec.3 <- seq(pop.min.3,pop.k.3,((pop.k.3-pop.min.3)/pop.step))

## Vector based on total population (regions combined)
num.vec <- (num.vec.1+num.vec.2+num.vec.3)/3


## Vector based on hatchling population (Region 1)
poph.k.1 <- (hsr.avg*seh*clutch*ef*0.25) * pop.k.1 * ssd[age.max]; poph.min.1 <- (hsr.avg*seh*clutch*ef*0.80) * pop.min.1 * st.sd1[age.max]; xmidh.1 <- (poph.k.1 + poph.min.1)/2; numh.vec.1 <- seq(poph.min.1,poph.k.1,((poph.k.1-poph.min.1)/pop.step)); 

## Vector based on hatchling population (Region 2)
poph.k.2 <- (hsr.avg*seh*clutch*ef*0.25) * pop.k.2 * ssd[age.max]; poph.min.2 <- (hsr.avg*seh*clutch*ef*0.80) * pop.min.2 * st.sd2[age.max]; xmidh.2 <- (poph.k.2 + poph.min.2)/2; numh.vec.2 <- seq(poph.min.2,poph.k.2,((poph.k.2-poph.min.2)/pop.step))

## Vector based on hatchling population (Region 3)
poph.k.3 <- (hsr.avg*seh*clutch*ef*0.25) * pop.k.3 * ssd[age.max]; poph.min.3 <- (hsr.avg*seh*clutch*ef*0.80) * pop.min.3 * st.sd3[age.max]; xmidh.3 <- (poph.k.3 + poph.min.3)/2; numh.vec.3 <- seq(poph.min.3,poph.k.3,((poph.k.3-poph.min.3)/pop.step))

## Vector based on hatchling population (regions combined)
numh.vec <- (numh.vec.1+numh.vec.2+numh.vec.3)/3


## Vector based on 8-12 year old population (Region 1)
poplg.k.1 <- sum((pop.k.1 * ssd[8:age.max]) * 2); poplg.min.1 <- sum((pop.min.1 * st.sd1[8:age.max]) * 2); xmidlg.1 <- (poplg.k.1 + poplg.min.1)/2; numlg.vec.1 <- seq(poplg.min.1,poplg.k.1,((poplg.k.1-poplg.min.1)/pop.step))

## Vector based on 8-12 year old population (Region 2)
poplg.k.2 <- sum((pop.k.2 * ssd[8:age.max]) * 2); poplg.min.2 <- sum((pop.min.2 * st.sd2[8:age.max]) * 2); xmidlg.2 <- (poplg.k.2 + poplg.min.2)/2; numlg.vec.2 <- seq(poplg.min.2,poplg.k.2,((poplg.k.2-poplg.min.2)/pop.step))

## Vector based on 8-12 year old population (Region 3)
poplg.k.3 <- sum((pop.k.3 * ssd[8:age.max]) * 2); poplg.min.3 <- sum((pop.min.3 * st.sd3[8:age.max]) * 2); xmidlg.3 <- (poplg.k.3 + poplg.min.3)/2; numlg.vec.3 <- seq(poplg.min.3,poplg.k.3,((poplg.k.3-poplg.min.3)/pop.step))

## Vector based on 8-12 year old population (Region 3)
numlg.vec <- (numlg.vec.1+numlg.vec.2+numlg.vec.3)/3


## Set up base logistic form
quant <- seq(0,1,(1/pop.step)) ## sets up quantile vector
logistic.dens <- dlogis(quant, location = 0, scale = 0.25, log = FALSE) ## estimates logistic density function with long left tail
log.dens <- (logistic.dens - min(logistic.dens)) / max(logistic.dens - min(logistic.dens))

## Calculate logistic function between min and max s0
s0.vec <- (seq(s0.min,s0.max,-((s0.min-s0.max)/pop.step))); range.s0.vec <- range(s0.vec)[2] - range(s0.vec)[1]; s0.vec.logis <- (range.s0.vec * log.dens) + min(s0.vec)

## Calculate logistic function between min and max s1
s1.vec <- (seq(s1.min,s1.max,-((s1.min-s1.max)/pop.step))); range.s1.vec <- range(s1.vec)[2] - range(s1.vec)[1]; s1.vec.logis <- (range.s1.vec * log.dens) + min(s1.vec)

## Calculate logistic function between min and max s2
s2.vec <- (seq(s2.min,s2.max,-((s2.min-s2.max)/pop.step))); range.s2.vec <- range(s2.vec)[2] - range(s2.vec)[1]; s2.vec.logis <- (range.s2.vec * log.dens) + min(s2.vec)

## Calculate logistic function between min and max s3
s3.vec <- (seq(s3.min,s3.max,-((s3.min-s3.max)/pop.step))); range.s3.vec <- range(s3.vec)[2] - range(s3.vec)[1]; s3.vec.logis <- (range.s3.vec * log.dens) + min(s3.vec)

## Calculate logistic function between min and max s1p
pfb.vec <- (seq(pfb.min,pfb.max,-((pfb.min-pfb.max)/pop.step))); range.pfb.vec <- range(pfb.vec)[2] - range(pfb.vec)[1]; pfb.vec.logis <- (range.pfb.vec * log.dens) + min(pfb.vec)

## Plot dd functions
row <- 3
col <- 3
par(mfrow=c(row,col))
plot(numh.vec.1,s0.vec.logis,xlab="N hatchlings",ylab="s0",xlim=range(-0.5:((max(numh.vec.1))+1)),ylim=range(c(0,1)),type='l')
title(main="Region 1")
plot(numh.vec.2,s0.vec.logis,xlab="N hatchlings",ylab="s0",xlim=range(-0.5:((max(numh.vec.2))+1)),ylim=range(c(0,1)),type='l')
title(main="Region 2")
plot(numh.vec.3,s0.vec.logis,xlab="N hatchlings",ylab="s0",xlim=range(-0.5:((max(numh.vec.3))+1)),ylim=range(c(0,1)),type='l')
title(main="Region 3")
plot(numlg.vec.1,s1.vec.logis,xlab="N >2.0 m",ylab="s1",xlim=range(-0.5:((max(numlg.vec.1))+1)),ylim=range(c(0,1)),type='l')
plot(numlg.vec.2,s1.vec.logis,xlab="N >2.0 m",ylab="s1",xlim=range(-0.5:((max(numlg.vec.2))+1)),ylim=range(c(0,1)),type='l')
plot(numlg.vec.3,s1.vec.logis,xlab="N >2.0 m",ylab="s1",xlim=range(-0.5:((max(numlg.vec.3))+1)),ylim=range(c(0,1)),type='l')
plot(num.vec.1,s2.vec.logis,xlab="N total",ylab="s2",xlim=range(-0.5:((max(num.vec.1))+1)),ylim=range(c(0,1)),type='l')
plot(num.vec.2,s2.vec.logis,xlab="N total",ylab="s2",xlim=range(-0.5:((max(num.vec.2))+1)),ylim=range(c(0,1)),type='l')
plot(num.vec.3,s2.vec.logis,xlab="N total",ylab="s2",xlim=range(-0.5:((max(num.vec.3))+1)),ylim=range(c(0,1)),type='l')
par(mfrow=c(1,2))

row <- 2
col <- 3
par(mfrow=c(row,col))
plot(num.vec.1,s3.vec.logis,xlab="N total",ylab="s3",xlim=range(-0.5:((max(num.vec.1))+1)),ylim=range(c(0,1)),type='l')
title(main="Region 1")
plot(num.vec.2,s3.vec.logis,xlab="N total",ylab="s3",xlim=range(-0.5:((max(num.vec.2))+1)),ylim=range(c(0,1)),type='l')
title(main="Region 2")
plot(num.vec.3,s3.vec.logis,xlab="N total",ylab="s3",xlim=range(-0.5:((max(num.vec.3))+1)),ylim=range(c(0,1)),type='l')
title(main="Region 3")
plot(num.vec.1,pfb.vec.logis,xlab="N total",ylab="prop f breeding",xlim=range(-0.5:((max(num.vec.1))+1)),ylim=range(c(0,1)),type='l')
plot(num.vec.2,pfb.vec.logis,xlab="N total",ylab="prop f breeding",xlim=range(-0.5:((max(num.vec.2))+1)),ylim=range(c(0,1)),type='l')
plot(num.vec.3,pfb.vec.logis,xlab="N total",ylab="prop f breeding",xlim=range(-0.5:((max(num.vec.3))+1)),ylim=range(c(0,1)),type='l')
par(mfrow=c(1,2))

## Combine dd functions across regions
row <- 2
col <- 3
par(mfrow=c(row,col))
plot(numh.vec,s0.vec.logis,xlab="N hatchlings",ylab="s0",xlim=range(-0.5:((max(numh.vec))+1)),ylim=range(c(0,1)),type='l')
plot(numlg.vec,s1.vec.logis,xlab="N >2.0 m",ylab="s1",xlim=range(-0.5:((max(numlg.vec))+1)),ylim=range(c(0,1)),type='l')
plot(num.vec,s2.vec.logis,xlab="N total",ylab="s2",xlim=range(-0.5:((max(num.vec))+1)),ylim=range(c(0,1)),type='l')
plot(num.vec,s3.vec.logis,xlab="N total",ylab="s3",xlim=range(-0.5:((max(num.vec))+1)),ylim=range(c(0,1)),type='l')
plot(num.vec,pfb.vec.logis,xlab="N total",ylab="prop f breeding",xlim=range(-0.5:((max(num.vec))+1)),ylim=range(c(0,1)),type='l')
par(mfrow=c(1,2))

## Create data frame
log.func <- data.frame(num.vec.1,num.vec.2,num.vec.3,numh.vec.1,numh.vec.2,numh.vec.3,numlg.vec.1,numlg.vec.2,numlg.vec.3,s0.vec.logis,s1.vec.logis,s2.vec.logis,s3.vec.logis,pfb.vec.logis)

## Estimate 4-parameter logistic model coefficients
s0.param.1 <- getInitial(s0.vec.logis ~ SSfpl(numh.vec.1, a.s0.1, b.s0.1, xmid.s0.1, scal.s0.1), data=log.func)
s0.param.2 <- getInitial(s0.vec.logis ~ SSfpl(numh.vec.2, a.s0.2, b.s0.2, xmid.s0.2, scal.s0.2), data=log.func)
s0.param.3 <- getInitial(s0.vec.logis ~ SSfpl(numh.vec.3, a.s0.3, b.s0.3, xmid.s0.3, scal.s0.3), data=log.func)
s0.param <- getInitial(s0.vec.logis ~ SSfpl(numh.vec, a.s0, b.s0, xmid.s0, scal.s0), data=log.func)

s1.param.1 <- getInitial(s1.vec.logis ~ SSfpl(numlg.vec.1, a.s1.1, b.s1.1, xmid.s1.1, scal.s1.1), data=log.func)
s1.param.2 <- getInitial(s1.vec.logis ~ SSfpl(numlg.vec.2, a.s1.2, b.s1.2, xmid.s1.2, scal.s1.2), data=log.func)
s1.param.3 <- getInitial(s1.vec.logis ~ SSfpl(numlg.vec.3, a.s1.3, b.s1.3, xmid.s1.3, scal.s1.3), data=log.func)
s1.param <- getInitial(s1.vec.logis ~ SSfpl(numlg.vec, a.s1, b.s1, xmid.s1, scal.s1), data=log.func)

s2.param.1 <- getInitial(s2.vec.logis ~ SSfpl(num.vec.1, a.s2.1, b.s2.1, xmid.s2.1, scal.s2.1), data=log.func)
s2.param.2 <- getInitial(s2.vec.logis ~ SSfpl(num.vec.2, a.s2.2, b.s2.2, xmid.s2.2, scal.s2.2), data=log.func)
s2.param.3 <- getInitial(s2.vec.logis ~ SSfpl(num.vec.3, a.s2.3, b.s2.3, xmid.s2.3, scal.s2.3), data=log.func)
s2.param <- getInitial(s2.vec.logis ~ SSfpl(num.vec, a.s2, b.s2, xmid.s2, scal.s2), data=log.func)

s3.param.1 <- getInitial(s3.vec.logis ~ SSfpl(num.vec.1, a.s3.1, b.s3.1, xmid.s3.1, scal.s3.1), data=log.func)
s3.param.2 <- getInitial(s3.vec.logis ~ SSfpl(num.vec.2, a.s3.2, b.s3.2, xmid.s3.2, scal.s3.2), data=log.func)
s3.param.3 <- getInitial(s3.vec.logis ~ SSfpl(num.vec.3, a.s3.3, b.s3.3, xmid.s3.3, scal.s3.3), data=log.func)
s3.param <- getInitial(s3.vec.logis ~ SSfpl(num.vec, a.s3, b.s3, xmid.s3, scal.s3), data=log.func)

pfb.param.1 <- getInitial(pfb.vec.logis ~ SSfpl(num.vec.1, a.pfb.1, b.pfb.1, xmid.pfb.1, scal.pfb.1), data=log.func)
pfb.param.2 <- getInitial(pfb.vec.logis ~ SSfpl(num.vec.2, a.pfb.2, b.pfb.2, xmid.pfb.2, scal.pfb.2), data=log.func)
pfb.param.3 <- getInitial(pfb.vec.logis ~ SSfpl(num.vec.3, a.pfb.3, b.pfb.3, xmid.pfb.3, scal.pfb.3), data=log.func)
pfb.param <- getInitial(pfb.vec.logis ~ SSfpl(num.vec, a.pfb, b.pfb, xmid.pfb, scal.pfb), data=log.func)

## ************************************************************************************************

## dd scenario selector
if (dd.select == 1) s0.param.1 <- s0.param else s0.param.1 <- s0.param.1; if (dd.select == 1) s0.param.2 <- s0.param else s0.param.2 <- s0.param.2; if (dd.select == 1) s0.param.3 <- s0.param else s0.param.3 <- s0.param.3
if (dd.select == 1) s1.param.1 <- s1.param else s1.param.1 <- s1.param.1; if (dd.select == 1) s1.param.2 <- s1.param else s1.param.2 <- s1.param.2; if (dd.select == 1) s1.param.3 <- s1.param else s1.param.3 <- s1.param.3
if (dd.select == 1) s2.param.1 <- s2.param else s2.param.1 <- s2.param.1; if (dd.select == 1) s2.param.2 <- s2.param else s2.param.2 <- s2.param.2; if (dd.select == 1) s2.param.3 <- s2.param else s2.param.3 <- s2.param.3
if (dd.select == 1) s3.param.1 <- s3.param else s3.param.1 <- s3.param.1; if (dd.select == 1) s3.param.2 <- s3.param else s3.param.2 <- s3.param.2; if (dd.select == 1) s3.param.3 <- s3.param else s3.param.3 <- s3.param.3
if (dd.select == 1) pfb.param.1 <- pfb.param else pfb.param.1 <- pfb.param.1; if (dd.select == 1) pfb.param.2 <- pfb.param else pfb.param.2 <- pfb.param.2; if (dd.select == 1) pfb.param.3 <- pfb.param else pfb.param.3 <- pfb.param.3

## Set storage vectors for mean and confidence intervals for population sizes
popd.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); popd.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); popd.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1); popd.perm <- matrix(0,nrow=perm,ncol=tlimit+1)

## Set storage vectors for mean and confidence intervals for age class vectors
n0fd.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n0fd.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n0fd.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n1fd.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n1fd.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n1fd.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n2fd.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n2fd.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n2fd.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n3fd.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n3fd.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n3fd.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n4fd.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n4fd.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n4fd.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n5fd.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n5fd.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n5fd.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n6fd.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n6fd.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n6fd.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n7fd.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n7fd.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n7fd.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n8fd.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n8fd.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n8fd.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n9fd.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n9fd.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n9fd.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n10fd.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n10fd.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n10fd.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n11fd.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n11fd.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n11fd.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n12fd.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n12fd.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n12fd.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n13fd.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n13fd.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n13fd.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);

n0md.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n0md.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n0md.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n1md.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n1md.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n1md.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n2md.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n2md.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n2md.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n3md.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n3md.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n3md.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n4md.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n4md.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n4md.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n5md.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n5md.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n5md.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n6md.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n6md.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n6md.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n7md.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n7md.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n7md.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n8md.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n8md.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n8md.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n9md.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n9md.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n9md.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n10md.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n10md.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n10md.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n11md.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n11md.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n11md.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n12md.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n12md.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n12md.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
n13md.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n13md.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); n13md.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);

eggs.1.perm <- matrix(0,nrow=perm,ncol=tlimit+1); eggs.2.perm <- matrix(0,nrow=perm,ncol=tlimit+1); eggs.3.perm <- matrix(0,nrow=perm,ncol=tlimit+1);
h.a1.1.fix.perm <- matrix(0,nrow=perm,ncol=tlimit+1); h.a1.2.fix.perm <- matrix(0,nrow=perm,ncol=tlimit+1); h.a1.3.fix.perm <- matrix(0,nrow=perm,ncol=tlimit+1)

## Stochastic catastrophe scenario modelled from 0.14G (Reed et al. 2004)
p.cat <- 0.147/rep.vals[2]
sev.cat.vec <- seq(0,100,10)
sev.cat.prob <- (2510*(0.9376^(sev.cat.vec)))/100

sev.cat.dat <- data.frame(sev.cat.vec,sev.cat.prob)
sev.cat.fit <- nls(sev.cat.vec ~ (a.cat*(sev.cat.prob^b.cat)),
		data = sev.cat.dat,
		start = list(a.cat = 50, b.cat = -0.2),
		trace = TRUE)
sum.sev.cat.fit <- summary(sev.cat.fit)
sev.cat.coeff <- as.numeric(sum.sev.cat.fit$parameters)

## Recent proportions per age class (mean 1996-1998)
################################################################################
## Region 1
## Adelaide/Mary/Wildman/W Alligator/S Alligator/E Alligator/Coopers/Murganella
################################################################################
end.sd1.tmp1 <- (c(0.2420,0.0375,0.0492,0.0510,0.0634,0.5568))
ssdd.1.tmp <- (ssd[6:age.max]/(sum(ssd[6:age.max])))
p6.ssd.1 <- (ssdd.1.tmp) * (end.sd1.tmp1[6])/2
end.sd1.tmp <- rep(0,age.max)
end.sd1.tmp[1:5] <- end.sd1.tmp1[1:5]/2
end.sd1.tmp[6:age.max] <- p6.ssd.1
end.sd1 <- c(end.sd1.tmp,end.sd1.tmp) ## Region 1 start vector

##########################################################################
## Region 2
## Victoria/Fitzmaurice/Moyle/Daly/Reynolds/Finniss/
## Hart/Rose/Roper/Towns/Nathan/McArthur/Johnson/Wearyan/Robinson/Calvert
##########################################################################
end.sd2.tmp1 <- (c(0.0501,0.0258,0.0388,0.0434,0.0712,0.7707))
ssdd.2.tmp <- (ssd[6:age.max]/(sum(ssd[6:age.max])))
p6.ssd.2 <- (ssdd.2.tmp) * (end.sd2.tmp1[6])/2
end.sd2.tmp <- rep(0,age.max)
end.sd2.tmp[1:5] <- end.sd2.tmp1[1:5]/2
end.sd2.tmp[6:age.max] <- p6.ssd.2
end.sd2 <- c(end.sd2.tmp,end.sd2.tmp) ## Region 3 start vector

###############################################################################################################
## Region 3
## King/Goomadeer/Liverpool/Tomkinson/Cadell/Blyth/Glyde/Goyder/Woolen/Kalarwoi/
## Buckingham/Warawurowoi/Kurala/Habgood/Darwarunga/Baralminar/Gobalpa/Goromuru/Cato/Peter John/Burungbirinung
###############################################################################################################
end.sd3.tmp1 <- (c(0.3374,0.0397,0.0441,0.0324,0.0468,0.4997))
ssdd.3.tmp <- (ssd[6:age.max]/(sum(ssd[6:age.max])))
p6.ssd.3 <- (ssdd.3.tmp) * (end.sd3.tmp1[6])/2
end.sd3.tmp <- rep(0,age.max)
end.sd3.tmp[1:5] <- end.sd3.tmp1[1:5]/2
end.sd3.tmp[6:age.max] <- p6.ssd.3
end.sd3 <- c(end.sd3.tmp,end.sd3.tmp) ## Region 3 start vector

## Choose which age distribution to use in start matrix
## pick initial vectors
if (sd.init == 1) sd1 <- (st.sd1) else sd1 <- (end.sd1)
if (sd.init == 1) sd2 <- (st.sd2) else sd2 <- (end.sd2)
if (sd.init == 1) sd3 <- (st.sd3) else sd3 <- (end.sd3)

for (z in 1:perm) {

## define start survival

        s0.d.1 <- as.numeric(SSfpl((pop.init.1*sd1[age.max]*seh*clutch*hsr.avg*ef*pfb.min),s0.param.1[1],s0.param.1[2],s0.param.1[3],s0.param.1[4]))
            if (s0.d.1 < s0.max) s0.d.1 <- s0.max else s0.d.1 <- s0.d.1
	    if (s0.d.1 > s0.min) s0.d.1 <- s0.min else s0.d.1 <- s0.d.1

        s0.d.2 <- as.numeric(SSfpl((pop.init.2*sd2[age.max]*seh*clutch*hsr.avg*ef*pfb.min),s0.param.2[1],s0.param.2[2],s0.param.2[3],s0.param.2[4]))
            if (s0.d.2 < s0.max) s0.d.2 <- s0.max else s0.d.2 <- s0.d.2
	    if (s0.d.2 > s0.min) s0.d.2 <- s0.min else s0.d.2 <- s0.d.2

        s0.d.3 <- as.numeric(SSfpl((pop.init.3*sd3[age.max]*seh*clutch*hsr.avg*ef*pfb.min),s0.param.3[1],s0.param.3[2],s0.param.3[3],s0.param.3[4]))
            if (s0.d.3 < s0.max) s0.d.3 <- s0.max else s0.d.3 <- s0.d.3
	    if (s0.d.3 > s0.min) s0.d.3 <- s0.min else s0.d.3 <- s0.d.3

        s1.d.1 <- as.numeric(SSfpl((sum(pop.init.1*sd1[8:age.max])*2),s1.param.1[1],s1.param.1[2],s1.param.1[3],s1.param.1[4]))
            if (s1.d.1 < s1.max) s1.d.1 <- s1.max else s1.d.1 <- s1.d.1
	    if (s1.d.1 > s1.min) s1.d.1 <- s1.min else s1.d.1 <- s1.d.1

        s1.d.2 <- as.numeric(SSfpl((sum(pop.init.1*sd2[8:age.max])*2),s1.param.2[1],s1.param.2[2],s1.param.2[3],s1.param.2[4]))
            if (s1.d.2 < s1.max) s1.d.2 <- s1.max else s1.d.2 <- s1.d.2
	    if (s1.d.2 > s1.min) s1.d.2 <- s1.min else s1.d.2 <- s1.d.2

        s1.d.3 <- as.numeric(SSfpl((sum(pop.init.1*sd3[8:age.max])*2),s1.param.3[1],s1.param.3[2],s1.param.3[3],s1.param.3[4]))
            if (s1.d.3 < s1.max) s1.d.3 <- s1.max else s1.d.3 <- s1.d.3
	    if (s1.d.3 > s1.min) s1.d.3 <- s1.min else s1.d.3 <- s1.d.3

        s2.d.1 <- as.numeric(SSfpl(pop.init.1,s2.param.1[1],s2.param.1[2],s2.param.1[3],s2.param.1[4]))
            if (s2.d.1 < s2.max) s2.d.1 <- s2.max else s2.d.1 <- s2.d.1
	    if (s2.d.1 > s2.min) s2.d.1 <- s2.min else s2.d.1 <- s2.d.1

        s2.d.2 <- as.numeric(SSfpl(pop.init.2,s2.param.2[1],s2.param.2[2],s2.param.2[3],s2.param.2[4]))
            if (s2.d.2 < s2.max) s2.d.2 <- s2.max else s2.d.2 <- s2.d.2
	    if (s2.d.2 > s2.min) s2.d.2 <- s2.min else s2.d.2 <- s2.d.2

        s2.d.3 <- as.numeric(SSfpl(pop.init.3,s2.param.3[1],s2.param.3[2],s2.param.3[3],s2.param.3[4]))
            if (s2.d.3 < s2.max) s2.d.3 <- s2.max else s2.d.3 <- s2.d.3
	    if (s2.d.3 > s2.min) s2.d.3 <- s2.min else s2.d.3 <- s2.d.3

        s3.d.1 <- as.numeric(SSfpl(pop.init.1,s3.param.1[1],s3.param.1[2],s3.param.1[3],s3.param.1[4]))
            if (s3.d.1 < s3.max) s3.d.1 <- s3.max else s3.d.1 <- s3.d.1
	    if (s3.d.1 > s3.min) s3.d.1 <- s3.min else s3.d.1 <- s3.d.1

        s3.d.2 <- as.numeric(SSfpl(pop.init.2,s3.param.2[1],s3.param.2[2],s3.param.2[3],s3.param.2[4]))
            if (s3.d.2 < s3.max) s3.d.2 <- s3.max else s3.d.2 <- s3.d.2
	    if (s3.d.2 > s3.min) s3.d.2 <- s3.min else s3.d.2 <- s3.d.2

        s3.d.3 <- as.numeric(SSfpl(pop.init.3,s3.param.3[1],s3.param.3[2],s3.param.3[3],s3.param.3[4]))
            if (s3.d.3 < s3.max) s3.d.3 <- s3.max else s3.d.3 <- s3.d.3
	    if (s3.d.3 > s3.min) s3.d.3 <- s3.min else s3.d.3 <- s3.d.3
        
## define start proportion females breeding
        pfb.d.1 <- as.numeric(SSfpl(pop.init.1,pfb.param.1[1],pfb.param.1[2],pfb.param.1[3],pfb.param.1[4]))
            if (pfb.d.1 < pfb.max) pfb.d.1 <- pfb.max else pfb.d.1 <- pfb.d.1
	    if (pfb.d.1 > pfb.min) pfb.d.1 <- pfb.min else pfb.d.1 <- pfb.d.1

        pfb.d.2 <- as.numeric(SSfpl(pop.init.2,pfb.param.2[1],pfb.param.2[2],pfb.param.2[3],pfb.param.2[4]))
            if (pfb.d.2 < pfb.max) pfb.d.2 <- pfb.max else pfb.d.2 <- pfb.d.2
	    if (pfb.d.2 > pfb.min) pfb.d.2 <- pfb.min else pfb.d.2 <- pfb.d.2

        pfb.d.3 <- as.numeric(SSfpl(pop.init.3,pfb.param.3[1],pfb.param.3[2],pfb.param.3[3],pfb.param.3[4]))
            if (pfb.d.3 < pfb.max) pfb.d.3 <- pfb.max else pfb.d.3 <- pfb.d.3
	    if (pfb.d.3 > pfb.min) pfb.d.3 <- pfb.min else pfb.d.3 <- pfb.d.3

## New survival vector
s.vec.d.1 <- sf.vec ## Region 1
s.vec.d.1[1] <- s0.d.1; s.vec.d.1[2] <- s1.d.1; s.vec.d.1[3] <- s2.d.1; s.vec.d.1[4] <- s3.d.1

s.vec.d.2 <- sf.vec ## Region 2
s.vec.d.2[1] <- s0.d.2; s.vec.d.2[2] <- s1.d.2; s.vec.d.2[3] <- s2.d.2; s.vec.d.2[4] <- s3.d.2

s.vec.d.3 <- sf.vec ## Region 3
s.vec.d.3[1] <- s0.d.3; s.vec.d.3[2] <- s1.d.3; s.vec.d.3[3] <- s2.d.3; s.vec.d.3[4] <- s3.d.3

## Stochastic flooding events - re-calculate seh
seh.rand <- runif(tlimit,min=seh.min,max=seh.max)

## New fecundity vector
##m.vec.d.1 <- c(0,0,0,0,0,0,0,0,0,pr10*seh*s.vec.d.1[1]*hsr.avg*clutch*ef*pfb.d.1,pr11*seh*s.vec.d.1[1]*hsr.avg*clutch*ef*pfb.d.1,seh*s.vec.d.1[1]*hsr.avg*clutch*ef*pfb.d.1) ## Region 1
##m.vec.d.2 <- c(0,0,0,0,0,0,0,0,0,pr10*seh*s.vec.d.2[1]*hsr.avg*clutch*ef*pfb.d.2,pr11*seh*s.vec.d.2[1]*hsr.avg*clutch*ef*pfb.d.2,seh*s.vec.d.2[1]*hsr.avg*clutch*ef*pfb.d.2) ## Region 2
##m.vec.d.3 <- c(0,0,0,0,0,0,0,0,0,pr10*seh*s.vec.d.3[1]*hsr.avg*clutch*ef*pfb.d.3,pr11*seh*s.vec.d.3[1]*hsr.avg*clutch*ef*pfb.d.3,seh*s.vec.d.3[1]*hsr.avg*clutch*ef*pfb.d.3) ## Region 3

## proportional proportional breeding females
prf.vec <- rf.vec / (max(rf.vec))

pfb.d.vec.1 <- prf.vec * pfb.d.1
pfb.d.vec.2 <- prf.vec * pfb.d.2
pfb.d.vec.3 <- prf.vec * pfb.d.3

m.vec.d.1 <- pfb.d.vec.1[1:age.max]*s.vec.d.1[1]*seh*clutch*hsr.avg
m.vec.d.2 <- pfb.d.vec.2[1:age.max]*s.vec.d.2[1]*seh*clutch*hsr.avg
m.vec.d.3 <- pfb.d.vec.3[1:age.max]*s.vec.d.3[1]*seh*clutch*hsr.avg

##total females in population
f.1 <- pop.min.1*sr; f.2 <- pop.min.2*sr; f.3 <- pop.min.3*sr

##total males in population
mal.1 <- pop.min.1*(1-sr); mal.2 <- pop.min.2*(1-sr); mal.3 <- pop.min.3*(1-sr)

## The normal matrices
ad.1 <- matrix(data<-0,nrow<-k,ncol<-k); ad.2 <- matrix(data<-0,nrow<-k,ncol<-k); ad.3 <- matrix(data<-0,nrow<-k,ncol<-k)

## Add survival vectors to matrix
## Region 1
diag(ad.1[2:age.max,1:(age.max-1)]) <- s.vec.d.1[2:age.max] ## adds s.vec to female quadrant
ad.1[age.max,age.max] <- s.vec.d.1[age.max+1] ## put final survival value in female age.max/ag.max cell
diag(ad.1[(age.max+2):k,(age.max+1):(k-1)]) <- s.vec.d.1[2:age.max] ## adds s.vec to male quadrant
ad.1[k,k] <- s.vec.d.1[age.max+1] ## put final survival value in male age.max/ag.max cell
ad.1[1,1:(age.max)] <- m.vec.d.1[1:(age.max)] ## adds female fecundity
ad.1[(age.max+1),1:(age.max)] <- m.vec.d.1[1:(age.max)] ## adds male fecundity

## Region 2
diag(ad.2[2:age.max,1:(age.max-1)]) <- s.vec.d.2[2:age.max] ## adds s.vec to female quadrant
ad.2[age.max,age.max] <- s.vec.d.2[age.max+1] ## put final survival value in female age.max/ag.max cell
diag(ad.2[(age.max+2):k,(age.max+1):(k-1)]) <- s.vec.d.2[2:age.max] ## adds s.vec to male quadrant
ad.2[k,k] <- s.vec.d.2[age.max+1] ## put final survival value in male age.max/ag.max cell
ad.2[1,1:(age.max)] <- m.vec.d.2[1:(age.max)] ## adds female fecundity
ad.2[(age.max+1),1:(age.max)] <- m.vec.d.2[1:(age.max)] ## adds male fecundity

## Region 3
diag(ad.3[2:age.max,1:(age.max-1)]) <- s.vec.d.3[2:age.max] ## adds s.vec to female quadrant
ad.3[age.max,age.max] <- s.vec.d.3[age.max+1] ## put final survival value in female age.max/ag.max cell
diag(ad.3[(age.max+2):k,(age.max+1):(k-1)]) <- s.vec.d.3[2:age.max] ## adds s.vec to male quadrant
ad.3[k,k] <- s.vec.d.3[age.max+1] ## put final survival value in male age.max/ag.max cell
ad.3[1,1:(age.max)] <- m.vec.d.3[1:(age.max)] ## adds female fecundity
ad.3[(age.max+1),1:(age.max)] <- m.vec.d.3[1:(age.max)] ## adds male fecundity

## eigenvalues & eigenvectors
dd.1 <- eigen(ad.1)$values; wd.1 <- eigen(ad.1)$vectors
dd.2 <- eigen(ad.2)$values; wd.2 <- eigen(ad.2)$vectors
dd.3 <- eigen(ad.3)$values; wd.3 <- eigen(ad.3)$vectors

## max lambda
maxlambdad.1 <- max.lambda(ad.1); maxlambdad.2 <- max.lambda(ad.2); maxlambdad.3 <- max.lambda(ad.3)

## max r
rd.1 <- max.r(ad.1); rd.2 <- max.r(ad.2); rd.3 <- max.r(ad.3)

## Stable stage distribution
ssdd.1 <- stable.stage.dist(ad.1); ssdd.2 <- stable.stage.dist(ad.2); ssdd.3 <- stable.stage.dist(ad.3)

## Plot start, end and stable sd together
## rowa <- 3
## cola <- 3
## par(mfrow=c(rowa,cola))
## plot(st.sd1[1:age.max],type="l",xlab="age class",ylab="proportion",ylim=range(c(0,0.25)))
## title(main="Region 1")
## plot(st.sd2[1:age.max],type="l",xlab="age class",ylab="proportion",ylim=range(c(0,0.25)))
## title(main="Region 2")
## plot(st.sd3[1:age.max],type="l",xlab="age class",ylab="proportion",ylim=range(c(0,0.25)))
## title(main="Region 3")
## plot(end.sd1[1:age.max],type="l",xlab="age class",ylab="proportion",ylim=range(c(0,0.25)))
## plot(end.sd2[1:age.max],type="l",xlab="age class",ylab="proportion",ylim=range(c(0,0.25)))
## plot(end.sd3[1:age.max],type="l",xlab="age class",ylab="proportion",ylim=range(c(0,0.25)))
## plot(ssdd.1[1:age.max],type="l",xlab="age class",ylab="proportion",ylim=range(c(0,0.25)))
## plot(ssdd.2[1:age.max],type="l",xlab="age class",ylab="proportion",ylim=range(c(0,0.25)))
## plot(ssdd.3[1:age.max],type="l",xlab="age class",ylab="proportion",ylim=range(c(0,0.25)))
## par(mfrow=c(1,2))

## pick initial vectors
if (sd.init == 1) nd.1 <- (st.sd1*Nd.init*nprop1) else nd.1 <- (end.sd1*Nd.init*nprop1)
if (sd.init == 1) nd.2 <- (st.sd2*Nd.init*nprop2) else nd.2 <- (end.sd2*Nd.init*nprop2)
if (sd.init == 1) nd.3 <- (st.sd3*Nd.init*nprop3) else nd.3 <- (end.sd3*Nd.init*nprop3)

## Use predicted end stage distribution if chosen (calculated from 1978-2004 prediction)
pend.sd1 <- c(0.185317428,0.048232706,0.033799522,0.024016873,0.018890421,0.015627349,0.013520960,0.012048097,0.011243499,0.010418121,0.010036887,0.009558343,0.204583673)
pend.sd2 <- c(0.186060407,0.047706299,0.033141949,0.023287064,0.018470048,0.015476684,0.013615699,0.012270972,0.011412802,0.010425666,0.010003002,0.009460751,0.207920786)
pend.sd3 <- c(0.17814605,0.05122357,0.03674009,0.02647407,0.02035241,0.01645651,0.01391335,0.01197071,0.01139597,0.01079933,0.01107524,0.01145758,0.18737123)

if ((sd.init == 2) & (sd.source == 1)) nd.1 <- (pend.sd1*Nd.init*nprop1) else nd.1 <- nd.1
if ((sd.init == 2) & (sd.source == 1)) nd.2 <- (pend.sd2*Nd.init*nprop2) else nd.2 <- nd.2
if ((sd.init == 2) & (sd.source == 1)) nd.3 <- (pend.sd3*Nd.init*nprop3) else nd.3 <- nd.3

##set population size year step vector
popd.vec.1 <- rep(0,tlimit+1); popd.vec.1[1] <- sum(nd.1)
popd.vec.2 <- rep(0,tlimit+1); popd.vec.2[1] <- sum(nd.2)
popd.vec.3 <- rep(0,tlimit+1); popd.vec.3[1] <- sum(nd.3)

##set storage vectors
## Female
## Region 1
n0fd.1.vec <- rep(0,tlimit+1)
n1fd.1.vec <- rep(0,tlimit+1); n1fd.1.vec[1] <- nd.1[1]
n2fd.1.vec <- rep(0,tlimit+1); n2fd.1.vec[1] <- nd.1[2]
n3fd.1.vec <- rep(0,tlimit+1); n3fd.1.vec[1] <- nd.1[3]
n4fd.1.vec <- rep(0,tlimit+1); n4fd.1.vec[1] <- nd.1[4]
n5fd.1.vec <- rep(0,tlimit+1); n5fd.1.vec[1] <- nd.1[5]
n6fd.1.vec <- rep(0,tlimit+1); n6fd.1.vec[1] <- nd.1[6]
n7fd.1.vec <- rep(0,tlimit+1); n7fd.1.vec[1] <- nd.1[7]
n8fd.1.vec <- rep(0,tlimit+1); n8fd.1.vec[1] <- nd.1[8]
n9fd.1.vec <- rep(0,tlimit+1); n9fd.1.vec[1] <- nd.1[9]
n10fd.1.vec <- rep(0,tlimit+1); n10fd.1.vec[1] <- nd.1[10]
n11fd.1.vec <- rep(0,tlimit+1); n11fd.1.vec[1] <- nd.1[11]
n12fd.1.vec <- rep(0,tlimit+1); n12fd.1.vec[1] <- nd.1[12]
n13fd.1.vec <- rep(0,tlimit+1); n13fd.1.vec[1] <- nd.1[13]

## Region 2
n0fd.2.vec <- rep(0,tlimit+1)
n1fd.2.vec <- rep(0,tlimit+1); n1fd.2.vec[1] <- nd.2[1]
n2fd.2.vec <- rep(0,tlimit+1); n2fd.2.vec[1] <- nd.2[2]
n3fd.2.vec <- rep(0,tlimit+1); n3fd.2.vec[1] <- nd.2[3]
n4fd.2.vec <- rep(0,tlimit+1); n4fd.2.vec[1] <- nd.2[4]
n5fd.2.vec <- rep(0,tlimit+1); n5fd.2.vec[1] <- nd.2[5]
n6fd.2.vec <- rep(0,tlimit+1); n6fd.2.vec[1] <- nd.2[6]
n7fd.2.vec <- rep(0,tlimit+1); n7fd.2.vec[1] <- nd.2[7]
n8fd.2.vec <- rep(0,tlimit+1); n8fd.2.vec[1] <- nd.2[8]
n9fd.2.vec <- rep(0,tlimit+1); n9fd.2.vec[1] <- nd.2[9]
n10fd.2.vec <- rep(0,tlimit+1); n10fd.2.vec[1] <- nd.2[10]
n11fd.2.vec <- rep(0,tlimit+1); n11fd.2.vec[1] <- nd.2[11]
n12fd.2.vec <- rep(0,tlimit+1); n12fd.2.vec[1] <- nd.2[12]
n13fd.2.vec <- rep(0,tlimit+1); n13fd.2.vec[1] <- nd.2[13]

## Region 3
n0fd.3.vec <- rep(0,tlimit+1)
n1fd.3.vec <- rep(0,tlimit+1); n1fd.3.vec[1] <- nd.3[1]
n2fd.3.vec <- rep(0,tlimit+1); n2fd.3.vec[1] <- nd.3[2]
n3fd.3.vec <- rep(0,tlimit+1); n3fd.3.vec[1] <- nd.3[3]
n4fd.3.vec <- rep(0,tlimit+1); n4fd.3.vec[1] <- nd.3[4]
n5fd.3.vec <- rep(0,tlimit+1); n5fd.3.vec[1] <- nd.3[5]
n6fd.3.vec <- rep(0,tlimit+1); n6fd.3.vec[1] <- nd.3[6]
n7fd.3.vec <- rep(0,tlimit+1); n7fd.3.vec[1] <- nd.3[7]
n8fd.3.vec <- rep(0,tlimit+1); n8fd.3.vec[1] <- nd.3[8]
n9fd.3.vec <- rep(0,tlimit+1); n9fd.3.vec[1] <- nd.3[9]
n10fd.3.vec <- rep(0,tlimit+1); n10fd.3.vec[1] <- nd.3[10]
n11fd.3.vec <- rep(0,tlimit+1); n11fd.3.vec[1] <- nd.3[11]
n12fd.3.vec <- rep(0,tlimit+1); n12fd.3.vec[1] <- nd.3[12]
n13fd.3.vec <- rep(0,tlimit+1); n13fd.3.vec[1] <- nd.3[13]

## Male
## Region 1
n0md.1.vec <- rep(0,tlimit+1)
n1md.1.vec <- rep(0,tlimit+1); n1md.1.vec[1] <- nd.1[age.max+1]
n2md.1.vec <- rep(0,tlimit+1); n2md.1.vec[1] <- nd.1[age.max+2]
n3md.1.vec <- rep(0,tlimit+1); n3md.1.vec[1] <- nd.1[age.max+3]
n4md.1.vec <- rep(0,tlimit+1); n4md.1.vec[1] <- nd.1[age.max+4]
n5md.1.vec <- rep(0,tlimit+1); n5md.1.vec[1] <- nd.1[age.max+5]
n6md.1.vec <- rep(0,tlimit+1); n6md.1.vec[1] <- nd.1[age.max+6]
n7md.1.vec <- rep(0,tlimit+1); n7md.1.vec[1] <- nd.1[age.max+7]
n8md.1.vec <- rep(0,tlimit+1); n8md.1.vec[1] <- nd.1[age.max+8]
n9md.1.vec <- rep(0,tlimit+1); n9md.1.vec[1] <- nd.1[age.max+9]
n10md.1.vec <- rep(0,tlimit+1); n10md.1.vec[1] <- nd.1[age.max+10]
n11md.1.vec <- rep(0,tlimit+1); n11md.1.vec[1] <- nd.1[age.max+11]
n12md.1.vec <- rep(0,tlimit+1); n12md.1.vec[1] <- nd.1[age.max+12]
n13md.1.vec <- rep(0,tlimit+1); n13md.1.vec[1] <- nd.1[age.max+13]

## Region 2
n0md.2.vec <- rep(0,tlimit+1)
n1md.2.vec <- rep(0,tlimit+1); n1md.2.vec[1] <- nd.2[age.max+1]
n2md.2.vec <- rep(0,tlimit+1); n2md.2.vec[1] <- nd.2[age.max+2]
n3md.2.vec <- rep(0,tlimit+1); n3md.2.vec[1] <- nd.2[age.max+3]
n4md.2.vec <- rep(0,tlimit+1); n4md.2.vec[1] <- nd.2[age.max+4]
n5md.2.vec <- rep(0,tlimit+1); n5md.2.vec[1] <- nd.2[age.max+5]
n6md.2.vec <- rep(0,tlimit+1); n6md.2.vec[1] <- nd.2[age.max+6]
n7md.2.vec <- rep(0,tlimit+1); n7md.2.vec[1] <- nd.2[age.max+7]
n8md.2.vec <- rep(0,tlimit+1); n8md.2.vec[1] <- nd.2[age.max+8]
n9md.2.vec <- rep(0,tlimit+1); n9md.2.vec[1] <- nd.2[age.max+9]
n10md.2.vec <- rep(0,tlimit+1); n10md.2.vec[1] <- nd.2[age.max+10]
n11md.2.vec <- rep(0,tlimit+1); n11md.2.vec[1] <- nd.2[age.max+11]
n12md.2.vec <- rep(0,tlimit+1); n12md.2.vec[1] <- nd.2[age.max+12]
n13md.2.vec <- rep(0,tlimit+1); n13md.2.vec[1] <- nd.2[age.max+13]

## Region 3
n0md.3.vec <- rep(0,tlimit+1)
n1md.3.vec <- rep(0,tlimit+1); n1md.3.vec[1] <- nd.3[age.max+1]
n2md.3.vec <- rep(0,tlimit+1); n2md.3.vec[1] <- nd.3[age.max+2]
n3md.3.vec <- rep(0,tlimit+1); n3md.3.vec[1] <- nd.3[age.max+3]
n4md.3.vec <- rep(0,tlimit+1); n4md.3.vec[1] <- nd.3[age.max+4]
n5md.3.vec <- rep(0,tlimit+1); n5md.3.vec[1] <- nd.3[age.max+5]
n6md.3.vec <- rep(0,tlimit+1); n6md.3.vec[1] <- nd.3[age.max+6]
n7md.3.vec <- rep(0,tlimit+1); n7md.3.vec[1] <- nd.3[age.max+7]
n8md.3.vec <- rep(0,tlimit+1); n8md.3.vec[1] <- nd.3[age.max+8]
n9md.3.vec <- rep(0,tlimit+1); n9md.3.vec[1] <- nd.3[age.max+9]
n10md.3.vec <- rep(0,tlimit+1); n10md.3.vec[1] <- nd.3[age.max+10]
n11md.3.vec <- rep(0,tlimit+1); n11md.3.vec[1] <- nd.3[age.max+11]
n12md.3.vec <- rep(0,tlimit+1); n12md.3.vec[1] <- nd.3[age.max+12]
n13md.3.vec <- rep(0,tlimit+1); n13md.3.vec[1] <- nd.3[age.max+13]

##set year step vector
yrd.vec <- rep(0,tlimit+1)
yrd.vec[1] <- 0

	for (j in 1:tlimit) {
    		yrd.vec[j+1] <- j
	}

## modify year vector for calendar year
yrd.vec <- yrd.vec + st.yr

## Egg harvest information
eh.yr.vec <- seq(1980,2002,1) ## egg harvests started in 1980
eh.vec <- c(135,2758,327,298,2320,3518,3737,4401,5300,6497,12010,9212,15298,12379,17322,19033,29044,21979,10812,13976,11987,15478,17536) ## harvested counts
ve.vec <- c(0,100,0,85,1354,2493,2236,2760,3410,3886,8859,5491,9919,8538,12881,13106,21872,15777,8076,8153,7016,8658,10605) ## viable eggs
pev.vec <- ve.vec/eh.vec ## proportion eggs viable
pev.mean <- mean(pev.vec[5:(length(pev.vec))]) ## mean proportion viable after initial years
pev.lo <- quantile((pev.vec[5:(length(pev.vec))]),probs=0.025)
pev.up <- quantile((pev.vec[5:(length(pev.vec))]),probs=0.975)

h.a1.1.vec <- eh.vec*pharv.1
h.a1.2.vec <- eh.vec*pharv.2
h.a1.3.vec <- eh.vec*pharv.3
h.a1.1.fix <- rep(0,tlimit+1); 	h.a1.2.fix <- rep(0,tlimit+1); 	h.a1.3.fix <- rep(0,tlimit+1)
h.a1.1.fix[3:(tlimit+1)] <- h.a1.1.vec[1:(tlimit-1)]
h.a1.2.fix[3:(tlimit+1)] <- h.a1.2.vec[1:(tlimit-1)]
h.a1.3.fix[3:(tlimit+1)] <- h.a1.3.vec[1:(tlimit-1)]

h.a1.sub.1 <- which(is.na(h.a1.1.fix) == "TRUE"); h.a1.sub.2 <- which(is.na(h.a1.2.fix) == "TRUE"); h.a1.sub.3 <- which(is.na(h.a1.3.fix) == "TRUE")
if (length(h.a1.sub.1) > 0) h.a1.1.fix[h.a1.sub.1] <- h.a1.1.fix[(min(h.a1.sub.1)-1)]
if (length(h.a1.sub.2) > 0) h.a1.2.fix[h.a1.sub.2] <- h.a1.2.fix[(min(h.a1.sub.2)-1)]
if (length(h.a1.sub.3) > 0) h.a1.3.fix[h.a1.sub.3] <- h.a1.3.fix[(min(h.a1.sub.3)-1)]

eggs.vec.1 <- rep(0,tlimit+1); eggs.vec.2 <- rep(0,tlimit+1); eggs.vec.3 <- rep(0,tlimit+1)

## Egg harvest plots
##row <- 2
##col <- 2
##par(mfrow=c(row,col))
##plot(eh.yr.vec,eh.vec,xlab="year",ylab="eggs",type="l")
##lines(eh.yr.vec,ve.vec,col="green")
##plot(eh.yr.vec,pev.vec,xlab="year",ylab="proportion viable eggs",type="l",col="green")
##abline(h=pev.mean,col="black")
##abline(h=pev.lo,col="red")
##abline(h=pev.up,col="red")
##par(mfrow=c(1,2))

##then iterate
	for (ii in 1:tlimit) {
		nd.1 <- ad.1 %*% nd.1; nd.2 <- ad.2 %*% nd.2; nd.3 <- ad.3 %*% nd.3

	## Set any negative values in nds to 1
	zero.sub.1 <- as.numeric(which(nd.1<1)); nd.1[zero.sub.1] <- 1
	zero.sub.2 <- as.numeric(which(nd.2<1)); nd.2[zero.sub.2] <- 1
	zero.sub.3 <- as.numeric(which(nd.3<1)); nd.3[zero.sub.3] <- 1

	## Set any NAs in the nd vectors to 1
	na.sub.1 <- which(is.na(nd.1)); na.sub.2 <- which(is.na(nd.2)); na.sub.3 <- which(is.na(nd.3))
	nd.1[na.sub.1] <- 1; nd.2[na.sub.2] <- 1; nd.3[na.sub.3] <- 1

	## Calculate severity of a potential catastrophe
	sev.cat.input <- runif(1,0,(max(sev.cat.prob)))
	sev.cat <- ((sev.cat.coeff[1]*(sev.cat.input^sev.cat.coeff[2])))/100

	if (sev.cat > 1) sev.cat <- 0.99

	## calculate whether a catastrophe occurs
	if (runif(1,0,1) <= p.cat) nd.1 <- nd.1*(1-sev.cat) else nd.1 <- nd.1
	if (runif(1,0,1) <= p.cat) nd.2 <- nd.2*(1-sev.cat) else nd.2 <- nd.2
	if (runif(1,0,1) <= p.cat) nd.3 <- nd.3*(1-sev.cat) else nd.3 <- nd.3

	## Set any negative values in nds to 1
	zero.sub.1 <- as.numeric(which(nd.1<1)); nd.1[zero.sub.1] <- 1
	zero.sub.2 <- as.numeric(which(nd.2<1)); nd.2[zero.sub.2] <- 1
	zero.sub.3 <- as.numeric(which(nd.3<1)); nd.3[zero.sub.3] <- 1

	## Set any NAs in the nd vectors to 1
	na.sub.1 <- which(is.na(nd.1)); na.sub.2 <- which(is.na(nd.2)); na.sub.3 <- which(is.na(nd.3))
	nd.1[na.sub.1] <- 1; nd.2[na.sub.2] <- 1; nd.3[na.sub.3] <- 1

	## Put values into storage vectors
	## Total population
	popd.vec.1[ii+1] <- sum(nd.1); popd.vec.2[ii+1] <- sum(nd.2); popd.vec.3[ii+1] <- sum(nd.3)

	## Females
	n0fd.1.vec[ii+1] <- (hsr.avg*seh*ef*clutch*pfb.d.1*(n12fd.1.vec[ii])) ## only calculates for 12-year-old reproductive output (hatchlings produced)
	n0fd.2.vec[ii+1] <- (hsr.avg*seh*ef*clutch*pfb.d.2*(n12fd.2.vec[ii])) ## only calculates for 12-year-old reproductive output (hatchlings produced)
	n0fd.3.vec[ii+1] <- (hsr.avg*seh*ef*clutch*pfb.d.3*(n12fd.3.vec[ii])) ## only calculates for 12-year-old reproductive output (hatchlings produced)
    	n1fd.1.vec[ii+1] <- nd.1[1]; n1fd.2.vec[ii+1] <- nd.2[1]; n1fd.3.vec[ii+1] <- nd.3[1]
    	n2fd.1.vec[ii+1] <- nd.1[2]; n2fd.2.vec[ii+1] <- nd.2[2]; n2fd.3.vec[ii+1] <- nd.3[2]
    	n3fd.1.vec[ii+1] <- nd.1[3]; n3fd.2.vec[ii+1] <- nd.2[3]; n3fd.3.vec[ii+1] <- nd.3[3]
    	n4fd.1.vec[ii+1] <- nd.1[4]; n4fd.2.vec[ii+1] <- nd.2[4]; n4fd.3.vec[ii+1] <- nd.3[4]
    	n5fd.1.vec[ii+1] <- nd.1[5]; n5fd.2.vec[ii+1] <- nd.2[5]; n5fd.3.vec[ii+1] <- nd.3[5]
    	n6fd.1.vec[ii+1] <- nd.1[6]; n6fd.2.vec[ii+1] <- nd.2[6]; n6fd.3.vec[ii+1] <- nd.3[6]
    	n7fd.1.vec[ii+1] <- nd.1[7]; n7fd.2.vec[ii+1] <- nd.2[7]; n7fd.3.vec[ii+1] <- nd.3[7]
    	n8fd.1.vec[ii+1] <- nd.1[8]; n8fd.2.vec[ii+1] <- nd.2[8]; n8fd.3.vec[ii+1] <- nd.3[8]
    	n9fd.1.vec[ii+1] <- nd.1[9]; n9fd.2.vec[ii+1] <- nd.2[9]; n9fd.3.vec[ii+1] <- nd.3[9]
    	n10fd.1.vec[ii+1] <- nd.1[10]; n10fd.2.vec[ii+1] <- nd.2[10]; n10fd.3.vec[ii+1] <- nd.3[10]
    	n11fd.1.vec[ii+1] <- nd.1[11]; n11fd.2.vec[ii+1] <- nd.2[11]; n11fd.3.vec[ii+1] <- nd.3[11]
    	n12fd.1.vec[ii+1] <- nd.1[12]; n12fd.2.vec[ii+1] <- nd.2[12]; n12fd.3.vec[ii+1] <- nd.3[12]
    	n13fd.1.vec[ii+1] <- nd.1[13]; n13fd.2.vec[ii+1] <- nd.2[13]; n13fd.3.vec[ii+1] <- nd.3[13]

	## Males
	n0md.1.vec[ii+1] <- (hsr.avg*seh*ef*clutch*pfb.d.1*(n12fd.1.vec[ii])) ## only calculates for 12-year-old reproductive output (hatchlings produced)
	n0md.2.vec[ii+1] <- (hsr.avg*seh*ef*clutch*pfb.d.2*(n12fd.2.vec[ii])) ## only calculates for 12-year-old reproductive output (hatchlings produced)
	n0md.3.vec[ii+1] <- (hsr.avg*seh*ef*clutch*pfb.d.3*(n12fd.3.vec[ii])) ## only calculates for 12-year-old reproductive output (hatchlings produced)
    	n1md.1.vec[ii+1] <- nd.1[age.max+1]; n1md.2.vec[ii+1] <- nd.2[age.max+1]; n1md.3.vec[ii+1] <- nd.3[age.max+1]
    	n2md.1.vec[ii+1] <- nd.1[age.max+2]; n2md.2.vec[ii+1] <- nd.2[age.max+2]; n2md.3.vec[ii+1] <- nd.3[age.max+2]
    	n3md.1.vec[ii+1] <- nd.1[age.max+3]; n3md.2.vec[ii+1] <- nd.2[age.max+3]; n3md.3.vec[ii+1] <- nd.3[age.max+3]
    	n4md.1.vec[ii+1] <- nd.1[age.max+4]; n4md.2.vec[ii+1] <- nd.2[age.max+4]; n4md.3.vec[ii+1] <- nd.3[age.max+4]
    	n5md.1.vec[ii+1] <- nd.1[age.max+5]; n5md.2.vec[ii+1] <- nd.2[age.max+5]; n5md.3.vec[ii+1] <- nd.3[age.max+5]
    	n6md.1.vec[ii+1] <- nd.1[age.max+6]; n6md.2.vec[ii+1] <- nd.2[age.max+6]; n6md.3.vec[ii+1] <- nd.3[age.max+6]
    	n7md.1.vec[ii+1] <- nd.1[age.max+7]; n7md.2.vec[ii+1] <- nd.2[age.max+7]; n7md.3.vec[ii+1] <- nd.3[age.max+7]
    	n8md.1.vec[ii+1] <- nd.1[age.max+8]; n8md.2.vec[ii+1] <- nd.2[age.max+8]; n8md.3.vec[ii+1] <- nd.3[age.max+8]
    	n9md.1.vec[ii+1] <- nd.1[age.max+9]; n9md.2.vec[ii+1] <- nd.2[age.max+9]; n9md.3.vec[ii+1] <- nd.3[age.max+9]
    	n10md.1.vec[ii+1] <- nd.1[age.max+10]; n10md.2.vec[ii+1] <- nd.2[age.max+10]; n10md.3.vec[ii+1] <- nd.3[age.max+10]
    	n11md.1.vec[ii+1] <- nd.1[age.max+11]; n11md.2.vec[ii+1] <- nd.2[age.max+11]; n11md.3.vec[ii+1] <- nd.3[age.max+11]
    	n12md.1.vec[ii+1] <- nd.1[age.max+12]; n12md.2.vec[ii+1] <- nd.2[age.max+12]; n12md.3.vec[ii+1] <- nd.3[age.max+12]
    	n13md.1.vec[ii+1] <- nd.1[age.max+13]; n13md.2.vec[ii+1] <- nd.2[age.max+13]; n13md.3.vec[ii+1] <- nd.3[age.max+13]
    
        ## Set negative density feedback function for the matrix

        ## Redefine survival probabilities
        s0.d.1 <- as.numeric(SSfpl((nd.1[age.max]*seh*clutch*hsr.avg*ef*pfb.min),s0.param.1[1],s0.param.1[2],s0.param.1[3],s0.param.1[4]))
            if (s0.d.1 < s0.max) s0.d.1 <- s0.max else s0.d.1 <- s0.d.1
	    if (s0.d.1 > s0.min) s0.d.1 <- s0.min else s0.d.1 <- s0.d.1

        s0.d.2 <- as.numeric(SSfpl((nd.2[age.max]*seh*clutch*hsr.avg*ef*pfb.min),s0.param.2[1],s0.param.2[2],s0.param.2[3],s0.param.2[4]))
            if (s0.d.2 < s0.max) s0.d.2 <- s0.max else s0.d.2 <- s0.d.2
	    if (s0.d.2 > s0.min) s0.d.2 <- s0.min else s0.d.2 <- s0.d.2

        s0.d.3 <- as.numeric(SSfpl((nd.3[age.max]*seh*clutch*hsr.avg*ef*pfb.min),s0.param.3[1],s0.param.3[2],s0.param.3[3],s0.param.3[4]))
            if (s0.d.3 < s0.max) s0.d.3 <- s0.max else s0.d.3 <- s0.d.3
	    if (s0.d.3 > s0.min) s0.d.3 <- s0.min else s0.d.3 <- s0.d.3

        s1.d.1 <- as.numeric(SSfpl((sum(nd.1[8:age.max])*2),s1.param.1[1],s1.param.1[2],s1.param.1[3],s1.param.1[4]))
            if (s1.d.1 < s1.max) s1.d.1 <- s1.max else s1.d.1 <- s1.d.1
	    if (s1.d.1 > s1.min) s1.d.1 <- s1.min else s1.d.1 <- s1.d.1

        s1.d.2 <- as.numeric(SSfpl((sum(nd.2[8:age.max])*2),s1.param.2[1],s1.param.2[2],s1.param.2[3],s1.param.2[4]))
            if (s1.d.2 < s1.max) s1.d.2 <- s1.max else s1.d.2 <- s1.d.2
	    if (s1.d.2 > s1.min) s1.d.2 <- s1.min else s1.d.2 <- s1.d.2

        s1.d.3 <- as.numeric(SSfpl((sum(nd.3[8:age.max])*2),s1.param.3[1],s1.param.3[2],s1.param.3[3],s1.param.3[4]))
            if (s1.d.3 < s1.max) s1.d.3 <- s1.max else s1.d.3 <- s1.d.3
	    if (s1.d.3 > s1.min) s1.d.3 <- s1.min else s1.d.3 <- s1.d.3

        s2.d.1 <- as.numeric(SSfpl(sum(nd.1),s2.param.1[1],s2.param.1[2],s2.param.1[3],s2.param.1[4]))
            if (s2.d.1 < s2.max) s2.d.1 <- s2.max else s2.d.1 <- s2.d.1
	    if (s2.d.1 > s2.min) s2.d.1 <- s2.min else s2.d.1 <- s2.d.1

        s2.d.2 <- as.numeric(SSfpl(sum(nd.2),s2.param.2[1],s2.param.2[2],s2.param.2[3],s2.param.2[4]))
            if (s2.d.2 < s2.max) s2.d.2 <- s2.max else s2.d.2 <- s2.d.2
	    if (s2.d.2 > s2.min) s2.d.2 <- s2.min else s2.d.2 <- s2.d.2

        s2.d.3 <- as.numeric(SSfpl(sum(nd.3),s2.param.3[1],s2.param.3[2],s2.param.3[3],s2.param.3[4]))
            if (s2.d.3 < s2.max) s2.d.3 <- s2.max else s2.d.3 <- s2.d.3
	    if (s2.d.3 > s2.min) s2.d.3 <- s2.min else s2.d.3 <- s2.d.3

        s3.d.1 <- as.numeric(SSfpl(sum(nd.1),s3.param.1[1],s3.param.1[2],s3.param.1[3],s3.param.1[4]))
            if (s3.d.1 < s3.max) s3.d.1 <- s3.max else s3.d.1 <- s3.d.1
	    if (s3.d.1 > s3.min) s3.d.1 <- s3.min else s3.d.1 <- s3.d.1

        s3.d.2 <- as.numeric(SSfpl(sum(nd.2),s3.param.2[1],s3.param.2[2],s3.param.2[3],s3.param.2[4]))
            if (s3.d.2 < s3.max) s3.d.2 <- s3.max else s3.d.2 <- s3.d.2
	    if (s3.d.2 > s3.min) s3.d.2 <- s3.min else s3.d.2 <- s3.d.2

        s3.d.3 <- as.numeric(SSfpl(sum(nd.3),s3.param.3[1],s3.param.3[2],s3.param.3[3],s3.param.3[4]))
            if (s3.d.3 < s3.max) s3.d.3 <- s3.max else s3.d.3 <- s3.d.3
	    if (s3.d.3 > s3.min) s3.d.3 <- s3.min else s3.d.3 <- s3.d.3
        
	## re-define proportion females breeding
        pfb.d.1 <- as.numeric(SSfpl(sum(nd.1),pfb.param.1[1],pfb.param.1[2],pfb.param.1[3],pfb.param.1[4]))
            if (pfb.d.1 < pfb.max) pfb.d.1 <- pfb.max else pfb.d.1 <- pfb.d.1
	    if (pfb.d.1 > pfb.min) pfb.d.1 <- pfb.min else pfb.d.1 <- pfb.d.1

        pfb.d.2 <- as.numeric(SSfpl(sum(nd.2),pfb.param.2[1],pfb.param.2[2],pfb.param.2[3],pfb.param.2[4]))
            if (pfb.d.2 < pfb.max) pfb.d.2 <- pfb.max else pfb.d.2 <- pfb.d.2
	    if (pfb.d.2 > pfb.min) pfb.d.2 <- pfb.min else pfb.d.2 <- pfb.d.2

        pfb.d.3 <- as.numeric(SSfpl(sum(nd.3),pfb.param.3[1],pfb.param.3[2],pfb.param.3[3],pfb.param.3[4]))
            if (pfb.d.3 < pfb.max) pfb.d.3 <- pfb.max else pfb.d.3 <- pfb.d.3
	    if (pfb.d.3 > pfb.min) pfb.d.3 <- pfb.min else pfb.d.3 <- pfb.d.3

	## New survival vector
	s.vec.d.1 <- sf.vec ## Region 1
	s.vec.d.1[1] <- s0.d.1; s.vec.d.1[2] <- s1.d.1; s.vec.d.1[3] <- s2.d.1; s.vec.d.1[4] <- s3.d.1
	s.vec.d.2 <- sf.vec ## Region 2
	s.vec.d.2[1] <- s0.d.2; s.vec.d.2[2] <- s1.d.2; s.vec.d.2[3] <- s2.d.2; s.vec.d.2[4] <- s3.d.2
	s.vec.d.3 <- sf.vec ## Region 3
	s.vec.d.3[1] <- s0.d.3; s.vec.d.3[2] <- s1.d.3; s.vec.d.3[3] <- s2.d.3; s.vec.d.3[4] <- s3.d.3

	## Stochastic seh (flooding loss)
	seh <- ifelse(stoch.select == 2, seh <- seh.rand[ii], seh <- 0.314)

	## New fecundity vector
	##m.vec.d.1 <- c(0,0,0,0,0,0,0,0,0,pr10*seh*s.vec.d.1[1]*hsr.avg*clutch*ef*pfb.d.1,pr11*seh*s.vec.d.1[1]*hsr.avg*clutch*ef*pfb.d.1,seh*s.vec.d.1[1]*hsr.avg*clutch*ef*pfb.d.1) ## Region 1
	##m.vec.d.2 <- c(0,0,0,0,0,0,0,0,0,pr10*seh*s.vec.d.2[1]*hsr.avg*clutch*ef*pfb.d.2,pr11*seh*s.vec.d.2[1]*hsr.avg*clutch*ef*pfb.d.2,seh*s.vec.d.2[1]*hsr.avg*clutch*ef*pfb.d.2) ## Region 2
	##m.vec.d.3 <- c(0,0,0,0,0,0,0,0,0,pr10*seh*s.vec.d.3[1]*hsr.avg*clutch*ef*pfb.d.3,pr11*seh*s.vec.d.3[1]*hsr.avg*clutch*ef*pfb.d.3,seh*s.vec.d.3[1]*hsr.avg*clutch*ef*pfb.d.3) ## Region 3

	pfb.d.vec.1 <- prf.vec * pfb.d.1
	pfb.d.vec.2 <- prf.vec * pfb.d.2
	pfb.d.vec.3 <- prf.vec * pfb.d.3

	m.vec.d.1 <- pfb.d.vec.1[1:age.max]*s.vec.d.1[1]*seh*clutch*hsr.avg
	m.vec.d.2 <- pfb.d.vec.2[1:age.max]*s.vec.d.2[1]*seh*clutch*hsr.avg
	m.vec.d.3 <- pfb.d.vec.3[1:age.max]*s.vec.d.3[1]*seh*clutch*hsr.avg

	## Adjust survival and fertility vectors for demographic stochasticity
	## Poisson-re-sampled fertility
	n.f.1 <- m.vec.d.1*nd.1[1:age.max]; n.f.2 <- m.vec.d.2*nd.2[1:age.max]; n.f.3 <- m.vec.d.3*nd.3[1:age.max]
	md.vec.1 <- sum(m.vec.d.1)*(rpois(age.max,n.f.1)/sum(n.f.1))
	md.vec.2 <- sum(m.vec.d.2)*(rpois(age.max,n.f.2)/sum(n.f.2))
	md.vec.3 <- sum(m.vec.d.3)*(rpois(age.max,n.f.3)/sum(n.f.3))

	## Binomial-re-sampled survival
	n.s0.1 <- round(((nd.1[13])*hsr.avg*clutch*seh*pfb.d.vec.1[13])+((nd.1[12])*hsr.avg*clutch*seh*pfb.d.vec.1[12])+((nd.1[11])*hsr.avg*clutch*seh*pfb.d.vec.1[11])+((nd.1[10])*hsr.avg*clutch*seh*pfb.d.vec.1[10])+((nd.1[9])*hsr.avg*clutch*seh*pfb.d.vec.1[9])+((nd.1[8])*hsr.avg*clutch*seh*pfb.d.vec.1[8]))
	n.s0.2 <- round(((nd.2[13])*hsr.avg*clutch*seh*pfb.d.vec.2[13])+((nd.2[12])*hsr.avg*clutch*seh*pfb.d.vec.2[12])+((nd.2[11])*hsr.avg*clutch*seh*pfb.d.vec.2[11])+((nd.2[10])*hsr.avg*clutch*seh*pfb.d.vec.2[10])+((nd.2[9])*hsr.avg*clutch*seh*pfb.d.vec.2[9])+((nd.2[8])*hsr.avg*clutch*seh*pfb.d.vec.2[8]))
	n.s0.3 <- round(((nd.3[13])*hsr.avg*clutch*seh*pfb.d.vec.3[13])+((nd.3[12])*hsr.avg*clutch*seh*pfb.d.vec.3[12])+((nd.3[11])*hsr.avg*clutch*seh*pfb.d.vec.3[11])+((nd.3[10])*hsr.avg*clutch*seh*pfb.d.vec.3[10])+((nd.3[9])*hsr.avg*clutch*seh*pfb.d.vec.3[9])+((nd.3[8])*hsr.avg*clutch*seh*pfb.d.vec.3[8]))
	s.vec.d.1[1] <- rbinom(1,n.s0.1,s.vec.d.1[1])/n.s0.1
	s.vec.d.2[1] <- rbinom(1,n.s0.2,s.vec.d.2[1])/n.s0.2
	s.vec.d.3[1] <- rbinom(1,n.s0.3,s.vec.d.3[1])/n.s0.3

	n.s1.1 <- round(nd.1[1]); n.s1.2 <- round(nd.2[1]); n.s1.3 <- round(nd.3[1])
	s.vec.d.1[2] <- rbinom(1,n.s1.1,s.vec.d.1[2])/n.s1.1
	s.vec.d.2[2] <- rbinom(1,n.s1.2,s.vec.d.2[2])/n.s1.2
	s.vec.d.3[2] <- rbinom(1,n.s1.3,s.vec.d.3[2])/n.s1.3

	n.s2.1 <- round(nd.1[2]); n.s2.2 <- round(nd.2[2]); n.s2.3 <- round(nd.3[2])
	s.vec.d.1[3] <- rbinom(1,n.s2.1,s.vec.d.1[3])/n.s2.1
	s.vec.d.2[3] <- rbinom(1,n.s2.2,s.vec.d.2[3])/n.s2.2
	s.vec.d.3[3] <- rbinom(1,n.s2.3,s.vec.d.3[3])/n.s2.3

	n.s3.1 <- round(nd.1[3]); n.s3.2 <- round(nd.2[3]); n.s3.3 <- round(nd.3[3])
	s.vec.d.1[4] <- rbinom(1,n.s3.1,s.vec.d.1[4])/n.s3.1
	s.vec.d.2[4] <- rbinom(1,n.s3.2,s.vec.d.2[4])/n.s3.2
	s.vec.d.3[4] <- rbinom(1,n.s3.3,s.vec.d.3[4])/n.s3.3


	## ****************************************************************
	## nd vectors updated for egg harvest
	## ****************************************************************

	## Re-adjust h.a1 vectors for time-series used; 1-year olds that would have resulted from harvest eggs
	if (sd.init == 1) h.a1.1 <- ef*seh*s.vec.d.1[1]*h.a1.1.fix[ii+1] else h.a1.1 <- ef*seh*s.vec.d.1[1]*(ptharv.1*((n8fd.1.vec[ii]*clutch*pfb.d.vec.1[8])+(n9fd.1.vec[ii]*clutch*pfb.d.vec.1[9])+(n10fd.1.vec[ii]*clutch*pfb.d.vec.1[10])+(n11fd.1.vec[ii]*clutch*pfb.d.vec.1[11])+(n12fd.1.vec[ii]*clutch*pfb.d.vec.1[12])+(n13fd.1.vec[ii]*clutch*pfb.d.vec.1[13])))
	if (sd.init == 1) h.a1.2 <- ef*seh*s.vec.d.2[1]*h.a1.2.fix[ii+1] else h.a1.2 <- ef*seh*s.vec.d.2[1]*(ptharv.2*((n8fd.2.vec[ii]*clutch*pfb.d.vec.1[8])+(n9fd.2.vec[ii]*clutch*pfb.d.vec.2[9])+(n10fd.2.vec[ii]*clutch*pfb.d.vec.2[10])+(n11fd.2.vec[ii]*clutch*pfb.d.vec.2[11])+(n12fd.2.vec[ii]*clutch*pfb.d.vec.2[12])+(n13fd.2.vec[ii]*clutch*pfb.d.vec.2[13])))
	if (sd.init == 1) h.a1.3 <- ef*seh*s.vec.d.3[1]*h.a1.3.fix[ii+1] else h.a1.3 <- ef*seh*s.vec.d.3[1]*(ptharv.3*((n8fd.3.vec[ii]*clutch*pfb.d.vec.1[8])+(n9fd.3.vec[ii]*clutch*pfb.d.vec.3[9])+(n10fd.3.vec[ii]*clutch*pfb.d.vec.3[10])+(n11fd.3.vec[ii]*clutch*pfb.d.vec.3[11])+(n12fd.3.vec[ii]*clutch*pfb.d.vec.3[12])+(n13fd.3.vec[ii]*clutch*pfb.d.vec.3[13])))

	## updated nd vectors
	if ((nd.1[1]) < (h.a1.1*hsr.avg)) nd.1[1] <- 0 else nd.1[1] <- (nd.1[1] - (h.a1.1*hsr.avg)) ## Adjust nd vector for would-be 1st-years
	if ((nd.1[age.max+1]) < (h.a1.2*(1-hsr.avg))) nd.1[age.max+1] <- 0 else nd.1[age.max+1] <- (nd.1[age.max+1] - (h.a1.1*(1-hsr.avg))) ## Adjust nd vector for would-be 1st-years
	if ((nd.2[1]) < (h.a1.2*hsr.avg)) nd.2[1] <- 0 else nd.2[1] <- (nd.2[1] - (h.a1.2*hsr.avg)) ## Adjust nd vector for would-be 1st-years
	if ((nd.2[age.max+1]) < (h.a1.2*(1-hsr.avg))) nd.2[age.max+1] <- 0 else nd.2[age.max+1] <- (nd.2[age.max+1] - (h.a1.2*(1-hsr.avg))) ## Adjust nd vector for would-be 1st-years
	if ((nd.3[1]) < (h.a1.3*hsr.avg)) nd.3[1] <- 0 else nd.3[1] <- (nd.3[1] - (h.a1.3*hsr.avg)) ## Adjust nd vector for would-be 1st-years
	if ((nd.3[age.max+1]) < (h.a1.3*(1-hsr.avg))) nd.3[age.max+1] <- 0 else nd.3[age.max+1] <- (nd.3[age.max+1] - (h.a1.3*(1-hsr.avg))) ## Adjust nd vector for would-be 1st-years

	## Set any negative values in nds to 1
	zero.sub.1 <- as.numeric(which(nd.1<1)); nd.1[zero.sub.1] <- 1
	zero.sub.2 <- as.numeric(which(nd.2<1)); nd.2[zero.sub.2] <- 1
	zero.sub.3 <- as.numeric(which(nd.3<1)); nd.3[zero.sub.3] <- 1

	## Egg production vectors (estimate proportion of total eggs harvested; retrospective analysis only)
	eggs.1 <- (pfb.d.vec.1[8]*n8fd.1.vec[ii]*clutch)+(pfb.d.vec.1[9]*n9fd.1.vec[ii]*clutch)+(pfb.d.vec.1[10]*n10fd.1.vec[ii]*clutch)+(pfb.d.vec.1[11]*n11fd.1.vec[ii]*clutch)+(pfb.d.vec.1[12]*n12fd.1.vec[ii]*clutch)+(pfb.d.vec.1[13]*n13fd.1.vec[ii]*clutch)
	eggs.2 <- (pfb.d.vec.2[8]*n8fd.2.vec[ii]*clutch)+(pfb.d.vec.2[9]*n9fd.2.vec[ii]*clutch)+(pfb.d.vec.2[10]*n10fd.2.vec[ii]*clutch)+(pfb.d.vec.2[11]*n11fd.2.vec[ii]*clutch)+(pfb.d.vec.2[12]*n12fd.2.vec[ii]*clutch)+(pfb.d.vec.2[13]*n13fd.1.vec[ii]*clutch)
	eggs.3 <- (pfb.d.vec.3[8]*n8fd.3.vec[ii]*clutch)+(pfb.d.vec.3[9]*n9fd.3.vec[ii]*clutch)+(pfb.d.vec.3[10]*n10fd.3.vec[ii]*clutch)+(pfb.d.vec.3[11]*n11fd.3.vec[ii]*clutch)+(pfb.d.vec.3[12]*n12fd.3.vec[ii]*clutch)+(pfb.d.vec.3[13]*n13fd.1.vec[ii]*clutch)

	eggs.vec.1[ii+1] <- eggs.1
	eggs.vec.2[ii+1] <- eggs.2
	eggs.vec.3[ii+1] <- eggs.3

	## ****************************************************************

	###################################################################################################
	## Historic harvest function
	## adjust population vector at each time step by uniform distribution of annual historical harvest
	###################################################################################################

		for (h in 1:k) {
			if (harv.hist.vec[k] > nd.1[k]) nd.1[k] <- nd.1[k] else nd.1[k] <- (nd.1[k] - harv.hist.vec[k])
			if (harv.hist.vec[k] > nd.2[k]) nd.1[k] <- nd.2[k] else nd.2[k] <- (nd.2[k] - harv.hist.vec[k])
			if (harv.hist.vec[k] > nd.3[k]) nd.1[k] <- nd.3[k] else nd.3[k] <- (nd.3[k] - harv.hist.vec[k])
		}

	## Set any negative values in nds to 1
	zero.sub.1 <- as.numeric(which(nd.1<0)); nd.1[zero.sub.1] <- 1
	zero.sub.2 <- as.numeric(which(nd.2<0)); nd.2[zero.sub.2] <- 1
	zero.sub.3 <- as.numeric(which(nd.3<0)); nd.3[zero.sub.3] <- 1

	## Kill functions
	## adjust kill rate by population size if indicated by user
	if (kill.vary == 0) kill.1 <- kill.1 else kill.1 <- sum(nd.1)*(kill.pc/100)
	if (kill.vary == 0) kill.2 <- kill.2 else kill.2 <- sum(nd.2)*(kill.pc/100)
	if (kill.vary == 0) kill.3 <- kill.3 else kill.3 <- sum(nd.3)*(kill.pc/100)

        ## *****************************************************************
        ## Kill adult males only (> 3 m) 13+ -yr olds
        ## *****************************************************************

	if (kill.type == 1) mkill.1 <- kill.1 else mkill.1 <- 0
	if (kill.type == 1) mkill.2 <- kill.2 else mkill.2 <- 0
	if (kill.type == 1) mkill.3 <- kill.3 else mkill.3 <- 0

        if ((nd.1[k]) < mkill.1) nd.1[k] <- 1 else nd.1[k] <- (nd.1[k] - mkill.1) ## Adjust nd vector for kill
        if ((nd.2[k]) < mkill.2) nd.2[k] <- 1 else nd.2[k] <- (nd.2[k] - mkill.2) ## Adjust nd vector for kill
        if ((nd.3[k]) < mkill.3) nd.3[k] <- 1 else nd.3[k] <- (nd.3[k] - mkill.3) ## Adjust nd vector for kill
        
        ## *****************************************************************

        ## *****************************************************************
        ## Kill animals > 2.5 metres (13+ -yr-old females; 9+ -year old males)
        ## *****************************************************************

	if (kill.type == 2) akill.1 <- kill.1 else akill.1 <- 0
	if (kill.type == 2) akill.2 <- kill.2 else akill.2 <- 0
	if (kill.type == 2) akill.3 <- kill.3 else akill.3 <- 0

	n.kill.f1 <- akill.1/2; n.kill.f2 <- akill.2/2; n.kill.f3 <- akill.3/2 ## uniform distribution of kill among males (9-13) & females (13)
        n.kill.m1 <- akill.1/2; n.kill.m2 <- akill.2/2; n.kill.m3 <- akill.3/2 ## uniform distribution of kill among males (9-13) & females (13)
        nm.kill.1 <- n.kill.m1/5; nm.kill.2 <- n.kill.m2/5; nm.kill.3 <- n.kill.m3/5 ## uniform distribution of kill among males (9-13)

        kill.all.vec.1 <- rep(0,k); kill.all.vec.2 <- rep(0,k); kill.all.vec.3 <- rep(0,k)
        kill.all.vec.1[age.max] <- n.kill.f1; kill.all.vec.2[age.max] <- n.kill.f2; kill.all.vec.3[age.max] <- n.kill.f3
	kill.all.vec.1[22:k] <- nm.kill.1; kill.all.vec.2[22:k] <- nm.kill.2; kill.all.vec.3[22:k] <- nm.kill.3

            if ((nd.1[age.max] + sum(nd.1[22:k])) < akill.1) { 
		nd.1[age.max] <- 1; nd.1[22:k] <- 1
	    }
		nd.1 <- (nd.1 - kill.all.vec.1)

            if ((nd.2[age.max] + sum(nd.2[22:k])) < akill.2) {
		nd.2[age.max] <- 1; nd.2[22:k] <- 1
	     }
		nd.2 <- (nd.2 - kill.all.vec.2)

            if ((nd.3[age.max] + sum(nd.3[22:k])) < akill.3) {
		nd.3[age.max] <- 1; nd.3[22:k] <- 1
	    }
		nd.3 <- (nd.3 - kill.all.vec.2)
        
        ## *****************************************************************

        ## *****************************************************************
        ## Kill constant number of juveniles & adults as per NT Croc MP
        ## *****************************************************************

	if (kill.type == 3) adkill.1 <- kill.num.ad/3 else adkill.1 <- 0
	if (kill.type == 3) adkill.2 <- kill.num.ad/3 else adkill.2 <- 0
	if (kill.type == 3) adkill.3 <- kill.num.ad/3 else adkill.3 <- 0

	if (kill.type == 3) jvkill.1 <- kill.num.juv/3 else jvkill.1 <- 0
	if (kill.type == 3) jvkill.2 <- kill.num.juv/3 else jvkill.2 <- 0
	if (kill.type == 3) jvkill.3 <- kill.num.juv/3 else jvkill.3 <- 0

	n.jv.kill.1 <- jvkill.1/(age.max-1); n.jv.kill.2 <- jvkill.2/(age.max-1); n.jv.kill.3 <- jvkill.3/(age.max-1);

      kill.all.vec.1 <- rep(0,k); kill.all.vec.2 <- rep(0,k); kill.all.vec.3 <- rep(0,k)
      kill.all.vec.1[age.max] <- adkill.1/2; kill.all.vec.2[age.max] <- adkill.2/2; kill.all.vec.3[age.max] <- adkill.3/2
      kill.all.vec.1[k] <- adkill.1/2; kill.all.vec.2[k] <- adkill.2/2; kill.all.vec.3[k] <- adkill.3/2
      kill.all.vec.1[1:(age.max-1)] <- n.jv.kill.1/2; kill.all.vec.2[1:(age.max-1)] <- n.jv.kill.2/2; kill.all.vec.3[1:(age.max-1)] <- n.jv.kill.3/2
      kill.all.vec.1[(age.max+1):(k-1)] <- n.jv.kill.1/2; kill.all.vec.2[(age.max+1):(k-1)] <- n.jv.kill.2/2; kill.all.vec.3[(age.max+1):(k-1)] <- n.jv.kill.3/2

	## Set any negative values in nds to 1
	zero.sub.1 <- as.numeric(which(nd.1<1)); nd.1[zero.sub.1] <- 1
	zero.sub.2 <- as.numeric(which(nd.2<1)); nd.2[zero.sub.2] <- 1
	zero.sub.3 <- as.numeric(which(nd.3<1)); nd.3[zero.sub.3] <- 1

        ## *****************************************************************

	## Fill matrix with new survival & fecundity vectors
	## Re-initiliase matrices
	ad.1 <- matrix(data<-0,nrow<-k,ncol<-k)
	ad.2 <- matrix(data<-0,nrow<-k,ncol<-k)
	ad.3 <- matrix(data<-0,nrow<-k,ncol<-k)

	## Region 1
	diag(ad.1[2:age.max,1:(age.max-1)]) <- s.vec.d.1[2:age.max] ## adds s.vec to female quadrant
	ad.1[age.max,age.max] <- s.vec.d.1[age.max+1] ## put final survival value in female age.max/ag.max cell
	diag(ad.1[(age.max+2):k,(age.max+1):(k-1)]) <- s.vec.d.1[2:age.max] ## adds s.vec to male quadrant
	ad.1[k,k] <- s.vec.d.1[age.max+1] ## put final survival value in male age.max/ag.max cell
	ad.1[1,1:(age.max)] <- m.vec.d.1[1:(age.max)] ## adds female fecundity
	ad.1[(age.max+1),1:(age.max)] <- m.vec.d.1[1:(age.max)] ## adds male fecundity

	## Region 2
	diag(ad.2[2:age.max,1:(age.max-1)]) <- s.vec.d.2[2:age.max] ## adds s.vec to female quadrant
	ad.2[age.max,age.max] <- s.vec.d.2[age.max+1] ## put final survival value in female age.max/ag.max cell
	diag(ad.2[(age.max+2):k,(age.max+1):(k-1)]) <- s.vec.d.2[2:age.max] ## adds s.vec to male quadrant
	ad.2[k,k] <- s.vec.d.2[age.max+1] ## put final survival value in male age.max/ag.max cell
	ad.2[1,1:(age.max)] <- m.vec.d.2[1:(age.max)] ## adds female fecundity
	ad.2[(age.max+1),1:(age.max)] <- m.vec.d.2[1:(age.max)] ## adds male fecundity

	## Region 3
	diag(ad.3[2:age.max,1:(age.max-1)]) <- s.vec.d.3[2:age.max] ## adds s.vec to female quadrant
	ad.3[age.max,age.max] <- s.vec.d.3[age.max+1] ## put final survival value in female age.max/ag.max cell
	diag(ad.3[(age.max+2):k,(age.max+1):(k-1)]) <- s.vec.d.3[2:age.max] ## adds s.vec to male quadrant
	ad.3[k,k] <- s.vec.d.3[age.max+1] ## put final survival value in male age.max/ag.max cell
	ad.3[1,1:(age.max)] <- m.vec.d.3[1:(age.max)] ## adds female fecundity
	ad.3[(age.max+1),1:(age.max)] <- m.vec.d.3[1:(age.max)] ## adds male fecundity

    
	} ## end ii loop

## Make iterative plots
row <- 1
col <- 1
par(mfrow=c(row,col))
plot(yrd.vec,popd.vec.1,xlab="year",ylab="N",ylim=range(c(0,(max(cbind(popd.vec.1,popd.vec.2,popd.vec.3))))),type="l")
lines(yrd.vec,popd.vec.2,col="red")
lines(yrd.vec,popd.vec.3,col="green")
par(mfrow=c(1,1))

## Calculate log vectors
log.popd.vec.1<-log10(popd.vec.1); log.popd.vec.2<-log10(popd.vec.2); log.popd.vec.3<-log10(popd.vec.3)

## Population sizes after tlimit years
n.final.1 <- popd.vec.1[tlimit+1]; n.final.2 <- popd.vec.2[tlimit+1]; n.final.3 <- popd.vec.3[tlimit+1]; n.final <- n.final.1 + n.final.2 + n.final.3

## Total population size over all Regions
popd.vec <- popd.vec.1 + popd.vec.2 + popd.vec.3

## store population vectors
popd.1.perm[z,] <- popd.vec.1; popd.2.perm[z,] <- popd.vec.2; popd.3.perm[z,] <- popd.vec.3; popd.perm[z,] <- popd.vec

## Store age-class vectors
n0fd.1.perm[z,] <- n0fd.1.vec; n0fd.2.perm[z,] <- n0fd.2.vec; n0fd.3.perm[z,] <- n0fd.3.vec
n1fd.1.perm[z,] <- n1fd.1.vec; n1fd.2.perm[z,] <- n1fd.2.vec; n1fd.3.perm[z,] <- n1fd.3.vec
n2fd.1.perm[z,] <- n2fd.1.vec; n2fd.2.perm[z,] <- n2fd.2.vec; n2fd.3.perm[z,] <- n2fd.3.vec
n3fd.1.perm[z,] <- n3fd.1.vec; n3fd.2.perm[z,] <- n3fd.2.vec; n3fd.3.perm[z,] <- n3fd.3.vec
n4fd.1.perm[z,] <- n4fd.1.vec; n4fd.2.perm[z,] <- n4fd.2.vec; n4fd.3.perm[z,] <- n4fd.3.vec
n5fd.1.perm[z,] <- n5fd.1.vec; n5fd.2.perm[z,] <- n5fd.2.vec; n5fd.3.perm[z,] <- n5fd.3.vec
n6fd.1.perm[z,] <- n6fd.1.vec; n6fd.2.perm[z,] <- n6fd.2.vec; n6fd.3.perm[z,] <- n6fd.3.vec
n7fd.1.perm[z,] <- n7fd.1.vec; n7fd.2.perm[z,] <- n7fd.2.vec; n7fd.3.perm[z,] <- n7fd.3.vec
n8fd.1.perm[z,] <- n8fd.1.vec; n8fd.2.perm[z,] <- n8fd.2.vec; n8fd.3.perm[z,] <- n8fd.3.vec
n9fd.1.perm[z,] <- n9fd.1.vec; n9fd.2.perm[z,] <- n9fd.2.vec; n9fd.3.perm[z,] <- n9fd.3.vec
n10fd.1.perm[z,] <- n10fd.1.vec; n10fd.2.perm[z,] <- n10fd.2.vec; n10fd.3.perm[z,] <- n10fd.3.vec
n11fd.1.perm[z,] <- n11fd.1.vec; n11fd.2.perm[z,] <- n11fd.2.vec; n11fd.3.perm[z,] <- n11fd.3.vec
n12fd.1.perm[z,] <- n12fd.1.vec; n12fd.2.perm[z,] <- n12fd.2.vec; n12fd.3.perm[z,] <- n12fd.3.vec
n13fd.1.perm[z,] <- n13fd.1.vec; n13fd.2.perm[z,] <- n13fd.2.vec; n13fd.3.perm[z,] <- n13fd.3.vec

n0md.1.perm[z,] <- n0md.1.vec; n0md.2.perm[z,] <- n0md.2.vec; n0md.3.perm[z,] <- n0md.3.vec
n1md.1.perm[z,] <- n1md.1.vec; n1md.2.perm[z,] <- n1md.2.vec; n1md.3.perm[z,] <- n1md.3.vec
n2md.1.perm[z,] <- n2md.1.vec; n2md.2.perm[z,] <- n2md.2.vec; n2md.3.perm[z,] <- n2md.3.vec
n3md.1.perm[z,] <- n3md.1.vec; n3md.2.perm[z,] <- n3md.2.vec; n3md.3.perm[z,] <- n3md.3.vec
n4md.1.perm[z,] <- n4md.1.vec; n4md.2.perm[z,] <- n4md.2.vec; n4md.3.perm[z,] <- n4md.3.vec
n5md.1.perm[z,] <- n5md.1.vec; n5md.2.perm[z,] <- n5md.2.vec; n5md.3.perm[z,] <- n5md.3.vec
n6md.1.perm[z,] <- n6md.1.vec; n6md.2.perm[z,] <- n6md.2.vec; n6md.3.perm[z,] <- n6md.3.vec
n7md.1.perm[z,] <- n7md.1.vec; n7md.2.perm[z,] <- n7md.2.vec; n7md.3.perm[z,] <- n7md.3.vec
n8md.1.perm[z,] <- n8md.1.vec; n8md.2.perm[z,] <- n8md.2.vec; n8md.3.perm[z,] <- n8md.3.vec
n9md.1.perm[z,] <- n9md.1.vec; n9md.2.perm[z,] <- n9md.2.vec; n9md.3.perm[z,] <- n9md.3.vec
n10md.1.perm[z,] <- n10md.1.vec; n10md.2.perm[z,] <- n10md.2.vec; n10md.3.perm[z,] <- n10md.3.vec
n11md.1.perm[z,] <- n11md.1.vec; n11md.2.perm[z,] <- n11md.2.vec; n11md.3.perm[z,] <- n11md.3.vec
n12md.1.perm[z,] <- n12md.1.vec; n12md.2.perm[z,] <- n12md.2.vec; n12md.3.perm[z,] <- n12md.3.vec
n13md.1.perm[z,] <- n13md.1.vec; n13md.2.perm[z,] <- n13md.2.vec; n13md.3.perm[z,] <- n13md.3.vec

## Store egg vectors
eggs.1.perm[z,] <- eggs.vec.1; eggs.2.perm[z,] <- eggs.vec.2; eggs.3.perm[z,] <- eggs.vec.3
h.a1.1.fix.perm[z,] <- h.a1.1.fix; h.a1.2.fix.perm[z,] <- h.a1.2.fix; h.a1.3.fix.perm[z,] <- h.a1.3.fix

print(z)

} ## end z loop

## Make mean end age distribution vectors (proportion females in each age class)
sdend1.1.perm <- n1fd.1.perm / popd.1.perm; sdend1.2.perm <- n1fd.2.perm / popd.2.perm; sdend1.3.perm <- n1fd.3.perm / popd.3.perm
sdend2.1.perm <- n2fd.1.perm / popd.1.perm; sdend2.2.perm <- n2fd.2.perm / popd.2.perm; sdend2.3.perm <- n2fd.3.perm / popd.3.perm
sdend3.1.perm <- n3fd.1.perm / popd.1.perm; sdend3.2.perm <- n3fd.2.perm / popd.2.perm; sdend3.3.perm <- n3fd.3.perm / popd.3.perm
sdend4.1.perm <- n4fd.1.perm / popd.1.perm; sdend4.2.perm <- n4fd.2.perm / popd.2.perm; sdend4.3.perm <- n4fd.3.perm / popd.3.perm
sdend5.1.perm <- n5fd.1.perm / popd.1.perm; sdend5.2.perm <- n5fd.2.perm / popd.2.perm; sdend5.3.perm <- n5fd.3.perm / popd.3.perm
sdend6.1.perm <- n6fd.1.perm / popd.1.perm; sdend6.2.perm <- n6fd.2.perm / popd.2.perm; sdend6.3.perm <- n6fd.3.perm / popd.3.perm
sdend7.1.perm <- n7fd.1.perm / popd.1.perm; sdend7.2.perm <- n7fd.2.perm / popd.2.perm; sdend7.3.perm <- n7fd.3.perm / popd.3.perm
sdend8.1.perm <- n8fd.1.perm / popd.1.perm; sdend8.2.perm <- n8fd.2.perm / popd.2.perm; sdend8.3.perm <- n8fd.3.perm / popd.3.perm
sdend9.1.perm <- n9fd.1.perm / popd.1.perm; sdend9.2.perm <- n9fd.2.perm / popd.2.perm; sdend9.3.perm <- n9fd.3.perm / popd.3.perm
sdend10.1.perm <- n10fd.1.perm / popd.1.perm; sdend10.2.perm <- n10fd.2.perm / popd.2.perm; sdend10.3.perm <- n10fd.3.perm / popd.3.perm
sdend11.1.perm <- n11fd.1.perm / popd.1.perm; sdend11.2.perm <- n11fd.2.perm / popd.2.perm; sdend11.3.perm <- n11fd.3.perm / popd.3.perm
sdend12.1.perm <- n12fd.1.perm / popd.1.perm; sdend12.2.perm <- n12fd.2.perm / popd.2.perm; sdend12.3.perm <- n12fd.3.perm / popd.3.perm
sdend13.1.perm <- n13fd.1.perm / popd.1.perm; sdend13.2.perm <- n13fd.2.perm / popd.2.perm; sdend13.3.perm <- n13fd.3.perm / popd.3.perm

## Add age matrices together to form 4 stages
nst1.1.perm <- n0fd.1.perm + n1fd.1.perm + n1md.1.perm
nst1.2.perm <- n0fd.2.perm + n1fd.2.perm + n1md.2.perm
nst1.3.perm <- n0fd.3.perm + n1fd.3.perm + n1md.3.perm

nst2.1.perm <- n2fd.1.perm + n2md.1.perm
nst2.2.perm <- n2fd.2.perm + n2md.2.perm
nst2.3.perm <- n2fd.3.perm + n2md.3.perm

nst3.1.perm <- n3fd.1.perm + n4fd.1.perm + n5fd.1.perm + n6fd.1.perm + n7fd.1.perm + n3md.1.perm + n4md.1.perm + n5md.1.perm
nst3.2.perm <- n3fd.2.perm + n4fd.2.perm + n5fd.2.perm + n6fd.2.perm + n7fd.2.perm + n3md.2.perm + n4md.2.perm + n5md.2.perm
nst3.3.perm <- n3fd.3.perm + n4fd.3.perm + n5fd.3.perm + n6fd.3.perm + n7fd.3.perm + n3md.3.perm + n4md.3.perm + n5md.3.perm

nst4.1.perm <- n8fd.1.perm + n9fd.1.perm + n10fd.1.perm + n11fd.1.perm + n12fd.1.perm + n13fd.1.perm + n6md.1.perm + n7md.1.perm + n8md.1.perm + n9md.1.perm + n10md.1.perm + n11md.1.perm + n12md.1.perm + n13md.1.perm
nst4.2.perm <- n8fd.2.perm + n9fd.2.perm + n10fd.2.perm + n11fd.2.perm + n12fd.2.perm + n13fd.2.perm + n6md.2.perm + n7md.2.perm + n8md.2.perm + n9md.2.perm + n10md.2.perm + n11md.2.perm + n12md.2.perm + n13md.2.perm
nst4.3.perm <- n8fd.3.perm + n9fd.3.perm + n10fd.3.perm + n11fd.3.perm + n12fd.3.perm + n13fd.3.perm + n6md.3.perm + n7md.3.perm + n8md.3.perm + n9md.3.perm + n10md.3.perm + n11md.3.perm + n12md.3.perm + n13md.3.perm

## Divide age matrices by population matrix
pst1.1.perm <- nst1.1.perm / popd.1.perm
pst1.2.perm <- nst1.2.perm / popd.2.perm
pst1.3.perm <- nst1.3.perm / popd.3.perm

pst2.1.perm <- nst2.1.perm / popd.1.perm
pst2.2.perm <- nst2.2.perm / popd.2.perm
pst2.3.perm <- nst2.3.perm / popd.3.perm

pst3.1.perm <- nst3.1.perm / popd.1.perm
pst3.2.perm <- nst3.2.perm / popd.2.perm
pst3.3.perm <- nst3.3.perm / popd.3.perm

pst4.1.perm <- nst4.1.perm / popd.1.perm
pst4.2.perm <- nst4.2.perm / popd.2.perm
pst4.3.perm <- nst4.3.perm / popd.3.perm

## Mean & confidence intervals over z permutations
popd.mean.1 <- rep(0,tlimit+1); popd.mean.2 <- rep(0,tlimit+1); popd.mean.3 <- rep(0,tlimit+1); popd.mean <- rep(0,tlimit+1)
popd.up.1 <- rep(0,tlimit+1); popd.up.2 <- rep(0,tlimit+1); popd.up.3 <- rep(0,tlimit+1); popd.up <- rep(0,tlimit+1)
popd.lo.1 <- rep(0,tlimit+1); popd.lo.2 <- rep(0,tlimit+1); popd.lo.3 <- rep(0,tlimit+1); popd.lo <- rep(0,tlimit+1)
nst1.mean.1 <- rep(0,tlimit+1); nst1.mean.2 <- rep(0,tlimit+1); nst1.mean.3 <- rep(0,tlimit+1)
nst2.mean.1 <- rep(0,tlimit+1); nst2.mean.2 <- rep(0,tlimit+1); nst2.mean.3 <- rep(0,tlimit+1)
nst3.mean.1 <- rep(0,tlimit+1); nst3.mean.2 <- rep(0,tlimit+1); nst3.mean.3 <- rep(0,tlimit+1)
nst4.mean.1 <- rep(0,tlimit+1); nst4.mean.2 <- rep(0,tlimit+1); nst4.mean.3 <- rep(0,tlimit+1)
nst1.lo.1 <- rep(0,tlimit+1); nst1.lo.2 <- rep(0,tlimit+1); nst1.lo.3 <- rep(0,tlimit+1)
nst2.lo.1 <- rep(0,tlimit+1); nst2.lo.2 <- rep(0,tlimit+1); nst2.lo.3 <- rep(0,tlimit+1)
nst3.lo.1 <- rep(0,tlimit+1); nst3.lo.2 <- rep(0,tlimit+1); nst3.lo.3 <- rep(0,tlimit+1)
nst4.lo.1 <- rep(0,tlimit+1); nst4.lo.2 <- rep(0,tlimit+1); nst4.lo.3 <- rep(0,tlimit+1)
nst1.up.1 <- rep(0,tlimit+1); nst1.up.2 <- rep(0,tlimit+1); nst1.up.3 <- rep(0,tlimit+1)
nst2.up.1 <- rep(0,tlimit+1); nst2.up.2 <- rep(0,tlimit+1); nst2.up.3 <- rep(0,tlimit+1)
nst3.up.1 <- rep(0,tlimit+1); nst3.up.2 <- rep(0,tlimit+1); nst3.up.3 <- rep(0,tlimit+1)
nst4.up.1 <- rep(0,tlimit+1); nst4.up.2 <- rep(0,tlimit+1); nst4.up.3 <- rep(0,tlimit+1)
pst1.mean.1 <- rep(0,tlimit+1); pst1.mean.2 <- rep(0,tlimit+1); pst1.mean.3 <- rep(0,tlimit+1)
pst2.mean.1 <- rep(0,tlimit+1); pst2.mean.2 <- rep(0,tlimit+1); pst2.mean.3 <- rep(0,tlimit+1)
pst3.mean.1 <- rep(0,tlimit+1); pst3.mean.2 <- rep(0,tlimit+1); pst3.mean.3 <- rep(0,tlimit+1)
pst4.mean.1 <- rep(0,tlimit+1); pst4.mean.2 <- rep(0,tlimit+1); pst4.mean.3 <- rep(0,tlimit+1)
pst1.lo.1 <- rep(0,tlimit+1); pst1.lo.2 <- rep(0,tlimit+1); pst1.lo.3 <- rep(0,tlimit+1)
pst2.lo.1 <- rep(0,tlimit+1); pst2.lo.2 <- rep(0,tlimit+1); pst2.lo.3 <- rep(0,tlimit+1)
pst3.lo.1 <- rep(0,tlimit+1); pst3.lo.2 <- rep(0,tlimit+1); pst3.lo.3 <- rep(0,tlimit+1)
pst4.lo.1 <- rep(0,tlimit+1); pst4.lo.2 <- rep(0,tlimit+1); pst4.lo.3 <- rep(0,tlimit+1)
pst1.up.1 <- rep(0,tlimit+1); pst1.up.2 <- rep(0,tlimit+1); pst1.up.3 <- rep(0,tlimit+1)
pst2.up.1 <- rep(0,tlimit+1); pst2.up.2 <- rep(0,tlimit+1); pst2.up.3 <- rep(0,tlimit+1)
pst3.up.1 <- rep(0,tlimit+1); pst3.up.2 <- rep(0,tlimit+1); pst3.up.3 <- rep(0,tlimit+1)
pst4.up.1 <- rep(0,tlimit+1); pst4.up.2 <- rep(0,tlimit+1); pst4.up.3 <- rep(0,tlimit+1)
sdend1.1 <- rep(0,tlimit+1); sdend1.2 <- rep(0,tlimit+1); sdend1.3 <- rep(0,tlimit+1)
sdend2.1 <- rep(0,tlimit+1); sdend2.2 <- rep(0,tlimit+1); sdend2.3 <- rep(0,tlimit+1)
sdend3.1 <- rep(0,tlimit+1); sdend3.2 <- rep(0,tlimit+1); sdend3.3 <- rep(0,tlimit+1)
sdend4.1 <- rep(0,tlimit+1); sdend4.2 <- rep(0,tlimit+1); sdend4.3 <- rep(0,tlimit+1)
sdend5.1 <- rep(0,tlimit+1); sdend5.2 <- rep(0,tlimit+1); sdend5.3 <- rep(0,tlimit+1)
sdend6.1 <- rep(0,tlimit+1); sdend6.2 <- rep(0,tlimit+1); sdend6.3 <- rep(0,tlimit+1)
sdend7.1 <- rep(0,tlimit+1); sdend7.2 <- rep(0,tlimit+1); sdend7.3 <- rep(0,tlimit+1)
sdend8.1 <- rep(0,tlimit+1); sdend8.2 <- rep(0,tlimit+1); sdend8.3 <- rep(0,tlimit+1)
sdend9.1 <- rep(0,tlimit+1); sdend9.2 <- rep(0,tlimit+1); sdend9.3 <- rep(0,tlimit+1)
sdend10.1 <- rep(0,tlimit+1); sdend10.2 <- rep(0,tlimit+1); sdend10.3 <- rep(0,tlimit+1)
sdend11.1 <- rep(0,tlimit+1); sdend11.2 <- rep(0,tlimit+1); sdend11.3 <- rep(0,tlimit+1)
sdend12.1 <- rep(0,tlimit+1); sdend12.2 <- rep(0,tlimit+1); sdend12.3 <- rep(0,tlimit+1)
sdend13.1 <- rep(0,tlimit+1); sdend13.2 <- rep(0,tlimit+1); sdend13.3 <- rep(0,tlimit+1)
eggs.mean.1 <- rep(0,tlimit+1); eggs.mean.2 <- rep(0,tlimit+1); eggs.mean.3 <- rep(0,tlimit+1)
h.a1.1.fix.mean <- rep(0,tlimit+1); h.a1.2.fix.mean <- rep(0,tlimit+1); h.a1.3.fix.mean <- rep(0,tlimit+1)

for (d in 1:(tlimit+1)) {
	popd.mean.1[d] <- mean(popd.1.perm[,d]); popd.mean.2[d] <- mean(popd.2.perm[,d]); popd.mean.3[d] <- mean(popd.3.perm[,d]); popd.mean[d] <- mean(popd.perm[,d])
	popd.up.1[d] <- quantile(popd.1.perm[,d],probs=0.975); popd.up.2[d] <- quantile(popd.2.perm[,d],probs=0.975); popd.up.3[d] <- quantile(popd.3.perm[,d],probs=0.975); popd.up[d] <- quantile(popd.perm[,d],probs=0.975)
	popd.lo.1[d] <- quantile(popd.1.perm[,d],probs=0.025); popd.lo.2[d] <- quantile(popd.2.perm[,d],probs=0.025); popd.lo.3[d] <- quantile(popd.3.perm[,d],probs=0.025); popd.lo[d] <- quantile(popd.perm[,d],probs=0.025)

	nst1.mean.1[d] <- mean(nst1.1.perm[,d]); nst1.mean.2[d] <- mean(nst1.2.perm[,d]); nst1.mean.3[d] <- mean(nst1.3.perm[,d])
	nst1.up.1[d] <- quantile(nst1.1.perm[,d],probs=0.975); nst1.up.2[d] <- quantile(nst1.2.perm[,d],probs=0.975); nst1.up.3[d] <- quantile(nst1.3.perm[,d],probs=0.975)
	nst1.lo.1[d] <- quantile(nst1.1.perm[,d],probs=0.025); nst1.lo.2[d] <- quantile(nst1.2.perm[,d],probs=0.025); nst1.lo.3[d] <- quantile(nst1.3.perm[,d],probs=0.025)

	nst2.mean.1[d] <- mean(nst2.1.perm[,d]); nst2.mean.2[d] <- mean(nst2.2.perm[,d]); nst2.mean.3[d] <- mean(nst2.3.perm[,d])
	nst2.up.1[d] <- quantile(nst2.1.perm[,d],probs=0.975); nst2.up.2[d] <- quantile(nst2.2.perm[,d],probs=0.975); nst2.up.3[d] <- quantile(nst2.3.perm[,d],probs=0.975)
	nst2.lo.1[d] <- quantile(nst2.1.perm[,d],probs=0.025); nst2.lo.2[d] <- quantile(nst2.2.perm[,d],probs=0.025); nst2.lo.3[d] <- quantile(nst2.3.perm[,d],probs=0.025)

	nst3.mean.1[d] <- mean(nst3.1.perm[,d]); nst3.mean.2[d] <- mean(nst3.2.perm[,d]); nst3.mean.3[d] <- mean(nst3.3.perm[,d])
	nst3.up.1[d] <- quantile(nst3.1.perm[,d],probs=0.975); nst3.up.2[d] <- quantile(nst3.2.perm[,d],probs=0.975); nst3.up.3[d] <- quantile(nst3.3.perm[,d],probs=0.975)
	nst3.lo.1[d] <- quantile(nst3.1.perm[,d],probs=0.025); nst3.lo.2[d] <- quantile(nst3.2.perm[,d],probs=0.025); nst3.lo.3[d] <- quantile(nst3.3.perm[,d],probs=0.025)

	nst4.mean.1[d] <- mean(nst4.1.perm[,d]); nst4.mean.2[d] <- mean(nst4.2.perm[,d]); nst4.mean.3[d] <- mean(nst4.3.perm[,d])
	nst4.up.1[d] <- quantile(nst4.1.perm[,d],probs=0.975); nst4.up.2[d] <- quantile(nst4.2.perm[,d],probs=0.975); nst4.up.3[d] <- quantile(nst4.3.perm[,d],probs=0.975)
	nst4.lo.1[d] <- quantile(nst4.1.perm[,d],probs=0.025); nst4.lo.2[d] <- quantile(nst4.2.perm[,d],probs=0.025); nst4.lo.3[d] <- quantile(nst4.3.perm[,d],probs=0.025)

	pst1.mean.1[d] <- mean(pst1.1.perm[,d]); pst1.mean.2[d] <- mean(pst1.2.perm[,d]); pst1.mean.3[d] <- mean(pst1.3.perm[,d])
	pst1.up.1[d] <- quantile(pst1.1.perm[,d],probs=0.975); pst1.up.2[d] <- quantile(pst1.2.perm[,d],probs=0.975); pst1.up.3[d] <- quantile(pst1.3.perm[,d],probs=0.975)
	pst1.lo.1[d] <- quantile(pst1.1.perm[,d],probs=0.025); pst1.lo.2[d] <- quantile(pst1.2.perm[,d],probs=0.025); pst1.lo.3[d] <- quantile(pst1.3.perm[,d],probs=0.025)

	pst2.mean.1[d] <- mean(pst2.1.perm[,d]); pst2.mean.2[d] <- mean(pst2.2.perm[,d]); pst2.mean.3[d] <- mean(pst2.3.perm[,d])
	pst2.up.1[d] <- quantile(pst2.1.perm[,d],probs=0.975); pst2.up.2[d] <- quantile(pst2.2.perm[,d],probs=0.975); pst2.up.3[d] <- quantile(pst2.3.perm[,d],probs=0.975)
	pst2.lo.1[d] <- quantile(pst2.1.perm[,d],probs=0.025); pst2.lo.2[d] <- quantile(pst2.2.perm[,d],probs=0.025); pst2.lo.3[d] <- quantile(pst2.3.perm[,d],probs=0.025)

	pst3.mean.1[d] <- mean(pst3.1.perm[,d]); pst3.mean.2[d] <- mean(pst3.2.perm[,d]); pst3.mean.3[d] <- mean(pst3.3.perm[,d])
	pst3.up.1[d] <- quantile(pst3.1.perm[,d],probs=0.975); pst3.up.2[d] <- quantile(pst3.2.perm[,d],probs=0.975); pst3.up.3[d] <- quantile(pst3.3.perm[,d],probs=0.975)
	pst3.lo.1[d] <- quantile(pst3.1.perm[,d],probs=0.025); pst3.lo.2[d] <- quantile(pst3.2.perm[,d],probs=0.025); pst3.lo.3[d] <- quantile(pst3.3.perm[,d],probs=0.025)

	pst4.mean.1[d] <- mean(pst4.1.perm[,d]); pst4.mean.2[d] <- mean(pst4.2.perm[,d]); pst4.mean.3[d] <- mean(pst4.3.perm[,d])
	pst4.up.1[d] <- quantile(pst4.1.perm[,d],probs=0.975); pst4.up.2[d] <- quantile(pst4.2.perm[,d],probs=0.975); pst4.up.3[d] <- quantile(pst4.3.perm[,d],probs=0.975)
	pst4.lo.1[d] <- quantile(pst4.1.perm[,d],probs=0.025); pst4.lo.2[d] <- quantile(pst4.2.perm[,d],probs=0.025); pst4.lo.3[d] <- quantile(pst4.3.perm[,d],probs=0.025)

	sdend1.1[d] <- mean(sdend1.1.perm[,d]); sdend1.2[d] <- mean(sdend1.2.perm[,d]); sdend1.3[d] <- mean(sdend1.3.perm[,d])
	sdend2.1[d] <- mean(sdend2.1.perm[,d]); sdend2.2[d] <- mean(sdend2.2.perm[,d]); sdend2.3[d] <- mean(sdend2.3.perm[,d])
	sdend3.1[d] <- mean(sdend3.1.perm[,d]); sdend3.2[d] <- mean(sdend3.2.perm[,d]); sdend3.3[d] <- mean(sdend3.3.perm[,d])
	sdend4.1[d] <- mean(sdend4.1.perm[,d]); sdend4.2[d] <- mean(sdend4.2.perm[,d]); sdend4.3[d] <- mean(sdend4.3.perm[,d])
	sdend5.1[d] <- mean(sdend5.1.perm[,d]); sdend5.2[d] <- mean(sdend5.2.perm[,d]); sdend5.3[d] <- mean(sdend5.3.perm[,d])
	sdend6.1[d] <- mean(sdend6.1.perm[,d]); sdend6.2[d] <- mean(sdend6.2.perm[,d]); sdend6.3[d] <- mean(sdend6.3.perm[,d])
	sdend7.1[d] <- mean(sdend7.1.perm[,d]); sdend7.2[d] <- mean(sdend7.2.perm[,d]); sdend7.3[d] <- mean(sdend7.3.perm[,d])
	sdend8.1[d] <- mean(sdend8.1.perm[,d]); sdend8.2[d] <- mean(sdend8.2.perm[,d]); sdend8.3[d] <- mean(sdend8.3.perm[,d])
	sdend9.1[d] <- mean(sdend9.1.perm[,d]); sdend9.2[d] <- mean(sdend9.2.perm[,d]); sdend9.3[d] <- mean(sdend9.3.perm[,d])
	sdend10.1[d] <- mean(sdend10.1.perm[,d]); sdend10.2[d] <- mean(sdend10.2.perm[,d]); sdend10.3[d] <- mean(sdend10.3.perm[,d])
	sdend11.1[d] <- mean(sdend11.1.perm[,d]); sdend11.2[d] <- mean(sdend11.2.perm[,d]); sdend11.3[d] <- mean(sdend11.3.perm[,d])
	sdend12.1[d] <- mean(sdend12.1.perm[,d]); sdend12.2[d] <- mean(sdend12.2.perm[,d]); sdend12.3[d] <- mean(sdend12.3.perm[,d])
	sdend13.1[d] <- mean(sdend13.1.perm[,d]); sdend13.2[d] <- mean(sdend13.2.perm[,d]); sdend13.3[d] <- mean(sdend13.3.perm[,d])

	eggs.mean.1 <- mean(eggs.1.perm[,d]); eggs.mean.2 <- mean(eggs.2.perm[,d]); eggs.mean.3 <- mean(eggs.3.perm[,d])
	h.a1.1.fix.mean <- mean(h.a1.1.fix.perm[,d]); h.a1.2.fix.mean <- mean(h.a1.2.fix.perm[,d]); h.a1.3.fix.mean <- mean(h.a1.3.fix.perm[,d])

} ## end d loop

if (sum(popd.up.1) == 0) popd.up.1 <- popd.mean.1 else popd.up.1 <- popd.up.1
if (sum(popd.up.2) == 0) popd.up.2 <- popd.mean.2 else popd.up.2 <- popd.up.2
if (sum(popd.up.3) == 0) popd.up.3 <- popd.mean.3 else popd.up.3 <- popd.up.3
if (sum(popd.lo.1) == 0) popd.lo.1 <- popd.mean.1 else popd.lo.1 <- popd.lo.1
if (sum(popd.lo.2) == 0) popd.lo.2 <- popd.mean.2 else popd.lo.2 <- popd.lo.2
if (sum(popd.lo.3) == 0) popd.lo.3 <- popd.mean.3 else popd.lo.3 <- popd.lo.3

## Make vectors of mean end distributions
sdend.1 <- c(sdend1.1[tlimit+1],sdend2.1[tlimit+1],sdend3.1[tlimit+1],sdend4.1[tlimit+1],sdend5.1[tlimit+1],sdend6.1[tlimit+1],sdend7.1[tlimit+1],sdend8.1[tlimit+1],sdend9.1[tlimit+1],sdend10.1[tlimit+1],sdend11.1[tlimit+1],sdend12.1[tlimit+1],sdend13.1[tlimit+1])
sdend.2 <- c(sdend1.2[tlimit+1],sdend2.2[tlimit+1],sdend3.2[tlimit+1],sdend4.2[tlimit+1],sdend5.2[tlimit+1],sdend6.2[tlimit+1],sdend7.2[tlimit+1],sdend8.2[tlimit+1],sdend9.2[tlimit+1],sdend10.2[tlimit+1],sdend11.2[tlimit+1],sdend12.2[tlimit+1],sdend13.2[tlimit+1])
sdend.3 <- c(sdend1.3[tlimit+1],sdend2.3[tlimit+1],sdend3.3[tlimit+1],sdend4.3[tlimit+1],sdend5.3[tlimit+1],sdend6.3[tlimit+1],sdend7.3[tlimit+1],sdend8.3[tlimit+1],sdend9.3[tlimit+1],sdend10.3[tlimit+1],sdend11.3[tlimit+1],sdend12.3[tlimit+1],sdend13.3[tlimit+1])

## Estimate proportion of total eggs harvested (retrospective analysis)
peggs.tot.1 <- h.a1.1.fix.mean / eggs.mean.1
peggs.tot.2 <- h.a1.2.fix.mean / eggs.mean.2
peggs.tot.3 <- h.a1.3.fix.mean / eggs.mean.3

## Projection plots
row <- 2
col <- 2
par(mfrow=c(row,col))
plot(yrd.vec,popd.mean.1,xlab="year",ylab="N",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0:(max(popd.mean.1)+(0.25*(max(popd.vec.1))))),type='l')
title(main = "Region 1")
##mtext(scenario, side=3, outer=TRUE, line = -19, adj=0.5, cex = 1.2, col="red")
##mtext("% killed per year = ", side=3, outer=TRUE, line = -20, adj=0.5, cex = 0.9)
##mtext(kill.pc, side=3, outer=TRUE, line = -21, adj=0.5, cex = 0.9)
##mtext(kill.vary.scenario, side=3, outer=TRUE, line = -22, adj=0.5, cex = 0.9)
lines(yrd.vec,popd.up.1,col="red")
lines(yrd.vec,popd.lo.1,col="red")

plot(yrd.vec,popd.mean.2,xlab="year",ylab="N",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0:(max(popd.mean.2)+(0.25*(max(popd.mean.2))))),type='l')
title(main = "Region 2")
lines(yrd.vec,popd.up.2,col="red")
lines(yrd.vec,popd.lo.2,col="red")

plot(yrd.vec,popd.mean.3,xlab="year",ylab="N",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0:(max(popd.mean.3)+(0.25*(max(popd.mean.3))))),type='l')
title(main = "Region 3")
lines(yrd.vec,popd.up.3,col="red")
lines(yrd.vec,popd.lo.3,col="red")

plot(yrd.vec,popd.mean,xlab="year",ylab="N",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0:(max(popd.mean)+(0.25*(max(popd.mean))))),type='l')
title(main = "Regions Combined")
lines(yrd.vec,popd.up,col="red")
lines(yrd.vec,popd.lo,col="red")

par(mfrow=c(1,2))



####################################################################################################
## Summary
####################################################################################################

##################################
maxlambdad.1 ## Maxlambda Region 1
maxlambdad.2 ## Maxlambda Region 2
maxlambdad.3 ## Maxlambda Region 3
##################################

###########################################
n.final.1 ## Final population size Region 1
n.final.2 ## Final population size Region 2
n.final.3 ## Final population size Region 3
###########################################

#################################################
n.final ## Total population size over all regions
#################################################

############################################################
popd.mean.1[tlimit+1] ## Final mean population size Region 1
popd.mean.2[tlimit+1] ## Final mean population size Region 2
popd.mean.3[tlimit+1] ## Final mean population size Region 3
############################################################

###################################################
popd.mean[tlimit+1] ## Mean total
popd.lo[tlimit+1] ## Lower 95 % confidence interval
popd.up[tlimit+1] ## Upper 95 % confidence interval
###################################################

##################################
kill.pc ## % killed per year
kill.vary.scenario
##################################

########################################
## Mean egg harvest proportions
## Only valid for retrospective analysis
peggs.tot.1
peggs.tot.2
peggs.tot.3
########################################

####################################################################################################
####################################################################################################

## Stages (by proportion)
## Region 1
row <- 2
col <- 2
par(mfrow=c(row,col))
plot(yrd.vec,pst1.mean.1,xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,0.55),type='l')
lines(yrd.vec,pst1.up.1,col="red")
lines(yrd.vec,pst1.lo.1,col="red")
title(main = "Stage 1")
##mtext("Region 1", side=3, outer=TRUE, line = -3, adj=0.5, cex = 1.5, col="red")
plot(yrd.vec,pst2.mean.1,xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,0.55),type='l')
lines(yrd.vec,pst2.up.1,col="red")
lines(yrd.vec,pst2.lo.1,col="red")
title(main = "Stage 2")
plot(yrd.vec,pst3.mean.1,xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,0.55),type='l')
lines(yrd.vec,pst3.up.1,col="red")
lines(yrd.vec,pst3.lo.1,col="red")
title(main = "Stage 3")
plot(yrd.vec,pst4.mean.1,xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,0.55),type='l')
lines(yrd.vec,pst4.up.1,col="red")
lines(yrd.vec,pst4.lo.1,col="red")
title(main = "Stage 4")
par(mfrow=c(1,2))

## Region 2
row <- 2
col <- 2
par(mfrow=c(row,col))
plot(yrd.vec,pst1.mean.2,xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,0.55),type='l')
lines(yrd.vec,pst1.up.2,col="red")
lines(yrd.vec,pst1.lo.2,col="red")
title(main = "Stage 1")
##mtext("Region 1", side=3, outer=TRUE, line = -3, adj=0.5, cex = 1.5, col="red")
plot(yrd.vec,pst2.mean.2,xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,0.55),type='l')
lines(yrd.vec,pst2.up.2,col="red")
lines(yrd.vec,pst2.lo.2,col="red")
title(main = "Stage 2")
plot(yrd.vec,pst3.mean.2,xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,0.55),type='l')
lines(yrd.vec,pst3.up.2,col="red")
lines(yrd.vec,pst3.lo.2,col="red")
title(main = "Stage 3")
plot(yrd.vec,pst4.mean.2,xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,0.55),type='l')
lines(yrd.vec,pst4.up.2,col="red")
lines(yrd.vec,pst4.lo.2,col="red")
title(main = "Stage 4")
par(mfrow=c(1,2))

## Region 3
row <- 2
col <- 2
par(mfrow=c(row,col))
plot(yrd.vec,pst1.mean.3,xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,0.55),type='l')
lines(yrd.vec,pst1.up.3,col="red")
lines(yrd.vec,pst1.lo.3,col="red")
title(main = "Stage 1")
##mtext("Region 1", side=3, outer=TRUE, line = -3, adj=0.5, cex = 1.5, col="red")
plot(yrd.vec,pst2.mean.3,xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,0.55),type='l')
lines(yrd.vec,pst2.up.3,col="red")
lines(yrd.vec,pst2.lo.3,col="red")
title(main = "Stage 2")
plot(yrd.vec,pst3.mean.3,xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,0.55),type='l')
lines(yrd.vec,pst3.up.3,col="red")
lines(yrd.vec,pst3.lo.3,col="red")
title(main = "Stage 3")
plot(yrd.vec,pst4.mean.3,xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,0.55),type='l')
lines(yrd.vec,pst4.up.3,col="red")
lines(yrd.vec,pst4.lo.3,col="red")
title(main = "Stage 4")
par(mfrow=c(1,2))

## Temporal change in age composition according to Parks data

obs.pst1.1 <- c(0.2094,0.2801,0.2187,0.2359,0.1427,0.4984,0.3828,0.2400,0.2766,0.3261,0.3502,0.2238,0.2569,0.4830,0.3686,0.3369,0.1875,0.1480,0.1450,0.3153,0.1907,0.2074)
obs.pst2.1 <- c(0.1718,0.1782,0.1162,0.1773,0.1365,0.0685,0.0781,0.1449,0.0626,0.0961,0.0765,0.1303,0.1501,0.0751,0.1337,0.0636,0.0643,0.0873,0.0863,0.0739,0.0532,0.0263)
obs.pst3.1 <- c(0.3595,0.2673,0.2631,0.2204,0.3164,0.1121,0.0742,0.1849,0.1715,0.1482,0.0989,0.1436,0.1486,0.0794,0.1314,0.0760,0.0913,0.1610,0.1133,0.1157,0.1306,0.1132)
obs.pst4.1 <- c(0.2593,0.2743,0.4020,0.3664,0.4045,0.3209,0.4648,0.4302,0.4892,0.4296,0.4744,0.5022,0.4444,0.3625,0.3662,0.5235,0.6568,0.6037,0.6554,0.4951,0.6255,0.6530)

obs.pst1.2 <- c(0.1010,0.0000,0.0000,0.0000,0.0000,0.0000,0.0558,0.0787,0.0943,0.0690,0.0637,0.0692,0.1044,0.0477,0.0835,0.0344,0.0665,0.0187,0.0524,0.0656,0.0239)
obs.pst2.2 <- c(0.1263,0.0000,0.0000,0.0000,0.0000,0.0000,0.0413,0.0208,0.0450,0.0548,0.0483,0.0576,0.0414,0.0197,0.0355,0.0362,0.0490,0.0272,0.0349,0.0273,0.0417)
obs.pst3.2 <- c(0.3704,0.0000,0.0000,0.0000,0.0000,0.0000,0.1286,0.0973,0.1072,0.1357,0.0989,0.0807,0.1162,0.0592,0.0924,0.0869,0.1261,0.0732,0.1309,0.1276,0.0417)
obs.pst4.2 <- c(0.4024,0.0000,0.0000,0.0000,0.0000,0.0000,0.7743,0.8031,0.7535,0.7405,0.7890,0.7925,0.7380,0.8733,0.7886,0.8425,0.7584,0.8808,0.7818,0.7795,0.8926)

obs.pst1.3 <- c(0.5581,0.4607,0.4136,0.3998,0.4812,0.5511,0.4815,0.4025,0.4345,0.4880,0.3574,0.3688,0.5378,0.2997,0.5695,0.5570,0.4984,0.3392,0.4988,0.3021,0.3290,0.2411,0.3947,0.3947,0.2107)
obs.pst2.3 <- c(0.3876,0.2216,0.2738,0.2802,0.1857,0.1187,0.1576,0.2017,0.1544,0.1489,0.1071,0.0865,0.0554,0.1020,0.1002,0.0565,0.0266,0.1079,0.0970,0.0934,0.0899,0.0985,0.0726,0.0716,0.0192)
obs.pst3.3 <- c(0.0543,0.2618,0.1655,0.2201,0.2273,0.1737,0.1911,0.1934,0.1980,0.1948,0.1520,0.1685,0.1203,0.1313,0.0863,0.0523,0.0747,0.1245,0.0938,0.1197,0.1373,0.1514,0.1014,0.0993,0.0383)
obs.pst4.3 <- c(0.0000,0.0559,0.1471,0.0999,0.1058,0.1565,0.1698,0.2024,0.2132,0.1683,0.3835,0.3761,0.2865,0.4669,0.2440,0.3341,0.4003,0.4284,0.3105,0.4848,0.4438,0.5091,0.4312,0.4345,0.7318)

## Stages (by proportion, with observed values)
## Region 1
row <- 2
col <- 2
par(mfrow=c(row,col))
plot(yrd.vec[2:tlimit+1],pst1.mean.1[2:tlimit],xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,(max(c(obs.pst1.1,obs.pst2.1,obs.pst3.1,obs.pst4.1)))),type='l')
lines(yrd.vec[2:tlimit+1],pst1.up.1[2:tlimit],col="red")
lines(yrd.vec[2:tlimit+1],pst1.lo.1[2:tlimit],col="red")
#points(yrd.vec[2:tlimit+1],obs.pst1.1[4:22],pch=19)
title(main = "Stage 1")
##mtext("Region 1", side=3, outer=TRUE, line = -3, adj=0.5, cex = 1.5, col="red")
plot(yrd.vec[2:tlimit+1],pst2.mean.1[2:tlimit],xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,(max(c(obs.pst1.1,obs.pst2.1,obs.pst3.1,obs.pst4.1)))),type='l')
lines(yrd.vec[2:tlimit+1],pst2.up.1[2:tlimit],col="red")
lines(yrd.vec[2:tlimit+1],pst2.lo.1[2:tlimit],col="red")
#points(yrd.vec[2:tlimit+1],obs.pst2.1[4:22],pch=19)
title(main = "Stage 2")
plot(yrd.vec[2:tlimit+1],pst3.mean.1[2:tlimit],xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,(max(c(obs.pst1.1,obs.pst2.1,obs.pst3.1,obs.pst4.1)))),type='l')
lines(yrd.vec[2:tlimit+1],pst3.up.1[2:tlimit],col="red")
lines(yrd.vec[2:tlimit+1],pst3.lo.1[2:tlimit],col="red")
#points(yrd.vec[2:tlimit+1],obs.pst3.1[4:22],pch=19)
title(main = "Stage 3")
plot(yrd.vec[2:tlimit+1],pst4.mean.1[2:tlimit],xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,(max(c(obs.pst1.1,obs.pst2.1,obs.pst3.1,obs.pst4.1)))),type='l')
lines(yrd.vec[2:tlimit+1],pst4.up.1[2:tlimit],col="red")
lines(yrd.vec[2:tlimit+1],pst4.lo.1[2:tlimit],col="red")
points(yrd.vec[2:tlimit+1],obs.pst4.1[4:22],pch=19)
#title(main = "Stage 4")
par(mfrow=c(1,2))

## Region 2
row <- 2
col <- 2
par(mfrow=c(row,col))
plot(yrd.vec[2:tlimit+1],pst1.mean.2[2:tlimit],xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,(max(c(obs.pst1.2,obs.pst2.2,obs.pst3.2,obs.pst4.2)))),type='l')
lines(yrd.vec[2:tlimit+1],pst1.up.2[2:tlimit],col="red")
lines(yrd.vec[2:tlimit+1],pst1.lo.2[2:tlimit],col="red")
points(yrd.vec[7:21],obs.pst1.2[7:21],pch=19)
title(main = "Stage 1")
##mtext("Region 1", side=3, outer=TRUE, line = -3, adj=0.5, cex = 1.5, col="red")
plot(yrd.vec[2:tlimit+1],pst2.mean.2[2:tlimit],xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,(max(c(obs.pst1.2,obs.pst2.2,obs.pst3.2,obs.pst4.2)))),type='l')
lines(yrd.vec[2:tlimit+1],pst2.up.2[2:tlimit],col="red")
lines(yrd.vec[2:tlimit+1],pst2.lo.2[2:tlimit],col="red")
points(yrd.vec[7:21],obs.pst2.2[7:21],pch=19)
title(main = "Stage 2")
plot(yrd.vec[2:tlimit+1],pst3.mean.2[2:tlimit],xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,(max(c(obs.pst1.2,obs.pst2.2,obs.pst3.2,obs.pst4.2)))),type='l')
lines(yrd.vec[2:tlimit+1],pst3.up.2[2:tlimit],col="red")
lines(yrd.vec[2:tlimit+1],pst3.lo.2[2:tlimit],col="red")
points(yrd.vec[7:21],obs.pst3.2[7:21],pch=19)
title(main = "Stage 3")
plot(yrd.vec[2:tlimit+1],pst4.mean.2[2:tlimit],xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,(max(c(obs.pst1.2,obs.pst2.2,obs.pst3.2,obs.pst4.2)))),type='l')
lines(yrd.vec[2:tlimit+1],pst4.up.2[2:tlimit],col="red")
lines(yrd.vec[2:tlimit+1],pst4.lo.2[2:tlimit],col="red")
points(yrd.vec[7:21],obs.pst4.2[7:21],pch=19)
title(main = "Stage 4")
par(mfrow=c(1,2))

## Region 3
row <- 2
col <- 2
par(mfrow=c(row,col))
plot(yrd.vec[2:tlimit+1],pst1.mean.3[2:tlimit],xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,(max(c(obs.pst1.3,obs.pst2.3,obs.pst3.3,obs.pst4.3)))),type='l')
lines(yrd.vec[2:tlimit+1],pst1.up.3[2:tlimit],col="red")
lines(yrd.vec[2:tlimit+1],pst1.lo.3[2:tlimit],col="red")
#points(yrd.vec[2:tlimit+1],obs.pst1.3[7:25],pch=19)
title(main = "Stage 1")
##mtext("Region 1", side=3, outer=TRUE, line = -3, adj=0.5, cex = 1.5, col="red")
plot(yrd.vec[2:tlimit+1],pst2.mean.3[2:tlimit],xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,(max(c(obs.pst1.3,obs.pst2.3,obs.pst3.3,obs.pst4.3)))),type='l')
lines(yrd.vec[2:tlimit+1],pst2.up.3[2:tlimit],col="red")
lines(yrd.vec[2:tlimit+1],pst2.lo.3[2:tlimit],col="red")
#points(yrd.vec[2:tlimit+1],obs.pst2.3[7:25],pch=19)
title(main = "Stage 2")
plot(yrd.vec[2:tlimit+1],pst3.mean.3[2:tlimit],xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,(max(c(obs.pst1.3,obs.pst2.3,obs.pst3.3,obs.pst4.3)))),type='l')
lines(yrd.vec[2:tlimit+1],pst3.up.3[2:tlimit],col="red")
lines(yrd.vec[2:tlimit+1],pst3.lo.3[2:tlimit],col="red")
#points(yrd.vec[2:tlimit+1],obs.pst3.3[7:25],pch=19)
title(main = "Stage 3")
plot(yrd.vec[2:tlimit+1],pst4.mean.3[2:tlimit],xlab="year",ylab="Proportion",xlim=range((st.yr-0.5):(st.yr+tlimit+1)),ylim=range(0,(max(c(obs.pst1.3,obs.pst2.3,obs.pst3.3,obs.pst4.3)))),type='l')
lines(yrd.vec[2:tlimit+1],pst4.up.3[2:tlimit],col="red")
lines(yrd.vec[2:tlimit+1],pst4.lo.3[2:tlimit],col="red")
#points(yrd.vec[2:tlimit+1],obs.pst4.3[7:25],pch=19)
title(main = "Stage 4")
par(mfrow=c(1,2))


## Check proportions per age class (mean 1996-1998)
################################################################################
## Region 1
## Adelaide/Mary/Wildman/W Alligator/S Alligator/E Alligator/Coopers/Murganella
################################################################################
rowa <- 3
cola <- 3
par(mfrow=c(rowa,cola))
plot(end.sd1[1:age.max],type="h",xlab="age class",ylab="observed proportion",ylim=range(c(0,max(c(end.sd1,sdend.1,ssdd.1)))))
plot(sdend.1[1:age.max],type="h",xlab="age class",ylab="predicted end proportion",ylim=range(c(0,max(c(end.sd1,sdend.1,ssdd.1)))))
title(main="Region 1")
plot(ssdd.1[1:age.max],type="h",xlab="age class",ylab="stable proportion",ylim=range(c(0,max(c(end.sd1,sdend.1,ssdd.1)))))

##########################################################################
## Region 2
## Victoria/Fitzmaurice/Moyle/Daly/Reynolds/Finniss/
## Hart/Rose/Roper/Towns/Nathan/McArthur/Johnson/Wearyan/Robinson/Calvert
##########################################################################
plot(end.sd2[1:age.max],type="h",xlab="age class",ylab="observed prop",ylim=range(c(0,max(c(end.sd2,sdend.2,ssdd.2)))))
plot(sdend.2[1:age.max],type="h",xlab="age class",ylab="predicted end prop",ylim=range(c(0,max(c(end.sd2,sdend.2,ssdd.2)))))
title(main="Region 2")
plot(ssdd.2[1:age.max],type="h",xlab="age class",ylab="stable prop",ylim=range(c(0,max(c(end.sd2,sdend.2,ssdd.2)))))

###############################################################################################################
## Region 3
## King/Goomadeer/Liverpool/Tomkinson/Cadell/Blyth/Glyde/Goyder/Woolen/Kalarwoi/
## Buckingham/Warawurowoi/Kurala/Habgood/Darwarunga/Baralminar/Gobalpa/Goromuru/Cato/Peter John/Burungbirinung
###############################################################################################################
plot(end.sd3[1:age.max],type="h",xlab="age class",ylab="observed prop",ylim=range(c(0,max(c(end.sd3,sdend.3,ssdd.3)))))
plot(sdend.3[1:age.max],type="h",xlab="age class",ylab="predicted end prop",ylim=range(c(0,max(c(end.sd3,sdend.3,ssdd.3)))))
title(main="Region 3")
plot(ssdd.3[1:age.max],type="h",xlab="age class",ylab="stable prop",ylim=range(c(0,max(c(end.sd3,sdend.3,ssdd.3)))))
par(mfrow=c(1,2))


## Replot proportions by stage
## Observed end proportions
ost1end.1 <- sum(end.sd1[1]); ost2end.1 <- sum(end.sd1[2:3]); ost3end.1 <- sum(end.sd1[4:12]); ost4end.1 <- sum(end.sd1[age.max])
ost1end.2 <- sum(end.sd2[1]); ost2end.2 <- sum(end.sd2[2:3]); ost3end.2 <- sum(end.sd2[4:12]); ost4end.2 <- sum(end.sd2[age.max])
ost1end.3 <- sum(end.sd3[1]); ost2end.3 <- sum(end.sd3[2:3]); ost3end.3 <- sum(end.sd3[4:12]); ost4end.3 <- sum(end.sd3[age.max])
ostend.1 <- c(ost1end.1,ost2end.1,ost3end.1,ost4end.1); ostend.2 <- c(ost1end.2,ost2end.2,ost3end.2,ost4end.2); ostend.3 <- c(ost1end.3,ost2end.3,ost3end.3,ost4end.3)

## Predicted end proportions
pst1end.1 <- pst1.1.perm[tlimit+1]; pst2end.1 <- pst2.1.perm[tlimit+1]; pst3end.1 <- pst3.1.perm[tlimit+1]; pst4end.1 <- pst4.1.perm[tlimit+1]
pst1end.2 <- pst1.2.perm[tlimit+1]; pst2end.2 <- pst2.2.perm[tlimit+1]; pst3end.2 <- pst3.2.perm[tlimit+1]; pst4end.2 <- pst4.2.perm[tlimit+1]
pst1end.3 <- pst1.3.perm[tlimit+1]; pst2end.3 <- pst2.3.perm[tlimit+1]; pst3end.3 <- pst3.3.perm[tlimit+1]; pst4end.3 <- pst4.3.perm[tlimit+1]
pstend.1 <- c(pst1end.1,pst2end.1,pst3end.1,pst4end.1); pstend.2 <- c(pst1end.2,pst2end.2,pst3end.2,pst4end.2); pstend.3 <- c(pst1end.3,pst2end.3,pst3end.3,pst4end.3)

## Stable end proportions
sst1end.1 <- sum(ssdd.1[1]); sst2end.1 <- sum(ssdd.1[2:3]); sst3end.1 <- sum(ssdd.1[4:12]); sst4end.1 <- sum(ssdd.1[age.max])
sst1end.2 <- sum(ssdd.2[1]); sst2end.2 <- sum(ssdd.2[2:3]); sst3end.2 <- sum(ssdd.2[4:12]); sst4end.2 <- sum(ssdd.2[age.max])
sst1end.3 <- sum(ssdd.3[1]); sst2end.3 <- sum(ssdd.3[2:3]); sst3end.3 <- sum(ssdd.3[4:12]); sst4end.3 <- sum(ssdd.3[age.max])
sstend.1 <- c(sst1end.1,sst2end.1,sst3end.1,sst4end.1); sstend.2 <- c(sst1end.2,sst2end.2,sst3end.2,sst4end.2); sstend.3 <- c(sst1end.3,sst2end.3,sst3end.3,sst4end.3)

################################################################################
## Region 1
## Adelaide/Mary/Wildman/W Alligator/S Alligator/E Alligator/Coopers/Murganella
################################################################################
rowa <- 3
cola <- 3
par(mfrow=c(rowa,cola))
plot(ostend.1,type="h",xlab="age class",ylab="observed proportion",ylim=range(c(0,max(c(ostend.1,pstend.1,sstend.1)))))
plot(pstend.1,type="h",xlab="age class",ylab="predicted end proportion",ylim=range(c(0,max(c(ostend.1,pstend.1,sstend.1)))))
title(main="Region 1")
plot(sstend.1,type="h",xlab="age class",ylab="stable proportion",ylim=range(c(0,max(c(ostend.1,pstend.1,sstend.1)))))

##########################################################################
## Region 2
## Victoria/Fitzmaurice/Moyle/Daly/Reynolds/Finniss/
## Hart/Rose/Roper/Towns/Nathan/McArthur/Johnson/Wearyan/Robinson/Calvert
##########################################################################
plot(ostend.2,type="h",xlab="age class",ylab="observed proportion",ylim=range(c(0,max(c(ostend.2,pstend.2,sstend.2)))))
plot(pstend.2,type="h",xlab="age class",ylab="predicted end proportion",ylim=range(c(0,max(c(ostend.2,pstend.2,sstend.2)))))
title(main="Region 2")
plot(sstend.2,type="h",xlab="age class",ylab="stable proportion",ylim=range(c(0,max(c(ostend.2,pstend.2,sstend.2)))))

###############################################################################################################
## Region 3
## King/Goomadeer/Liverpool/Tomkinson/Cadell/Blyth/Glyde/Goyder/Woolen/Kalarwoi/
## Buckingham/Warawurowoi/Kurala/Habgood/Darwarunga/Baralminar/Gobalpa/Goromuru/Cato/Peter John/Burungbirinung
###############################################################################################################
plot(ostend.3,type="h",xlab="age class",ylab="observed proportion",ylim=range(c(0,max(c(ostend.3,pstend.3,sstend.3)))))
plot(pstend.3,type="h",xlab="age class",ylab="predicted end proportion",ylim=range(c(0,max(c(ostend.3,pstend.3,sstend.3)))))
title(main="Region 3")
plot(sstend.3,type="h",xlab="age class",ylab="stable proportion",ylim=range(c(0,max(c(ostend.3,pstend.3,sstend.3)))))
par(mfrow=c(1,2))


## Replot proportions by age (1-6 years)
## Observed end proportions
oa6end.1 <- end.sd1.tmp1/2; oa6end.2 <- end.sd2.tmp1/2; oa6end.3 <- end.sd3.tmp1/2

## Predicted end proportions
pa6end.1 <- sdend.1[1:6]; pa6end.1[6] <- sum(sdend.1[6:age.max]); pa6end.2 <- sdend.2[1:6]; pa6end.2[6] <- sum(sdend.2[6:age.max]); pa6end.3 <- sdend.3[1:6]; pa6end.3[6] <- sum(sdend.3[6:age.max])

## Stable end proportions
sa6end.1 <- ssdd.1[1:6]; sa6end.1[6] <- sum(ssdd.1[6:age.max]); sa6end.2 <- ssdd.2[1:6]; sa6end.2[6] <- sum(ssdd.2[6:age.max]); sa6end.3 <- ssdd.3[1:6]; sa6end.3[6] <- sum(ssdd.3[6:age.max])

################################################################################
## Region 1
## Adelaide/Mary/Wildman/W Alligator/S Alligator/E Alligator/Coopers/Murganella
################################################################################
rowa <- 3
cola <- 3
par(mfrow=c(rowa,cola))
plot(oa6end.1,type="h",xlab="age class",ylab="observed proportion",ylim=range(c(0,max(c(oa6end.1,pa6end.1,sa6end.1)))))
plot(pa6end.1,type="h",xlab="age class",ylab="predicted end proportion",ylim=range(c(0,max(c(oa6end.1,pa6end.1,sa6end.1)))))
title(main="Region 1")
plot(sa6end.1,type="h",xlab="age class",ylab="stable proportion",ylim=range(c(0,max(c(oa6end.1,pa6end.1,sa6end.1)))))

##########################################################################
## Region 2
## Victoria/Fitzmaurice/Moyle/Daly/Reynolds/Finniss/
## Hart/Rose/Roper/Towns/Nathan/McArthur/Johnson/Wearyan/Robinson/Calvert
##########################################################################
plot(oa6end.2,type="h",xlab="age class",ylab="observed proportion",ylim=range(c(0,max(c(oa6end.2,pa6end.2,sa6end.2)))))
plot(pa6end.2,type="h",xlab="age class",ylab="predicted end proportion",ylim=range(c(0,max(c(oa6end.2,pa6end.2,sa6end.2)))))
title(main="Region 2")
plot(sa6end.2,type="h",xlab="age class",ylab="stable proportion",ylim=range(c(0,max(c(oa6end.2,pa6end.2,sa6end.2)))))

###############################################################################################################
## Region 3
## King/Goomadeer/Liverpool/Tomkinson/Cadell/Blyth/Glyde/Goyder/Woolen/Kalarwoi/
## Buckingham/Warawurowoi/Kurala/Habgood/Darwarunga/Baralminar/Gobalpa/Goromuru/Cato/Peter John/Burungbirinung
###############################################################################################################
plot(oa6end.3,type="h",xlab="age class",ylab="observed proportion",ylim=range(c(0,max(c(oa6end.3,pa6end.3,sa6end.3)))))
plot(pa6end.3,type="h",xlab="age class",ylab="predicted end proportion",ylim=range(c(0,max(c(oa6end.3,pa6end.3,sa6end.3)))))
title(main="Region 3")
plot(sa6end.3,type="h",xlab="age class",ylab="stable proportion",ylim=range(c(0,max(c(oa6end.3,pa6end.3,sa6end.3)))))
par(mfrow=c(1,2))



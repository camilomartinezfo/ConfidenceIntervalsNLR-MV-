library(optimr)

data <- read.csv2("abarcodata.csv")
data$dD <- as.numeric(data$dD)
data$D <- as.numeric(data$D)
data$Drelativa <- as.numeric(data$Drelativa)
data$D <- data$D/10
data$dD <- data$dD/10
data$Drelativa <- (data$dD)/data$D
pos <- which(data$D>150)
data <- data[-c(pos),]
pos2 <- which(data$Drelativa==1)
data <- data[-c(pos2),]


## Maximum Likelihood Estimators
loglik.bertalanffy = function(parameters, year, D)
{
      logL =  tryCatch(
            {
                  A = parameters["A"]
                  k = parameters["k"]
                  c = parameters["c"]
                  sigma =  parameters["sigma"]
                  model = A*((1-exp(-k*year))^(c))
                  logL <- sum(dnorm(D, model, sigma, log = TRUE)) 
            },
            error = function(cond)
            {
                  logL = -999999
            })
      print(logL)
      return(-logL)
}

parameters = c(A =  140, k = 0.007, c = 1.50, sigma = 1)

modelo1.ML = with(subset(data,!is.na(D) & !is.na(year)),
                  optimr(par = parameters,
                         fn  = loglik.bertalanffy,
                         year   = year,
                         D   = D,
                         control = list(maxiter=1000),
                         method = "Nelder-Mead",
                         hessian = TRUE))

FIM <- solve(modelo1.ML$hessian)

se <- sqrt(diag(FIM))

t <- modelo1.ML$par/se

pval <- 2*(1-pt(abs(t),length(data$year)-4))

# Displaying results
results <- as.matrix(cbind(modelo1.ML$par,se,t,pval))
colnames (results) <- c("parameter","se","t","p")
rownames (results) <- c("A","k","c","sigma")
print(results,digits=5)



## IC: delta method 
fgh2 <- deriv(D ~ A*((1-exp(-k*year))^c), c("A", "k", "c"), 
              function(A, k, c, year){} ) 
range(data$year)

x.new <- seq(1, 200, by=0.1)
beta2.est <- c(modelo1.ML$par)[-4]
f.new <- fgh2(beta2.est[1],beta2.est[2], beta2.est[3], x.new)

g.new <- attr(f.new,"gradient")
VCOV <- FIM[-4,-4]
rownames(VCOV) <- c("A","k","c")
colnames(VCOV) <- c("A","k","c")
V.beta2 <- VCOV
GS=rowSums((g.new%*%V.beta2)*g.new)

s.t <- qt(0.975,44)
deltaf <- sqrt(GS)*s.t
df.delta <- data.frame(x=x.new, f=f.new, lwr.conf=f.new-deltaf, upr.conf=f.new+deltaf)
head(df.delta)
tail(df.delta)

sigma2.est <- modelo1.ML$par[4]
deltay <- sqrt(GS + sigma2.est^2)*s.t


ci.lines<-function(){
      
      yv <- f.new
      ci<-deltaf
      uyv<-yv+ci
      lyv<-yv-ci
      lines(x.new,uyv,lty=3,lwd=2,col="black")
      lines(x.new,lyv,lty=3,lwd=2,col="black")
}

par(mfrow=c(1,1),mai=c(0.9,0.9,0.5,0.5),mgp=c(2.0,0.6,0), cex=1.2)
plot(data$year,data$D,xlim=c(0,200),ylim=c(0,150),pch=16,bg="white",
     col="darkgray",cex=1,lwd=1.9,xlab="Age (year)",ylab="Diameter (cm)",cex.lab=1.2)
ci.lines()
curve(beta2.est[1]*((1-exp(-beta2.est[2]*x))^(beta2.est[3])), add = TRUE, col = "black", lwd =2)


# Approach Bates & Watts (1988, p. 59)
#-------------------------------------------------------------------------------

# Matrix of regressors
X = matrix(cbind(rep(1,48),data$year),48,2)
X

# Vector of observed response values
y <- matrix(data$D)
y

# Gauss-Newton Method

abarco = function(A, k, c, year){
   D = A*((1-exp(-k*year))^c)
   return(as.matrix(D))}

#  Optimal parametric estimators
theta = matrix(c(modelo1.ML$par[1],modelo1.ML$par[2],modelo1.ML$par[3]),3,1)
pred = round(abarco(data$year,theta[1],theta[2],theta[3]),3)
pred

# Table 2.1 from Bates & Watts (1988), p. 41
table = as.data.frame(cbind(rep(1:length(data$year)), X[,2], y, pred))
names(table) = c("Obs", "x", "y", "predicted_0")

# Calculate the residuals and add them to the table
table$residual_0 = table$y - table$predicted

# Derived matrix
derivative = function(year,A,k,c){
   d1 = round(A*((1-exp(-k*year))^c),4)
   d2 = round(A*((1-exp(-k*year))^c),2)
   return(as.matrix(cbind(d1,d2)))
}

# The derived matrix is calculated for all observations
# and the initial estimators.

V_0 = derivative(data$year,theta[1],theta[2],theta[3])
V_0

table$der_theta1_0 = V_0[,1]
table$der_theta2_0 = V_0[,2]

# Calculation of the residual sum of squares for the initial iteration.
rss0 = sum(table$residual_0^2)

# QR decomposition of the derived matrix V^0
QR =qr(V_0)
QR$qr

# Recovering of matrices Q and R from a QR object.
# In this case we obtain Q1 and R1
Q1 = qr.Q(QR)
Q1
R1 = qr.R(QR)
R1

# Bands for expected responses at year = 
y0.100 = round(abarco(theta[1],theta[2],theta[3],50),5)
y0.100

# Derivate vector
v = derivative(50,theta[1],theta[2],theta[3])
v

# v^T*(R1)^-1
vTinvR = v%*%solve(R1)
vTinvR

# Norm of the vector vTinvR
norm1 = sum(vTinvR*vTinvR)^0.5
norm1

# Trazar bandas para la respuesta esperada
s = sqrt(rss0/44)
EFE = qf(0.95,4,44)
P = 4

xnew <- seq(1,200,0.1) #range
ynew = theta[1]*((1-exp(-theta[2]*xnew))^theta[3])
par(mfrow=c(1,1),mai=c(0.9,0.9,0.5,0.5),mgp=c(2.0,0.6,0), cex=1.2)
plot(data$year, data$D, pch=16, col = "darkgray",
     xlab = "Age (year)",
     ylab = "Diameter (cm)",
     xlim = c(1,210),
     ylim = c(0,150), cex.lab =1.2)
lines(xnew,ynew,lwd=2)
lines(xnew,ynew+s*norm1*sqrt(P*EFE),lwd=2,lty=3,col="black")
lines(xnew,ynew-s*norm1*sqrt(P*EFE),lwd=2,lty=3,col="black")


## IC: Monte Carlo method

FittedML <- function(x){
      diametro = beta2.est[1]*((1-exp(-beta2.est[2]*x))^(beta2.est[3]))
      return(diametro)
}

predictNLS_MV <- function(EXPR, modelo, var.pred, newdata, level = 0.95, 
                          nsim = 10000, ...)
{
      require(MASS, quietly = TRUE)
      
      ## coefficients
      COEF <- modelo1.ML$par[-4]
      
      ## get variance-covariance matrix
      FIM <- solve(modelo$hessian)
      VCOV <- FIM[-4,-4]
      rownames(VCOV) <- c("A","k","c")
      colnames(VCOV) <- c("A","k","c")
      
      ## augment variance-covariance matrix for 'mvrnorm' 
      ## by adding a column/row for 'error in x'
      NCOL <- ncol(VCOV)
      ADD1 <- c(rep(0, NCOL))
      ADD1 <- matrix(ADD1, ncol = 1)
      colnames(ADD1) <- var.pred
      VCOV <- cbind(VCOV, ADD1)
      ADD2 <- c(rep(0, NCOL + 1))
      ADD2 <- matrix(ADD2, nrow = 1)
      rownames(ADD2) <- var.pred
      VCOV <- rbind(VCOV, ADD2)
      
      NR <- nrow(newdata)
      varPLACE <- ncol(VCOV)
      outMAT <- NULL  
      
      for (i in 1:NR) {
            
            ## get predictor values and optional errors
            predVAL <- newdata[i, 1]
            predERROR <- 0
            names(predVAL) <- var.pred
            names(predERROR) <- var.pred
            
            ## create mean vector for 'mvrnorm'
            MU <- c(COEF, predVAL)
            
            ## create variance-covariance matrix for 'mvrnorm'
            ## by putting error^2 in lower-right position of VCOV
            newVCOV <- VCOV
            newVCOV[varPLACE, varPLACE] <- predERROR^2
            
            ## create MC simulation matrix
            simMAT <- mvrnorm(n = nsim, mu = MU, Sigma = newVCOV, 
                              empirical = TRUE)
            
            ## evaluate expression on rows of simMAT
            EVAL <- eval(EXPR, envir = as.data.frame(simMAT))
            
            ## collect statistics
            PRED <- data.frame(predVAL)
            colnames(PRED) <- var.pred
            FITTED <- FittedML(data.frame(PRED))
            MEAN.sim <- mean(EVAL, na.rm = TRUE)
            SD.sim <- sd(EVAL, na.rm = TRUE)
            QUANT <- quantile(EVAL, c((1 - level)/2, level + (1 - level)/2),
                              na.rm=TRUE)
            outMAT <- rbind(outMAT, c(FITTED, MEAN.sim, SD.sim))
      }
      
      colnames(outMAT) <- c("fit", "mean", "sd")
      rownames(outMAT) <- NULL
      
      return(outMAT)  
}


fitted <- predictNLS_MV(EXPR=expression(A*((1-exp(-k*year))^c)), modelo=modelo1.ML, 
                           var.pred = "year", newdata = data.frame(x.new))
fitted <- as.data.frame(fitted)

Predicted <- c()
for(i in 1:1991){
      Predicted[i] <- fitted$fit[[i]]
}

Desv <- c()
for(i in 1:1991){
      Desv[i] <- fitted$sd[[i]]
}

plot(data$year,data$D,xlim=c(0,200),ylim=c(0,150),pch=16,bg="white",
     col="darkgray",cex=1,lwd=1.9,xlab="Age (year)",ylab="Diameter (cm)",cex.lab=1.2)
curve(beta2.est[1]*((1-exp(-beta2.est[2]*x))^(beta2.est[3])), add = TRUE, col = "black", lwd =2)
valort <- qt(0.975,45)
lines(x.new,Predicted+Desv*valort,lwd=2,lty=3)
lines(x.new,Predicted-Desv*valort,lwd=2,lty=3)


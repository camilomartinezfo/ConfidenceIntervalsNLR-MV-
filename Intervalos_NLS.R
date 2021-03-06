
# reading data.
data = read.csv("Puromycin.csv")


################################################################################
################ Least Squares Estimators ##############################
################################################################################

nonlinearmod = nls(Velocity ~ beta1*Conc/(beta2 + Conc), data = data, 
                   start = list(beta1=205, beta2=0.08),
                   trace = TRUE)
summary(nonlinearmod)


# Confidence intervals for the mean response.

# Delta method

#-------------------------------------------------------------------------------

fgh2 <- deriv(Velocity ~ beta1*Conc/(beta2 + Conc), c("beta1", "beta2"), 
              function(beta1, beta2, Conc){} ) 

x.new <- seq(0.01, 1.2, by=0.01)
beta2.est <- coef(nonlinearmod)
f.new <- fgh2(beta2.est[1],beta2.est[2], x.new)

g.new <- attr(f.new,"gradient")
V.beta2 <- vcov(nonlinearmod)
GS=rowSums((g.new%*%V.beta2)*g.new)

t <- qt(0.975,10)
deltaf <- sqrt(GS)*t

ci.lines<-function(){
  
  yv <- f.new
  ci<-deltaf
  uyv<-yv+ci
  lyv<-yv-ci
  lines(x.new,uyv, lty=3, lwd=2, col = "black")
  lines(x.new,lyv, lty=3, lwd=2, col = "black")
    }

par(mfrow=c(1,1),mai=c(0.9,0.9,0.5,0.5),mgp=c(2.0,0.6,0),cex=1.2)
plot(data$Conc,data$Velocity,pch=16,main="", cex.lab=1.2, 
     xlab = "Concentration (ppm)", 
     ylab= expression(Velocity ~ (counts/min^2)), xlim = c(0, 1.2),
     ylim = c(50,220), col = "darkgray", las = 1)
x <- seq(0,1.2,0.01)
curve(212.7*x/(0.06412+ x), add = TRUE, col = "black", lwd =2)
ci.lines()


# Bates & Watts Approach (1988, p. 59)
#-------------------------------------------------------------------------------

# Matrix containing regressors
X = matrix(cbind(rep(1,12),data$Conc),12,2)
X

# Vector containing response values
y <- matrix(data$Velocity)
y

# Gauss-Newton method
# The expectation function eta(theta) is approximated by
# eta(theta) = eta(theta^0) + V^0(theta-theta^0), as a result
# of first order Taylor series
# Residuals zeta(theta) are approximated by
# zeta(theta) = zeta^0 - V^0*delta, with delta = theta-theta^0.
# A Gauss increment delta^0 is calculates to minimize
# the approximate residual sum of squares ||zeta^0-V^0*delta||^2

puromycin = function(Conc,theta1,theta2){
  Vel = theta1*Conc/(theta2+Conc)
  return(as.matrix(Vel))
}

#  Initial theta = (205,0.08)^T
theta = matrix(c(205.0,0.08),2,1)
pred = round(puromycin(data$Conc,theta[1],theta[2]),2)
pred

# Creating a dataframe to reproduce Table 2.1 from 
# Bates & Watts (1988), p. 41
table = as.data.frame(cbind(rep(1:length(data$Conc)), X[,2], y, pred))
names(table) = c("Obs", "x", "y", "predicted_0")

# Calculating residuals and appending them to the 
# table dataframe 
table$residual_0 = table$y - table$predicted

# Derivative matrix
derivative = function(Conc,theta1,theta2){
  d1 = round(Conc/(theta2+Conc),4)
  d2 = round(-theta1*Conc/(theta2+Conc)^2,2)
  return(as.matrix(cbind(d1,d2)))
}

# Calculating derivative matrix  (V^0) for all observations 
# and initial theta
V_0 = derivative(data$Conc,theta[1],theta[2])
V_0

table$der_theta1_0 = V_0[,1]
table$der_theta2_0 = V_0[,2]

# Calculating residual sum of squares for initial iteration
rss0 = sum(table$residual_0^2)

# QR decomposition of the derivative matrix V^0, with  
# the result given in a compact form
QR =qr(V_0)
QR$qr

# Reconstructing Q and R matrices from a QR object.
# In this case we obtain Q1 and R1
Q1 = qr.Q(QR)
Q1
R1 = qr.R(QR)
R1

# Checking that V_0=Q1*R1
V_0_alt = Q1%*%R1
V_0_alt


# Defining a residual vector from the corresponding column
# in table dataframe
z_0 = matrix(table$residual)

# Calculating w1=Q1^T*z_0
w1 = t(Q1)%*%z_0

# We know that R1*delta_0=w1. Then delta_0=inverse(R1)*w1
delta_0 = round(solve(R1)%*%w1,3)
delta_0

#  Updating theta = (213.03,0.063)^T
theta1 = theta + delta_0
table$predicted_1 = round(puromycin(data$Conc,theta1[1],theta1[2]),2)

# Calculating new residuals 
table$residual_1 = table$y - table$predicted_1

# Calculating residual sum of squares for new iteration
rss1 = sum(table$residual_1^2)

# rss1 (1205) is smaller than rss (3155). So we move to theta1
# and perform another iteration.
# The process is repeated until there is  change in the 
# residual sum of squares.
# Convergence is reached at theta = (212.7,0.0641)
theta_opt = matrix(c(212.7,0.0641),2,1)

# Optimal predicted values
table$predicted_opt = round(puromycin(data$Conc,theta_opt[1],theta_opt[2]),2)

# Calculating derivative matrix  (V^0) for all observations 
# and optimal theta
V_0_opt = derivative(data$Conc,theta_opt[1],theta_opt[2])
V_0_opt

table$der_theta1_opt = V_0_opt[,1]
table$der_theta2_opt = V_0_opt[,2]

# Calculating optimal residuals 
table$residual_opt = table$y - table$predicted_opt

# Calculating residual sum of squares for initial iteration
rss_opt = sum(table$residual_opt^2)

# QR decomposition of the derivative matrix V^0, with  
# the result given in a compact form
QR_opt =qr(V_0_opt)
QR_opt$qr

# Reconstructing Q and R matrices from a QR object.
# In this case we obtain Q1 and R1
Q1_opt = qr.Q(QR_opt)
Q1_opt
R1_opt = qr.R(QR_opt)
R1_opt

# Bands for the expected responses at Conc=0.4
y0.4 = round(puromycin(0.4,theta_opt[1],theta_opt[2]),1)
y0.4

# Derivative vector
v = derivative(0.4,theta_opt[1],theta_opt[2])
v

# v^T*inverse(R1)
vTinvR = v%*%solve(R1_opt)
vTinvR

# Norm  of the vector vTinvR
norm1 = sum(vTinvR*vTinvR)^0.5
norm1

# Plotting bands for  expected response
s = sqrt(rss_opt/10)
EFE = qf(0.95,2,10)
P = 2

xnew <- seq(0,1.2,0.01) #range
ynew = theta_opt[1]*xnew/(theta_opt[2] + xnew) 
par(mfrow=c(1,1),mai=c(0.9,0.9,0.5,0.5),mgp=c(2.0,0.6,0), cex=1.2)
plot(data$Conc, data$Velocity, pch=16, col = "darkgray", las =1,
     xlab = "Concentration (ppm)",
     ylab = expression(Velocity ~ (counts/min^2)),
     xlim = c(0,1.2),
     ylim = c(50,220), cex.lab =1.2)
lines(xnew,ynew,lwd=2)
lines(xnew,ynew+s*norm1*sqrt(P*EFE),lwd=2,lty=3,col="black")
lines(xnew,ynew-s*norm1*sqrt(P*EFE),lwd=2,lty=3,col="black")


# Monte Carlo simulation (Spiess, 2013)
#-------------------------------------------------------------------------------

predictNLS_MC <- function(object, var.pred, newdata, level = 0.95, 
                          nsim = 10000, ...)
{
      require(MASS, quietly = TRUE)
      
      ## get right-hand side of formula
      RHS <- as.list(object$call$formula)[[3]]
      EXPR <- as.expression(RHS)
      
      ## coefficients
      COEF <- coef(object)
      
      ## get variance-covariance matrix
      VCOV <- vcov(object)
      
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
            predVAL <- newdata[i, 1.2]
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
            FITTED <- predict(object, newdata = data.frame(PRED))
            MEAN.sim <- mean(EVAL, na.rm = TRUE)
            SD.sim <- sd(EVAL, na.rm = TRUE)
            QUANT <- quantile(EVAL, c((1 - level)/2, level + (1 - level)/2),
                              na.rm=TRUE)
            RES <- c(FITTED, MEAN.sim, SD.sim, QUANT[1], QUANT[2])
            outMAT <- rbind(outMAT, RES)
      }
      
      colnames(outMAT) <- c("fit", "mean", "sd", names(QUANT[1]), 
                            names(QUANT[2]))
      rownames(outMAT) <- NULL
      
      return(outMAT)  
}

Ajustados <- predictNLS_MC(object=nonlinearmod, var.pred = "Conc", 
                           newdata = data.frame(x.new))
Ajustados <- as.data.frame(Ajustados)
plot(data$Conc, data$Velocity, pch=16, col = "darkgray", las =1, cex.lab=1.2,
     xlab = "Concentration (ppm)",
     ylab = expression(Velocity ~ (counts/min^2)),ylim=c(50,220), xlim = c(0,1.2) )
lines(seq(0.01,1.2,0.01),Ajustados$fit,lwd=2)
lines(seq(0.01,1.2,0.01),Ajustados$`2.5%`,lwd=2, lty=3)
lines(seq(0.01,1.2,0.01),Ajustados$`97.5%`,lwd=2, lty=3)


# Graphic representation of the intervals.

plot(data$Conc,data$Velocity,pch=16,main="", cex.lab=1.2, 
     xlab = "Concentration (ppm)", 
     ylab= expression(Velocity ~ (counts/min^2)), xlim = c(0, 1.2),
     ylim = c(50,220), col = "darkgray")
x <- seq(0,1.2,0.01)
curve(212.7*x/(0.06412+ x), add = TRUE, col = "black", lwd =2)
ci.lines<-function(){
  
  yv <- f.new
  ci<-deltaf
  uyv<-yv+ci
  lyv<-yv-ci
  lines(x.new,uyv, lty=5, lwd=2, col = "coral")
  lines(x.new,lyv, lty=5, lwd=2, col = "coral")
}
ci.lines()
lines(xnew,ynew+s*norm1*sqrt(P*EFE),lwd=2,lty=1, col = "yellow3" )
lines(xnew,ynew-s*norm1*sqrt(P*EFE),lwd=2,lty=1, col = "yellow3" )
lines(seq(0.01,1.2,0.01),Ajustados$`2.5%`,lwd=2, lty=5, col = "turquoise4")
lines(seq(0.01,1.2,0.01),Ajustados$`97.5%`,lwd=2, lty=5, col = "turquoise4")
legend(x = "bottomright", c("Delta", "Bates & Watts", "Monte Carlo"),
       lty =c(5,1,5), lwd = c(2,2,2), 
       col = c("coral", "yellow3", "turquoise4"), cex = 0.7,
       title = "Method")




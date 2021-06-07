
# lectura de datos.
data = read.csv("Puromycin.csv")


################################################################################
################ Estimadores de Máxima verosimilitud ###########################
################################################################################

nonlinearmod = nls(Velocity ~ beta1*Conc/(beta2 + Conc), data = data, 
                   start = list(beta1=205, beta2=0.08),
                   trace = TRUE)
summary(nonlinearmod)

# Intervalos de confianza para la respuesta media.


# Método delta
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
  lines(x.new,uyv, lty=3, lwd=2)
  lines(x.new,lyv, lty=3, lwd=2)
    }

plot(data$Conc,data$Velocity,pch=20,main="", cex.lab=1.2, 
     xlab = "Concentración (ppm)", 
     ylab= expression(Velocity ~ (counts/min^2)), xlim = c(0, 1.2),
     ylim = c(50,220), col = "red")
x <- seq(0,1.2,0.01)
curve(212.7*x/(0.06412+ x), add = TRUE, col = "black", lwd =2)
ci.lines()


# Aproximación Bates & Watts (1988, p. 59)
#-------------------------------------------------------------------------------

theta1 = 2.127e+02
theta2 = 6.412e-02 
xnew <- seq(min(data$Conc),1.2,0.01) 
ynew = theta1*xnew/(theta2 + xnew) 
par(mfrow=c(1,1),mai=c(0.9,0.9,0.5,0.5),mgp=c(2.0,0.6,0),cex=1.2)
plot(data$Conc, data$Velocity, pch=20, col = "red", las =1,
     xlab = "Concentration (ppm)",
     ylab = expression(Velocity ~ (counts/min^2)), ylim = c(50,220),
     xlim = c(0,1.2))
lines(xnew,ynew,lwd=2)
lines(xnew,ynew+summary(nonlinearmod)$sigma,lwd=2,lty=3)
lines(xnew,ynew-summary(nonlinearmod)$sigma,lwd=2,lty=3)


# Simulación de Monte Carlo
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
plot(data$Conc, data$Velocity, pch=20, col = "red", las =1, cex.lab=1.2,
     xlab = "Concentration (ppm)",
     ylab = expression(Velocity ~ (counts/min^2)),ylim=c(45,220))
lines(seq(0.01,1.2,0.01),Ajustados$fit,lwd=2)
lines(seq(0.01,1.2,0.01),Ajustados$`2.5%`,lwd=2, lty=3)
lines(seq(0.01,1.2,0.01),Ajustados$`97.5%`,lwd=2, lty=3)



################################################################################
################ Estimadores de Máxima verosimilitud ###########################
################################################################################

# Estimador de máxima verosimilitud

library(optimr)

loglik.Puromycin= function(parameters, C, V)
{
  logL =  tryCatch(
    {
      beta1 = parameters["b1"]
      beta2 = parameters["b2"] 
      sigma =  parameters["sigma"]
      model =  beta1*C/(beta2 + C)
      logL <- sum(dnorm(V, model, sigma, log = TRUE)) 
    },
    error = function(cond)
    {
      logL = -999999
    })
  print(logL)
  return(-logL)
}

parameters = c(b1 = 212.7, b2 =0.06412, sigma = 1)
modelML = with(subset(data,!is.na(Velocity) & !is.na(Conc)),
               optimr(par = parameters,
                      fn  = loglik.Puromycin,
                      C   = Conc,
                      V   = Velocity,
                      control = list(maxiter=1000),
                      method = "Nelder-Mead",
                      hessian = TRUE))

# Matriz de información de fisher

FIM <- solve(modelML$hessian)

vcov_ml_R <- FIM[-3,-3] 

# Desviación estándar

se <- sqrt(diag(FIM))

# Cuántil t
t <- modelML$par/se

# Valor P
pval <- 2*(1-pt(abs(t),length(data$Conc)-3))

results <- as.matrix(cbind(modelML$par,se,t,pval))
colnames (results) <- c("parameter","se","t","p")
rownames (results) <- c("beta1","beta2","sigma")
print(results,digits=5)

# Estimadores paramétricos.
betas <- modelML$par


## Intervalos de confianza para la respuesta media.


# Método delta
#-------------------------------------------------------------------------------

fgh2 <- deriv(Velocity ~ b1*Conc/(b2 + Conc), c("b1", "b2"), 
              function(b1,b2,Conc){} ) 

x.new <- seq(0.01, 1.2, by=0.01)
f.new <- fgh2(222.50792836,0.07393871, x.new)

g.new <- attr(f.new,"gradient")
V.beta2 <- vcov_ml_R
GS=rowSums((g.new%*%V.beta2)*g.new)

s.t <- qt(0.957, 9) # Sigma también cuenta cómo parámetro
deltaf <- sqrt(GS)*s.t

ci.lines<-function(){
  
  yv <- f.new
  ci<-deltaf
  uyv<-yv+ci
  lyv<-yv-ci
  lines(x.new,uyv,lty=3, lwd=2)
  lines(x.new,lyv,lty=3, lwd=2)
}

# Intervalos de confianza para la respuesta media
plot(data$Conc,data$Velocity,pch=20,main="", cex.lab=1.2, 
     xlab = "Concentración (ppm)", 
     ylab= expression(Velocity ~ (counts/min^2)), xlim = c(0, 1.2),
     ylim = c(50,220), col = "red")
x <- seq(0,1.2,0.01)
curve(222.50792836*x/(0.07393871 + x), add = TRUE, col = "black", lwd=2)
ci.lines()


# Aproximación Bates & Watts (1988, p. 59)
#-------------------------------------------------------------------------------

theta1 = 222.50792836
theta2 = 0.07393871
xnew <- seq(min(data$Conc),1.2,0.01) 
ynew = theta1*xnew/(theta2 + xnew) 
par(mfrow=c(1,1),mai=c(0.9,0.9,0.5,0.5),mgp=c(2.0,0.6,0),cex=1.2)
plot(data$Conc, data$Velocity, pch=20, col = "red", las =1,
     xlab = "Concentration (ppm)",
     ylab = expression(Velocity ~ (counts/min^2)), ylim = c(50,220),
     xlim = c(0,1.2))
lines(xnew,ynew,lwd=2)
lines(xnew,ynew+summary(nonlinearmod)$sigma,lwd=2,lty=3)
lines(xnew,ynew-summary(nonlinearmod)$sigma,lwd=2,lty=3)


# Simulación de Monte Carlo 
#-------------------------------------------------------------------------------

FittedML <- function(Conc){
      Velocidad = modelML$par[1]*Conc/(modelML$par[2] + Conc)
      return(Velocidad)
}

predictNLS_MV <- function(EXPR, modelo, var.pred, newdata, level = 0.95, 
                          nsim = 10000, ...)
{
      require(MASS, quietly = TRUE)
      
      ## coefficients
      COEF <- modelo$par[-3]
      
      ## get variance-covariance matrix
      FIM <- solve(modelo$hessian)
      VCOV <- FIM[-3,-3] 
      
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
            RES <- c(FITTED, MEAN.sim, SD.sim, QUANT[1], QUANT[2])
            outMAT <- rbind(outMAT, RES)
      }
      
      colnames(outMAT) <- c("fit", "mean", "sd", names(QUANT[1]), 
                            names(QUANT[2]))
      rownames(outMAT) <- NULL
      
      return(outMAT)  
}

Ajustados <- predictNLS_MV(EXPR=expression(b1*Conc/(b2 + Conc)), modelo=modelML, 
                           var.pred = "Conc", newdata = data.frame(x.new))
Ajustados <- as.data.frame(Ajustados)
plot(data$Conc, data$Velocity, pch=20, col = "red", las =1,
     xlab = "Concentration (ppm)",
     ylab = expression(Velocity ~ (counts/min^2)),ylim=c(45,220))
lines(seq(0.01,1.2,0.01),Ajustados$fit,lwd=2)
lines(seq(0.01,1.2,0.01),Ajustados$`2.5%`,lwd=2, lty=3)
lines(seq(0.01,1.2,0.01),Ajustados$`97.5%`,lwd=2, lty=3)





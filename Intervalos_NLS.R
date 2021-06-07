
# lectura de datos.
data = read.csv("Puromycin.csv")


################################################################################
################ Estimadores de Mínimos Cuadrados ##############################
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
  lines(x.new,uyv, lty=3, lwd=2, col = "black")
  lines(x.new,lyv, lty=3, lwd=2, col = "black")
    }

plot(data$Conc,data$Velocity,pch=16,main="", cex.lab=1.2, 
     xlab = "Concentración (ppm)", 
     ylab= expression(Velocity ~ (counts/min^2)), xlim = c(0, 1.2),
     ylim = c(50,220), col = "darkgray")
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
plot(data$Conc, data$Velocity, pch=16, col = "darkgray", las =1,
     xlab = "Concentration (ppm)",
     ylab = expression(Velocity ~ (counts/min^2)), ylim = c(50,220),
     xlim = c(0,1.2))
lines(xnew,ynew,lwd=2)
lines(xnew,ynew+summary(nonlinearmod)$sigma,lwd=2,lty=3, col = "black")
lines(xnew,ynew-summary(nonlinearmod)$sigma,lwd=2,lty=3, col = "black")


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
plot(data$Conc, data$Velocity, pch=16, col = "darkgray", las =1, cex.lab=1.2,
     xlab = "Concentration (ppm)",
     ylab = expression(Velocity ~ (counts/min^2)),ylim=c(45,220))
lines(seq(0.01,1.2,0.01),Ajustados$fit,lwd=2)
lines(seq(0.01,1.2,0.01),Ajustados$`2.5%`,lwd=2, lty=3)
lines(seq(0.01,1.2,0.01),Ajustados$`97.5%`,lwd=2, lty=3)


# Representación gráfica de los intervalos.

plot(data$Conc,data$Velocity,pch=16,main="", cex.lab=1.2, 
     xlab = "Concentración (ppm)", 
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
lines(xnew,ynew+summary(nonlinearmod)$sigma,lwd=2,lty=1, col = "yellow3" )
lines(xnew,ynew-summary(nonlinearmod)$sigma,lwd=2,lty=1, col = "yellow3" )
lines(seq(0.01,1.2,0.01),Ajustados$`2.5%`,lwd=2, lty=5, col = "turquoise4")
lines(seq(0.01,1.2,0.01),Ajustados$`97.5%`,lwd=2, lty=5, col = "turquoise4")
legend(x = "bottomright", c("Delta", "Bates & Watts", "Monte Carlo"),
       lty =c(5,1,5), lwd = c(2,2,2), 
       col = c("coral", "yellow3", "turquoise4"), cex = 0.7,
       title = "Método")




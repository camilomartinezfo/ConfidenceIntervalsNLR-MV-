## Simulación Monte Carlo para modelo con mínimos cuadrados

predictNLS <- function(object, var.pred, newdata, level = 0.95, nsim = 10000, ...)
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
            simMAT <- mvrnorm(n = nsim, mu = MU, Sigma = newVCOV, empirical = TRUE)
            
            ## evaluate expression on rows of simMAT
            EVAL <- eval(EXPR, envir = as.data.frame(simMAT))
            
            ## collect statistics
            PRED <- data.frame(predVAL)
            colnames(PRED) <- var.pred
            FITTED <- predict(object, newdata = data.frame(PRED))
            MEAN.sim <- mean(EVAL, na.rm = TRUE)
            SD.sim <- sd(EVAL, na.rm = TRUE)
            QUANT <- quantile(EVAL, c((1 - level)/2, level + (1 - level)/2),na.rm=TRUE)
            RES <- c(FITTED, MEAN.sim, SD.sim, QUANT[1], QUANT[2])
            outMAT <- rbind(outMAT, RES)
      }
      
      colnames(outMAT) <- c("fit", "mean", "sd", names(QUANT[1]), names(QUANT[2]))
      rownames(outMAT) <- NULL
      
      return(outMAT)  
}

Ajustados <- predictNLS(object=nonlinearmod, var.pred = "Conc", newdata = data.frame(x.new))
Ajustados <- as.data.frame(Ajustados)
plot(data$Conc, data$Velocity, pch=20, col = "red", las =1,
     xlab = "Concentration (ppm)",
     ylab = expression(Velocity ~ (counts/min^2)),ylim=c(45,220))
lines(seq(0.01,1.2,0.01),Ajustados$fit,lwd=2)
lines(seq(0.01,1.2,0.01),Ajustados$`2.5%`,lwd=2, lty=3)
lines(seq(0.01,1.2,0.01),Ajustados$`97.5%`,lwd=2, lty=3)



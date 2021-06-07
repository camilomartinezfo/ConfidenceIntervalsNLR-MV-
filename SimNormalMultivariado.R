nonlinearmod


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
      colnames(ADD1) <- names(var.pred)
      VCOV <- cbind(VCOV, ADD1)
      ADD2 <- c(rep(0, NCOL + 1))
      ADD2 <- matrix(ADD2, nrow = 1)
      rownames(ADD2) <- names(var.pred)
      VCOV <- rbind(VCOV, ADD2) 
      
      NR <- nrow(newdata)
      varPLACE <- ncol(VCOV)
      outMAT <- NULL  
      
      for (i in 1:NR) {
            
            ## get predictor values and optional errors
            predVAL <- newdata[i, 1]
            predERROR <- 0
            names(predVAL) <- names(var.pred)  
            names(predERROR) <- names(var.pred) 
            
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
            colnames(PRED) <- names(var.pred)   
            FITTED <- predict(object, newdata = data.frame(PRED))
            MEAN.sim <- mean(EVAL, na.rm = TRUE)
            SD.sim <- sd(EVAL, na.rm = TRUE)
            QUANT <- quantile(EVAL, c((1 - level)/2, level + (1 - level)/2),na.rm=TRUE)
            RES <- c(FITTED, MEAN.sim, SD.sim, QUANT[1], QUANT[2])
            outMAT <- rbind(outMAT, RES)
      }
      
      colnames(outMAT) <- c("fit", "mean", "sd", "median", "mad", names(QUANT[1]), names(QUANT[2]))
      rownames(outMAT) <- NULL
      
      return(outMAT)  
}

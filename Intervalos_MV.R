# lectura de datos.
data = read.csv("Puromycin.csv")


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
f.new <- fgh2(modelML$par[1],modelML$par[2], x.new)

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
      lines(x.new,uyv,lty=3, lwd=2, col="black")
      lines(x.new,lyv,lty=3, lwd=2, col="black")
}

# Intervalos de confianza para la respuesta media
par(mfrow=c(1,1),mai=c(0.9,0.9,0.5,0.5),mgp=c(2.0,0.6,0),cex=1.2)
plot(data$Conc,data$Velocity,pch=16,main="", cex.lab=1.2, 
     xlab = "Concentración (ppm)", 
     ylab= expression(Velocity ~ (counts/min^2)), xlim = c(0, 1.2),
     ylim = c(50,220), col = "darkgray")
x <- seq(0,1.2,0.01)
curve(modelML$par[1]*x/(modelML$par[2] + x), add = TRUE, col = "black", lwd=2)
ci.lines()


# Aproximación Bates & Watts (1988, p. 59)
#-------------------------------------------------------------------------------

# Matriz de regresores
X = matrix(cbind(rep(1,12),data$Conc),12,2)
X

# Vector de los valores de respuesta observados
y <- matrix(data$Velocity)
y

# Método de Gauss-Newton 

puromycin = function(Conc,theta1,theta2){
   Vel = theta1*Conc/(theta2+Conc)
   return(as.matrix(Vel))
}

#  Teta optimo = los de max verosimilitud
theta = matrix(c(modelML$par[1],modelML$par[2]),2,1)
pred = round(puromycin(data$Conc,theta[1],theta[2]),2)
pred

# Tabla 2.1 de Bates & Watts (1988), p. 41
table = as.data.frame(cbind(rep(1:length(data$Conc)), X[,2], y, pred))
names(table) = c("Obs", "x", "y", "predicted_0")

# Calcular los residuos y agregarlos a la tabla
table$residual_0 = table$y - table$predicted

# Matriz derivada
derivative = function(Conc,theta1,theta2){
   d1 = round(Conc/(theta2+Conc),4)
   d2 = round(-theta1*Conc/(theta2+Conc)^2,2)
   return(as.matrix(cbind(d1,d2)))
}

# Se calcula la matriz derivada para todas las observaciones 
# y el teta inicial.

V_0 = derivative(data$Conc,theta[1],theta[2])
V_0

table$der_theta1_0 = V_0[,1]
table$der_theta2_0 = V_0[,2]

# Cálculo de la suma de cuadrados residual para la iteración inicial.
rss0 = sum(table$residual_0^2)

# Descomposición QR de la matriz derivada V^0, con
# el resultado dado en forma compacta.
QR =qr(V_0)
QR$qr

# Reconstrucción de matrices Q y R a partir de un objeto QR.
# En este caso obtenemos Q1 y R1
Q1 = qr.Q(QR)
Q1
R1 = qr.R(QR)
R1

# Bandas para las respuestas esperadas en Conc = 0.4
y0.4 = round(puromycin(0.4,theta[1],theta[2]),1)
y0.4

# Vector derivado
v = derivative(0.4,theta[1],theta[2])
v

# v^T*inversa(R1)
vTinvR = v%*%solve(R1)
vTinvR

# Norma del vector v en vR/magnitud del vector.
norm1 = sum(vTinvR*vTinvR)^0.5
norm1

# Trazar bandas para la respuesta esperada
s = sqrt(rss0/10)
EFE = qf(0.95,2,10)
P = 2

xnew <- seq(0,1.2,0.01) #range
ynew = theta[1]*xnew/(theta[2] + xnew) 
par(mfrow=c(1,1),mai=c(0.9,0.9,0.5,0.5),mgp=c(2.0,0.6,0), cex=1.2)
plot(data$Conc, data$Velocity, pch=16, col = "darkgray", las =1,
     xlab = "Concentration (ppm)",
     ylab = expression(Velocity ~ (counts/min^2)),
     xlim = c(0,1.2),
     ylim = c(50,220), cex.lab =1.2)
lines(xnew,ynew,lwd=2)
lines(xnew,ynew+s*norm1*sqrt(P*EFE),lwd=2,lty=3,col="black")
lines(xnew,ynew-s*norm1*sqrt(P*EFE),lwd=2,lty=3,col="black")



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
plot(data$Conc, data$Velocity, pch=16, col = "darkgray", las =1,cex.lab=1.2,
     xlab = "Concentration (ppm)",
     ylab = expression(Velocity ~ (counts/min^2)),ylim=c(45,220))
lines(seq(0.01,1.2,0.01),Ajustados$fit,lwd=2)
lines(seq(0.01,1.2,0.01),Ajustados$`2.5%`,lwd=2, lty=3)
lines(seq(0.01,1.2,0.01),Ajustados$`97.5%`,lwd=2, lty=3)




# Representación gráfica de los intervalos.

ci.lines<-function(){
      
      yv <- f.new
      ci<-deltaf
      uyv<-yv+ci
      lyv<-yv-ci
      lines(x.new,uyv,lty=5, lwd=2, col="coral")
      lines(x.new,lyv,lty=5, lwd=2, col="coral")
}


plot(data$Conc,data$Velocity,pch=16,main="", cex.lab=1.2, 
     xlab = "Concentración (ppm)", 
     ylab= expression(Velocity ~ (counts/min^2)), xlim = c(0, 1.2),
     ylim = c(50,220), col = "darkgray")
x <- seq(0,1.2,0.01)
curve(modelML$par[1]*x/(modelML$par[2] + x), add = TRUE, col = "black", lwd=2)
ci.lines()

lines(xnew,ynew+s*norm1*sqrt(P*EFE),lwd=2,lty=1,col="yellow3")
lines(xnew,ynew-s*norm1*sqrt(P*EFE),lwd=2,lty=1,col="yellow3")

lines(seq(0.01,1.2,0.01),Ajustados$`2.5%`,lwd=2, lty=5,col="turquoise4")
lines(seq(0.01,1.2,0.01),Ajustados$`97.5%`,lwd=2, lty=5,col="turquoise4")

legend(x="bottomright",c("Delta","Bates & Watts","Monte Carlo"),lty=c(5,1,5),
       lwd=c(2,2,2), col=c("coral","yellow3","turquoise4"),cex=0.7,title="Método")


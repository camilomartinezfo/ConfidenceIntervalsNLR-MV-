
# lectura de datos.
data = read.csv("Puromycin.csv")

# Estimador de mínimos cuadrados ordinarios no lineales.

nonlinearmod = nls(Velocity ~ beta1*Conc/(beta2 + Conc), data = data, 
                   start = list(beta1=205, beta2=0.08),
                   trace = TRUE)
summary(nonlinearmod)

# Intervalos de confianza para la respuesta media.

# Método delta

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

ci.lines<-function(nonlinearmod){
  
  yv <- f.new
  ci<-deltaf
  uyv<-yv+ci
  lyv<-yv-ci
  lines(x.new,uyv,lty=3)
  lines(x.new,lyv,lty=3)
    }

plot(data$Conc,data$Velocity,pch=20,main="", cex.lab=1.5,
     cex.main=1.5, xlab = "Concentración (ppm)", 
     ylab= expression(Velocity ~ (counts/min^2)), xlim = c(0, 1.2),
     ylim = c(50,210), col = "red")
x <- seq(0,1.2,0.01)
curve(212.7*x/(0.06412+ x), add = TRUE, col = "black", lwd =2)
ci.lines(nonlinearmod)

# Aproximación Bates & Watts (1988, p. 59)

theta1 = 2.127e+02
theta2 = 6.412e-02 
xnew <- seq(min(data$Conc),max(data$Conc),0.01) 
ynew = theta1*xnew/(theta2 + xnew) 
par(mfrow=c(1,1),mai=c(0.9,0.9,0.5,0.5),mgp=c(2.0,0.6,0),cex=1.2)
plot(data$Conc, data$Velocity, pch=20, col = "red", las =1,
     xlab = "Concentration (ppm)",
     ylab = expression(Velocity ~ (counts/min^2)))
lines(xnew,ynew,lwd=2)
lines(xnew,ynew+summary(nonlinearmod)$sigma,lwd=2,lty=3)
lines(xnew,ynew-summary(nonlinearmod)$sigma,lwd=2,lty=3)

# Simulación de Monte Carlo



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

# We then maximize this function using optim.
parameters = c(b1 = 212.7, b2 =0.06412, sigma = 1)
modelML = with(subset(data,!is.na(Velocity) & !is.na(Conc)),
               optimr(par = parameters,
                      fn  = loglik.Puromycin,
                      C   = Conc,
                      V   = Velocity,
                      control = list(maxiter=1000),
                      method = "Nelder-Mead",
                      hessian = TRUE))


FIM <- solve(modelML$hessian)


se <- sqrt(diag(FIM))

t <- modelML$par/se

pval <- 2*(1-pt(abs(t),length(my_data$dbh.cm)-4))

# Displaying results
results <- as.matrix(cbind(modelML$par,se,t,pval))
colnames (results) <- c("parameter","se","t","p")
rownames (results) <- c("beta1","beta2","sigma")
print(results,digits=5)

# Storing estimated coefficients
etas <- modelML$par
etas



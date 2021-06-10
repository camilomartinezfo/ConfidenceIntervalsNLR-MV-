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


## Maxima verosimilitud
loglik.sweetgum = function(parameters, year, D)
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
                         fn  = loglik.sweetgum,
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

s.t <- qt(0.975,45)
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


plot(data$year,data$D,xlim=c(0,200),ylim=c(0,150),pch=16,bg="white",
     col="darkgray",cex=1,lwd=1.9,xlab="Age (years)",ylab="Diameter (cm)",cex.lab=1.2)
ci.lines()






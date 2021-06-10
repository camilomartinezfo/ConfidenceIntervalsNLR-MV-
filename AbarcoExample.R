library(optimr)

datos <- read.csv2("abarcodata.csv")
datos$dD <- as.numeric(datos$dD)
datos$D <- as.numeric(datos$D)
datos$Drelativa <- as.numeric(datos$Drelativa)
datos$D <- datos$D/10
datos$dD <- datos$dD/10
datos$Drelativa <- (datos$dD)/datos$D
pos <- which(datos$D>150)
datos <- datos[-c(pos),]
pos2 <- which(datos$Drelativa==1)
datos <- datos[-c(pos2),]


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

modelo1.ML = with(subset(datos,!is.na(D) & !is.na(year)),
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

pval <- 2*(1-pt(abs(t),length(datos$year)-4))

# Displaying results
results <- as.matrix(cbind(modelo1.ML$par,se,t,pval))
colnames (results) <- c("parameter","se","t","p")
rownames (results) <- c("A","k","c","sigma")
print(results,digits=5)



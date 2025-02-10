rm(list=ls()) 

# sets working directories:
setwd("C:\\Users\\alici\\OneDrive\\Escritorio\\R - Bayesian statistics\\r_proy\\Code\\CaseStudies\\BART")

library("coda")
library("R2jags")





ntrials <- 90 # Number of trials for the BART
p = 0.15

Data   <- list(
  matrix (data = as.numeric (as.matrix (read.table ("BillSober.txt"))[-1,]), 90, 8),
  matrix (data = as.numeric (as.matrix (read.table ("GeorgeSober.txt"))[-1,]), 90, 8)
)

head(Data)

npers <- as.numeric(length(Data))	
cash  <- npumps <- matrix (, npers, ntrials)
d     <- array ( , c(npers, ntrials, 30)) # Data in binary format

for (i in 1:npers){
  
  cash[i,]   <- (Data[[i]][,7]!=0)*1	# Cash or burst?
  npumps[i,] <- Data[[i]][,6]				  # Nr. of pumps
  
  
  for (j in 1:ntrials){
    if (npumps[i,j]>0) {d[i, j, 1:npumps[i,j]] <- rep(0, npumps[i,j])}
    if (cash[i,j]==1) {d[i, j, (npumps[i,j]+1)] <- 1}
  }
}

options <- cash + npumps  # Nr. of decision possibilities

mod_string = "model{
  
  # Priors
  mug ~ dunif(0,10)
  sigmag ~ dunif(0,10)
  lambdag = 1/(sigmag^2)

  beta ~ dunif(0,10)

  
  # Choice Data
  for (i in 1:npers){
    gplus[i] ~ dnorm(mug,lambdag)T(0,)
    omega[i] = -gplus[i]/log(1-p)
      
      for (j in 1:ntrials){
        for (k in 1:options[i,j]){
          theta[i,j,k] = 1-(1/(1+max(-15,min(15,exp(beta*(k-omega[i]))))))
          d[i,j,k] ~ dbern(theta[i,j,k])
        }
      }
    }
}"

set.seed(113)

#list(npers=npers, p=p, ntrials=ntrials, options=options, d=d) 

data_jags <- list("npers", "p", "ntrials", "options", "d") # to be passed on to JAGS

myinits <- list(
  list(mug = 1.2, sigmag = 0.1, mub = 0.8, sigmab = 0.8),
  list(mug = 1.5, sigmag = 0.2, mub = 1.0, sigmab = 1.2))



# parameters to be monitored:
parameters <- c("beta","mug", "sigmag", "gplus", "mub", "sigmab", "omega")


samples_gb <- jags(data_jags, inits=myinits, parameters, model.file = "modelo_gb.txt", 
                n.chains=2, n.iter=5000, n.burnin=2000, n.thin=1)

autojags(samples_gb)
traceplot(samples_gb)
samples_gb



gplus.b <- samples_gb$BUGSoutput$sims.list$gplus[,1]
gplus.g <- samples_gb$BUGSoutput$sims.list$gplus[,2]

omega.b  <- samples_gb$BUGSoutput$sims.list$omega[,1]
omega.g  <- samples_gb$BUGSoutput$sims.list$omega[,2]

beta.b  <- samples_gb$BUGSoutput$sims.list$beta[,1]
beta.g  <- samples_gb$BUGSoutput$sims.list$beta[,2]

dev.off()


par (cex.main = 2.5, cex.lab = 1.5, cex.axis = 1.5, mar = c(4, 5, 3, 1), las = 1)
layout (matrix (1:8, 2, 4, byrow = T))

hist(npumps[1,], xlab=" ", main = "B: # pumps", breaks=c(0:max(npumps[1,])), xlim = range(c(0:9)), col="lightblue", axes=F)
axis(2)
axis(1, at=seq(.5,8.5,by=1), labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9"))

plot(density(gplus.b), xlab = expression (gamma^'+'),
     main = expression (paste ("B: Posterior ", gamma^'+')), xlim = c(0.5,2.0), bty = 'n')

plot (density(omega.b), xlab = expression (omega),
      main = expression (paste ("B: Posterior ", omega)), xlim = c(2 ,12), bty = 'n')
plot (density(beta.b), xlab = expression (beta),
      main = expression (paste ("BG: Posterior ", beta)), xlim = c(0.2, 1.2), bty = 'n')



hist(npumps[2,], xlab=" ", main = "G: # pumps", breaks=c(0:max(npumps[2,])), xlim = range(c(0:9)), col="lightblue", axes=F)
axis(2)
axis(1, at=seq(.5,8.5,by=1), labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9"))

plot(density(gplus.g), xlab = expression (gamma^'+'),
     main = expression (paste ("G: Posterior ", gamma^'+')), xlim = c(0.5,2.0), bty = 'n')

plot (density(omega.g), xlab = expression (omega),
      main = expression (paste ("G: Posterior ", omega)), xlim = c(2 ,12), bty = 'n')
plot (density(beta.g), xlab = expression (beta),
      main = expression (paste ("BG: Posterior ", beta)), xlim = c(0.1, 1.2), bty = 'n')

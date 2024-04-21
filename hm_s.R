rm(list=ls()) 

# sets working directories:
setwd("C:\\Users\\alici\\OneDrive\\Escritorio\\R - Bayesian statistics\\r_proy\\Code\\CaseStudies\\BART")

library("coda")
library("R2jags")


ntrials <- 180 # Number of trials for the BART
p = 0.15

Data   <- matrix (data = as.numeric (as.matrix (read.table ("Bill_George_Sober.txt"))[-1,]), 180, 9)

d     <- array ( , c(ntrials, 30)) # Data in binary format

cash  <- (Data[,7]!=0)*1	# Cash or burst?
npumps <- Data[,6]					  # Nr. of pumps
options <- cash + npumps 

for (j in 1:ntrials)
{
  if (npumps[j]>0) {d[j, 1:npumps[j]] <- rep (0, npumps[j])}
  if (cash[j]==1) {d[j, (npumps[j]+1)] <- 1}
}


set.seed(113)


data_jags <- list("p", "ntrials", "options", "d") # to be passed on to JAGS

myinits <-	list(
  list(gplus = 1.2, beta = 0.5),
  list(gplus = 1, beta = 0.3)
  )

# parameters to be monitored:
parameters <- c("beta","gplus", "omega")

samples_s <- jags(data_jags, inits=myinits, parameters, model.file = "modelo_simple.txt", 
                  n.chains=2, n.iter=7000, n.burnin=3500, n.thin=1)


autojags(samples_s)
traceplot(samples_s)

samples_s


gplus.m <- samples_s$BUGSoutput$sims.list$gplus
beta.m  <- samples_s$BUGSoutput$sims.list$beta
omega.m  <- samples_s$BUGSoutput$sims.list$omega


dev.off()


par (cex.main = 2.5, cex.lab = 1.5, cex.axis = 1.5, mar = c(3, 5, 5, 1), las = 1)
layout (matrix (1:4, 1, 4, byrow = T))

hist(npumps, xlab=" ", main = "# pumps", breaks=c(0:max(npumps)), xlim = range(c(0:9)), col="lightblue", axes=F)
axis(2)
axis(1, at=seq(.5,8.5,by=1), labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9"))

plot(density(gplus.m), xlab = expression (gamma^'+'),
     main = expression (paste ("Posterior ", gamma^'+')), xlim = c(0.5,2.0), bty = 'n')

plot (density(omega.m), xlab = expression (omega),
      main = expression (paste ("Posterior ", omega)), xlim = c(2,12), bty = 'n')

plot (density(beta.m), xlab = expression (beta),
      main = expression (paste ("Posterior ", beta)), xlim = c(0.1,1.2), bty = 'n')

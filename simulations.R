#######
library("splines2")
library("quadprog")
library("doParallel")
#library("diptest")
#library("moments")
#library("multimode")
#library(doSNOW)
#library("benchden") #claw distribution
library(doRNG)
#######


n <- 200 #sample size
nsim <- 200 #simulation size
B <- 400 #bootstrap sample size for each sample
alpha <- 0.05 
eps_values <- c(0.001, 0.01, 0.1, 1) #candidate epsilon

ncores <- max(1, parallel::detectCores() - 1) #number of cores registered for parallel
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

results <- vector("list", length(eps_values))

for (i in seq_along(eps_values)) {
  eps <- eps_values[i]
  t1 <- Sys.time()
  sim_result <- foreach(
    r = 1:nsim,
    .combine = "rbind",
    .packages = c("splines2", "quadprog", "moments",
                  "multimode", "benchden"),
    .options.RNG = 123
  ) %dorng% {
    z <- rbinom(n, 1, 0.6)
    x <- rep(0, n)
    x[z == 1] <- rnorm(sum(z == 1), 0, 1)
    x[z == 0] <- rnorm(sum(z == 0), 3.5, 1)
    ans <- bmodetest(x,B = B,cv = TRUE,parallel = FALSE,eps = eps)
    c(pvalue = ans$pvalue, reject = as.numeric(ans$pvalue < alpha))
  }
  elapsed <- difftime(Sys.time(), t1, units = "secs")
  results[[i]] <- data.frame(
    eps = eps,
    power = mean(sim_result[, "reject"]),
    elapsed_seconds = as.numeric(elapsed)
  )
  cat(paste("eps =", eps,"completed\n"))
}
parallel::stopCluster(cl)
power_results <- do.call(rbind, results)
power_results
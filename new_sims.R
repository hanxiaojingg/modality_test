library("splines2")
library("quadprog")
library("doParallel")
library("diptest")
library("moments")
library("multimode")
library(doSNOW)
library("benchden")

set.seed(1234)
n=200
uu=runif(n*1000);x=rnorm(n*1000,0,1);x[uu>0.8]=rnorm(sum(uu>0.8),0,3);
x[uu<0.4]=rnorm(sum(uu<0.4),0,1)
x=matrix(x,1000,n)
out <- matrix(0,1000,6)
for(i in 1:1000){
  out[i,1]=bmodetest(x[i,],B=1000,eps1=5,eps2=2,od=2/7)$pvalue
}

n=200
d = seq(0,5,0.2)
d=2
u = runif(n)
x = rnorm(n,0,1)
x[u<0.4] = rnorm(sum(u<0.4),d,1)
hist(x,breaks = 20)
out1 = modetest(x, method = "ACR")
out2 = modetest(x, method = "SI")
out3 = modetest(x, method = "HY", lowsup=quantile(x,.01),uppsup=quantile(x,.99))
out4 = modetest(x, method = "FM")
out5 = modetest(x, method = "HH")
out6 = modetest(x, method = "CH")
out1$p.value
out2$p.value
out3$p.value
out4$p.value
out5$p.value
out6$p.value

n=200
reps = 10
d = seq(0,5,1)
out = matrix(0,6,length(d))
for(i in 1:length(d)){
  set.seed(1234)
  u = runif(n*reps)
  x = rnorm(n*reps,0,1)
  x[u<0.4] = rnorm(sum(u<0.4),d[i],1)
  x = matrix(x,reps,n)
  pval = matrix(0,6,reps)
  for(j in 1:reps){
    pval[1,j] = modetest(x[j,], method = "ACR")$p.value
    pval[2,j] = modetest(x[j,], method = "SI")$p.value
    pval[3,j] = modetest(x[j,], method = "HY", lowsup=quantile(x[j,],.01),uppsup=quantile(x[j,],.99))$p.value
    pval[4,j] = modetest(x[j,], method = "FM")$p.value
    pval[5,j] = modetest(x[j,], method = "HH")$p.value
    pval[6,j] = modetest(x[j,], method = "CH")$p.value
  }
  out[,i] = apply(pval, 1, function(x){sum(x<0.05)/reps})
}


n = 200
reps = 10
d = seq(0, 5, 1)
out = matrix(0, nrow = 6, ncol = length(d))

for (i in seq_along(d)) {
  set.seed(1234)
  u = runif(n * reps)
  x = rnorm(n * reps, 0, 1)
  x[u < 0.4] = rnorm(sum(u < 0.4), d[i], 1)
  x = matrix(x, reps, n)
  pval = matrix(0, nrow = 6, ncol = reps)
  hy_bounds = t(apply(x, 1, quantile, probs = c(0.01, 0.99)))
  
  for (j in seq_len(reps)) {
    xj = x[j, ]
    pval[1, j] = modetest(xj, method = "ACR")$p.value
    pval[2, j] = modetest(xj, method = "SI")$p.value
    pval[3, j] = modetest(
      xj,
      method = "HY",
      lowsup = hy_bounds[j, 1],
      uppsup = hy_bounds[j, 2]
    )$p.value
    pval[4, j] = modetest(xj, method = "FM")$p.value
    pval[5, j] = modetest(xj, method = "HH")$p.value
    pval[6, j] = modetest(xj, method = "CH")$p.value
  }
  out[, i] = rowMeans(pval < 0.05)
}





######parallel#######
######norm_mix1######
#####################

### step 1 ###
library(foreach)
library(doParallel)

n = 200
n = 800
reps = 500
d=c(0,1,seq(2,4,0.1))
dist_names = c("mixnorm")   # later extend this
methods = c("ACR", "SI", "HY", "FM", "HH", "CH")

scenarios = expand.grid(
  dist = dist_names,
  d_index = seq_along(d),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

### step 2 ###
generate_x = function(dist_name, di, n, reps) {
  if (dist_name == "mixnorm") {
    set.seed(1234)
    u = runif(n * reps)
    x = rnorm(n * reps, 0, 1)
    x[u < 0.4] = rnorm(sum(u < 0.4), di, 1)
    x = matrix(x, reps, n)
    return(x)
  }
  
  stop("Unknown distribution")
}

### step 3 ###
x_list = vector("list", nrow(scenarios))

for (s in seq_len(nrow(scenarios))) {
  dist_name = scenarios$dist[s]
  di = d[scenarios$d_index[s]]
  
  x_list[[s]] = generate_x(dist_name, di, n, reps)
}

### step 4 ###
run_methods = function(x, methods) {
  reps = nrow(x)
  pval = matrix(0, nrow = length(methods), ncol = reps)
  
  hy_bounds = t(apply(x, 1, quantile, probs = c(0.01, 0.99)))
  
  for (j in seq_len(reps)) {
    xj = x[j, ]
    
    pval[1, j] = modetest(xj, method = "ACR")$p.value
    pval[2, j] = modetest(xj, method = "SI")$p.value
    pval[3, j] = modetest(
      xj,
      method = "HY",
      lowsup = hy_bounds[j, 1],
      uppsup = hy_bounds[j, 2]
    )$p.value
    pval[4, j] = modetest(xj, method = "FM")$p.value
    pval[5, j] = modetest(xj, method = "HH")$p.value
    pval[6, j] = modetest(xj, method = "CH")$p.value
  }
  
  rowMeans(pval < 0.05)
}

### step 5 ###
format_secs <- function(x) {
  x <- round(x)
  hh <- x %/% 3600
  mm <- (x %% 3600) %/% 60
  ss <- x %% 60
  sprintf("%02d:%02d:%02d", hh, mm, ss)
}

ncores = max(1, parallel::detectCores() - 1)

# outfile = "" makes worker messages appear in the console
cl = parallel::makeCluster(ncores, outfile = "")
doSNOW::registerDoSNOW(cl)

# load package on workers
parallel::clusterEvalQ(cl, library(multimode))

t1 = Sys.time()
n_tasks = nrow(scenarios)

cat("Simulation queue:\n")
print(data.frame(
  task = seq_len(n_tasks),
  dist = scenarios$dist,
  d = d[scenarios$d_index]
))

pb <- txtProgressBar(min = 0, max = n_tasks, style = 3)

progress <- function(n) {
  elapsed = as.numeric(difftime(Sys.time(), t1, units = "secs"))
  eta = if (n > 0) elapsed / n * (n_tasks - n) else NA_real_
  
  cat(sprintf(
    "\n[%s] completed %d/%d | elapsed: %s | rough ETA left: %s\n",
    format(Sys.time(), "%H:%M:%S"),
    n, n_tasks,
    format_secs(elapsed),
    if (is.na(eta)) "--:--:--" else format_secs(eta)
  ))
  
  setTxtProgressBar(pb, n)
}

opts <- list(
  progress = progress,
  preschedule = FALSE
)

res_list = foreach(
  s = seq_len(nrow(scenarios)),
  .packages = "multimode",
  .options.snow = opts,
  .inorder = FALSE
) %dopar% {
  
  dist_now = scenarios$dist[s]
  d_now = d[scenarios$d_index[s]]
  
  cat(sprintf("[%s] START dist=%s, d=%s\n",
              format(Sys.time(), "%H:%M:%S"),
              dist_now, d_now))
  flush.console()
  
  out_s = run_methods(x_list[[s]], methods)
  
  cat(sprintf("[%s] DONE  dist=%s, d=%s\n",
              format(Sys.time(), "%H:%M:%S"),
              dist_now, d_now))
  flush.console()
  
  list(
    dist = dist_now,
    d = d_now,
    out = out_s
  )
}

close(pb)
parallel::stopCluster(cl)

t2 = Sys.time()
t2 - t1

### step 6 ###
result_df = do.call(rbind, lapply(res_list, function(z) {
  data.frame(
    dist = z$dist,
    d = z$d,
    method = methods,
    power = z$out
  )
}))
out1 = result_df

saveRDS(out1, file = "out_norm_800.rds")

out1 = readRDS("out_norm_800.rds")

######parallel#######
######norm_mix2######
#####################

### step 1 ###
library(foreach)
library(doParallel)

n = 200
reps = 200
#d = seq(0, 5, 1)
d = c(0,1,seq(2,4.6,0.1))
#d = d[1:18]
#d=c(0,1,seq(2,3.9,0.1))
dist_names = c("mixnorm2")   # later extend this
methods = c("SI", "HY", "FM", "HH", "CH")

scenarios = expand.grid(
  dist = dist_names,
  d_index = seq_along(d),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

### step 2 ###
generate_x = function(dist_name, di, n, reps) {
  if (dist_name == "mixnorm2") {
    set.seed(1234)
    u = runif(n * reps)
    x = rnorm(n * reps, 0, 1)
    x[u < 0.4] = rnorm(sum(u < 0.4), di, 1)
    x[u > 0.8] = rnorm(sum(u > 0.8), 0, 3)
    x = matrix(x, reps, n)
    return(x)
  }
  
  stop("Unknown distribution")
}

### step 3 ###
x_list = vector("list", nrow(scenarios))

for (s in seq_len(nrow(scenarios))) {
  dist_name = scenarios$dist[s]
  di = d[scenarios$d_index[s]]
  
  x_list[[s]] = generate_x(dist_name, di, n, reps)
}

### step 4 ###
run_methods = function(x, methods) {
  reps = nrow(x)
  pval = matrix(0, nrow = length(methods), ncol = reps)
  
  hy_bounds = t(apply(x, 1, quantile, probs = c(0.01, 0.99)))
  
  for (j in seq_len(reps)) {
    xj = x[j, ]
    
    pval[1, j] = modetest(xj, method = "SI")$p.value
    pval[2, j] = modetest(
      xj,
      method = "HY",
      lowsup = hy_bounds[j, 1],
      uppsup = hy_bounds[j, 2]
    )$p.value
    pval[3, j] = modetest(xj, method = "FM")$p.value
    pval[4, j] = modetest(xj, method = "HH")$p.value
    pval[5, j] = modetest(xj, method = "CH")$p.value
  }
  
  rowMeans(pval < 0.05)
}

### step 5 ###
format_secs <- function(x) {
  x <- round(x)
  hh <- x %/% 3600
  mm <- (x %% 3600) %/% 60
  ss <- x %% 60
  sprintf("%02d:%02d:%02d", hh, mm, ss)
}

ncores = max(1, parallel::detectCores() - 1)

# outfile = "" makes worker messages appear in the console
cl = parallel::makeCluster(ncores, outfile = "")
doSNOW::registerDoSNOW(cl)

# load package on workers
parallel::clusterEvalQ(cl, library(multimode))

t1 = Sys.time()
n_tasks = nrow(scenarios)

cat("Simulation queue:\n")
print(data.frame(
  task = seq_len(n_tasks),
  dist = scenarios$dist,
  d = d[scenarios$d_index]
))

pb <- txtProgressBar(min = 0, max = n_tasks, style = 3)

progress <- function(n) {
  elapsed = as.numeric(difftime(Sys.time(), t1, units = "secs"))
  eta = if (n > 0) elapsed / n * (n_tasks - n) else NA_real_
  
  cat(sprintf(
    "\n[%s] completed %d/%d | elapsed: %s | rough ETA left: %s\n",
    format(Sys.time(), "%H:%M:%S"),
    n, n_tasks,
    format_secs(elapsed),
    if (is.na(eta)) "--:--:--" else format_secs(eta)
  ))
  
  setTxtProgressBar(pb, n)
}

opts <- list(
  progress = progress,
  preschedule = FALSE
)

res_list = foreach(
  s = seq_len(nrow(scenarios)),
  .packages = "multimode",
  .options.snow = opts,
  .inorder = FALSE
) %dopar% {
  
  dist_now = scenarios$dist[s]
  d_now = d[scenarios$d_index[s]]
  
  cat(sprintf("[%s] START dist=%s, d=%s\n",
              format(Sys.time(), "%H:%M:%S"),
              dist_now, d_now))
  flush.console()
  
  out_s = run_methods(x_list[[s]], methods)
  
  cat(sprintf("[%s] DONE  dist=%s, d=%s\n",
              format(Sys.time(), "%H:%M:%S"),
              dist_now, d_now))
  flush.console()
  
  list(
    dist = dist_now,
    d = d_now,
    out = out_s
  )
}

close(pb)
parallel::stopCluster(cl)

t2 = Sys.time()
t2 - t1

### step 6 ###
result_df = do.call(rbind, lapply(res_list, function(z) {
  data.frame(
    dist = z$dist,
    d = z$d,
    method = methods,
    power = z$out
  )
}))
out2 = result_df
# ### save data
# saveRDS(
#   list(
#     scenarios = scenarios,
#     d = d,
#     x_list = x_list
#   ),
#   file = "sim_norm_200.rds"
# )
# 
# sim_data = readRDS("sim_data.rds")
# x_list = sim_data$x_list

saveRDS(out2, file = "out2_norm_800.rds")

##################
######chi_mix######
#####################

### step 1 ###
n = 800
reps = 500
d = 10:30
d = 10:24
dist_names = c("mixchi")   # later extend this
methods = c("SI", "HY", "FM", "HH", "CH")

scenarios = expand.grid(
  dist = dist_names,
  d_index = seq_along(d),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

### step 2 ###
generate_x = function(dist_name, di, n, reps) {
  if (dist_name == "mixchi") {
    set.seed(1234)
    u = runif(n * reps)
    x = rchisq(n * reps, 5)
    x[u < 0.5] = rchisq(sum(u < 0.5), di)
    x = matrix(x, reps, n)
    return(x)
  }
  
  stop("Unknown distribution")
}

### step 3 ###
x_list = vector("list", nrow(scenarios))

for (s in seq_len(nrow(scenarios))) {
  dist_name = scenarios$dist[s]
  di = d[scenarios$d_index[s]]
  
  x_list[[s]] = generate_x(dist_name, di, n, reps)
}

### step 4 ###
run_methods = function(x, methods) {
  reps = nrow(x)
  pval = matrix(0, nrow = length(methods), ncol = reps)
  
  hy_bounds = t(apply(x, 1, quantile, probs = c(0.01, 0.99)))
  
  for (j in seq_len(reps)) {
    xj = x[j, ]
    
    pval[1, j] = modetest(xj, method = "SI")$p.value
    pval[2, j] = modetest(
      xj,
      method = "HY",
      lowsup = hy_bounds[j, 1],
      uppsup = hy_bounds[j, 2]
    )$p.value
    pval[3, j] = modetest(xj, method = "FM")$p.value
    pval[4, j] = modetest(xj, method = "HH")$p.value
    pval[5, j] = modetest(xj, method = "CH")$p.value
  }
  
  rowMeans(pval < 0.05)
}

### step 5 ###
format_secs <- function(x) {
  x <- round(x)
  hh <- x %/% 3600
  mm <- (x %% 3600) %/% 60
  ss <- x %% 60
  sprintf("%02d:%02d:%02d", hh, mm, ss)
}

ncores = max(1, parallel::detectCores() - 1)

# outfile = "" makes worker messages appear in the console
cl = parallel::makeCluster(ncores, outfile = "")
doSNOW::registerDoSNOW(cl)

# load package on workers
parallel::clusterEvalQ(cl, library(multimode))

t1 = Sys.time()
n_tasks = nrow(scenarios)

cat("Simulation queue:\n")
print(data.frame(
  task = seq_len(n_tasks),
  dist = scenarios$dist,
  d = d[scenarios$d_index]
))

pb <- txtProgressBar(min = 0, max = n_tasks, style = 3)

progress <- function(n) {
  elapsed = as.numeric(difftime(Sys.time(), t1, units = "secs"))
  eta = if (n > 0) elapsed / n * (n_tasks - n) else NA_real_
  
  cat(sprintf(
    "\n[%s] completed %d/%d | elapsed: %s | rough ETA left: %s\n",
    format(Sys.time(), "%H:%M:%S"),
    n, n_tasks,
    format_secs(elapsed),
    if (is.na(eta)) "--:--:--" else format_secs(eta)
  ))
  
  setTxtProgressBar(pb, n)
}

opts <- list(
  progress = progress,
  preschedule = FALSE
)

res_list = foreach(
  s = seq_len(nrow(scenarios)),
  .packages = "multimode",
  .options.snow = opts,
  .inorder = FALSE
) %dopar% {
  
  dist_now = scenarios$dist[s]
  d_now = d[scenarios$d_index[s]]
  
  cat(sprintf("[%s] START dist=%s, d=%s\n",
              format(Sys.time(), "%H:%M:%S"),
              dist_now, d_now))
  flush.console()
  
  out_s = run_methods(x_list[[s]], methods)
  
  cat(sprintf("[%s] DONE  dist=%s, d=%s\n",
              format(Sys.time(), "%H:%M:%S"),
              dist_now, d_now))
  flush.console()
  
  list(
    dist = dist_now,
    d = d_now,
    out = out_s
  )
}

close(pb)
parallel::stopCluster(cl)

t2 = Sys.time()
t2 - t1

### step 6 ###
result_df = do.call(rbind, lapply(res_list, function(z) {
  data.frame(
    dist = z$dist,
    d = z$d,
    method = methods,
    power = z$out
  )
}))
out3 = result_df

saveRDS(out3, file = "out3_chisq_800.rds")

 


############claw 200 ####
registerDoParallel(9)
t1 = Sys.time()
n = 200
reps = 1000
set.seed(1234)
x = rberdev(n*reps, dnum=23)
x = matrix(x, reps, n)
pval = matrix(0,8,reps)

registerDoParallel(9)
t1 = Sys.time()
for(j in 301:500){
  xj = x[j, ]
  pval[1, j] = modetest(xj, method = "SI")$p.value
  pval[2, j] = modetest(
    xj,
    method = "HY",
    lowsup = quantile(xj,0.01),
    uppsup = quantile(xj,0.99)
  )$p.value
  pval[3, j] = modetest(xj, method = "FM")$p.value
  pval[4, j] = modetest(xj, method = "HH")$p.value
  pval[5, j] = modetest(xj, method = "CH")$p.value
  pval[6, j] = modetest(xj, method = "ACR")$p.value
  pval[7, j] = bmodetest(xj,B=1000,eps1=5,eps2=2)$pvalue
  pval[8, j] = dip.test(xj)$p.value
  print(j)
}
claw_200 = rowMeans(pval[,1:500] < 0.05)
t2 = Sys.time()
t2-t1
stopImplicitCluster()

saveRDS(claw_200, file = "out_claw_200.rds")

############claw 800##############
registerDoParallel(9)
t1 = Sys.time()
n = 800
reps = 1000
set.seed(1234)
x = rberdev(n*reps, dnum=23)
x = matrix(x, reps, n)
pval = matrix(0,8,reps)

registerDoParallel(9)
t1 = Sys.time()
for(j in 143:300){
  xj = x[j, ]
  pval[1, j] = modetest(xj, method = "SI")$p.value
  pval[2, j] = modetest(
    xj,
    method = "HY",
    lowsup = quantile(xj,0.01),
    uppsup = quantile(xj,0.99)
  )$p.value
  pval[3, j] = modetest(xj, method = "FM")$p.value
  pval[4, j] = modetest(xj, method = "HH")$p.value
  pval[5, j] = modetest(xj, method = "CH")$p.value
  pval[6, j] = modetest(xj, method = "ACR")$p.value
  pval[7, j] = bmodetest(xj,B=1000,eps1=5,eps2=2)$pvalue
  pval[8, j] = dip.test(xj)$p.value
  print(j)
}
claw_800 = rowMeans(pval[,1:300] < 0.05)
t2 = Sys.time()
t2-t1
stopImplicitCluster()

saveRDS(pval, file = "out_claw_800.rds")

############claw 800##############smaller lambda
registerDoParallel(9)
t1 = Sys.time()
n = 800
reps = 1000
set.seed(1234)
x = rberdev(n*reps, dnum=23)
x = matrix(x, reps, n)
pval = matrix(0,8,reps)

registerDoParallel(9)
t1 = Sys.time()
for(j in 1:200){
  xj = x[j, ]
  pval[1, j] = modetest(xj, method = "SI")$p.value
  pval[2, j] = modetest(
    xj,
    method = "HY",
    lowsup = quantile(xj,0.01),
    uppsup = quantile(xj,0.99)
  )$p.value
  pval[3, j] = modetest(xj, method = "FM")$p.value
  pval[4, j] = modetest(xj, method = "HH")$p.value
  pval[5, j] = modetest(xj, method = "CH")$p.value
  pval[6, j] = modetest(xj, method = "ACR")$p.value
  pval[7, j] = bmodetest(xj,B=1000,eps1=5,eps2=2)$pvalue
  pval[8, j] = dip.test(xj)$p.value
  print(j)
}
claw_800 = rowMeans(pval[,1:88] < 0.05)
t2 = Sys.time()
t2-t1
stopImplicitCluster()

saveRDS(pval, file = "out2_claw_800.rds")

############## heavy tail 200 ###########
n = 200
reps = 1000
set.seed(1234)
x = rcauchy(n*reps)
x = matrix(x, reps, n)
pval = matrix(0,8,reps)

registerDoParallel(9)
t1 = Sys.time()
for(j in 101:200){
  xj = x[j, ]
  pval[1, j] = modetest(xj, method = "SI")$p.value
  pval[2, j] = modetest(
    xj,
    method = "HY",
    lowsup = quantile(xj,0.01),
    uppsup = quantile(xj,0.99)
  )$p.value
  pval[3, j] = modetest(xj, method = "FM")$p.value
  pval[4, j] = modetest(xj, method = "HH")$p.value
  pval[5, j] = modetest(xj, method = "CH")$p.value
  pval[6, j] = modetest(xj, method = "ACR")$p.value
  pval[7, j] = bmodetest(xj,B=1000,eps1=5,eps2=2)$pvalue
  pval[8, j] = dip.test(xj)$p.value
  print(j)
}
cauchy_200 = rowMeans(pval[,1:200] < 0.05)
t2 = Sys.time()
t2-t1
stopImplicitCluster()

saveRDS(cauchy_200, file = "out_cauchy_200.rds")

############## heavy tail 800 ###########
n = 800
reps = 1000
set.seed(1234)
x = rcauchy(n*reps)
x = matrix(x, reps, n)
pval = matrix(0,8,reps)

registerDoParallel(9)
t1 = Sys.time()
for(j in 101:200){
  xj = x[j, ]
  pval[1, j] = modetest(xj, method = "SI")$p.value
  pval[2, j] = modetest(
    xj,
    method = "HY",
    lowsup = quantile(xj,0.01),
    uppsup = quantile(xj,0.99)
  )$p.value
  pval[3, j] = modetest(xj, method = "FM")$p.value
  pval[4, j] = modetest(xj, method = "HH")$p.value
  pval[5, j] = modetest(xj, method = "CH")$p.value
  pval[6, j] = modetest(xj, method = "ACR")$p.value
  pval[7, j] = bmodetest(xj,B=1000,eps1=5,eps2=2)$pvalue
  pval[8, j] = dip.test(xj)$p.value
  print(j)
}
cauchy_800 = rowMeans(pval[,1:200] < 0.05)
t2 = Sys.time()
t2-t1
stopImplicitCluster()

saveRDS(cauchy_800, file = "out_cauchy_800.rds")


############## heavy tail 800 ###########smaller lambda
n = 800
reps = 1000
set.seed(1234)
x = rcauchy(n*reps)
x = matrix(x, reps, n)
pval = matrix(0,8,reps)

registerDoParallel(9)
t1 = Sys.time()
for(j in 348:500){
  xj = x[j, ]
  pval[1, j] = modetest(xj, method = "SI")$p.value
  pval[2, j] = modetest(
    xj,
    method = "HY",
    lowsup = quantile(xj,0.01),
    uppsup = quantile(xj,0.99)
  )$p.value
  pval[3, j] = modetest(xj, method = "FM")$p.value
  pval[4, j] = modetest(xj, method = "HH")$p.value
  pval[5, j] = modetest(xj, method = "CH")$p.value
  pval[6, j] = modetest(xj, method = "ACR")$p.value
  pval[7, j] = bmodetest(xj,B=1000,eps1=5,eps2=2,v=10)$pvalue
  pval[8, j] = dip.test(xj)$p.value
  print(j)
}
cauchy_800 = rowMeans(pval[,1:447] < 0.05)
t2 = Sys.time()
t2-t1
stopImplicitCluster()

saveRDS(cauchy_800, file = "out2_cauchy_800.rds")


###########test lambda on claw
registerDoParallel(9)
t1 = Sys.time()
n = 800
reps = 1000
set.seed(1234)
x = rberdev(n*reps, dnum=23)
x = matrix(x, reps, n)
pval = matrix(0,8,reps)

registerDoParallel(9)
t1 = Sys.time()
for(j in 1:100){
  xj = x[j, ]
  pval[1, j] = bmodetest(xj,B=1000,eps1=5,eps2=2,v=2)$pvalue
  pval[2, j] = bmodetest(xj,B=1000,eps1=5,eps2=2,v=5)$pvalue
  pval[3, j] = bmodetest(xj,B=1000,eps1=5,eps2=2,v=10)$pvalue
  pval[4, j] = bmodetest(xj,B=1000,eps1=5,eps2=2,v=50)$pvalue
  pval[5, j] = bmodetest(xj,B=1000,eps1=5,eps2=2,v=100)$pvalue
  pval[6, j] = bmodetest(xj,B=1000,eps1=5,eps2=2,v=500)$pvalue
  pval[7, j] = bmodetest(xj,B=1000,eps1=5,eps2=2,v=1000)$pvalue
  print(j)
}
claw_800 = rowMeans(pval[,1:100] < 0.05)
t2 = Sys.time()
t2-t1
stopImplicitCluster()
pval[1:7,1:20]
saveRDS(pval[1:7,1:20], file = "lam_claw_800.rds")







# Visualize
hist(cauchy_samples, breaks = 100, xlim = c(-20, 20), 
     main = "Cauchy Distribution", col = "skyblue")

# Install and load the actuar package
#library(actuar)

# Generate 1000 samples from a Pareto distribution
# shape = 3 (determines tail heaviness), scale = 1
pareto_samples <- rpareto(1000, shape = 3, scale = 1)

# Visualize
hist(pareto_samples, breaks = 50, xlim = c(0, 10), 
     main = "Pareto Distribution", col = "salmon")
# Compare Pareto (heavy) vs Normal (thin)
x <- sort(rpareto(1000, shape = 2, scale = 1))
y <- sort(rnorm(1000, mean = 0, sd = 1))

plot(log(x), log(1 - (1:1000)/1000), col="red", 
     main="Log-Log Survival Plot (Pareto vs Normal)", 
     xlab="Log(x)", ylab="Log(Survival Probability)")
points(log(y), log(1 - (1:1000)/1000), col="blue")
legend("topright", legend=c("Pareto", "Normal"), col=c("red", "blue"), pch=1)





#############sim1 norm mix
par(mfrow=c(1,2), mar=c(4,4,4,4))
## n = 200, normal mixture1
d=c(0,1,seq(2,4.8,0.1))
dip = c(0,0,0,0.002,0.002,0.004,0.006,0.007,0.01,0.014,0.032,0.041,0.064,0.093,0.128,
        0.169,0.237,0.316,0.407,0.504,0.586,0.673,0.747,0.821,0.882,0.917,0.944,0.967,0.98,0.989,0.994)
ker = c(0.047,0.042,0.046,0.042,0.048,0.047,0.056,0.068,0.088,0.12,0.143,0.188,0.232,0.301,0.374,
        0.472,0.564,0.643,0.712,0.789,0.86,0.916,0.938,0.965,0.984,0.99,0.998,1,1,1,1)
spl = c(0.005,0.002,0.014,0.014,0.016,0.028,0.039,0.051,0.079,0.117,0.166,0.229,0.324,0.423,0.526,
        0.635,0.742,0.824,0.892,0.941,0.968,0.979,0.99,0.998,1,1,1,1,1,1,1)
plot(d,spl,type="l",ylab='proportion rejection',main = "0.6N(0,1)+0.4N(d,1), n=200",lwd=2,lty=1,col='darkred')
points(d,spl,pch=15,cex=0.6,col='darkred')
lines(d,ker,lwd=2.5,lty=4,col='steelblue')
points(d,ker,pch=4,cex=0.6,col='steelblue')
lines(d,dip,lwd=2.5,lty=3,col='darkgreen')
points(d,dip,pch=3,cex=0.6,col='darkgreen')
abline(h=0.05,lty=9)
abline(v=2.442,lty=9)
#legend("topleft",c('splines','ACR','diptest'),cex=1,lwd=2,pch=c(15,4,3),lty=c(2,4,3),col=c('darkred','steelblue','darkgreen'))
out1=readRDS("out_norm_200.rds")
methods = c("SI", "HY", "FM", "HH", "CH")
method_cols = c("black", "orange", "purple", "brown", "magenta", "cyan4")
for (k in seq_along(methods)) {
  meth = methods[k]
  tmp = out1[out1$method == meth, ]
  tmp = tmp%>%arrange(d)
  lines(tmp$d, tmp$power, lwd = 2, lty = 5, col = method_cols[k])
  points(tmp$d, tmp$power, pch = 5, cex = 0.6, col = method_cols[k])
}

legend("topleft",
       legend = c("splines", "ACR", "diptest", methods),
       lwd = 2,
       pch = c(15, 4, 3, rep(5, 6)),
       lty = c(1, 4, 3, rep(5, 6)),
       col = c("darkred", "steelblue", "darkgreen", method_cols))



d=c(0,1,seq(2,4,0.1))
dip=c(0,0,0,0,0,0,0,0.005,
      0.009,0.02,0.033,0.065,0.145,0.256,0.402,0.572,0.743,0.875,0.955,0.984,0.998,0.999,1)
ker=c(0.046,0.037,0.04,0.042,0.038,0.045,0.056,0.078,
      0.107,0.161,0.249,0.387,0.544,0.678,0.829,0.928,0.974,0.994,0.998,0.999,1,1,1)
spl=c(0.003,0.003,0.004,0.004,0.004,0.019,0.042,0.085,0.169,0.285,0.451,0.643,0.8,0.916,0.969,0.992,
      0.999,1,1,1,1,1,1)
plot(d,spl,type="l",ylab='proportion rejection',main = "0.6N(0,1)+0.4N(d,1), n=800",lwd=2.5,lty=1,col='darkred')
points(d,spl,pch=15,cex=0.6,col='darkred')
lines(d,ker,lwd=2.5,lty=4,col='steelblue')
points(d,ker,pch=4,cex=0.6,col='steelblue')
lines(d,dip,lwd=2.5,lty=3,col='darkgreen')
points(d,dip,pch=3,cex=0.6,col='darkgreen')
abline(h=0.05,lty=9)
abline(v=2.442,lty=9)
out1=readRDS("out_norm_800.rds")
methods = c("SI", "HY", "FM", "HH", "CH")
method_cols = c("black", "orange", "purple", "brown", "magenta", "cyan4")
for (k in seq_along(methods)) {
  meth = methods[k]
  tmp = out1[out1$method == meth, ]
  tmp = tmp%>%arrange(d)
  lines(tmp$d, tmp$power, lwd = 2, lty = 5, col = method_cols[k])
  points(tmp$d, tmp$power, pch = 5, cex = 0.6, col = method_cols[k])
}

legend("topleft",
       legend = c("splines", "ACR", "diptest", methods),
       lwd = 2,
       pch = c(15, 4, 3, rep(5, 6)),
       lty = c(1, 4, 3, rep(5, 6)),
       col = c("darkred", "steelblue", "darkgreen", method_cols))







#######################sim2 nor_mix2
par(mfrow=c(1,2), mar=c(4,4,4,4))
## n = 200, normal mixture2
d=c(0,1,seq(2,4.6,0.1))
spl=c(0.011,0.014,0.02,0.028,0.033,0.049,0.058,0.085,0.111,0.157,0.209,0.274,0.347,0.43,0.532,0.606,
      0.704,0.793,0.856,0.893,0.931,0.952,0.967,0.982,0.989,0.995,0.998,0.999,1)
dip=c(0,0,0,0,0,0.002,0.01,0.017,0.019,0.029,0.045,0.066,0.096,0.136,0.181,0.239,0.31,0.378,
      0.458,0.547,0.646,0.719,0.776,0.836,0.875,0.914,0.938,0.964,0.976)
ker=c(0.044,0.047,0.079,0.082,0.084,0.106,0.125,0.145,0.17,0.221,
      0.27,0.317,0.385,0.46,0.511,0.588,0.659,0.724,0.801,0.849,0.889,0.92,0.944,0.961,0.977,0.988,0.996,0.997,0.999)
plot(d,spl,type="l",ylab='proportion rejection',main = "0.4N(0,1)+0.4N(d,1)+0.2N(0,9), n=200",lwd=2.5,lty=1,ylim=c(0,1),col='darkred')
points(d,spl,pch=15,cex=0.6,col='darkred')
lines(d,ker,lwd=2.5,lty=4,col='steelblue')
points(d,ker,pch=4,cex=0.6,col='steelblue')
lines(d,dip,lwd=2.5,lty=3,col='darkgreen')
points(d,dip,pch=3,cex=0.6,col='darkgreen')
abline(h=0.05,lty=9)
abline(v=2,lty=9)
out1=readRDS("out2_norm_200.rds")
methods = c("SI", "HY", "FM", "HH", "CH")
method_cols = c("black", "orange", "purple", "brown", "magenta", "cyan4")
for (k in seq_along(methods)) {
  meth = methods[k]
  tmp = out1[out1$method == meth, ]
  tmp = tmp%>%arrange(d)
  lines(tmp$d, tmp$power, lwd = 2, lty = 5, col = method_cols[k])
  points(tmp$d, tmp$power, pch = 5, cex = 0.6, col = method_cols[k])
}

legend("topleft",
       legend = c("splines", "ACR", "diptest", methods),
       lwd = 2,
       pch = c(15, 4, 3, rep(5, 6)),
       lty = c(1, 4, 3, rep(5, 6)),
       col = c("darkred", "steelblue", "darkgreen", method_cols))




d=c(0,seq(1,1.8,0.2),seq(1.9,3.9,0.1))
d=c(0,1,seq(2,3.9,0.1))
spl=c(0.009,0.007,0.021,0.035,0.054,0.105,0.149,0.27,0.4,0.538,0.695,0.827,0.909,0.97,0.99,
      0.997,0.999,1,1,1,1,1)
dip=c(0.001,0,0.001,0.001,0.002,0.005,0.009,0.019,0.041,
      0.066,0.137,0.228,0.365,0.529,0.698,0.823,0.915,0.966,0.989,0.998,0.999,1)
ker=c(0.05,0.049,0.103,0.121,0.142,0.175,0.238,0.284,0.355,
      0.478,0.583,0.714,0.816,0.907,0.962,0.981,0.994,0.998,0.999,1,1,1)
plot(d,spl,type="l",ylab='proportion rejection',main = "0.4N(0,1)+0.4N(d,1)+0.2N(0,9), n=800",lwd=2.5,lty=1,col='darkred')
points(d,spl,pch=15,cex=0.6,col='darkred')
lines(d,ker,lwd=2.5,lty=4,col='steelblue')
points(d,ker,pch=4,cex=0.6,col='steelblue')
lines(d,dip,lwd=2.5,lty=3,col='darkgreen')
points(d,dip,pch=3,cex=0.6,col='darkgreen')
abline(h=0.05,lty=9)
abline(v=2,lty=9)
legend("topleft",c('splines','kernel','diptest'),cex=1,lwd=2,pch=c(15,4,3),lty=c(2,4,3),col=c('darkred','steelblue','darkgreen'))


out1=readRDS("out2_norm_800.rds")
methods = c("SI", "HY", "FM", "HH", "CH")
method_cols = c("black", "orange", "purple", "brown", "magenta", "cyan4")
for (k in seq_along(methods)) {
  meth = methods[k]
  tmp = out1[out1$method == meth, ]
  tmp = tmp%>%arrange(d)
  lines(tmp$d, tmp$power, lwd = 2, lty = 5, col = method_cols[k])
  points(tmp$d, tmp$power, pch = 5, cex = 0.6, col = method_cols[k])
}

legend("topleft",
       legend = c("splines", "ACR", "diptest", methods),
       lwd = 2,
       pch = c(15, 4, 3, rep(5, 6)),
       lty = c(1, 4, 3, rep(5, 6)),
       col = c("darkred", "steelblue", "darkgreen", method_cols))


##########################sim3 ch



par(mfrow=c(1,2), mar=c(4,4,4,4))
d=c(8:30)
d=10:30
spl=c(0.015,0.016,0.017,0.028,0.036,00.057,0.085,0.171,0.25,0.378,0.536,0.656,0.762,
      0.867,0.924,0.968,0.984,0.992,0.997,0.999,1)
ker=c(0.066,0.069,0.067,0.089,0.084,0.088,0.116,0.161,0.204,0.294,
      0.422,0.52,0.642,0.782,0.866,0.922,0.963,0.977,0.991,0.997,0.998)
dip=c(0.002,0.003,0.002,0.006,0.003,0.004,0.017,0.029,0.061,0.088,
      0.125,00.224,0.353,0.469,0.641,0.73,0.831,0.906,0.943,0.967,0.988)
plot(d,spl,type="l",ylab='proportion rejection',main = "0.5chisq(5)+0.5chisq(d), n=200",lwd=2,lty=1,col='darkred')
points(d,spl,pch=15,cex=0.6,col='darkred')
lines(d,ker,lwd=2.5,lty=4,col='steelblue')
points(d,ker,pch=4,cex=0.6,col='steelblue')
lines(d,dip,lwd=2.5,lty=3,col='darkgreen')
points(d,dip,pch=3,cex=0.6,col='darkgreen')
abline(h=0.05,lty=9)
abline(v=15.25,lty=9)
#legend("topleft",c('splines','kernel','diptest'),cex=1,lwd=2,pch=c(15,4,3),lty=c(2,4,3),col=c('darkred','steelblue','darkgreen'))
out1=readRDS("out3_chisq_200.rds")
methods = c("SI", "HY", "FM", "HH", "CH")
method_cols = c("black", "orange", "purple", "brown", "magenta", "cyan4")
for (k in seq_along(methods)) {
  meth = methods[k]
  tmp = out1[out1$method == meth, ]
  tmp = tmp%>%arrange(d)
  lines(tmp$d, tmp$power, lwd = 2, lty = 5, col = method_cols[k])
  points(tmp$d, tmp$power, pch = 5, cex = 0.6, col = method_cols[k])
}

legend("topleft",
       legend = c("splines", "ACR", "diptest", methods),
       lwd = 2,
       pch = c(15, 4, 3, rep(5, 6)),
       lty = c(1, 4, 3, rep(5, 6)),
       col = c("darkred", "steelblue", "darkgreen", method_cols))





d=c(10:24)
spl=c(0.009,0.009,0.008,0.014,0.038,0.058,0.141,0.323,0.604,0.84,0.967,0.994,0.999,1,1)
dip=c(0,0,0,0.002,0.002,0.003,0.006,0.021,0.08,0.242,0.466,0.738,0.914,0.974,0.998)
ker=c(0.059,0.076,0.063,0.056,0.07,0.077,0.104,
      0.22,0.427,0.667,0.856,0.955,0.992,0.998,1)
plot(d,spl,type="l",ylab='proportion rejection',main = "0.5chisq(5)+0.5chisq(d), n=800",lwd=2.5,lty=1,col='darkred')
points(d,spl,pch=15,cex=0.6,col='darkred')
lines(d,ker,lwd=2.5,lty=4,col='steelblue')
points(d,ker,pch=4,cex=0.6,col='steelblue')
lines(d,dip,lwd=2.5,lty=3,col='darkgreen')
points(d,dip,pch=3,cex=0.6,col='darkgreen')
abline(h=0.05,lty=9)
abline(v=15.25,lty=9)
#legend("topleft",c('splines','kernel','diptest'),cex=1,lwd=2,pch=c(15,4,3),lty=c(2,4,3),col=c('darkred','steelblue','darkgreen'))
out1=readRDS("out3_chisq_800.rds")
methods = c("SI", "HY", "FM", "HH", "CH")
method_cols = c("black", "orange", "purple", "brown", "magenta", "cyan4")
for (k in seq_along(methods)) {
  meth = methods[k]
  tmp = out1[out1$method == meth, ]
  tmp = tmp%>%arrange(d)
  lines(tmp$d, tmp$power, lwd = 2, lty = 5, col = method_cols[k])
  points(tmp$d, tmp$power, pch = 5, cex = 0.6, col = method_cols[k])
}

legend("topleft",
       legend = c("splines", "ACR", "diptest", methods),
       lwd = 2,
       pch = c(15, 4, 3, rep(5, 6)),
       lty = c(1, 4, 3, rep(5, 6)),
       col = c("darkred", "steelblue", "darkgreen", method_cols))




claw_200 = readRDS("out_claw_200.rds")
claw_800 = readRDS("out_claw_800.rds")
claw_800 = rowMeans(claw_800[,1:300] < 0.05)
cauchy_200 = readRDS("out_cauchy_200.rds")
cauchy_800 = readRDS("out_cauchy_800.rds")
method = c("SI","HY","FM","HH","CH","ACR","SPLINES","DIP")
power = rbind(claw_200,claw_800,cauchy_200,cauchy_800)
colnames(power) = method



lam_claw_800 = readRDS("lam_claw_800.rds")
rownames(lam_claw_800) = c('lamda/2','lambda/5','lambda/10','lambda/50','lambda/100','lambda/500','lambda/1000')
rowMeans(lam_claw_800<0.05)




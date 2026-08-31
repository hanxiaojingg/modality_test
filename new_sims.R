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


###########test lambda on claw

n = 200
reps = 1000
set.seed(1234)
x = rberdev(n*reps, dnum=23)
x = matrix(x, reps, n)
pval = matrix(0,4,reps)

registerDoParallel(9)
t1 = Sys.time()
for(j in 1:100){
  xj = x[j, ]
  pval[1, j] = bmodetest(xj,B=500,lam=1)$pvalue
  pval[2, j] = bmodetest(xj,B=500,lam=0.1)$pvalue
  pval[3, j] = bmodetest(xj,B=500,lam=0.01)$pvalue
  pval[4, j] = bmodetest(xj,B=500,lam=0.001)$pvalue
  print(j)
}
rowMeans(pval[,1:100] < 0.05)
t2 = Sys.time()
t2-t1
stopImplicitCluster()

saveRDS(pval, file = "lam_claw_200.rds")


#################################################
############## cross validation #################
#################################################
bmodetest_cv <- function(y,lower = NULL, upper = NULL,B=1000,eps1=.5,eps2=.2,eps=0.01,cv=TRUE){
  nraw=length(y)
  q1=quantile(y,.01)
  q2=quantile(y,.99)
  rng=q2-q1
  if(is.null(lower)){s1=as.numeric(q1-.4*rng)}
  if(is.null(upper)){s2=as.numeric(q2+.4*rng)}
  s1=as.numeric(q1-.4*rng)
  s2=as.numeric(q2+.4*rng)
  capk=round(nraw^(1/7)*12)
  qy=as.numeric(quantile(y,1:capk/(capk+1)))
  ####add more knots on both sides
  if(s1<qy[1]){
    k1=min(floor(log((qy[1]-s1)/3/(qy[2]-qy[1])+1, 3/2) -1 ), round(capk/8))
    if(k1<2){
      qy = c(qy,s1)
      if((qy[1]-s1)/(qy[2]-qy[1])>=2) qy = c(qy,qy[1]-(qy[1]-s1)/2)
    }else{
      q1 = ((qy[1]-s1)/(qy[2]-qy[1]))^(1/(k1+1))
      qy = c(qy, qy[1] - (qy[2]-qy[1])*q1^(1:(k1)), s1)
    }
  }
  if(s2>qy[capk]){
    k2=min(floor(log((s2-qy[capk])/3/(qy[capk]-qy[capk-1])+1, 3/2) - 1), round(capk/8))
    if(k2<2){
      qy=c(qy,s2)
      if((s2-qy[capk])/(qy[capk]-qy[capk-1])>=2) qy = c(qy,(s2-qy[capk])/2+qy[capk])
    }else{
      q2 = ((s2-qy[capk])/(qy[capk]-qy[capk-1]))^(1/(k2+1))
      qy = c(qy, (qy[capk]-qy[capk-1])*q2^(1:(k2))+qy[capk], s2)
    }
  }
  qy = sort(qy)
  ####add more knots between the two modes
  d=min(diff(qy))
  dd=floor(diff(qy)/min(diff(qy)))
  m1=min(which(dd==1))
  m2=max(which(dd==1))
  if(sum(dd[m1:m2]>2)>0){ for(i in which(dd[m1:m2]>2)){qy=c(qy,(qy[m1+i-1]+qy[m1+i])/2)}
  }
  kn=sort(qy)
  
  s1=min(kn)
  s2=max(kn)
  yraw = y
  y = y[y>=s1 & y<=s2]
  n = length(y)
  if(!cv & is.null(lam)){
    K = kurtosis(y)
    if(K<2){
      lam = 10^2*n^(-1/7)
    } else if(K>2 & K<5){
      lam = 10^(4-K)*n^(-1/7)
    } else if(K>5 & K<9){
      lam = 10^(3/2-K/2)*n^(-1/7)
    } else{
      lam = 10^(-3)*n^(-1/7)
    }
  }
  if(is.null(eps)){
    eps=ifelse(abs(skewness(y))>0.7,eps2,eps1)
  }else{eps1=eps2=eps}
  m=length(kn)+1
  bspl=bSpline(y,degree=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  yp=0:4000/4000*(s2-s1)+s1
  s1=min(s1,min(yp))
  s2=max(s2,max(yp))
  bp=bSpline(yp,degree=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  slopes=bSpline(kn,degree=2,derivs=1,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  D2=bSpline(kn,degree=2,derivs=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  d=(s2-s1)/m
  ###  make all basis functions integrate to one
  dp=yp[2]-yp[1]
  avec=1:m
  for(i in 1:m){avec[i]=sum(bp[,i])*dp}
  for(i in 1:m){
    bspl[,i]=bspl[,i]/avec[i]
    bp[,i]=bp[,i]/avec[i]
    slopes[,i]=slopes[,i]/avec[i]
    D2[,i]=D2[,i]/avec[i]
  }
  av1=avec
  avec=rep(1,m)
  hmat=matrix(nrow=m,ncol=m)
  cvec=1:m
  for(i in 1:m){
    for(j in i:m){
      pr=bp[,i]*bp[,j]
      hmat[i,j]=sum(pr)*dp
      hmat[j,i]=hmat[i,j]
    }
    cvec[i]=sum(bspl[,i])/n
  }
  b0=rep(1/m,m)
  
  wmat=matrix(0,nrow=m,ncol=m-1)
  for(i in 1:(m-1)){wmat[i,i]=-1;wmat[i+1,i]=1}
  
  D1=matrix(0,m-2,m-1)
  for(i in 1:(m-2)){D1[i,i]=-1;D1[i,i+1]=1}
  D=D1%*%D2*d^(5/2)
  ##  get unimodal
  ans2=bmfit(y,kn,hmat,cvec,slopes,b0,wmat,D,bspl,lam,cv=cv)
  ans1=umfit(y,kn,hmat,cvec,slopes,b0,wmat,D,bspl,lam = ans2$lam,eps=eps, cv=FALSE)
  t1=as.numeric((ans1$crit-ans2$crit))
  #t2=as.numeric((ans1$crit-ans2$crit)/abs(ans1$crit))
  fhat2=round(bp%*%ans2$bhat,10)
  dfhat2=diff(fhat2)
  outtb=NULL
  if(sum(dfhat2>0)==0 |sum(dfhat2<0)==0){
    pvalue=2
  }else{
    md=max(which(dfhat2>0))
    if (!is.unsorted(fhat2[1:(md-1)]) ) {
      pvalue=2
    } else {
      cdf1=bp%*%ans1$bhat
      for(i in 2:4001){
        cdf1[i]=cdf1[i-1]+cdf1[i]
      }
      cdf1=cdf1-min(cdf1)
      cdf1=cdf1/cdf1[4001]
      outtb <- foreach(t=1:B,.combine = 'rbind') %dopar%{
        yb=sapply(1:n,function(o){u=runif(1);id=min(which(u<cdf1));alp=(cdf1[id]-u)/(cdf1[id]-cdf1[id-1]);alp*yp[id-1]+(1-alp)*yp[id]})
        modet(yb,kn,hmat,slopes,b0,wmat,D,bspl,av1,bp,lam=ans2$lam,eps=eps)
      } 
      pvalue=sum(outtb>t1)/B
    }
  }
  ans=new.env()
  ans$yp=yp
  ans$fhat1=bp%*%ans1$bhat
  ans$fhat2=bp%*%ans2$bhat
  ans$statistic=t1
  ans$tb=outtb
  ans$pvalue=pvalue
  ans$lam=c(ans1$lam,ans2$lam)
  ans$kn=kn
  ans$crit=c(ans1$crit,ans2$crit)
  ans$kurtosis=kurtosis(y)
  ans$skewness=skewness(y)
  ans$truncated=(nraw!=n)
  ans
}
##############
umfit=function(y,kn,hmat,cvec,slopes,b0,wmat,D,bspl,lam=NULL,eps,cv=cv){	
  n <- length(y)
  m=length(kn)+1
  k1=max(min(which(kn>quantile(y,.2)))-1,3)
  k2=min(max(which(kn<quantile(y,.8)))+1,m-4)         
  amatl1=list()
  for(k in k1:k2){
    amat=matrix(0,nrow=m+1,ncol=m)
    amat[1:k,]=slopes[1:k,]
    amat[(k+1):(m-1),]=-slopes[(k+1):(m-1),]
    amat[m,1]=1
    amat[m+1,m]=1
    epsvec=c(rep(0,k1-2),rep(eps/n^(2/7)/diff(range(kn))^2,k-k1+1),0,0,rep(eps/n^(2/7)/diff(range(kn))^2,k2-k+1),rep(0,m-k2-3),0,0)
    amatl1[[k-k1+1]]=list(amat=amat, epsvec=epsvec)
  }
  ## cross validation
  if(cv | is.null(lam)){
    zvec=t(wmat)%*%(cvec-hmat%*%b0-0.1*n^(-1/7)*t(D)%*%D%*%b0)
    qmat=t(wmat)%*%(hmat+0.1*n^(-1/7)*t(D)%*%D)%*%wmat 
    crit<-lapply(amatl1, function(x){ans <- quadprog::solve.QP(qmat,zvec,t(x[[1]]%*%wmat),x[[2]]-x[[1]]%*%b0);ans$value})
    amat1 <- amatl1[[which.min(crit)]][[1]]
    epsvec <- amatl1[[which.min(crit)]][[2]]
    ## number of folds
    L = 10
    lam_set = c(0.0001,0.001,0.01,0.1,1)
    fold_assignments <- sample(1:L, n, replace = TRUE)
    crit_cv = NULL
    for( lamt in lam_set){
      err = rep(0,L)
      for (l in 1:L) {
        train_data <- y[fold_assignments != l]
        val_data   <- y[fold_assignments == l]
        bspl_train = bspl[fold_assignments != l,]
        bspl_val =  bspl[fold_assignments == l,]
        cvec_train=cvec_val=1:m
        for(i in 1:m){
          cvec_train[i]=sum(bspl_train[,i])/n
          cvec_val[i]=sum(bspl_val[,i])/n
        }
        zvec=t(wmat)%*%(cvec_train-hmat%*%b0-lamt*n^(-1/7)*t(D)%*%D%*%b0)
        qmat=t(wmat)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%wmat 
        ans1=solve.QP(qmat,zvec,t(amat1%*%wmat),epsvec-amat1%*%b0)
        alphahat1=ans1$solution
        bhat1=wmat%*%alphahat1+b0
        err[l]=t(bhat1)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%bhat1-2*sum(cvec_val*bhat1)
      }
      crit_cv = c(crit_cv,mean(err))
    }
    lamt=lam_set[which(crit_cv==min(crit_cv))]
    zvec=t(wmat)%*%(cvec-hmat%*%b0-lamt*n^(-1/7)*t(D)%*%D%*%b0)
    qmat=t(wmat)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%wmat 
    # crit<-lapply(amatl1, function(x){ans <- quadprog::solve.QP(qmat,zvec,t(x[[1]]%*%wmat),x[[2]]-x[[1]]%*%b0);ans$value})
    # amat1 <- amatl1[[which.min(crit)]][[1]]
    # epsvec <- amatl1[[which.min(crit)]][[2]]
    ans1=solve.QP(qmat,zvec,t(amat1%*%wmat),epsvec-amat1%*%b0)
    alphahat1=ans1$solution
    bhat1=wmat%*%alphahat1+b0
    cr1=t(bhat1)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%bhat1-2*sum(cvec*bhat1)
  }else{
    lamt=lam
    zvec=t(wmat)%*%(cvec-hmat%*%b0-lamt*n^(-1/7)*t(D)%*%D%*%b0)
    qmat=t(wmat)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%wmat 
    crit<-lapply(amatl1, function(x){ans <- quadprog::solve.QP(qmat,zvec,t(x[[1]]%*%wmat),x[[2]]-x[[1]]%*%b0);ans$value})
    amat1 <- amatl1[[which.min(crit)]][[1]]
    epsvec <- amatl1[[which.min(crit)]][[2]]
    ans1=solve.QP(qmat,zvec,t(amat1%*%wmat),epsvec-amat1%*%b0)
    alphahat1=ans1$solution
    bhat1=wmat%*%alphahat1+b0
    cr1=t(bhat1)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%bhat1-2*sum(cvec*bhat1) 
  }
  
  ans1=new.env()
  ans1$bhat=bhat1
  ans1$lam=lamt
  ans1$crit=cr1
  ans1
}
############################################################################
bmfit=function(y,kn,hmat,cvec,slopes,b0,wmat,D,bspl,lam=NULL,cv=cv){
  n=length(y)
  m=dim(slopes)[2]
  m1=max(min(which(kn>quantile(y,.1)))-1,2)
  m2=min(max(which(kn<quantile(y,.9)))+1,m-2)
  if((m2-m1)<5){m1=max(2,m1-1);m2=min(m2+1,m-2)}
  trips=matrix(0,nrow=choose((m2-m1+1),3),ncol=3)
  nr=0   
  for(i in m1:(m2-4)){
    for(j in (i+2):(m2-2)){
      for(k in (j+2):(m2)){
        nr=nr+1
        trips[nr,]=c(i,j,k)
      }
    }
  }
  trips=trips[1:nr,]
  amatl2=list()
  for(i in 1:nr){
    amat=matrix(0,m+1,m)
    amat[1:(m-1),1:m]=slopes
    amat[(trips[i,1]+1):trips[i,2],]=-slopes[(trips[i,1]+1):trips[i,2],]
    amat[(trips[i,3]):(m-1),]=-slopes[(trips[i,3]):(m-1),]
    amat[m,1]=1
    amat[m+1,m]=1
    amatl2[[i]]=amat
  }
  if(cv | is.null(lam)){
    zvec=t(wmat)%*%(cvec-hmat%*%b0-0.1*n^(-1/7)*t(D)%*%D%*%b0)
    qmat=t(wmat)%*%(hmat+0.1*n^(-1/7)*t(D)%*%D)%*%wmat 
    crit<-lapply(amatl2, function(x){ans <- quadprog::solve.QP(qmat,zvec,t(x%*%wmat),-x%*%b0);ans$value})
    amat2 <- amatl2[[which.min(crit)]]
    L = 10
    lam_set = c(0.00001,0.0001,0.001,0.01,0.1,1)
    fold_assignments <- sample(1:L, n, replace = TRUE)
    crit_cv = NULL
    for( lamt in lam_set){
      err = rep(0,L)
      for (l in 1:L) {
        train_data <- y[fold_assignments != l]
        val_data   <- y[fold_assignments == l]
        bspl_train = bspl[fold_assignments != l,]
        bspl_val =  bspl[fold_assignments == l,]
        cvec_train=cvec_val=1:m
        for(i in 1:m){
          cvec_train[i]=sum(bspl_train[,i])/n
          cvec_val[i]=sum(bspl_val[,i])/n
        }
        zvec=t(wmat)%*%(cvec_train-hmat%*%b0-lamt*n^(-1/7)*t(D)%*%D%*%b0)
        qmat=t(wmat)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%wmat 
        ans2=solve.QP(qmat,zvec,t(amat2%*%wmat),-amat2%*%b0)
        alphahat2=ans2$solution
        bhat2=wmat%*%alphahat2+b0
        err[l]=t(bhat2)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%bhat2-2*sum(cvec_val*bhat2)
      }
      crit_cv = c(crit_cv,mean(err))
    }
    lamt=lam_set[which(crit_cv==min(crit_cv))]
    zvec=t(wmat)%*%(cvec-hmat%*%b0-lamt*n^(-1/7)*t(D)%*%D%*%b0)
    qmat=t(wmat)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%wmat 
    ans2=solve.QP(qmat,zvec,t(amat2%*%wmat),-amat2%*%b0)
    alphahat2=ans2$solution
    bhat2=wmat%*%alphahat2+b0
    cr2=t(bhat2)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%bhat2-2*sum(cvec*bhat2)
  }else{
    lamt=lam
    zvec=t(wmat)%*%(cvec-hmat%*%b0-lamt*n^(-1/7)*t(D)%*%D%*%b0)
    qmat=t(wmat)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%wmat 
    crit<-lapply(amatl2, function(x){ans <- quadprog::solve.QP(qmat,zvec,t(x%*%wmat),-x%*%b0);ans$value})
    amat2 <- amatl2[[which.min(crit)]]
    ans2=solve.QP(qmat,zvec,t(amat2%*%wmat),-amat2%*%b0)
    alphahat2=ans2$solution
    bhat2=wmat%*%alphahat2+b0
    cr2=t(bhat2)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%bhat2-2*sum(cvec*bhat2)
  }
  
  ans2=new.env()
  ans2$bhat=bhat2
  ans2$lam=lamt
  ans2$crit=cr2
  ans2
}
##################################
modet <- function(yb,kn,hmat,slopes,b0,wmat,D,bspl,av1,bp,lam,eps){
  n=length(yb)
  m=dim(slopes)[2]
  capk=length(kn) 
  s1=min(kn)
  s2=max(kn)
  bb=bSpline(yb,degree=2,knots=kn[2:(capk-1)],Boundary.knots=c(s1,s2),intercept=TRUE)
  m=dim(bb)[2]
  for(i in 1:m){bb[,i]=bb[,i]/av1[i]}
  cb=1:m
  for(i in 1:m){
    cb[i]=sum(bb[,i])/n
  }
  fit1=umfit(yb,kn,hmat,cb,slopes,b0,wmat,D,bspl,lam,eps=eps, cv=FALSE)
  fit2=bmfit(yb,kn,hmat,cb,slopes,b0,wmat,D,bspl,lam, cv=FALSE)
  fhat2=round(bp%*%fit2$bhat,10)
  dfhat2=diff(fhat2)
  if(sum(dfhat2>0)==0 |sum(dfhat2<0)==0){
    return (-1e+3)
  }else{md=max(which(dfhat2>0))
  if (!is.unsorted(fhat2[1:(md-1)]) ){
    return (-1e+3)
  }else{return(as.numeric(fit1$crit-fit2$crit))}
  }
}



# define a grid: lambda = 0.01, 0.1, 1
# partition data split the data into K roughtly equal-sized folds F1,...,FK
# outer loop: for each lambda
# inner loop: iterate over k=1,...,K
#             training set: Train = Data - Fk
#             find estimates b using train data and current lambda
#             Err_{lmabda,k} = Crit
# CV(lambda) = (Err1+...+ErrK)/K
# find the optimal lambda that has the smallest CV

y = rnorm(200,0,1)




n = 200
### generate sample from claw distribution
y = benchden::rberdev(n, dnum=23)/sd(y)
#hist(y,breaks=20)
registerDoParallel(8)
#ans=bmodetest(y, B=100) ### set B=1 to skip sampling and get the fits
ans_cv=bmodetest_cv(y,B=100)
stopImplicitCluster()
#par(mfrow=c(1,2))
hist(y,freq=FALSE,breaks=30)
lines(ans$yp,ans$fhat1,col=2)
lines(ans$yp,ans$fhat2,col=3)
hist(y,freq=FALSE,breaks=30)
lines(ans_cv$yp,ans_cv$fhat1,col=2)
lines(ans_cv$yp,ans_cv$fhat2,col=3)
ans$lam
ans_cv$lam
ans$crit
ans_cv$crit
ans$pvalue
ans_cv$pvalue




n = 200
ui = runif(n)
y = c(rnorm(sum(ui<0.4),0,1),rnorm(sum(ui>=0.4),4,1))
y = y/sd(y)
registerDoParallel(8)
ans=bmodetest(y, B=100) ### set B=1 to skip sampling and get the fits
ans_cv=bmodetest_cv(y,B=100)
stopImplicitCluster()
#par(mfrow=c(1,2))
hist(y,freq=FALSE,breaks=30)
lines(ans$yp,ans$fhat1,col=2)
lines(ans$yp,ans$fhat2,col=3)
hist(y,freq=FALSE,breaks=30)
lines(ans_cv$yp,ans_cv$fhat1,col=2)
lines(ans_cv$yp,ans_cv$fhat2,col=3)
ans$lam
ans_cv$lam
ans$crit
ans_cv$crit



n = 200
ui = runif(n)
y = c(rnorm(sum(ui<0.4),0,1),rnorm(sum(ui>=0.4),3,1))
y = y/sd(y)
#hist(y,breaks=20)
ans=bmodetest(y, B=1) ### set B=1 to skip sampling and get the fits
ans_cv=bmodetest_cv(y,B=1)
#par(mfrow=c(1,2))
hist(y,freq=FALSE,breaks=30)
lines(ans$yp,ans$fhat1,col=2)
lines(ans$yp,ans$fhat2,col=3)
hist(y,freq=FALSE,breaks=30)
lines(ans_cv$yp,ans_cv$fhat1,col=2)
lines(ans_cv$yp,ans_cv$fhat2,col=3)
ans$lam
ans_cv$lam
ans$crit
ans_cv$crit



n = 200
ui = runif(n)
y = c(rnorm(sum(ui<0.4),0,1),rnorm(sum(ui>=0.4),2,1))
y = y/sd(y)
#hist(y,breaks=20)
ans=bmodetest(y, B=1) ### set B=1 to skip sampling and get the fits
ans_cv=bmodetest_cv(y,B=1)
#par(mfrow=c(1,2))
hist(y,freq=FALSE,breaks=30)
lines(ans$yp,ans$fhat1,col=2)
lines(ans$yp,ans$fhat2,col=3)
hist(y,freq=FALSE,breaks=30)
lines(ans_cv$yp,ans_cv$fhat1,col=2)
lines(ans_cv$yp,ans_cv$fhat2,col=3)
ans$lam
ans_cv$lam
ans$crit
ans_cv$crit

par(mfrow=c(1,2))



set.seed(1234)
pvals = rep(100,100)
lambdas = rep(0,100)
n=200
registerDoParallel(9)
for(reps in 1:pvals){
  ### generate sample from claw distribution
  y = benchden::rberdev(n, dnum=23)/sd(y)
  ans_cv=bmodetest_cv(y,B=2,cv=FALSE)
  hist(y,freq=FALSE,breaks=30)
  lines(ans_cv$yp,ans_cv$fhat1,col=2)
  lines(ans_cv$yp,ans_cv$fhat2,col=3)
  lambdas[reps] = ans_cv$lam
  pvals[reps] = ans_cv$pvalue
}
stopImplicitCluster()
sum(pvals<=0.05) ##0.34


#################################################
############ mode on the endpoint ###############
######## use sample range as the support#########
#################################################
bmodetest <- function(y,lower = NULL, upper = NULL,B=1000,lam=NULL,eps1=.5,eps2=.2,eps=0.01,cv=TRUE,parallel=FALSE){
  n = length(y)
  y = sort(y)
  if(is.null(lower)){
    if(n > 5){
      s1 = min(y) - (y[5] - y[1])/4
    } else {
      s1 = min(y)
    }
  } else {
    s1 = lower
  }
  
  if(is.null(upper)){
    if(n > 5){
      s2 = y[n] + (y[n] - y[n-4])/4
    } else {
      s2 = max(y)
    }
  } else {
    s2 = upper
  }
  capk=round(n^(1/7)*12)
  qy=as.numeric(quantile(y,1:capk/(capk+1)))
  ####add more knots on both sides
  if(s1<qy[1]){
    k1=min(floor(log((qy[1]-s1)/3/(qy[2]-qy[1])+1, 3/2) -1 ), round(capk/8))
    if(k1<2){
      qy = c(qy,s1)
      if((qy[1]-s1)/(qy[2]-qy[1])>=2) qy = c(qy,qy[1]-(qy[1]-s1)/2)
    }else{
      q1 = ((qy[1]-s1)/(qy[2]-qy[1]))^(1/(k1+1))
      qy = c(qy, qy[1] - (qy[2]-qy[1])*q1^(1:(k1)), s1)
    }
  }
  if(s2>qy[capk]){
    k2=min(floor(log((s2-qy[capk])/3/(qy[capk]-qy[capk-1])+1, 3/2) - 1), round(capk/8))
    if(k2<2){
      qy=c(qy,s2)
      if((s2-qy[capk])/(qy[capk]-qy[capk-1])>=2) qy = c(qy,(s2-qy[capk])/2+qy[capk])
    }else{
      q2 = ((s2-qy[capk])/(qy[capk]-qy[capk-1]))^(1/(k2+1))
      qy = c(qy, (qy[capk]-qy[capk-1])*q2^(1:(k2))+qy[capk], s2)
    }
  }
  qy = sort(qy)
  ####add more knots between the two modes
  d=min(diff(qy))
  dd=floor(diff(qy)/min(diff(qy)))
  m1=min(which(dd==1))
  m2=max(which(dd==1))
  if(sum(dd[m1:m2]>2)>0){ for(i in which(dd[m1:m2]>2)){qy=c(qy,(qy[m1+i-1]+qy[m1+i])/2)}
  }
  kn=sort(qy)

  # if(!cv & is.null(lam)){
  #   K = kurtosis(y)
  #   if(K<2){
  #     lam = 10^2*n^(-1/7)
  #   } else if(K>2 & K<5){
  #     lam = 10^(4-K)*n^(-1/7)
  #   } else if(K>5 & K<9){
  #     lam = 10^(3/2-K/2)*n^(-1/7)
  #   } else{
  #     lam = 10^(-3)*n^(-1/7)
  #   }
  # }
  if(is.null(eps)){
    eps=ifelse(abs(skewness(y))>0.7,eps2,eps1)
  }else{eps1=eps2=eps}
  
  m=length(kn)+1
  bspl=bSpline(y,degree=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  yp=0:4000/4000*(s2-s1)+s1
  s1=min(s1,min(yp))
  s2=max(s2,max(yp))
  bp=bSpline(yp,degree=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  slopes=bSpline(kn,degree=2,derivs=1,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  D2=bSpline(kn,degree=2,derivs=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  d=(s2-s1)/m
  ###  make all basis functions integrate to one
  dp=yp[2]-yp[1]
  avec=1:m
  for(i in 1:m){avec[i]=sum(bp[,i])*dp}
  for(i in 1:m){
    bspl[,i]=bspl[,i]/avec[i]
    bp[,i]=bp[,i]/avec[i]
    slopes[,i]=slopes[,i]/avec[i]
    D2[,i]=D2[,i]/avec[i]
  }
  av1=avec
  avec=rep(1,m)
  hmat=matrix(nrow=m,ncol=m)
  cvec=1:m
  for(i in 1:m){
    for(j in i:m){
      pr=bp[,i]*bp[,j]
      hmat[i,j]=sum(pr)*dp
      hmat[j,i]=hmat[i,j]
    }
    cvec[i]=sum(bspl[,i])/n
  }
  b0=rep(1/m,m)
  
  wmat=matrix(0,nrow=m,ncol=m-1)
  for(i in 1:(m-1)){wmat[i,i]=-1;wmat[i+1,i]=1}
  
  D1=matrix(0,m-2,m-1)
  for(i in 1:(m-2)){D1[i,i]=-1;D1[i,i+1]=1}
  D=D1%*%D2*d^(5/2)
  ##  get unimodal
  ans2=bmfit(y,kn,hmat,cvec,slopes,b0,wmat,D,bspl,lam,cv=cv)
  ans1=umfit(y,kn,hmat,cvec,slopes,b0,wmat,D,bspl,lam = ans2$lam,eps=eps, cv=FALSE)
  t1=as.numeric((ans1$crit-ans2$crit))
  #t2=as.numeric((ans1$crit-ans2$crit)/abs(ans1$crit))
  fhat2=round(bp%*%ans2$bhat,10)
  dfhat2=diff(fhat2)
  outtb=NULL
  if(sum(dfhat2>0)==0 |sum(dfhat2<0)==0){
    pvalue=2
  }else{
    md=max(which(dfhat2>0))
    if (!is.unsorted(fhat2[1:(md-1)]) ) {
      pvalue=2
    } else {
      cdf1=bp%*%ans1$bhat
      for(i in 2:4001){
        cdf1[i]=cdf1[i-1]+cdf1[i]
      }
      cdf1=cdf1-min(cdf1)
      cdf1=cdf1/cdf1[4001]
      one_boot <- function(t) {
        yb=sapply(1:n,function(o){u=runif(1);id=min(which(u<cdf1));alp=(cdf1[id]-u)/(cdf1[id]-cdf1[id-1]);alp*yp[id-1]+(1-alp)*yp[id]})
        modet(yb,kn,hmat,slopes,b0,wmat,D,bspl,av1,bp,lam=ans2$lam,eps=eps)
      }
      if(parallel){
        outtb <- foreach(t=1:B,.combine='c') %dopar% {
          one_boot(t)
        }
      } else {
        outtb <- sapply(1:B, one_boot)
      }
      pvalue=sum(outtb>t1)/B
    }
  }
  ans=new.env()
  ans$yp=yp
  ans$fhat1=bp%*%ans1$bhat
  ans$fhat2=bp%*%ans2$bhat
  ans$statistic=t1
  ans$tb=outtb
  ans$pvalue=pvalue
  ans$lam=c(ans1$lam,ans2$lam)
  ans$kn=kn
  ans$crit=c(ans1$crit,ans2$crit)
  ans$kurtosis=kurtosis(y)
  ans$skewness=skewness(y)
  ans
}
##############
umfit=function(y,kn,hmat,cvec,slopes,b0,wmat,D,bspl,lam=NULL,eps,cv=cv){	
  n <- length(y)
  m=length(kn)+1
  amatl1=list()
  ## mode between t_k and t_{k+1} where k=1,...,m-2
  for(k in 1:(m-2)){
    amat=matrix(0,nrow=m+1,ncol=m)
    amat[1:k,]=slopes[1:k,]
    amat[(k+1):(m-1),]=-slopes[(k+1):(m-1),]
    amat[m,1]=1
    amat[m+1,m]=1
    epsvec=c(rep(eps/n^(2/7)/diff(range(kn))^2,k-1),0,0,rep(eps/n^(2/7)/diff(range(kn))^2,m-k-2),0,0)
    amatl1[[k]]=list(amat=amat, epsvec=epsvec)
  }
  amatl1[[m-1]]=list(amat = rbind(-slopes,c(rep(0,m-1),1)), epsvec = rep(0,m))
  amatl1[[m]]=list(amat = rbind(slopes,c(1, rep(0,m-1))), epsvec = rep(0,m))
  ## cross validation
  if(cv | is.null(lam)){
    zvec=t(wmat)%*%(cvec-hmat%*%b0-0.1*n^(-1/7)*t(D)%*%D%*%b0)
    qmat=t(wmat)%*%(hmat+0.1*n^(-1/7)*t(D)%*%D)%*%wmat 
    crit<-lapply(amatl1, function(x){ans <- quadprog::solve.QP(qmat,zvec,t(x[[1]]%*%wmat),x[[2]]-x[[1]]%*%b0);ans$value})
    amat1 <- amatl1[[which.min(crit)]][[1]]
    epsvec <- amatl1[[which.min(crit)]][[2]]
    ## number of folds
    L = 10
    lam_set = c(0.0001,0.001,0.01,0.1,1,10)
    fold_assignments <- sample(rep(1:L, length.out=n))
    crit_cv = NULL
    for( lamt in lam_set){
      err = rep(0,L)
      for (l in 1:L) {
        train_data <- y[fold_assignments != l]
        val_data   <- y[fold_assignments == l]
        bspl_train = bspl[fold_assignments != l,]
        bspl_val =  bspl[fold_assignments == l,]
        cvec_train=cvec_val=1:m
        n_train = nrow(bspl_train)
        n_val   = nrow(bspl_val)
        for(i in 1:m){
          cvec_train[i] = sum(bspl_train[,i]) / n_train
          cvec_val[i]   = sum(bspl_val[,i]) / n_val
        }
        zvec=t(wmat)%*%(cvec_train-hmat%*%b0-lamt*n_train^(-1/7)*t(D)%*%D%*%b0)
        qmat=t(wmat)%*%(hmat+lamt*n_train^(-1/7)*t(D)%*%D)%*%wmat 
        ans1=solve.QP(qmat,zvec,t(amat1%*%wmat),epsvec-amat1%*%b0)
        alphahat1=ans1$solution
        bhat1=wmat%*%alphahat1+b0
        err[l]=t(bhat1)%*%hmat%*%bhat1-2*sum(cvec_val*bhat1)
      }
      crit_cv = c(crit_cv,mean(err))
    }
    lamt=lam_set[which(crit_cv==min(crit_cv))]
    zvec=t(wmat)%*%(cvec-hmat%*%b0-lamt*n^(-1/7)*t(D)%*%D%*%b0)
    qmat=t(wmat)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%wmat 
    # crit<-lapply(amatl1, function(x){ans <- quadprog::solve.QP(qmat,zvec,t(x[[1]]%*%wmat),x[[2]]-x[[1]]%*%b0);ans$value})
    # amat1 <- amatl1[[which.min(crit)]][[1]]
    # epsvec <- amatl1[[which.min(crit)]][[2]]
    ans1=solve.QP(qmat,zvec,t(amat1%*%wmat),epsvec-amat1%*%b0)
    alphahat1=ans1$solution
    bhat1=wmat%*%alphahat1+b0
    cr1=t(bhat1)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%bhat1-2*sum(cvec*bhat1)
  }else{
    lamt=lam
    zvec=t(wmat)%*%(cvec-hmat%*%b0-lamt*n^(-1/7)*t(D)%*%D%*%b0)
    qmat=t(wmat)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%wmat 
    crit<-lapply(amatl1, function(x){ans <- quadprog::solve.QP(qmat,zvec,t(x[[1]]%*%wmat),x[[2]]-x[[1]]%*%b0);ans$value})
    amat1 <- amatl1[[which.min(crit)]][[1]]
    epsvec <- amatl1[[which.min(crit)]][[2]]
    ans1=solve.QP(qmat,zvec,t(amat1%*%wmat),epsvec-amat1%*%b0)
    alphahat1=ans1$solution
    bhat1=wmat%*%alphahat1+b0
    cr1=t(bhat1)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%bhat1-2*sum(cvec*bhat1) 
  }
  
  ans1=new.env()
  ans1$bhat=bhat1
  ans1$lam=lamt
  ans1$crit=cr1
  ans1
}
############################################################################
bmfit=function(y,kn,hmat,cvec,slopes,b0,wmat,D,bspl,lam=NULL,cv=cv){
  n=length(y)
  m=dim(slopes)[2]
  ## triplets (i,j,k) where i<j<k 
  trips=matrix(0,nrow=choose(m,3),ncol=3)
  nr=0
  amatl2=list()
  ## modes and antimodes between t_i and t_{i+1}, t_j and t_{j+1}, t_k and t_{k+1}, k<=m-2
  for(i in 1:(m-6)){
    for(j in (i+2):(m-4)){
      for(k in (j+2):(m-2)){
        nr=nr+1
        trips[nr,]=c(i,j,k)
        amat=matrix(0,m+1,m)
        amat[1:(m-1),1:m]=slopes
        amat[(i+1):j,]=-slopes[(i+1):j,]
        amat[k:(m-1),]=-slopes[k:(m-1),]
        amat[m,1]=1
        amat[m+1,m]=1
        amatl2[[nr]]=amat
        }
    }
  }
  ## modes at the lower endpoint, i.e. i=1 and mode is at t_1
  for(j in 3:(m-4)){
    for(k in (j+2):(m-2)){
      nr=nr+1
      trips[nr,]=c(1,j,k)
      amat=matrix(0,m,m)
      amat[1:(m-1),1:m]=slopes
      amat[1:j,]=-slopes[1:j,]
      amat[(k+1):(m-1),]=-slopes[(k+1):(m-1),]
      amat[m,m]=1
      amatl2[[nr]]=amat
    }
  }
  ## modes at upper endpoint, i.e. k=m-2 and mode is at t_{m-1}
  for(i in 1:(m-6)){
    for(j in (i+2):(m-4)){
      nr=nr+1
      trips[nr,]=c(i,j,m-2)
      amat=matrix(0,m,m)
      amat[1:(m-1),1:m]=slopes
      amat[(i+1):j,]=-slopes[(i+1):j,]
      amat[m,1]=1
      amatl2[[nr]]=amat
    }
  }
  ## both modes at the endpoints
  for(j in 3:(m-4)){
    nr=nr+1
    trips[nr,]=c(1,j,m-2)
    amat=matrix(0,m-1,m)
    amat[1:(m-1),1:m]=slopes
    amat[1:j,]=-slopes[1:j,]
    amatl2[[nr]]=amat
  }
  trips = trips[1:nr,,drop=FALSE]
  if(cv | is.null(lam)){
    zvec=t(wmat)%*%(cvec-hmat%*%b0-0.1*n^(-1/7)*t(D)%*%D%*%b0)
    qmat=t(wmat)%*%(hmat+0.1*n^(-1/7)*t(D)%*%D)%*%wmat 
    crit<-lapply(amatl2, function(x){ans <- quadprog::solve.QP(qmat,zvec,t(x%*%wmat),-x%*%b0);ans$value})
    amat2 <- amatl2[[which.min(crit)]]
    L = 10
    lam_set = c(0.00001,0.0001,0.001,0.01,0.1,1,10,100)
    fold_assignments <- sample(rep(1:L, length.out=n))
    crit_cv = NULL
    for( lamt in lam_set){
      err = rep(0,L)
      for (l in 1:L) {
        train_data <- y[fold_assignments != l]
        val_data   <- y[fold_assignments == l]
        bspl_train = bspl[fold_assignments != l,]
        bspl_val =  bspl[fold_assignments == l,]
        cvec_train=cvec_val=1:m
        n_train = nrow(bspl_train)
        n_val   = nrow(bspl_val)
        for(i in 1:m){
          cvec_train[i] = sum(bspl_train[,i]) / n_train
          cvec_val[i]   = sum(bspl_val[,i]) / n_val
        }
        zvec=t(wmat)%*%(cvec_train-hmat%*%b0-lamt*n_train^(-1/7)*t(D)%*%D%*%b0)
        qmat=t(wmat)%*%(hmat+lamt*n_train^(-1/7)*t(D)%*%D)%*%wmat 
        ans2=solve.QP(qmat,zvec,t(amat2%*%wmat),-amat2%*%b0)
        alphahat2=ans2$solution
        bhat2=wmat%*%alphahat2+b0
        err[l]=t(bhat2)%*%hmat%*%bhat2-2*sum(cvec_val*bhat2)
      }
      crit_cv = c(crit_cv,mean(err))
    }
    lamt=lam_set[which(crit_cv==min(crit_cv))]
    zvec=t(wmat)%*%(cvec-hmat%*%b0-lamt*n^(-1/7)*t(D)%*%D%*%b0)
    qmat=t(wmat)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%wmat 
    ans2=solve.QP(qmat,zvec,t(amat2%*%wmat),-amat2%*%b0)
    alphahat2=ans2$solution
    bhat2=wmat%*%alphahat2+b0
    cr2=t(bhat2)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%bhat2-2*sum(cvec*bhat2)
  }else{
    lamt=lam
    zvec=t(wmat)%*%(cvec-hmat%*%b0-lamt*n^(-1/7)*t(D)%*%D%*%b0)
    qmat=t(wmat)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%wmat 
    crit<-lapply(amatl2, function(x){ans <- quadprog::solve.QP(qmat,zvec,t(x%*%wmat),-x%*%b0);ans$value})
    amat2 <- amatl2[[which.min(crit)]]
    ans2=solve.QP(qmat,zvec,t(amat2%*%wmat),-amat2%*%b0)
    alphahat2=ans2$solution
    bhat2=wmat%*%alphahat2+b0
    cr2=t(bhat2)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%bhat2-2*sum(cvec*bhat2)
  }
  
  ans2=new.env()
  ans2$bhat=bhat2
  ans2$lam=lamt
  ans2$crit=cr2
  ans2
}
##################################
modet <- function(yb,kn,hmat,slopes,b0,wmat,D,bspl,av1,bp,lam,eps){
  n=length(yb)
  m=dim(slopes)[2]
  capk=length(kn) 
  s1=min(kn)
  s2=max(kn)
  bb=bSpline(yb,degree=2,knots=kn[2:(capk-1)],Boundary.knots=c(s1,s2),intercept=TRUE)
  m=dim(bb)[2]
  for(i in 1:m){bb[,i]=bb[,i]/av1[i]}
  cb=1:m
  for(i in 1:m){
    cb[i]=sum(bb[,i])/n
  }
  fit1=umfit(yb,kn,hmat,cb,slopes,b0,wmat,D,bspl,lam,eps=eps, cv=FALSE)
  fit2=bmfit(yb,kn,hmat,cb,slopes,b0,wmat,D,bspl,lam, cv=FALSE)
  fhat2=round(bp%*%fit2$bhat,10)
  dfhat2=diff(fhat2)
  if(sum(dfhat2>0)==0 |sum(dfhat2<0)==0){
    return (-1e+3)
  }else{md=max(which(dfhat2>0))
  if (!is.unsorted(fhat2[1:(md-1)]) ){
    return (-1e+3)
  }else{return(as.numeric(fit1$crit-fit2$crit))}
  }
}


n = 100
y = benchden::rberdev(n, dnum=23)/sd(y)
ans=bmodetest(y,B=1)
hist(y,freq=FALSE,breaks=30)
lines(ans$yp,ans$fhat1,col=2)
lines(ans$yp,ans$fhat2,col=3)
ans$lam
ans$pvalue



y = abs(rnorm(n,0,1))
y=y/sd(y)
registerDoParallel(9)
ans=bmodetest(y,B=100)
stopImplicitCluster()
hist(y,freq=FALSE,breaks=30)
lines(ans$yp,ans$fhat1,col=2)
lines(ans$yp,ans$fhat2,col=3)
rug(ans$kn)
ans$pvalue
ans$lam

##### sims #####
##### claw #####
set.seed(1234)
pvals = matrix(100,7,400)
lambdas = rep(0,400)
n=200
registerDoParallel(9)
for(j in 1:400){
  ### generate sample from claw distribution
  y = benchden::rberdev(n, dnum=23)
  y = y/sd(y)
  ans = bmodetest(y,B=200)
  lambdas[j] = ans$lam
  pvals[1, j] = ans$pvalue
  pvals[2, j] = modetest(y, method = "SI")$p.value
  pvals[3, j] = modetest(
    y,
    method = "HY",
    lowsup = quantile(y,0.01),
    uppsup = quantile(y,0.99)
  )$p.value
  pvals[4, j] = modetest(y, method = "FM")$p.value
  pvals[5, j] = modetest(y, method = "HH")$p.value
  pvals[6, j] = modetest(y, method = "CH")$p.value
  pvals[7, j] = modetest(y, method = "ACR")$p.value
  print(j)
}
stopImplicitCluster()
rowMeans(pvals < 0.05) 
#0.3250 0.0025 0.1200 0.1400 0.0725 0.6675 0.5475


##### sims #####
##### monotone #####
set.seed(1234)
pvals = matrix(100,7,400)
lambdas = rep(0,400)
n=200
registerDoParallel(9)
for(j in 1:400){
  ### generate sample from claw distribution
  y = abs(rnorm(n,0,1))
  y = y/sd(y)
  ans = bmodetest(y,B=200)
  lambdas[j] = ans$lam
  pvals[1, j] = ans$pvalue
  pvals[2, j] = modetest(y, method = "SI")$p.value
  pvals[3, j] = modetest(
    y,
    method = "HY",
    lowsup = quantile(y,0.01),
    uppsup = quantile(y,0.99)
  )$p.value
  pvals[4, j] = modetest(y, method = "FM")$p.value
  pvals[5, j] = modetest(y, method = "HH")$p.value
  pvals[6, j] = modetest(y, method = "CH")$p.value
  pvals[7, j] = modetest(y, method = "ACR")$p.value
  print(j)
}
stopImplicitCluster()
rowMeans(pvals < 0.05) 
## 0.0175 0.0000 0.0650 0.0875 0.0000 0.0950 0.0650



##### sims #####
##### bathtub #####
set.seed(1234)
pvals = matrix(100,7,400)
lambdas = rep(0,400)
n=200
registerDoParallel(9)
for(j in 1:400){
  ### generate sample from claw distribution
  u = runif(n)
  y = abs(rnorm(n,0,1))
  y[u<0.4] = - abs(rnorm(sum(u<0.4),0,1)) + 4
  #hist(y)
  y = y/sd(y)
  ans = bmodetest(y,B=200)
  lambdas[j] = ans$lam
  pvals[1, j] = ans$pvalue
  pvals[2, j] = modetest(y, method = "SI")$p.value
  pvals[3, j] = modetest(
    y,
    method = "HY",
    lowsup = quantile(y,0.01),
    uppsup = quantile(y,0.99)
  )$p.value
  pvals[4, j] = modetest(y, method = "FM")$p.value
  pvals[5, j] = modetest(y, method = "HH")$p.value
  pvals[6, j] = modetest(y, method = "CH")$p.value
  pvals[7, j] = modetest(y, method = "ACR")$p.value
  print(j)
}
stopImplicitCluster()
rowMeans(pvals < 0.05) 
##



hist(y,freq=FALSE,breaks=30)
lines(ans$yp,ans$fhat1,col=2)
lines(ans$yp,ans$fhat2,col=3)
ans$pvalue



###############################################
######### diagnose: lambda vs amat2 ###########
###############################################
findamat2 <- function(y,lower = NULL, upper = NULL){
  n = length(y)
  y = sort(y)
  if(is.null(lower)){
    if(n > 5){
      s1 = min(y) - (y[5] - y[1])/4
    } else {
      s1 = min(y)
    }
  } else {
    s1 = lower
  }
  
  if(is.null(upper)){
    if(n > 5){
      s2 = y[n] + (y[n] - y[n-4])/4
    } else {
      s2 = max(y)
    }
  } else {
    s2 = upper
  }
  capk=round(n^(1/7)*12)
  qy=as.numeric(quantile(y,1:capk/(capk+1)))
  ####add more knots on both sides
  if(s1<qy[1]){
    k1=min(floor(log((qy[1]-s1)/3/(qy[2]-qy[1])+1, 3/2) -1 ), round(capk/8))
    if(k1<2){
      qy = c(qy,s1)
      if((qy[1]-s1)/(qy[2]-qy[1])>=2) qy = c(qy,qy[1]-(qy[1]-s1)/2)
    }else{
      q1 = ((qy[1]-s1)/(qy[2]-qy[1]))^(1/(k1+1))
      qy = c(qy, qy[1] - (qy[2]-qy[1])*q1^(1:(k1)), s1)
    }
  }
  if(s2>qy[capk]){
    k2=min(floor(log((s2-qy[capk])/3/(qy[capk]-qy[capk-1])+1, 3/2) - 1), round(capk/8))
    if(k2<2){
      qy=c(qy,s2)
      if((s2-qy[capk])/(qy[capk]-qy[capk-1])>=2) qy = c(qy,(s2-qy[capk])/2+qy[capk])
    }else{
      q2 = ((s2-qy[capk])/(qy[capk]-qy[capk-1]))^(1/(k2+1))
      qy = c(qy, (qy[capk]-qy[capk-1])*q2^(1:(k2))+qy[capk], s2)
    }
  }
  qy = sort(qy)
  ####add more knots between the two modes
  d=min(diff(qy))
  dd=floor(diff(qy)/min(diff(qy)))
  m1=min(which(dd==1))
  m2=max(which(dd==1))
  if(sum(dd[m1:m2]>2)>0){ for(i in which(dd[m1:m2]>2)){qy=c(qy,(qy[m1+i-1]+qy[m1+i])/2)}
  }
  kn=sort(qy)
  
  m=length(kn)+1
  bspl=bSpline(y,degree=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  yp=0:4000/4000*(s2-s1)+s1
  s1=min(s1,min(yp))
  s2=max(s2,max(yp))
  bp=bSpline(yp,degree=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  slopes=bSpline(kn,degree=2,derivs=1,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  D2=bSpline(kn,degree=2,derivs=2,knots=kn[2:(m-2)],Boundary.knots=c(s1,s2),intercept=TRUE)
  d=(s2-s1)/m
  ###  make all basis functions integrate to one
  dp=yp[2]-yp[1]
  avec=1:m
  for(i in 1:m){avec[i]=sum(bp[,i])*dp}
  for(i in 1:m){
    bspl[,i]=bspl[,i]/avec[i]
    bp[,i]=bp[,i]/avec[i]
    slopes[,i]=slopes[,i]/avec[i]
    D2[,i]=D2[,i]/avec[i]
  }
  av1=avec
  avec=rep(1,m)
  hmat=matrix(nrow=m,ncol=m)
  cvec=1:m
  for(i in 1:m){
    for(j in i:m){
      pr=bp[,i]*bp[,j]
      hmat[i,j]=sum(pr)*dp
      hmat[j,i]=hmat[i,j]
    }
    cvec[i]=sum(bspl[,i])/n
  }
  b0=rep(1/m,m)
  
  wmat=matrix(0,nrow=m,ncol=m-1)
  for(i in 1:(m-1)){wmat[i,i]=-1;wmat[i+1,i]=1}
  
  D1=matrix(0,m-2,m-1)
  for(i in 1:(m-2)){D1[i,i]=-1;D1[i,i+1]=1}
  D=D1%*%D2*d^(5/2)

  ## triplets (i,j,k) where i<j<k 
  trips=matrix(0,nrow=choose(m,3)*2,ncol=3)
  nr=0
  amatl2=list()
  ## modes and antimodes between t_i and t_{i+1}, t_j and t_{j+1}, t_k and t_{k+1}, k<=m-2
  for(i in 1:(m-6)){
    for(j in (i+2):(m-4)){
      for(k in (j+2):(m-2)){
        nr=nr+1
        trips[nr,]=c(i,j,k)
        amat=matrix(0,m+1,m)
        amat[1:(m-1),1:m]=slopes
        amat[(i+1):j,]=-slopes[(i+1):j,]
        amat[k:(m-1),]=-slopes[k:(m-1),]
        amat[m,1]=1
        amat[m+1,m]=1
        amatl2[[nr]]=amat
      }
    }
  }
  ## modes at the lower endpoint, i.e. i=1 and mode is at t_1
  for(j in 3:(m-4)){
    for(k in (j+2):(m-2)){
      nr=nr+1
      trips[nr,]=c(1000,j,k)
      amat=matrix(0,m,m)
      amat[1:(m-1),1:m]=slopes
      amat[1:j,]=-slopes[1:j,]
      amat[(k+1):(m-1),]=-slopes[(k+1):(m-1),]
      amat[m,m]=1
      amatl2[[nr]]=amat
    }
  }
  ## modes at upper endpoint, i.e. k=m-2 and mode is at t_{m-1}
  for(i in 1:(m-6)){
    for(j in (i+2):(m-4)){
      nr=nr+1
      trips[nr,]=c(i,j,2000)
      amat=matrix(0,m,m)
      amat[1:(m-1),1:m]=slopes
      amat[(i+1):j,]=-slopes[(i+1):j,]
      amat[m,1]=1
      amatl2[[nr]]=amat
    }
  }
  ## both modes at the endpoints
  for(j in 3:(m-4)){
    nr=nr+1
    trips[nr,]=c(1000,j,2000)
    amat=matrix(0,m-1,m)
    amat[1:(m-1),1:m]=slopes
    amat[1:j,]=-slopes[1:j,]
    amatl2[[nr]]=amat
  }
  trips = trips[1:nr,,drop=FALSE]
  lam_set = c(0.00001,0.0001,0.001,0.01,0.1,1,10,100)
  
  amatindex = numeric(length(lam_set))
  mincrit = numeric(length(lam_set))
  secondcrit = numeric(length(lam_set))
  
  for(ll in seq_along(lam_set)){
    lamt = lam_set[ll]
    zvec = t(wmat)%*%(cvec-hmat%*%b0-lamt*n^(-1/7)*t(D)%*%D%*%b0)
    qmat = t(wmat)%*%(hmat+lamt*n^(-1/7)*t(D)%*%D)%*%wmat
    crit = sapply(amatl2, function(x){ans = quadprog::solve.QP(qmat,zvec,t(x%*%wmat),-x%*%b0)
      ans$value})
    ord = order(crit)
    amatindex[ll] = ord[1]
    mincrit[ll] = crit[ord[1]]
    secondcrit[ll] = crit[ord[2]]
  }
  out = data.frame(
    lam = lam_set,
    effective_lam = lam_set*n^(-1/7),
    amatindex = amatindex,
    i = trips[amatindex,1],
    j = trips[amatindex,2],
    k = trips[amatindex,3],
    mincrit = round(mincrit,5),
    secondcrit = round(secondcrit,5)
  )
  return(out)
}
findamat2(y)
##################################
set.seed(123)
n <- 200
R <- 2
## generators
dist_list <- list(
  # normal = function(n) {rnorm(n)},
  # half_normal = function(n) {abs(rnorm(n))},
  # exponential = function(n) {rexp(n)},
  # gamma = function(n) {rgamma(n, shape = 2)},
  # t3 = function(n) {rt(n, df = 3)},
  # bimodal = function(n) {z <- rbinom(n, 1, 0.5)
  #   rnorm(n, mean = ifelse(z == 1, -1.5, 1.5), sd = 0.7)
  # },
  claw = function(n) {benchden::rberdev(n, dnum = 23)},
  bathtub = function(n) {
    u <- runif(n)
    y <- abs(rnorm(n))
    id <- u < 0.5
    y[id] <- -abs(rnorm(sum(id))) + 4
    y
  },
  normal_mix_1 = function(n) {
    z <- rbinom(n, 1, 0.5)
    rnorm(n,
          mean = ifelse(z == 1, -1.5, 1.5),
          sd = 0.7)
  },
  ## symmetric, modes closer together
  normal_mix_2 = function(n) {
    z <- rbinom(n, 1, 0.5)
    rnorm(n,
          mean = ifelse(z == 1, -1, 1),
          sd = 0.7)
  },
  ## unequal mixture weights
  normal_mix_3 = function(n) {
    z <- rbinom(n, 1, 0.3)
    rnorm(n,
          mean = ifelse(z == 1, -2, 1),
          sd = 0.6)
  },
  ## different variances
  normal_mix_4 = function(n) {
    z <- rbinom(n, 1, 0.5)
    y <- numeric(n)
    y[z == 1] <- rnorm(sum(z == 1), -2, 0.4)
    y[z == 0] <- rnorm(sum(z == 0),  1, 1)
    y
  },
  chisq_mix_1 = function(n) {
    z <- rbinom(n, 1, 0.5)
    y <- numeric(n)
    y[z == 1] <- rchisq(sum(z == 1), df = 2)
    y[z == 0] <- rchisq(sum(z == 0), df = 5) + 6
    y
  },
  
  ## unequal weights
  chisq_mix_2 = function(n) {
    z <- rbinom(n, 1, 0.3)
    y <- numeric(n)
    y[z == 1] <- rchisq(sum(z == 1), df = 2)
    y[z == 0] <- rchisq(sum(z == 0), df = 4) + 7
    y
  },
  
  ## components with different skewness
  chisq_mix_3 = function(n) {
    z <- rbinom(n, 1, 0.5)
    y <- numeric(n)
    y[z == 1] <- rchisq(sum(z == 1), df = 1)
    y[z == 0] <- rchisq(sum(z == 0), df = 8) + 5
    y
  },
  ## normal + gamma
  normal_gamma_mix = function(n) {
    z <- rbinom(n, 1, 0.5)
    y <- numeric(n)
    y[z == 1] <- rnorm(sum(z == 1), 0, 0.6)
    y[z == 0] <- rgamma(sum(z == 0), shape = 3, scale = 0.5) + 3
    y
  },
  
  ## gamma mixture
  gamma_mix = function(n) {
    z <- rbinom(n, 1, 0.5)
    y <- numeric(n)
    y[z == 1] <- rgamma(sum(z == 1), shape = 2, scale = 0.5)
    y[z == 0] <- rgamma(sum(z == 0), shape = 5, scale = 0.4) + 4
    y
  }
)
all_results <- list()
for(dist_name in names(dist_list)) {
  cat("Running:", dist_name, "\n")
  dist_results <- list()
  for(r in 1:R) {
    y <- dist_list[[dist_name]](n)
    y <- y / sd(y)
    temp <- findamat2(y)
    temp$distribution <- dist_name
    temp$replicate <- r
    dist_results[[r]] <- temp
  }
  all_results[[dist_name]] <- do.call(rbind, dist_results)
}
results <- do.call(rbind, all_results)
rownames(results) <- NULL
results



amat_change <- aggregate(
  amatindex ~ distribution + replicate,
  data = results,
  FUN = function(x) length(unique(x))
)

names(amat_change)[3] <- "n_unique_amat"

aggregate(
  n_unique_amat ~ distribution,
  data = amat_change,
  FUN = function(x) {
    c(
      mean = mean(x),
      median = median(x),
      max = max(x)
    )
  }
)



## %% Set up  ----
setwd("/Volumes/SSD Hans/Github/EmpiricalIO")

library(haven)
library(tidyverse)
library(Hmisc)

# Helper packages
library(dplyr)     # for data wrangling
library(ggplot2)   # for awesome plotting

# Modeling packages
library(earth)     # for fitting MARS models
library(caret)     # for automating the tuning process

# Model interpretability packages
library(vip)       # for variable importance
library(pdp)       # for variable relationships
library(optimx)
library(plm)

## Uploading data ----
df <- read_dta("PS2/col_prod_data_clean.dta")
summary(df)
# write.csv(df, "df.csv")

# df$int <- rowSums(df[,c("rmats","renergy","rserv")])
# df$logint <- log(df$int)
df$select <- df$Irgd*df$Il*df$Irk*df$Irii

table(df$select, is.na(df$logrgd))
table(df$select, is.na(df$logrk))
## %% ACF estimation ----


obj <- function(theta, phi_hat, tdf) {
  #form ω_it( β_k,β_l )= ϕ_hat - β_k*k_it-β_ll*l_it
  k=tdf$logrk
  l=tdf$logl
  y=tdf$logrgd
  omega = phi_hat-theta[1]*k -theta[2]*l
  omegaep = y -theta[1]*k - theta[2]*l
  auxreg <- plm(omegaep~lag(omega),
                data=cbind(tdf,omega,omegaep),
                index=c("plant","year"))
  xiep <- residuals(auxreg)
  Z <- cbind(k,l)
  W <- diag(2)
  G = t(t(Z)%*%xiep)%*%W%*%(t(Z)%*%xiep)
  return(drop(G))
}

# First Stage ----
# I will use MARS to approximate yit=Φ_t(k_it,l_it,m_it)+ε_it
# Need to estimate for each industry

inds <- unique(df$sic3)
fit <- list()

for (i in 1:length(inds)) {
    #Select industry
    tdf <- data %>%
      filter(sic3==inds[i],
             select == 1)
    
    #Run MARS
    phi <- earth(
      y=tdf$logrgd,
      x=subset(tdf, select = c("logl","logrk","logrii")),
      degree = 3
    )
    # Get predictions
    phi_hat <- fitted(phi)
    
    ## Second stage --
    # init <- coef(lm(logrgd~logrk+logl+logrii, df))[c(2,3)]
    init <- c(0,0)
    fit[[i]] <- optimx(init,function(x)obj(x, phi_hat, tdf),
                  control = list(all.methods=TRUE))
}


fit[[1]]

coefs <- lapply(fit, function(x)coef(x[6,]))
lapply(fit, summary)
# omega <- df["logrgd"]-as.matrix(df[c("logrk","logl")])%*%t(coef(fit[[1]][6,]))
omega <- lapply(fit,
                function(x)df["logrgd"]-
                  as.matrix(df[c("logrk","logl")])%*%t(coef(x[6,])))
coef(fit[[1]])[6,]

sapply(fit, function(x)coef(x[6,]))

omega <- get_omega(df, coefs, inds)


summary(omega)
hist(omega)
lapply(omega, hist)

# Estimating naive results for comparison
naive <- lapply(inds, function(x)lm(logrgd~logrk+logl+logrii, df %>% filter(sic3==x)))

## Get SE Bootstrapping ----

## Bootstrap weird because of panel ----
table(df$plant)

# set.seed(33363)
# nbs=100
# draws <- list()
# for (i in 1:nbs) {
#     draws[[i]] <- sapply(inds, function(x){
#       t <- df %>% filter(sic3==x)
#       r <- runif(length(t), min=1,max = length(t))
#       as.integer(r)
#       }
#     )
# }
# df[draws[[1]][,1],] %>% count(plant)

obj_bs <- function(theta,phi_hat,tdf) {
  #form ω_it( β_k,β_l )= ϕ_hat - β_k*k_it-β_ll*l_it
  k=tdf$logrk
  l=tdf$logl
  y=tdf$logrgd
  omega = phi_hat-theta[1]*k -theta[2]*l
  omegaep = y -theta[1]*k - theta[2]*l
  auxreg <- plm(omegaep~lag(omega),
                data=cbind(tdf,omega,omegaep),
                index=c("year"))
  xiep <- residuals(auxreg)
  Z <- cbind(k,l)
  W <- diag(2)
  G = t(t(Z)%*%xiep)%*%W%*%(t(Z)%*%xiep)
  return(drop(G))
}

acf_bs <- function(inds,data,id){
    # data wrangling
    tdf <- data[id,] %>% 
      filter(sic3==inds) %>% #Select industry 
      filter(select == 1) #get rid of missing values
    
    #Run MARS
    phi <- earth(
      y=tdf$logrgd,
      x=subset(tdf, select = c("logl","logrk","logrii")),
      degree = 3
    )
    # Get predictions
    phi_hat <- fitted(phi)
    
    ## Second stage --
    init <- coef(lm(logrgd~logrk+logl+logrii, df))[c(2,3)]
    # init <- c(0,0)
    fit <- optimx(init,function(x)obj_bs(x,phi_hat,tdf),
                       method = "nlminb"
                       # control = list(all.methods=TRUE)
                       )
    return(coef(fit))
}

acf_bs(inds[1],df,1:74987)

library(boot)

boot(data = df, statistic = acf_bs, R=100, inds=inds[1], 
     strata = df$year, ncpus = 4)

bst <- list()
for (i in 1:length(inds)) {
  bst[[i]] <- boot(data = df %>% 
                     select(logrgd, logrk,logrk,
                            logl,logrii, plant, 
                            year, sic3, select), 
                   statistic = acf_bs, R=100, inds=inds[i], 
                   strata = df$year, ncpus = 4)
}

## Jackknife SE ----
library(resample)
obj_j <- function(theta,phi_hat,tdf) {
  #form ω_it( β_k,β_l )= ϕ_hat - β_k*k_it-β_ll*l_it
  k=tdf$logrk
  l=tdf$logl
  y=tdf$logrgd
  omega = phi_hat-theta[1]*k -theta[2]*l
  omegaep = y -theta[1]*k - theta[2]*l
  auxreg <- plm(omegaep~lag(omega),
                data=cbind(tdf,omega,omegaep),
                index=c("plant","year"))
  xiep <- residuals(auxreg)
  Z <- cbind(k,l)
  W <- diag(2)
  G = t(t(Z)%*%xiep)%*%W%*%(t(Z)%*%xiep)
  return(drop(G))
}

acf_j <- function(inds,data){
  # data wrangling
  tdf <- data %>% 
    filter(sic3==inds) %>% #Select industry 
    filter(select == 1) #get rid of missing values
  
  #Run MARS
  phi <- earth(
    y=tdf$logrgd,
    x=subset(tdf, select = c("logl","logrk","logrii")),
    degree = 3
  )
  # Get predictions
  phi_hat <- fitted(phi)
  
  ## Second stage --
  init <- coef(lm(logrgd~logrk+logl+logrii, df))[c(2,3)]
  # init <- c(0,0)
  fit <- optimx(init,function(x)obj_j(x,phi_hat,tdf),
                method = "nlminb"
                # control = list(all.methods=TRUE)
  )
  return(coef(fit))
}

jackknife(data = df %>% filter(sic3==inds[1]), 
          statistic =function(x)acf_j(inds[1],x) )

set.seed(55555)
my_jack<- function(data=df,R=100) {
  print("Starting Jackknife")
  coefs <- list()

  for (i in 1:length(inds)) {
    flush.console()
    #Get plants per industry
    plants <- unique(df %>% 
                       filter(sic3==inds[i]) %>% 
                       select(plant) %>% 
                       pull)
    s <- sample(plants, R)
    tmat <- list() 
    #Jackknife by blocks, removing one firm per iteration
    for(j in 1:R){
      #Remove one firm 
      tdf <- df %>% filter(plant!=s[j])
      #Estimate 
      tmat[[j]] <- acf_j(inds[i],tdf)
      if(is.integer(i/10)){print(paste0("Simulation ",i," done"))}
    }
    coefs[[i]] <- tmat
    print(paste0("Done for industry ", inds[i]))
  }
  return(coefs)
}

coefs_bs <- my_jack()
matrix(unlist(coefs_bs[1]), nrow = 100, byrow = TRUE)
coefs_bs <- lapply(coefs_bs,function(x) matrix(unlist(x), nrow = 100, byrow = TRUE))

## Results 1.a ----
# results <- tibble(Variables = c("Log K","Log L"), 
#                   `Coefficient Estimates (311)` = c(2*coef(fit[[1]])[6,1]-mean(coefs_bs[[1]][,1]),
#                                                     2*coef(fit[[1]])[6,2]-mean(coefs_bs[[1]][,2])),
#                   `Std. E. (311)` = std.error(coefs_bs[[1]]),
#                   `Coefficient Estimates (381)` = c(2*coef(fit[[2]])[6,1]-mean(coefs_bs[[2]][,1]),
#                                                     2*coef(fit[[2]])[6,2]-mean(coefs_bs[[2]][,2])),
#                   `Std. E. (381)` = std.error(coefs_bs[[2]]),
#                   `Coefficient Estimates (321)` = c(2*coef(fit[[3]])[6,1]-mean(coefs_bs[[3]][,1]),
#                                                     2*coef(fit[[3]])[6,2]-mean(coefs_bs[[3]][,2])),
#                   `Std. E. (321)` = std.error(coefs_bs[[3]]),
#                   `Coefficient Estimates (322)` = c(2*coef(fit[[4]])[6,1]-mean(coefs_bs[[4]][,1]),
#                                                     2*coef(fit[[4]])[6,2]-mean(coefs_bs[[4]][,2])),
#                   `Std. E. (322)` = std.error(coefs_bs[[4]]))
# 
# results
res <- get_bs_results(coefs_bs,coefs,inds,TRUE)
## 1 b ----

# sapply(omega,stats)
# 
# o_r <- sapply(omega,function(x){
#   stats<- quantile(x, probs = c(0, 0.25, 0.5, .75,1), na.rm=TRUE)
#   se <- std.error(x, na.rm = TRUE)
#   n <- sum(!is.na(x))
#   out <- rbind(stats[[1]],stats[[2]], stats[[3]], stats[[4]], stats[[5]],se,n)
#   return(out)
#   }
# )
# 
# row.names(o_r) <- c("Min","Q1","Mean","Q3","Max","SE","Obs")
# colnames(o_r) <- paste0("Industry ",inds)
# o_r

o_r <- sapply(omega,my_stats)
row.names(o_r) <- c("Min","Q1","Mean","Q3","Max","SE","Obs")
colnames(o_r) <- paste0("Industry ",inds)

dens <- lapply(omega, function(x)density(x[,1], na.rm=TRUE))
plot(NA,xlim=range(sapply(dens, "[", "x")), 
     ylim=range(sapply(dens, "[", "y")), 
     ann=FALSE)
mapply(lines, dens, col=1:4, 
       MoreArgs = list(xlab = "omega", lwd = 4, ylab="Density"))
legend("topright", legend=paste0("Industry: ",inds), fill=1:4)
title(main="Productivity distributions")


## 2 ----
#Drop all observation with 0 investment

coefs_n_inv <- lapply(inds,function(x)acf_j(x, df %>% filter(Irinv==1)))
coefs_bs_n_inv <- my_jack(data = df %>% filter(Irinv==1))

#Get coefficient estimates

get_bs_results <- function(bscfs=coefs_bs_n_inv, 
                           cfs=coefs_n_inv, 
                           ind=inds,
                           is_mat=FALSE){
  ifelse(is_mat, 
         cfs_bs <- bscfs,
         cfs_bs <- lapply(bscfs,function(x){ 
           matrix(unlist(x), nrow = length(x), byrow = TRUE)})
  )
  bs_mean_k <- sapply(cfs_bs,function(x)mean(x[,1]))
  bs_mean_l <- sapply(cfs_bs,function(x)mean(x[,2]))
  # bs_k <- mapply(function(x,y)2*x[1]-y,cfs,bs_mean_k)
  # bs_l <- mapply(function(x,y)2*x[2]-y,cfs,bs_mean_l)
  coefs <- matrix(unlist(cfs), nrow = length(cfs), byrow = TRUE)
  se <- sapply(cfs_bs,std.error)
  t <- cbind(coefs[,1],se[,1],coefs[,2],se[,2])
  tt <- t(t)
  row.names(tt) <- c("Log K","S.E.","Log L","S.E.")
  colnames(tt) <- paste0("Industry ",ind)
  print(tt)
  return(tt)
}

res_n_inv <- get_bs_results()

# Get productivity statistics

get_omega <- function(data=df, coefs=coefs_n_inv, ind=inds){
  mapply(function(x,z){
  y <- data %>% filter(sic3==x) %>% select(logrgd) %>% pull
  x <- data %>% filter(sic3==x) %>% select(logrk,logl) %>% as.matrix()
  omega <- y - x%*%t(z)
  return(omega)},ind,coefs)
}

# get_omega_bs <- function(data=df, bscfs=coefs_bs_n_inv, ind=inds){
# 
#   mapply(function(x,z){
#     y <- data %>% filter(sic3==x) %>% select(logrgd) %>% pull
#     x <- data %>% filter(sic3==x) %>% select(logrk,logl) %>% as.matrix()
#     omega <- y - x%*%t(z)
#     return(omega)},ind,coefs)
# }

o_n_inv <- get_omega(df %>% filter(Irinv==1), coefs_n_inv, inds)

my_stats <- function(x, ind=inds){
  stats<- quantile(x, probs = c(0, 0.25, 0.5, .75,1), na.rm=TRUE)
  se <- std.error(x, na.rm = TRUE)
  n <- sum(!is.na(x))
  out <- rbind(stats[[1]],stats[[2]], stats[[3]], stats[[4]], stats[[5]],se,n)
  
  # colnames(out) <- ind
  print(out)
  return(out)
}

o_sum <- sapply(o_n_inv,my_stats)
row.names(o_sum) <- c("Min","Q1","Mean","Q3","Max","SE","Obs")
colnames(o_sum) <- paste0("Industry ",inds)


# Compare with previous productivity estimates

dens2 <- lapply(o_n_inv,function(x)density(x, na.rm=TRUE))


par(mfrow=c(2,2))
for (i in 1:4) {
  plot(NA,xlim=range(c(dens[[i]]$x,dens2[[i]]$x)), 
       ylim=range(c(dens[[i]]$y,dens2[[i]]$y)), 
       ann=FALSE)
  mapply(lines, list(dens[[i]],dens2[[i]]), lty=1:2,
         MoreArgs = list(xlab = "omega", lwd = 4, col=i, ylab="Density"))
  legend("topright", legend=c("Full samle","Non-zero invesment"), 
         lty = 1:2, col = i, lwd = 2)
  title(main=paste0("Productivity distributions: Industry ",inds[i]))
}

save(o_sum,res_n_inv,res,inds,results, file = "prod.RData")
save.image()

# 3 ----
## Splitting by

df %>% 
  group_by(sic3, d_exp) %>% 
  summarise(n=n())

df %>% 
  group_by(sic3, d_imp) %>% 
  summarise(n=n())

df <- df %>% 
  group_by(sic3, year) %>% 
  mutate(avg_size=mean(logrk, na.rm =TRUE), big=logrk>avg_size) %>% 
  ungroup()

df %>% 
  group_by(sic3, big) %>% 
  summarise(n=n())

## Estimating 
### Exports
coefs_exp <- lapply(inds,function(x)acf_j(x, df %>% filter(d_exp==1)))
coefs_bs_exp <- my_jack(data = df %>% filter(d_exp==1))

coefs_n_exp <- lapply(inds,function(x)acf_j(x, df %>% filter(d_exp==0)))
coefs_bs_n_exp <- my_jack(data = df %>% filter(d_exp==0))
#### Imports
coefs_imp <- lapply(inds,function(x)acf_j(x, df %>% filter(d_imp==1)))
coefs_bs_imp <- my_jack(data = df %>% filter(d_imp==1))

coefs_n_imp <- lapply(inds,function(x)acf_j(x, df %>% filter(d_imp==0)))
coefs_bs_n_imp <- my_jack(data = df %>% filter(d_imp==0))

## Getting results
res_exp <- get_bs_results(coefs_bs_exp,coefs_exp,inds, FALSE)
res_imp <- get_bs_results(coefs_bs_imp,coefs_imp,inds, FALSE)

res_n_exp <- get_bs_results(coefs_bs_n_exp,coefs_n_exp,inds, FALSE)
res_n_imp <- get_bs_results(coefs_bs_n_imp,coefs_n_imp,inds, FALSE)

### Getting productivity

o_exp <- get_omega(df %>% filter(d_exp==1), coefs_exp, inds)
o_imp <- get_omega(df %>% filter(d_imp==1), coefs_imp, inds)

o_n_exp <- get_omega(df %>% filter(d_exp==0), coefs_n_exp, inds)
o_n_imp <- get_omega(df %>% filter(d_imp==0), coefs_n_imp, inds)

# Summary statistics of the productivity
o_sum_exp <- sapply(o_exp,my_stats)
row.names(o_sum_exp) <- c("Min","Q1","Mean","Q3","Max","SE","Obs")
colnames(o_sum_exp) <- paste0("Industry ",inds)

o_sum_imp <- sapply(o_imp,my_stats)
row.names(o_sum_imp) <- c("Min","Q1","Mean","Q3","Max","SE","Obs")
colnames(o_sum_imp) <- paste0("Industry ",inds)

o_sum_n_exp <- sapply(o_n_exp,my_stats)
row.names(o_sum_n_exp) <- c("Min","Q1","Mean","Q3","Max","SE","Obs")
colnames(o_sum_n_exp) <- paste0("Industry ",inds)

o_sum_n_imp <- sapply(o_n_imp,my_stats)
row.names(o_sum_n_imp) <- c("Min","Q1","Mean","Q3","Max","SE","Obs")
colnames(o_sum_n_imp) <- paste0("Industry ",inds)

## Compare densities

dens_exp <- lapply(o_exp, function(x)density(x, na.rm=TRUE))
dens_imp <- lapply(o_imp, function(x)density(x, na.rm=TRUE))

dens_n_exp <- lapply(o_n_exp, function(x)density(x, na.rm=TRUE))
dens_n_imp <- lapply(o_n_imp, function(x)density(x, na.rm=TRUE))

par(mfrow=c(2,2))
for (i in 1:4) {
  plot(NA,xlim=range(c(dens_n_exp[[i]]$x,dens_exp[[i]]$x)), 
       ylim=range(c(dens_n_exp[[i]]$y,dens_exp[[i]]$y)), 
       ann=FALSE)
  mapply(lines, list(dens_n_exp[[i]],dens_exp[[i]]), lty=1:2,
         MoreArgs = list(xlab = "omega", lwd = 4, col=i, ylab="Density"))
  legend("topright", legend=c("Non-Exporting","Exporting"), 
         lty = 1:2, col = i, lwd = 2)
  title(main=paste0("Productivity distributions: Industry ",inds[i]))
}

par(mfrow=c(2,2))
for (i in 1:4) {
  plot(NA,xlim=range(c(dens_n_imp[[i]]$x,dens_imp[[i]]$x)), 
       ylim=range(c(dens_n_imp[[i]]$y,dens_imp[[i]]$y)), 
       ann=FALSE)
  mapply(lines, list(dens_n_imp[[i]],dens_imp[[i]]), lty=1:2,
         MoreArgs = list(xlab = "omega", lwd = 4, col=i, ylab="Density"))
  legend("topright", legend=c("Non-Importing","Importing"), 
         lty = 1:2, col = i, lwd = 2)
  title(main=paste0("Productivity distributions: Industry ",inds[i]))
}

##Saving results to access in rmarkdown
save(list= ls(pattern = "^(o_sum|res|o_r)"), file = "prod.RData")

set.seed(123)

library(foreach)
library(matlib)

#number of products
J <- 100
#number of individuals
N <- 100000
#################################Q1####################################
#True coefficient
beta <-2
alpha <- 1
eta <- 1.5

#Observed product characteristics
x_j <- rnorm(J)
#Unobserved product characteristics
xi_j <- rnorm(J)
#perisonal shock J*N
epsilon_ij <- matrix(rgev(N*J),nrow=N)

#Observed cost shifter
w_j <- rnorm(J)
#Unobserved cost shifter
omega_j <- rnorm(J)

p_j <- matrix(eta*x_j+xi_j+w_j+omega_j)
delta_j <- beta*x_j-alpha*p_j+xi_j

#True utility
u_ij<- matrix(100,N,100)
foreach (j = 1:J) %do% {
  u_ij[,j] = delta_j[j]+epsilon_ij[,j]
    
}

#market share
#s_j<- matrix(0,J,1)
#s_j <- exp(delta_j)/(1+sum(exp(delta_j)))

#set mean utility of the outside good

epsilon_0 <- matrix(rgev(N*1),nrow=N)
#u_0 <- log(0.25*sum(exp(delta_j)))
u_0 <-log(0.25*sum(exp(delta_j)))
s_i0<- u_0+epsilon_0 #~20%
u_full_ij <- data.frame(cbind(s_i0,u_ij))
Choice <- matrix(colnames(u_full_ij)[apply(u_full_ij,1,which.max)])

s_j <- matrix(0,J+1,1)
foreach (j = 1:(J+1)) %do% {
  #cat("j ", j)
  #cat("Before s_j ", s_j[j])
  s_j[j] = sum(str_count(Choice, paste(paste("X",toString(j),sep=""),"\\b",sep="")))
  #cat("After s_j ", s_j[j])
  
}
#a)
s_j
s_j=s_j/sum(s_j)

#b)
temp <- matrix(0,(J+1),1)
temp <- log(s_j)-log(s_j[1])
delta_j_hat <- matrix(100,J,1)
delta_j_hat <-matrix(temp[2:(J+1),])
lm_delta = lm(delta_j_hat~x_j+p_j)

summary(lm_delta)

#c)
XX_j <- cbind(x_j,-p_j)
W_j <- cbind(x_j,w_j)
phi <- t(W_j) %*% W_j
theta <- inv(t(XX_j)%*%W_j%*%inv(phi)%*%t(W_j)%*%XX_j)%*%t(XX_j)%*%W_j%*%inv(phi)%*%t(W_j)%*%delta_j_hat
theta

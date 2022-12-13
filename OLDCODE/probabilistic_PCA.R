library("mvtnorm")
library("reshape2")
library("ggplot2")
library("dplyr")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# set the commonly used color
blue <- cbPalette[6]
red <- cbPalette[7]
green <- cbPalette[4]
orange <- cbPalette[2]
grey <- cbPalette[1]
yellow <- cbPalette[5]
purple <- cbPalette[8]
skyblue <- cbPalette[3]

# simulate with 100 diseases, 10 components and 10,000 individuals
D <- 100
M <- 2
N <- 10000

# define PCs
W0 <- matrix(c(rep(1, D/20), rep(0,D*19/20), rep(0,D*19/20), rep(1, D/20)) + rnorm(D*M, sd = 0.1), nrow = D, ncol = M) 
W1 <- matrix(c(rep(-1/20, D/20), rep(1/20,D/20),rep(0,D*9/10), rep(0,D*9/10), rep(1/20,D/20), rep(-1/20, D/20)) +  rnorm(D*M, sd = 0.01), nrow = D, ncol = M) 

# plot the lodings
melted_W0 <- melt(W0)
ggplot(melted_W0, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  scale_fill_gradient(name = "Number of \n contacts",low="white", high="red") +
  labs(x="PCs", y="Disease", title="W0") +
  theme(legend.position = "left", panel.background=element_blank())

melted_W1 <- melt(W1)
ggplot(melted_W1, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  scale_fill_gradient2(name = "Number of \n contacts",low="blue",mid="white", high="red",midpoint = 0) +
  labs(x="PCs", y="Disease",  title="W1") +
  theme(legend.position = "left", panel.background=element_blank())

# simulate PC components
Z <- rmvnorm(N, rep(0,M))

#################################################################################
# first attemp: generate disease age first then disease conditioning on age
#################################################################################

# generate T; generate age from pdf 0.5 + x, then scale to 30 years
num_ds <- rpois(N,10) + 1 # sample number of diseases for each individual
ds_age <- sapply(num_ds, function(x) ((runif(num_dis[x]) * 2 + 0.25)^.5 - .5)*20 )

# choose which disease to expose
T0 <- list() # output: T0 = 1 when there is a disease; 
T1 <- list() # T1 is the age of onset
for(n in 1:N){
  idx <- rep(NA, length(ds_age[[n]]))
  for(ds in 1:length(ds_age[[n]])){
    idx[ds] <- which.max( (W0 + W1 * ds_age[[n]][ds]) %*% Z[n,] ) 
  }
  t0 <- rep(0, D)
  t0[idx] <- 1
  t1 <- rep(0, D)
  t1[idx] <- ds_age[[n]]
  
  T0[[n]] <- t0
  T1[[n]] <- t1
}

#################################################################################
# second attemp: generate disease from a threshold model 
#################################################################################
generate_liability_disease <- function(S0, Drift, Variance, censored = 40){ # N: number of individauls 
  PATH <- matrix(nrow = length(S0), ncol = censored)
  V_matrix <- matrix(.2, censored, censored) + (1-.2) * diag(censored)
  Steps <- sapply(1:length(S0), function(x) rmvnorm(n = 1, mean = rep(Drift[x], censored),sigma =  V_matrix * Variance[x])) %>% t
  PATH <- apply(Steps, 1,  cumsum) %>% t
  PATH <- PATH + S0
  event <- apply(PATH, 1, function(x) ifelse(length(which(x > 0)) > 0, which(x > 0)[1], censored+1) )
  return(list(PATH, event))
}

T0 <- list() # output: T0 = 1 when there is a disease; 
T1 <- list() # T1 is the age of onset
N_sim <- 5
Trajectories <- list()
for(n in 1:N_sim){
  L <- matrix(nrow = D, ncol = 30)# save the trajectories  
  S0 <- rnorm(D, mean = 0, sd = 0.1)
    for(d in 1:D){
      Drift <- rep(NA, 30)
      for(t in 1:30){
        Drift[t] <- (W0[d,] + t * W1[d,] ) %*% Z[n,] 
        L[d,] <- generate_liability_disease(S0[d], Drift, 1, censored = 30)[[1]]
      }
    }
  Trajectories[[n]] <- L
}

paths <- Trajectories[[1]]

df_path <- data_frame(age = 41:70, path1 = paths[1,],path2 = paths[2,] ,path3 =paths[3,])

ggplot(df_path, aes(x = age)) + 
  theme(panel.background=element_blank()) + 
  geom_line(aes( y = path1, color = "Homogeneous wild type"),size = 1) +
  geom_line(aes( y = path2, color = "Heterogeneous"),size = 1) +
  geom_line(aes( y = path3, color = "Homogeneous risk alleles"),size = 1) +
  geom_hline(yintercept=10, linetype="dashed", color = red, size = 1) +
  scale_colour_manual(name="Genotype",values=c("Homogeneous wild type" = green, "Heterogeneous" = blue, "Homogeneous risk alleles" = red))



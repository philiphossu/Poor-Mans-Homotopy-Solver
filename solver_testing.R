library(mpoly); library(tidyverse)

# First Setup  ===================================================

# Creating a variables list so we know how many variables are used in target system
vars <- c("x1","x2","x3","x4")

# Enter target system as an mpoly object
target <- mp(c("x1^3*x2-x1^2+x2-2","x1^4*x2-x2^3+2*x2-127"))
target <- mp(c("x2-(x1^2)", "x2-(2 - x1^2)"))
target <- mp(c("x2-(x1+1)","x2-x1^2"))
target <- mp(c("x2-(x1+1)","x2-x1^2","x3-x1^2"))

# Create empty data frame
df <- tibble(x=rep(0,length(target)),y=rep(0,length(target)),z=rep(0,length(target)))
colnames(df) <- c("degree of nrow'th poly","startsys poly","startsys roots")

# Create & Solve Start System  ==================================

# Fill first column with total degree of each polynomial corresponding to the row #
df[1] <- mapply(totaldeg,target)

# Write pretty start systems in column 2
i <- 1;
while(i <= nrow(df)){
  df$`startsys poly`[i] <- paste0(vars[i],'^',as.character(df$`degree of nrow'th poly`[i]),"-1")
  i <- i + 1
}

# Setup matrix to hold information for start systems (for use in polyroot)
for_start_solving <- matrix(0, nrow=max(df$`degree of nrow'th poly`)+1, ncol=length(target))
for_start_solving[1,] <- -1

# Fill in support entries for polyroot to use in solving
i <- 1
while(i <= ncol(for_start_solving)){
  for_start_solving[df$`degree of nrow'th poly`[i]+1,i] <- 1
  i <- i + 1
}

# Solve start systems, store results in df
i <- 1
while(i <= ncol(for_start_solving)){
  df$`startsys roots`[i] <- list(polyroot(for_start_solving[,i]))
  i <- i + 1
}

# df

# Make combos of solutions, store in start_solutions and label
start_solutions <- as.tibble(expand.grid(df$`startsys roots`))
colnames(start_solutions) <- c(vars[1:length(start_solutions)])

# Setup Homotopy, Target Sys, Start Sys  ===============================

# Convert target system to a R function using mpoly
g <- as.function(target, varorder = vars)

# Convert start system to mpoly object then to R function
f <- as.function(mp(c(df$`startsys poly`)), varorder = vars)

# Define Homotopy
h <- function(t) function(v) (1-t)*g(v) + t*f(v)

# Define t-vals
t_vals <- rev(seq(0, 1, 0.05))[-1]
start_solutions[as.character(c(t_vals))] <- NA

# Add column in start solutions tibble to allow for intermediate solution recording
# start_solutions$tracking <- rep(list(0),nrow(start_solutions))

# Testing  =========================================================

# Calculate derivatives using mpoly

g  <- function(x) c(x[1]^2 + x[2]^2 - 1, x[1]*x[2] - 1)
dg <- function(x) matrix(c(
  2*x[1], 2*x[2],
  x[2],   x[1]
), nrow = 2, ncol = 2, byrow = TRUE)

g_mp <- mp(c("x1^2 + x2^2 - 1", "x1*x2 - 1"))

(test <- lapply(g_mp,deriv,var=vars[1:nrow(df)]))
# test <- as.mpoly(as.tibble(lapply(g_mp,deriv,var=vars[1:nrow(df)])))
test
as.mpoly(as.character(test[[1]]))
mp(test[[1]][1])

make_mp <- function(.){
  mp(.[1],.[2],.[3],.[4])
}

# as.function(c(deriv(g_mp[[1]],var=vars[1:nrow(df)])))
# test <- matrix(test,nrow=2,ncol=2)

# test <- mp(unlist(test, recursive = TRUE))
# as.function(test[,1])

testFun <- function(x) matrix(c(lapply(test,as.function)),nrow=2,ncol=2)

# CMD+SHIFT+C to remove comments

# i <- 1;
# for(i in 1:nrow(start_solutions)){
#   init <- as.vector(unname(unlist(start_solutions[i,])))[1:nrow(df)]
#   init_re <- init
#   # init_re <- Re(init)
#   # init_im <- Im(init)
#   #j <- 1
#   for (t in t_vals) {
#     o_re <- optim(
#       par = init_re,
#       fn = function(x) h(t)(x) %>% abs () %>% sum()
#     )
#     #o_im <- optim(
#     #  par = init_im,
#     #  fn = function(x) h(t)(x) %>% abs () %>% sum()
#     #)
#     init_re <- o_re$par
#     # init_im <- o_im$par
#     #start_solutions[[i,j+2]] <- list(o$par)
#     #j <- j + 1
#   }
#   #start_solutions[[i,j+2]] <- list(c(init[1:nrow(df)]))
#   #start_solutions$'0'[i] <- list(complex(real=init_re, imaginary=init_im))
#   start_solutions$'0'[i] <- list(c(init[1:nrow(df)]))
#
# }


# View Solutions  =====================================================
start_solutions$'0'
unique(start_solutions$'0')

start_solutions
start_solutions$'0.9'

# Check vs Bertini  ===================================================

library(bertini)

# target <- mp(c("x1^3*x2-x1^2+x2-2","x1^4*x2-x2^3+2*x2-127"))

# lhs <- mp(c("x2-(x1+1)","x2-x1^2"))

rhs <- mp(c("0","0"))

# poly_solve(lhs, rhs, c("x1","x2"))

poly_solve(target, rhs, c("x1","x2"))







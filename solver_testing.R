library(mpoly); library(tidyverse); library(kumerical); library(bertini)

# Define Functions  ===================================================

# Generate roots of unity for solving start systems
gen_roots_unity <- function(n){
  unity_sols <- c()
  for(k in 1:n){
    unity_sols <- c(unity_sols,complex(1,real=cos((2*pi*k)/n),imaginary=sin((2*pi*k)/n)))
  }
  unity_sols
}

# Path tracker which utilizes kumerical's multinewton function for
track_paths <- function(init,t_vals){
  for(i in t_vals){
    init <- multinewton(h(i),h_d(i),init)$root
  }
  init
}

# Generate Jacobian based on Dr. Kahle's function
jacobian <- function(poly_system, varorder = vars(poly_system)) {
  separate_mpolys <- lapply(poly_system, deriv, var = varorder)
  gradients <- lapply(separate_mpolys, as.function, varorder = varorder, silent = TRUE)
  J <- function(.) lapply(gradients, function(f) f(.))
  function(v) do.call(rbind, J(v))
}

# Define Variables, Setup  ================================================

# Create empty data frame for holding start system info
df <- tibble(x=rep(0,length(target)),y=rep(0,length(target)),z=rep(0,length(target)))
colnames(df) <- c("degree of nrow'th poly","startsys poly","startsys roots")

# Creating a variables list so we know how many variables are used in target system
vars <- c("x1","x2","x3","x4")

# Defining random complex gamma for the homotopy equation
gamma <- complex(1, runif(1), runif(1))

# Enter target system as an mpoly object
  # target <- mp(c("x1^3*x2-x1^2+x2-2","x1^4*x2-x2^3+2*x2-127"))
  # target <- mp(c("x2-(x1^2)", "x2-(2 - x1^2)"))
target <- mp(c("x1^2 + x2^2 - 1", "x1*x2 - 1"), varorder = vars)

# Fill first column of df with total degree of each polynomial corresponding to the row #
df[1] <- mapply(totaldeg,target)

# Write pretty start systems in column 2
i <- 1;
while(i <= nrow(df)){
  df$`startsys poly`[i] <- paste0(vars[i],'^',as.character(df$`degree of nrow'th poly`[i]),"-1")
  i <- i + 1
}

# Define t-vals for continuation
t_vals <- rev(seq(0, 1, 0.01))[-1]

# Solving Start System  ================================================

# Solve start systems, store results in df
i <- 1
while(i <= ncol(for_start_solving)){
  df$`startsys roots`[i] <- list(gen_roots_unity(df$`degree of nrow'th poly`[i]))
  i <- i + 1
}

# Make combos of solutions, store in start_solutions and label
start_solutions <- expand.grid(df$`startsys roots`)
colnames(start_solutions) <- c(vars[1:length(start_solutions)])

# Setup Homotopy, Target, Start, Derivs as Functions  ===================

# Convert target system to a R function using mpoly
g <- as.function(target, varorder = vars(target))
g_d <- jacobian(target)

# Convert start system to mpoly object then to R function
start_as_mp <- mp(c(df$`startsys poly`))
f <- as.function(start_as_mp, varorder = vars(start_as_mp))
f_gam <- function(.) gamma*f(.)

f_d <- jacobian(start_as_mp)
f_d_gam <- function(.) gamma*f_d(.)

# Define Homotopy
h <- function(t) function(v) (1-t)*g(v) + t*f_gam(v)
h_d <- function(t) function(v) (1-t)*g_d(v) + t*f_d_gam(v)

# Solving  ==============================================================

(solutions <- apply(start_solutions, 1, track_paths, t_vals))

# Check Solutions, Bertini  =============================================

# Should be ~0 for all entries if solutions were successful
apply(solutions, 1, g)

# Should be the same as the solutions we found
rhs <- mp(c("0","0"))
poly_solve(target, rhs, c("x1","x2"))







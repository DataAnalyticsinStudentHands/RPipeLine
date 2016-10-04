# Script to create data frames including variables of different
# type to test the genetic algorithm

# Number of columns per variable
cont <- 5     #continuous
natu <- 2     #natural integers
bool <- 3     #boolean
set.seed(1)
# The matrix X is initially created to store the exposures
X <- matrix(data = 0, nrow = 1000, ncol = cont + natu + bool)

# Continuous exposures have a uniform distribution over a user defined interval
for (i in 1:cont) {
  X[,i] <- runif(nrow(X), 0, 10)
}

# Natural exposures have a rounded uniform distribution over a user defined interval
for (i in (cont+1):(cont+natu)) {
  X[,i] <- round(runif(nrow(X), 0, 20))
}

# Boolean exposures have a bernoulli distribution
for (i in (cont+natu+1):(cont+natu+bool)) {
  X[,i] <- rbinom(nrow(X), 1, 0.4)
}

# Embed a rule: the logical array "event" will be true if a
# certain condition happens in the correspondent row in X

# In this case, we impose a baseline probability, and a higher
# probability of being a case if X(1) >= 6.0 and X(10) == TRUE
event <- matrix(data = 0, nrow = 1000, ncol = 1)
Pexp <- 0.5
Pbase <- 0.1

for (i in 1:nrow(event)) {
  if (X[i,1] > 6.0 && X[i,10] == 1)
    event[i,1] <- runif(1) < Pexp
  else event[i,1] <- runif(1) < Pbase
}

# Convert to data frame and save
X <- as.data.frame(X)
X[,(cont+1):(cont+natu)] <- apply(X[,(cont+1):(cont+natu)], 2, function(x) as.integer(x))
X[,(cont+natu+1):(cont+natu+bool)] <- apply(X[,(cont+natu+1):(cont+natu+bool)], 2, function(x) as.logical(x))
event <- as.data.frame(event)
colnames(event) <- c("Event")
X <- cbind(event, X)

save(X, file = "embedded_V1V10_6TRUE_05.RData")

# Script to create data frame to test the genetic algorithm

# The matrix X is initially created to store the exposures
X <- matrix(data = 0, nrow = 1000, ncol = 10)

# Each exposure has a uniform distribution over a user defined interval
for (i in 1:ncol(X)) {
  X[,i] <- runif(nrow(X), 0, 10)
}

# Embed a rule: the logical array "event" will be true if a
# certain condition happens in the correspondent row in X

# In this case, we impose a baseline probability, and a higher
# probability of being a case if X(1) >= 6.0
event <- matrix(data = 0, nrow = 1000, ncol = 1)
Pexp <- 0.5
Pbase <- 0.1

for (i in 1:nrow(event)) {
  if (X[i,1] < 6.0)
    event[i,1] <- runif(1) < Pbase
  else event[i,1] <- runif(1) < Pexp
}

# Convert to data frame and save
X <- as.data.frame(X)
event <- as.data.frame(event)
colnames(event) <- c("Event")
X <- cbind(event, X)

save(X, file = "embedded_V1_6_05.RData")

library(deSolve)

# time sequence 
time <- seq(0, 91, by = 1)

# parameters: a named vector
parameters <- c(r = 2, p = 1e+11, K = 1e+12)

# initial condition: a named vector
state <- c(y = 0, R = 1e+12)

# R function to calculate the value of the derivatives at each time value
# Use the names of the variables as defined in the vectors above
lotkaVolterra <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    dy = p
    dR = -p + r * ((K-R)/K) * R
    return(list(c(dy, dR)))
  })
}
## Integration with 'ode'
out <- ode(y = state, times = time, func = lotkaVolterra, parms = parameters)

## Plotting
out.df = as.data.frame(out) # required by ggplot: data object must be a data frame



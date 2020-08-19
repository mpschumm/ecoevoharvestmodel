# For obtaining structural mass
structural_mass <- function(l) {
  return(c_1*(l)^3)
}

# Calculate the instrinsic length-dependent E/S ratio that the fish is trying to attain
intrinsic_lambda <- function(l, l_bar, r, lambda_min, lambda_max) {
  return(lambda_min + (lambda_max-lambda_min)*(exp(r*(l-l_bar)))/(1+exp(r*(l-l_bar))))
}

# Energy intake for the fish per year, dependent on fish structural biomass
size_dependent_energy_intake <- function(S, I_0, I_1) {
  return(I_0*(S)^(I_1))
}

# Calculate the cost of the year of maintaining current structural and reversible biomass
maintenance_cost <- function(S, R, c_s, c_r) {
  return(c_s*S + c_r*R)
}

# Function for calculating a new, reduced level of reversible biomass, if net energy intake is negative
negative_I_net <- function(R, I_net, e_r) {
  if (I_net < 0) {
    return(R + (I_net/e_r))}
}

# Calculate the amount of energy needed to raise the R/S ratio to the length-dependent maximum
I_R_raising_lambda <- function(lambda_a, S, R) {
  return((lambda_a*S - R))
}

# Calculating proportion of any remaining energy to be converted to structural mass, to maintain R/S
remaining_energy_allocation <- function(e_s, e_r, lambda_a) {
  return((lambda_a*e_s)/(lambda_a*e_s + e_r))
}
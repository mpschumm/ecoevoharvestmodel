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

# Obtain R and S for a cohort prior to its reproduction (or failure to reproduce) that year
get_pre_reproductive_size <- function(R_past, l_past, a) {
  S_past = structural_mass(l_past)
  lambda_a = intrinsic_lambda(l_past, l_bar, r, lambda_min, lambda_max)
  I_net = size_dependent_energy_intake(S_past, I_0, I_1) - maintenance_cost(S_past, R_past, c_s, c_r)
  if (I_net < 0) {
    R <- max(0,negative_I_net(R, I_net, e_r))
  }
  if (I_net >= 0 & lambda_a > (R_past/S_past)) {
    # Gives amount of grams needed to raise R/S to lambda_l - not taking efficiency into account yet
    I_R = I_R_raising_lambda(lambda_a, S_past, R_past)
    if ((I_R/e_r - I_net) >=0) {
      R_past = I_net*e_r + R_past
      I_net = 0
    }
    if ((I_R/e_r - I_net) <0) {
      R_past = I_R + R_past
      I_net = I_net - I_R/e_r
    }
  }
  if (I_net >= 0) {
    R = R_past + remaining_energy_allocation(e_s, e_r, lambda_a)*e_r*I_net
    S = S_past + (1-remaining_energy_allocation(e_s, e_r, lambda_a))*e_s*I_net
  }
  l = (S/c_1)^(1/3)
  output <- c(R,l)
  return(output)
}
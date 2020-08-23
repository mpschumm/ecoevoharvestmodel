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
maintenance_cost <- function(S, E, c_S, c_E) {
  return(c_S*S + c_E*E)
}

# Function for calculating a new, reduced level of reversible biomass, if net energy intake is negative
negative_I_net <- function(E, I_net, e_E) {
  if (I_net < 0) {
    return(E + (I_net/e_E))}
}

# Calculate the amount of energy needed to raise the R/S ratio to the length-dependent maximum
I_E_raising_lambda <- function(lambda_a, S, E) {
  return((lambda_a*S - E))
}

# Calculating proportion of any remaining energy to be converted to structural mass, to maintain R/S
remaining_energy_allocation <- function(e_S, e_E, lambda_a) {
  return((lambda_a*e_S)/(lambda_a*e_S + e_E))
}

# Obtain R and S for a cohort prior to its reproduction (or failure to reproduce) that year
get_pre_reproductive_size <- function(E_past, l_past, a) {
  S_past = structural_mass(l_past)
  lambda_a = intrinsic_lambda(l_past, l_bar, r, lambda_min, lambda_max)
  I_net = size_dependent_energy_intake(S_past, I_0, I_1) - maintenance_cost(S_past, E_past, c_S, c_E)
  if (I_net < 0) {
    E <- max(0,negative_I_net(E, I_net, e_E))
  }
  if (I_net >= 0 & lambda_a > (E_past/S_past)) {
    # Gives amount of grams needed to raise E/S to lambda_l - not taking efficiency into account yet
    I_E = I_E_raising_lambda(lambda_a, S_past, E_past)
    if ((I_E/e_E - I_net) >=0) {
      E_past = I_net*e_E + E_past
      I_net = 0
    }
    if ((I_E/e_E - I_net) <0) {
      E_past = I_E + E_past
      I_net = I_net - I_E/e_E
    }
  }
  if (I_net >= 0) {
    E = E_past + remaining_energy_allocation(e_S, e_E, lambda_a)*e_E*I_net
    S = S_past + (1-remaining_energy_allocation(e_S, e_E, lambda_a))*e_S*I_net
  }
  l = (S/c_1)^(1/3)
  output <- c(E,l)
  return(output)
}
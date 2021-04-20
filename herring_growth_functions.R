# GROWTH FUNCTIONS

# For obtaining structural mass
structural_mass <- function(l) {
  return(c_1*(l)^c_2)
}

# Energy intake for the fish per year, dependent on fish structural biomass
size_dependent_energy_intake <- function(S) {
  # Numeric values from Ecological Archives E093-075-A1
  # Attack rate
  AR =  AR_1 * (((AR_2)*S/M_opt) * exp(1-((AR_2)*S/M_opt)) )^(AR_3)
  # Digestion
  H = H_1*((AR_2*S)^(H_2))
  return( timescale*((AR*AR_4)/(AR*AR_4*H+1))*K_e )
}

# Calculate the intrinsic length-dependent E/S ratio that the fish is trying to attain
intrinsic_lambda <- function(l, l_bar, r, lambda_min, lambda_max) {
  return(lambda_min + (lambda_max-lambda_min)*(exp(r*(l-l_bar)))/(1+exp(r*(l-l_bar))))
}

# Calculate the cost of the year of maintaining current structural and reversible biomass
maintenance_cost <- function(S, E, c_S, c_E) {
  return(c_S*S + c_E*E)
}

# Function for calculating a new, reduced level of reversible biomass, if net energy intake is negative
E_for_maintenance <- function(E, p_net, e_E) {
  E[E < 0] <- E[E < 0] + p_net*(1/e_E)
  return(E)
}

# Calculate the amount of energy needed to raise the R/S ratio to the length-dependent maximum
p_E_raising_lambda <- function(lambda_l, S, E) {
  return((lambda_l*S - E))
}

# Calculating proportion of any remaining energy to be converted to reversible mass, to maintain R/S
remaining_energy_allocation_E <- function(e_S, e_E, lambda_l) {
  return((lambda_l*e_S)/(lambda_l*e_S + e_E))
}

# Same as above, but for structural mass
remaining_energy_allocation_S <- function(e_S, e_E, lambda_l) {
  return(1-(lambda_l*e_S)/(lambda_l*e_S + e_E))
}

get_pre_reproductive_size <- function(E_past, l_past, resource, new_mu_exp) {
  S_past = structural_mass(l_past)
  S = S_past
  E = E_past
  lambda_l = intrinsic_lambda(l_past, l_bar, r, lambda_min, lambda_max)
  consumed_resource <- sum(size_dependent_energy_intake(S_past)*new_mu_exp*(resource/(K_half+resource)))
  p_net = size_dependent_energy_intake(S_past)*(resource/(K_half+resource)) - maintenance_cost(S_past, E_past, c_S, c_E)
  E = E_for_maintenance(E, p_net, e_E)
  p_net[p_net < 0] <- 0
  p_E <- p_E_raising_lambda(lambda_l, S, E)
  surplus <- (p_E*(1/e_E)) - p_net
  surplus <- surplus*(-1)
  E[surplus <= 0] <- p_net[surplus <= 0]*e_E + E[surplus <= 0]
  p_net[surplus <= 0] <- 0
  E[surplus > 0] <- p_E[surplus > 0] + E[surplus > 0]
  p_net[surplus > 0] <- p_net[surplus > 0] - p_E[surplus > 0]*(1/e_E)
  E[surplus > 0] = E[surplus > 0] + e_E*p_net[surplus > 0]*remaining_energy_allocation_E(e_S, e_E, lambda_l)[surplus > 0]
  S[surplus > 0] = S[surplus > 0] + e_S*p_net[surplus > 0]*remaining_energy_allocation_S(e_S, e_E, lambda_l)[surplus > 0]
  l = (S*(1/c_1))^(1/3)
  output <- list(E,l,consumed_resource)
  return(output)
}
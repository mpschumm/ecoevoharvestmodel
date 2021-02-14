library(pbapply)
library(truncnorm)

### COD MODEL

# Set longevity of the fish species, runtime of the model, and number of distinct genotypes
longevity_days = 9125
# set runtime to any arbitrary desired value, in units of days, that is a multiple of 365
runtime_days = 150*365
genotypes = 3
# Setting the length of the timestep in days - set to any fraction of 365
timescale= (365)/4
# Setting run time and longevity to units of time steps
runtime = round(runtime_days/timescale)
longevity = round(longevity_days/timescale)
# How long to run model for convergence before addition of genetic diversity
# in units of days
add_genotypes = 33440

# Set parameters
# Need to add density dependence
# The steepness of the increase in the ratio lambda (reversible over structural mass)
r = 6
# The length at which the increasing function of lambda reaches its inflection point
l_bar = 0.2
# Scaling parameter for obtaining structural mass from standard length
c_1 = 5787
c_2 = 3
# Multiplier for obtaining overall yearly energy intake as a function of structural mass
p_0 = 0.1*timescale
# Exponent for obtaining overall yearly energy intake as a function of structural mass
p_1 = 2/3
# Maximum ratio lambda (reversible over structural mass) a fish can attain over its life
lambda_max = 1.3
# Minimum ratio lambda a fish can attain over its life
lambda_min = 0
# Multiplier for obtaining cost of maintenance of structural mass
c_S = 0.003*timescale
# Multiplier for obtaining cost of maintenance of reversible mass
c_E = 0.0003*timescale
# Efficiency with which energy is converted into reversible mass
e_E = 0.9
# Efficiency with which energy is converted to structural mass
e_S = (1/3)
# Multiplier and exponent for the cost of reproduction (in grams)
r_0 = 6
r_1 = 0.6
# Effectively, reproductive investment (fraction of reversible mass devoted to reproduction)
w= 0.674
# Instantaneous base mortality rate, derived from daily instantaneous probability of survival for M=-0.2, multiplied to -1 to make positive
m = log(0.999452^timescale)*(-1)
# Instantaneous fishing mortality
Fishing = 0
# Ricker
alpha = 8.44e-09
b = 4.43e-15
# Parameter setting size-selection ogive steepness
q_f = 0.2
# Minimum catch size for the species (in meters)
mincatchsize = 0.5
# Total population biomass
b_t=0
# Increased natural mortality rate at E/S = 0 or S = 0 
# For roughly 0.5 year-based instantaneous mortality at ~1 yrs., Koster et al 2003
m_p_max=log(0.9986^timescale)*(-1)
m_c_max=log(0.95^timescale)*(-1)
# Egg packing constant
packing = 4.45e6
# Parameter for exponential decline in condition-dependent mortality with improved body condition
z_c=7
# Parameter for exponential decline in predation mortality with increased body length
z_p=8
# Resource density
RD = 9000000
# Consumption with resource density half-saturation constant
K_RD = 1000000
# Variability in resource
sd_RD = 0
# How much energy intake is affected by changes in resources (1=no change)
resource = 1

# Initializing the matrices for the first time point
timepoints <- vector(mode = "list", length = runtime)
timepoint1 <- list(matrix(data=0, nrow=longevity, ncol=genotypes), 
                   matrix(data=0, nrow=longevity, ncol=genotypes),
                   matrix(data=0, nrow=longevity, ncol=genotypes),
                   array(data=0, dim=c(longevity, genotypes, genotypes)))
timepoint1[[1]][1,] <- 3.938
timepoint1[[2]][1,] <- 0.1102
timepoint1[[3]][1,] <- 100
timepoint1[[4]][1,,] <- diag(genotypes)
timepoints[[1]] <- timepoint1

# Iterate through time steps
mapply( function(x) {
  timepoints[[x]] <<- model_run(timepoints[[x-1]], x)
}, 2:runtime)

# Function for an iteration of the model
model_run <-function(list, x) {
  # New E and l values for each cohort-genotype combo
  new_E_l <- get_pre_reproductive_size(list[[1]], list[[2]])
  # New population sizes for each cohort-genotype combo, after mortality in prior to breeding
  new_mu_exp <- mortality(list[[3]], new_E_l[[1]], new_E_l[[2]])
  # Add oldest and second-oldest cohorts
  new_mu_exp[nrow(new_mu_exp)-1,] <- new_mu_exp[nrow(new_mu_exp),] + new_mu_exp[nrow(new_mu_exp)-1,]
  # Empty out first row for young of year, for the mu_exp, E and l matrices
  new_mu_exp <- advance_age(new_mu_exp)
  new_E_l[[1]] <- advance_age(new_E_l[[1]])
  new_E_l[[2]] <- advance_age(new_E_l[[2]])
  # Empty out first row for young of year, for the mu matrix
  mu <- list[[4]]
  mu[2:nrow(mu),,] <- mu[1:(nrow(mu)-1),,]
  mu[1,,] <- 0
  # Reproduction is seasonal
  if (x%%round(365/timescale) == 0) {
    # Get the amount of eggs produced and the amounts produced by each cohort-genotype combo, and a new E after reversible energy is depleted by mating
    fecundity_object <- fecundity(new_E_l[[1]], new_E_l[[2]], new_mu_exp)
    new_E_l[[1]] <- fecundity_object[[3]]
    # Update first row of mu_exp with new individuals with expression based on past genotype
    inheritance_object <- inheritance(fecundity_object[[2]], mu, new_mu_exp, Ricker(fecundity_object[[1]]))
    # Updating mu_exp matrix
    new_mu_exp <- inheritance_object[[4]]
    # Give the YOY their sizes and lengths, if there are any
    if (sum(new_mu_exp[1,]) > 0) {
      new_E_l[[1]][1,] <- 3.938
      new_E_l[[2]][1,] <- 0.1102
    }
    # Updating mu matrix
    expressed <- as.vector(inheritance_object[[2]])
    genotype <- as.vector(inheritance_object[[3]])
    mu <- inheritance_object[[1]]
    new_cohort_genotypes<-mapply(function(k) scroll_through_phenotypes(k, expressed, genotype), seq(1, genotypes))
    mu[1,,] <- new_cohort_genotypes
  }
  # Return updated matrices for this time point (2 through runtime)
  return(list(new_E_l[[1]], new_E_l[[2]], new_mu_exp, mu))
}

# Function for advancing age
advance_age <- function(matrix) {
  matrix[2:(nrow(matrix)),] <- matrix[1:(nrow(matrix)-1),]
  matrix[1,] <- 0
  return(matrix)
}

### GROWTH FUNCTIONS

# For obtaining structural mass
structural_mass <- function(l) {
  return(c_1*(l)^c_2)
}

# Calculate the intrinsic length-dependent E/S ratio that the fish is trying to attain
intrinsic_lambda <- function(l, l_bar, r, lambda_min, lambda_max) {
  return(lambda_min + (lambda_max-lambda_min)*(exp(r*(l-l_bar)))/(1+exp(r*(l-l_bar))))
}

# Energy intake for the fish per year, dependent on fish structural biomass
size_dependent_energy_intake <- function(S, p_0, p_1) {
  return(resource*p_0*(S)^(p_1))
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

get_pre_reproductive_size <- function(E_past, l_past) {
  S_past = structural_mass(l_past)
  S = S_past
  E = E_past
  lambda_l = intrinsic_lambda(l_past, l_bar, r, lambda_min, lambda_max)
  p_net = size_dependent_energy_intake(S_past, p_0, p_1) - maintenance_cost(S_past, E_past, c_S, c_E)
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
  output <- list(E,l)
  return(output)
}

# Fishing mortality is size selective
selective_mortality <- function(length) {
  1 / (1 + exp(-q_f*(length - mincatchsize)*100))
}

# Mortality drops as E/S increases, as S increases
mortality <- function(mu_exp, E, l) {
  lambda <- (E)/structural_mass(l)
  lambda[is.na(lambda)] <- 0
  mu_exp <- mu_exp*exp( -m + Fishing*selective_mortality(l) - m_c_max*exp(-z_c*lambda) - (m_p_max-m)*exp(-z_p*l))
  return(mu_exp)
}

fecundity <- function(E, l, mu_exp) {
  # Calculate potential reproductive contribution in terms of biomass, with cost to reproduce substracted
  repro_biomass<-((w*E-(r_0*structural_mass(l))^(r_1)))*mu_exp
  repro_biomass_individual<-((w*E-(r_0*structural_mass(l))^(r_1)))
  # Zero out negative reproductive biomass - these individuals don't spawn
  repro_biomass[repro_biomass < 0] <- 0 
  repro_biomass_individual[repro_biomass_individual < 0] <- 0 
  E <- E - repro_biomass_individual
  total_fecundity <- sum(repro_biomass)
  # Convert from g. to kg.
  total_fecundity = total_fecundity/1000
  # Convert from kg. to fecundity
  total_fecundity = total_fecundity*(packing)
  return(list(total_fecundity, repro_biomass, E))
}

biomass_t <- function(E, l, mu_exp) {
  masses <- E + structural_mass(l)
  biomass <- sum(mu_exp*masses)
  return(biomass)
}

Ricker <- function(egg_production_total) {
  # Ricker for YOY
  N_0 = alpha * egg_production_total * exp(-b*egg_production_total)
  # Converting from thousands recruits to recruits, dividing by 2 because Ricker makes both male and female fish
  N_0 = N_0 * 1000
  N_0 = N_0/2
  return(N_0)
}

inheritance <- function(repro_biomass_cohorts, mu, mu_exp, N_0) {
  multiplier<-array(c(rep(repro_biomass_cohorts,ncol(repro_biomass_cohorts))), dim=c(dim(repro_biomass_cohorts),ncol(repro_biomass_cohorts)))
  genotype_contributions <- mu*multiplier
  genotype_contributions_vec <- apply(genotype_contributions, 3, sum)
  genotype_contributions_vec <- genotype_contributions_vec/sum(genotype_contributions_vec)
  genotype_contributions_vec[is.na(genotype_contributions_vec)] <- 0
  genotype_contributions_vec <- genotype_contributions_vec*N_0
  genotype_contributions_vec <- round(genotype_contributions_vec/10)
  expression_list <- mapply(function(z,y) round_special(rtruncnorm(y,1,ncol(repro_biomass_cohorts),mean=z,sd=1)), seq(1, ncol(repro_biomass_cohorts)), genotype_contributions_vec)
  genotype_list <- mapply(function(z,y) rep(z,y), seq(1, ncol(repro_biomass_cohorts)), genotype_contributions_vec)
  expression_vec <- tabulate(as.vector(expression_list))*10
  mu_exp[1,] <- expression_vec
  return(list(mu, expression_list, genotype_list, mu_exp))
}

scroll_through_phenotypes <- function(k, expressed, genotype) {
  row<-mapply(function(y) length(expressed[which(expressed==y & genotype==k)])/length(expressed[which(expressed==y)]), seq(1,genotypes))
  row[is.na(row)] <- 0
  return(row)
}

round_special <- function(value) {
  if (is.null(value)) {
    return(0)
  }
  if (!(is.null(value))) {
    return(round(value))
  }
}




vector <- mode("numeric", len)






### HERRING MODEL

# Set longevity of the fish species, runtime of the model, and number of distinct genotypes
longevity_days = 9125
# set runtime to any arbitrary desired value, in units of days, that is a multiple of 365
runtime_days = 150*365
genotypes = 2
# Setting the length of the timestep in days - set to any fraction of 365
timescale= (365)/4
# Setting run time and longevity to units of time steps
runtime = round(runtime_days/timescale)
longevity = round(longevity_days/timescale)
# How long to run model for convergence before addition of genetic diversity
# in units of days
add_genotypes = 33440

# Set parameters
# Need to add density dependence
# The steepness of the increase in the ratio lambda (reversible over structural mass)
r = 6
# The length at which the increasing function of lambda reaches its inflection point
l_bar = 0.2
# Scaling parameter for obtaining structural mass from standard length
c_1 = 5735
c_2 = 3.125
# Multiplier for obtaining overall yearly energy intake as a function of structural mass
p_0 = 0.1*timescale
# Exponent for obtaining overall yearly energy intake as a function of structural mass
p_1 = 2/3
# Maximum ratio lambda (reversible over structural mass) a fish can attain over its life
lambda_max = 1.3
# Minimum ratio lambda a fish can attain over its life
lambda_min = 0
# Multiplier for obtaining cost of maintenance of structural mass
c_S = 0.03*timescale
# Multiplier for obtaining cost of maintenance of reversible mass
c_E = 0.03*timescale
# Efficiency with which energy is converted into reversible mass
e_E = 0.9
# Efficiency with which energy is converted to structural mass
e_S = (1/3)
# Multiplier and exponent for the cost of reproduction (in grams)
r_0 = 6
r_1 = 0.6
# Effectively, reproductive investment (fraction of reversible mass devoted to reproduction)
w= 0.674
# Instantaneous base mortality rate, derived from daily instantaneous probability of survival for M=-0.2, multiplied to -1 to make positive
m = log(0.999452^timescale)*(-1)
# Instantaneous fishing mortality
Fishing = 0
# Ricker
alpha = 8.44e-09
b = 4.43e-15
# Parameter setting size-selection ogive steepness
q_f = 0.2
# Minimum catch size for the species (in meters)
mincatchsize = 0.5
# Total population biomass
b_t=0
# Increased natural mortality rate at E/S = 0 or S = 0 
# For roughly 0.5 year-based instantaneous mortality at ~1 yrs., Koster et al 2003
m_p_max=log(0.9986^timescale)*(-1)
m_c_max=log(0.95^timescale)*(-1)
# Egg packing constant
packing = 4.45e6
# Parameter for exponential decline in condition-dependent mortality with improved body condition
z_c=7
# Parameter for exponential decline in predation mortality with increased body length
z_p=8
# Resource density
RD = 9000000
# Consumption with resource density half-saturation constant
K_RD = 1000000
# Variability in resource
sd_RD = 0
# How much energy intake is affected by changes in resources (1=no change)
resource = 1
# Additional herring parameters
# Attack rate
AR_1=0.2
AR_2=0.4
AR_3=150
AR_4=0.5
# Digestion
H_1=4.8
H_2=-0.74
# Food conversion efficiency
K_e=0.5
# Optimal size for consuming part of resource
M_opt= 50
# Exponent for maintenance
c_SE_exp = 0.8

### GROWTH FUNCTIONS

# For obtaining structural mass
structural_mass <- function(l) {
  return(c_1*(l)^c_2)
}

# Calculate the intrinsic length-dependent E/S ratio that the fish is trying to attain
intrinsic_lambda <- function(l, l_bar, r, lambda_min, lambda_max) {
  return(lambda_min + (lambda_max-lambda_min)*(exp(r*(l-l_bar)))/(1+exp(r*(l-l_bar))))
}

# Energy intake for the fish per year, dependent on fish structural biomass
size_dependent_energy_intake <- function(S) {
  # Numeric values from Ecological Archives E093-075-A1
  # Attack rate
  AR = AR_1*((1.7*S)^AR_2) + AR_3*((S/M_opt)*exp(1-(S/M_opt)))^(AR_4)
  # Digestion
  H = H_1*((1.7*S)^(H_2))
  return(timescale*((AR*resource)/(AR*resource*H+1))*K_e)
}

# Calculate the cost of the year of maintaining current structural and reversible biomass
maintenance_cost <- function(S, E, c_S, c_E) {
  return(c_S*(S+E)^(c_SE_exp))
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

get_pre_reproductive_size <- function(E_past, l_past) {
  S_past = structural_mass(l_past)
  S = S_past
  E = E_past
  lambda_l = intrinsic_lambda(l_past, l_bar, r, lambda_min, lambda_max)
  p_net = size_dependent_energy_intake(S_past, p_0, p_1) - maintenance_cost(S_past, E_past, c_S, c_E)
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
  output <- list(E,l)
  return(output)
}

# Fishing mortality is size selective
selective_mortality <- function(length) {
  1 / (1 + exp(-q_f*(length - mincatchsize)*100))
}

# Mortality drops as E/S increases, as S increases
mortality <- function(mu_exp, E, l) {
  lambda <- (E)/structural_mass(l)
  mu_exp <- mu_exp*exp( -m + Fishing*selective_mortality(l) - m_c_max*exp(-z_c*lambda) - (m_p_max-m)*exp(-z_p*l))
  return(mu_exp)
}

# Shifting cohort ages forward
mu_exp[nrow(mu_exp),] <-mu_exp[nrow(mu_exp),] + mu_exp[nrow(mu_exp)-1,]
mu_exp <- as.matrix(rbind(mu_exp[1:nrow(mu_exp)-2,], mu_exp[nrow(mu_exp),]))

fecundity <- function(E, l) {
  # Calculate potential reproductive contribution in terms of biomass, with cost to reproduce substracted
  repro_biomass<-((w*E-(r_0*structural_mass(l))^(r_1))>0)
  # Zero out negative reproductive biomass - these individuals don't spawn
  repro_biomass[repro_biomass < 0] <- 0 
  total_fecundity <- sum(repro_biomass)
  relative_contributions <- repro_biomass*(1/total_fecundity)
  # Convert from g. to kg.
  total_fecundity = total_fecundity/1000
  # Convert from kg. to fecundity
  total_fecundity = total_fecundity*(packing)
  return(list(total_fecundity, repro_biomass))
}

biomass_t <- function(E, l, mu_exp) {
  masses <- E + structural_mass(l)
  biomass <- sum(mu_exp*masses)
  return(biomass)
}

Ricker <- function(egg_production_total) {
  # Ricker for YOY
  N_0 = alpha * egg_production_total * exp(-b*egg_production_total)
  # Converting from thousands recruits to recruits, dividing by 2 because Ricker makes both male and female fish
  N_0 = N_0 * 1000
  N_0 = N_0/2
  return(N_0)
}

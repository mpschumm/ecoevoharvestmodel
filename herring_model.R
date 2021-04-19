results_inflection <- vector(mode="numeric",length=10)
results_LAM <- vector(mode="numeric",length=10)
results_AAM <- vector(mode="numeric",length=10)
results_popn <- vector(mode = "list", length = 10)
results_biomass <- vector(mode = "list", length = 10)

for (results_counter in 1:10) {
  
  # Dependent packages
  library(pbapply)
  library(truncnorm)
  
  ### COD MODEL
  
  # run files containing functions for growth, reproduction and mortality processes
  
  source("mortality_functions.R")
  
  source("reproduction_functions.R")
  
  source("growth_functions.R")
  
  # Set longevity of the fish species, runtime of the model, and number of distinct genotypes
  longevity_days = 5475
  # set runtime to any arbitrary desired value, in units of days, that is a multiple of 365
  runtime_days =  400*91.25 # 150*365 # 20*365
  genotypes = 20
  # Setting the length of the timestep in days - set to any fraction of 365
  timescale= (365)/4
  # Setting run time and longevity to units of time steps
  runtime = round(runtime_days/timescale)
  longevity = round(longevity_days/timescale)
  # Whether fish breeds once a year (has a breeding season), or breeds continously throughout the year except for some period of time
  breeding_season = TRUE
  # How long to run model for convergence before addition of genetic diversity
  # in units of days
  add_genotypes = 21900
  
  # Set parameters
  # The steepness of the increase in the ratio lambda (reversible over structural mass)
  r = 6
  # The length at which the increasing function of lambda reaches its inflection point
  l_bar = 0
  # Scaling parameter for obtaining structural mass from standard length
  c_1 = 5735
  c_2 = 3.125
  # Maximum ratio lambda (reversible over structural mass) a fish can attain over its life
  lambda_max = 1.13
  # Minimum ratio lambda a fish can attain over its life
  lambda_min = 0
  # Multiplier for obtaining cost of maintenance of structural mass
  c_S = 0.03 # *timescale
  # Multiplier for obtaining cost of maintenance of reversible mass
  c_E = 0.03 # *timescale
  # Efficiency with which energy is converted into reversible mass
  e_E = 0.5
  # Efficiency with which energy is converted to structural mass
  e_S = 0.5
  # Multiplier and exponent for the cost of reproduction (in grams)
  r_0 = 0.7
  r_1 = 1
  # Effectively, reproductive investment (fraction of reversible mass devoted to reproduction)
  w= 0.674
  # Instantaneous base mortality rate, derived from daily instantaneous probability of survival for M=-0.2, multiplied to -1 to make positive
  m = log(0.999452)*(-1)
  # Instantaneous fishing mortality
  Fishing = 0 
  # Ricker
  alpha = 7.5e-4
  b = 2.84e-6
  # Parameter setting size-selection ogive steepness
  q_f = 0.2
  # Minimum catch size for the species (in meters)
  mincatchsize = 0.5
  # Total population biomass
  b_t=0
  # Increased natural mortality rate at E/S = 0 or S = 0 
  # For roughly 0.5 year-based instantaneous mortality at ~1 yrs., Koster et al 2003
  m_p_max=log(0.9986)*(-1)
  m_c_max=log(0.95)*(-1)
  # Egg packing constant
  packing = 3333333
  # Parameter for exponential decline in condition-dependent mortality with improved body condition
  z_c=7
  # Parameter for exponential decline in predation mortality with increased body length
  z_p=8
  # Amount of resource (initial)
  resource <- 1e+11
  # Resource carrying capacity
  resource_K <- 1e+11
  # Resource intrinsic growth rate
  resource_r <- 2.4
  # Resource variability
  r_SD = 20
  # Consumer functional response half-saturation constant
  K_half <- 1e+10
  # Additional herring parameters
  # Attack rate
  AR_1=150
  AR_2=1.7
  AR_3=0.5
  AR_4 = 25
  # Digestion
  H_1=4.8
  H_2=-0.74
  # Food conversion efficiency
  K_e=0.5
  # Optimal size for consuming part of resource
  M_opt= 50
  
  # Initializing the matrices for the first time point
  timepoints <- vector(mode = "list", length = runtime)
  # Each timepoint consists of a matrix known as E of reversible masses, organized in rows for ages and columns for genotypes; a matrix of the same organization, l, for lengths; a matrix of the same organization for numbers of individuals, mu_exp (for "expressed"); and an array where each z-level is a genotype and the cells contain the proportion of each age-phenotype combo that has that genotype
  timepoint1 <- list(matrix(data=0, nrow=longevity, ncol=genotypes), 
                     matrix(data=0, nrow=longevity, ncol=genotypes),
                     matrix(data=0, nrow=longevity, ncol=genotypes),
                     array(data=0, dim=c(longevity, genotypes, genotypes)))
  # YOY individuals sizes and lengths - the model begins with a set number of YOY individuals
  timepoint1[[1]][1,] <- 1.4
  timepoint1[[2]][1,] <- 0.057
  # Setting initial number of individuals
  timepoint1[[3]][1,] <- 100
  # Genotypes initially correspond perfectly to phenotypes
  timepoint1[[4]][1,,] <- diag(genotypes)
  timepoints[[1]] <- timepoint1
  # Initialize resources vector
  resources <- vector(mode="numeric", length=runtime)
  # Initialize resource consumption vector
  resource_consumption <- vector(mode="numeric", length=runtime)
  
  # Iterate through time steps
  pbmapply( function(x) {
    timepoints[[x]] <<- model_run(timepoints[[x-1]], x)
  }, 2:runtime)
  
  phenotype_averages <- matrix(nrow=3, ncol=genotypes)
  
  for (tracker_counter in 1:genotypes) {
    mature<-(w*as.vector(timepoints[[runtime]][[1]][,tracker_counter])-as.vector((r_0*(structural_mass(timepoints[[runtime]][[2]][,tracker_counter]))^(r_1))))
    phenotype_averages[1,tracker_counter] <- which(mature==min(mature[mature>0]))
    phenotype_averages[2,tracker_counter] <- timepoints[[runtime]][[2]][phenotype_averages[1,tracker_counter],tracker_counter]
    phenotype_averages[3,tracker_counter] <- timepoints[[runtime]][[3]][phenotype_averages[1,tracker_counter],tracker_counter]
  }
  
  phenotype_averages[3,] <- phenotype_averages[3,]/sum(phenotype_averages[3,])
  
  results_inflection[results_counter] <- sum(c(1:genotypes)*(timepoints[[runtime]][[3]][1,]/sum(timepoints[[runtime]][[3]][1,])))
  results_AAM[results_counter] <- sum(phenotype_averages[3,]*phenotype_averages[1,])
  results_LAM[results_counter] <- sum(phenotype_averages[3,]*phenotype_averages[2,])
  
  results_popn[[results_counter]]<-mapply(function(x) sum(timepoints[[x]][[3]]), seq(40,400, by=4))
  results_biomass[[results_counter]]<-mapply(function(x) sum((timepoints[[x]][[1]] + timepoints[[x]][[2]])*timepoints[[x]][[3]] ), seq(40,400, by=4))
  
}

# Function for an iteration of the model
model_run <-function(list, x) {
  # Adding fishing after a certain time point
  if (x==round(add_genotypes/timescale)) {
    Fishing <<- Fishing_add
  }
  # Adding differences between genotypes
  if (x==round(add_genotypes/timescale)) {
    r <<- seq(0, 12, by = (12-0)/(genotypes-1))
    r<<-matrix(rep(r,each=longevity),nrow=longevity)
  }
  new_E_l <<- list(list[[1]], list[[2]])
  new_mu_exp <<- list[[3]]
  # Generate vectors of random values, with some degree of temporal correlation phi
  times <- 1:round(timescale)
  mat <- as.matrix(dist(times, diag = T, upper= T))
  phi = 0.11
  AR_mat <- as.matrix(phi^mat)
  L <- t(chol(AR_mat))
  ar_var <- L%*%matrix(rnorm(round(timescale)*round(timescale),0,r_SD), ncol= round(timescale))*sqrt(1-phi^2)
  random_values <- ar_var[,1]
  mapply ( function(growing) {
    # New E and l values for each cohort-genotype combo
    new_E_l <<- get_pre_reproductive_size(new_E_l[[1]], new_E_l[[2]], resource, new_mu_exp)
    # Compute a random value for calculating r
    random_value <- random_values[growing]
    resource <<- resource + resource*((resource_K-resource)/resource_K)*resource_r*(exp(random_value)/(1+exp(random_value))) - new_E_l[[3]]
    # New population sizes for each cohort-genotype combo, after mortality in prior to breeding
    new_mu_exp <<- mortality(new_mu_exp, new_E_l[[1]], new_E_l[[2]])
  }, 1:round(timescale))
  # Record resource level
  resources[x] <<- resource
  # Record consumption level
  # Resource consumption is only for last day within period of days in 'timescale' (i.e., last day of season if timescale = 365/4)
  resource_consumption[x] <<- new_E_l[[3]]
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
  if ((breeding_season==T) & x%%round(365/timescale ) == 0 & x!=runtime) {
    # Get the amount of eggs produced and the amounts produced by each cohort-genotype combo, and a new E after reversible energy is depleted by mating
    fecundity_object <- fecundity(new_E_l[[1]], new_E_l[[2]], new_mu_exp)
    new_E_l[[1]] <- fecundity_object[[3]]
    # Update first row of mu_exp with new individuals with expression based on past genotype
    inheritance_object <- inheritance(fecundity_object[[2]], mu, new_mu_exp, Ricker(fecundity_object[[1]]))
    # Updating mu_exp matrix
    new_mu_exp <- inheritance_object[[4]]
    # Give the YOY their sizes and lengths, if there are any
    if (sum(new_mu_exp[1,]) > 0) {
      new_E_l[[1]][1,] <- 1.4
      new_E_l[[2]][1,] <- 0.057
    }
    # Create vectors with length=number of new individuals, where each value is the individual's expressed phenotype (for "expressed") and their genotype (for "genotype")
    mu <- inheritance_object[[1]]
    expressed <- inheritance_object[[2]]
    genotype <- inheritance_object[[3]]
    # length(expressed) should == length(genotype)
    # Filling in the genotypes of the new recruits
    # The scroll_through_genotypes function will create a row for each z-level in the new mu 3D array, where the row has the proportions of new recruits from each phenotype (column) that are in the genotype that this z-level of the array corresponds to
    # See "reproduction_function" file for more info in comments
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








# ^
# v

# Huss et al. 2012 herring model
# http://esapubs.org/archive///////ecol/E093/075/appendix-A.htm
# file:///Users/mschumm/O_lab/ecoevoharvestmodel/non_code/Huss_etal_AppendixA.htm

# key "parameteriz"

# genotype1r = 4
# genotype2r = 8
# genotype1l_bar = 0.2
# genotype2l_bar = 0.2
# genotype1w = 0.674
# genotype2w = 0.674

# genotype1r = 6
# genotype2r = 6
# genotype1l_bar = 0.2
# genotype2l_bar = 0.2
# genotype1w = 0.7
# genotype2w = 0.3

genotype1r = 6
genotype2r = 6
genotype1l_bar = 0.2
genotype2l_bar = 0
genotype1w = 0.674
genotype2w = 0.674

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

initialize_pop_matrices <- function(longevity, runtime, genotype) {
  # Initialize matrix of columns (years) and rows (age classes)
  # Multiplied by 3 because we'll be storing cohorts' reproductive biomass and lengths 
  age_dist=matrix(0L, nrow=longevity*3, ncol=runtime+1)
  # Add (arbitrarily) 1000000 individuals to age=10 to start with
  age_dist[round(3650/timescale), 1] = 1000000
  return(age_dist)
}

# Generate a matrix for each genotype
for (g in 1:genotypes) {
  assign(paste("age_dist_",g,sep=""),initialize_pop_matrices(longevity, runtime, g),envir = .GlobalEnv)
}

# Set parameters
# Need to add density dependence
# The steepness of the increase in the ratio lambda (reversible over structural mass)
r = 6
# The length at which the increasing function of lambda reaches its inflection point
l_bar = 0.2
# Scaling parameter for obtaining structural mass from standard length
c_1 = 5735
# Where is c2?
# Make an alternate of each of the 4 function files for herring
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
# Is maintenance really done?
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
AR_1=150
AR_2=1.7
AR_3=0.5
AR_4 = 25
# Digestion
H_1=4.8
H_2=-0.74
# Food conversion efficiency
K_e=0.5
# Optimal size for consuming part of resource
M_opt= 50
# Exponent for maintenance
c_SE_exp = 0.8

# GROWTH FUNCTIONS

# For obtaining structural mass
structural_mass <- function(l) {
  return(c_1*(l)^3.125)
}

# Calculate the instrinsic length-dependent E/S ratio that the fish is trying to attain
intrinsic_lambda <- function(l, l_bar, r, lambda_min, lambda_max) {
  return(lambda_min + (lambda_max-lambda_min)*(exp(r*(l-l_bar)))/(1+exp(r*(l-l_bar))))
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

# Calculate the cost of the year of maintaining current structural and reversible biomass
maintenance_cost <- function(S, E, c_S, c_E) {
  return(c_S*(S+E)^(c_SE_exp))
}

# Function for calculating a new, reduced level of reversible biomass, if net energy intake is negative
negative_p_net <- function(E, p_net, e_E) {
  if (p_net < 0) {
    return(E + (p_net/e_E))}
}

# Calculate the amount of energy needed to raise the R/S ratio to the length-dependent maximum
p_E_raising_lambda <- function(lambda_a, S, E) {
  return((lambda_a*S - E))
}

# Calculating proportion of any remaining energy to be converted to structural mass, to maintain R/S
remaining_energy_allocation <- function(e_S, e_E, lambda_a) {
  return((lambda_a*e_S)/(lambda_a*e_S + e_E))
}

# Obtain R and S for a cohort prior to its reproduction (or failure to reproduce) that year
get_pre_reproductive_size <- function(E_past, l_past, a) {
  S_past = structural_mass(l_past)
  S = S_past
  lambda_a = intrinsic_lambda(l_past, l_bar, r, lambda_min, lambda_max)
  p_net = size_dependent_energy_intake(S_past, p_0, p_1) - maintenance_cost(S_past, E_past, c_S, c_E)
  if (p_net < 0) {
    E <- max(0,negative_p_net(E_past, p_net, e_E))
  }
  if (p_net >= 0 & lambda_a > (E_past/S_past)) {
    # Gives amount of grams needed to raise E/S to lambda_l - not taking efficiency into account yet
    p_E = p_E_raising_lambda(lambda_a, S_past, E_past)
    if ((p_E/e_E - p_net) >=0) {
      E_past = p_net*e_E + E_past
      p_net = 0
    }
    if ((p_E/e_E - p_net) <0) {
      E_past = p_E + E_past
      p_net = p_net - p_E/e_E
    }
  }
  E = E_past
  if (p_net >= 0) {
    E = E_past + remaining_energy_allocation(e_S, e_E, lambda_a)*e_E*p_net
    S = S_past + (1-remaining_energy_allocation(e_S, e_E, lambda_a))*e_S*p_net
  }
  l = (S/c_1)^(1/3.125)
  output <- c(E,l)
  return(output)
}

# REPRODUCTION FUNCTIONS

# Calculate spawning stock biomass
fecundity <- function(ages, t, age_dist, g) {
  # Parameterizing genotypes
  if (t >= round(add_genotypes/timescale)) {
    w=0.674
    # Genotype parameterization
    if (g==1) {
      w=genotype1w
    }
    if (g==2) {
      w=genotype2w
    }
  }
  fecundity = 0
  for (a in ages) {
    if ((w*age_dist[a+(longevity),t]-(r_0*structural_mass(age_dist[a+(2*longevity),t]))^(r_1))>0) {
      fecundity = fecundity + (w*age_dist[a+(longevity),t]-(r_0*structural_mass(age_dist[a+(2*longevity),t]))^(r_1))*age_dist[a,t]
      age_dist[a+(longevity),t] = age_dist[a+(longevity),t] - w*age_dist[a+(longevity),t]
    }
  }
  # Convert from g. to kg.
  fecundity = fecundity/1000
  # Convert from kg. to fecundity
  fecundity = fecundity*(packing)
  assign(paste("age_dist_",g,sep=""),age_dist, envir = .GlobalEnv) 
  return(fecundity)
}

# Calculate total population biomass
biomass_t <- function(ages, t) {
  for (a in ages) {
    b_t = b_t + age_dist[a, t-1]*(structural_mass(age_dist[a+(2*longevity),t-1]))
  }
  return(b_t)
}

# Density dependent recruitment for this year
recruitment <- function(age, ages, t) {
  egg_production_total = 0
  for (g in 1:genotypes) {
    egg_production_total = egg_production_total + fecundity(ages, t, get(paste("age_dist_",g,sep="")),g)
  }
  relative_recruitment = vector(mode="numeric", length=genotypes)
  for (g in 1:genotypes) {
    relative_recruitment[g] = fecundity(ages, t, get(paste("age_dist_",g,sep="")),g)/egg_production_total
  }
  if (age==0) {
    N_0 = alpha * egg_production_total * exp(-b*egg_production_total)
    # Converting from millions recruits to recruits, dividing by 2 because Ricker makes both male and female fish
    N_0 = N_0 * 1000000
    N_0 = N_0/2
  }
  for (g in 1:genotypes) {
    age_dist = get(paste("age_dist_",g,sep=""))
    age_dist[1,t] = N_0*relative_recruitment[g]
    assign(paste("age_dist_",g,sep=""),age_dist, envir = .GlobalEnv) 
  }
}

# Add new recruits, and change the reproducing cohorts' amount of reversible mass
reproduce <- function(ages, t) {
  if (t>1) {
    recruitment(0, ages, t)
  }
}

# MORTALITY FUNCTIONS

# Fishing mortality is size selective
selective_mortality <- function(length) {
  1 / (1 + exp(-q_f*(length - mincatchsize)*100))
}

# Mortality drops as E/S increases, as S increases
mortality <- function(age) {
  cond <- 100*(((age_dist[age-1+(longevity),t-1]))+structural_mass(age_dist[age-1+(2*longevity),t-1]))/((100*age_dist[age-1+(2*longevity),t-1])^3)
  length = age_dist[age-1+(2*longevity),t-1]
  if (age < longevity) {
    N_t = age_dist[age-1,t-1]*exp( -m + -Fishing*selective_mortality(length) - m_c_max*exp(-z_c*cond) - (m_p_max-m)*exp(-z_p*length) )
  }
  if (age == longevity) {
    N_t = age_dist[age-1,t-1]*exp( -m + -Fishing*selective_mortality(length) - m_c_max*exp(-z_c*cond) - (m_p_max-m)*exp(-z_p*length))
    N_t = age_dist[age,t-1]*exp( -m + -Fishing*selective_mortality(length) - m_c_max*exp(-z_c*cond) - (m_p_max-m)*exp(-z_p*length))
  }
  return(N_t)
}

# Initial masses and lengths for age 0 cod from table 2.2.1 from Ianelli et al., and (for the lambda ratio) from Audzijonyte and Richards
# Calculate lengths (and reversible masses) at each age for the model with no density dependence
lengths_at_t1 <- function(genotypes) {
  t=1
  for (g in 1:genotypes) {
    age_dist = get(paste("age_dist_",g,sep=""))
    age_dist[1+(longevity), 1] = 3.938
    age_dist[1+2*(longevity), 1] = 0.1102
    for (i in 2:longevity) {
      new_E_l = get_pre_reproductive_size(age_dist[i-1+(longevity), 1],age_dist[i-1+2*(longevity), 1], i-1)
      age_dist[i+(longevity), 1] = new_E_l[1]
      age_dist[i+(2*longevity), 1] = new_E_l[2]
      if ((w*age_dist[i+(longevity),t]-(r_0*structural_mass(age_dist[i+(2*longevity),t]))^(r_1))>0)
      {
        age_dist[i+(longevity),t] = age_dist[i+(longevity),t] - w*age_dist[i+(longevity),t]*(timescale/365) 
      }
    }
    assign(paste("age_dist_",g,sep=""),age_dist, envir = .GlobalEnv) 
  }
}  

# Initialize lengths and reversible masses for ages
lengths_at_t1(genotypes)

# Set runtime so is extra time left in end of model, to account for recruitment lag
runtime <- (runtime - (round(365/timescale)))

# Run the model
for (t in 2:runtime) {
  random_multiplier <- max(rnorm(1,1,sd=sd_RD),0)
  resource = (random_multiplier*RD/(K_RD+ random_multiplier*RD))  
  b_t=0
  if (t == round(add_genotypes/timescale)) {
    Fishing = log(exp(-0)^(1/(365/timescale)))*(-1)
    sd_RD = 0
  }
  for (genotype in 1:genotypes) {
    age_dist = get(paste("age_dist_",genotype,sep=""))
    b_t = biomass_t(1:longevity, t)
    assign(paste("age_dist_",genotype,sep=""),age_dist, envir = .GlobalEnv) 
  }
  for (genotype in 1:genotypes) {
    
    # Parameterizing genotypes
    if (t >= round(add_genotypes/timescale)) {
      r=6
      # Genotype parameterization
      if (genotype==1) {
        # a_bar = 36.5
        r=genotype1r
        l_bar = genotype1l_bar
      }
      if (genotype==2) {
        # a_bar = 365*2
        r=genotype2r
        l_bar = genotype2l_bar
      }
    }
    
    age_dist = get(paste("age_dist_",genotype,sep=""))
    age_dist[1+(longevity), t] =  3.938
    age_dist[1+2*(longevity), t] = 0.1102
    for (age in 2:longevity) {
      new_E_l = get_pre_reproductive_size(age_dist[age-1+(longevity), t-1],age_dist[age-1+2*(longevity), t-1], age-1)
      age_dist[age+(longevity), t] = new_E_l[1]
      age_dist[age+(2*longevity), t] = new_E_l[2]
      
      # Mortality
      for (age in 2:longevity) {
        age_dist[age, t] <- mortality(age)
      }
    }
    assign(paste("age_dist_",genotype,sep=""),age_dist, envir = .GlobalEnv) 
  }
  if ((t%%round(365/timescale) == 0)) {
    reproduce(1:longevity, t) 
  }
}

### Analysis code

approach <- vector(mode="numeric", length=runtime)
for (i in 1:runtime) {
  for (j in 1:genotypes) {
    approach[i] <- approach[i]+ sum(get(paste("age_dist_", j, sep=""))[1:longevity,i])
  } 
}
# Plot population dynamics
plot(approach[seq(round(add_genotypes/timescale),length(approach), by=4)],lty=1, pch=19, type="l", ylab="Abundance",xlab="Time (units of 3 month intervals)", ylim=c(0,max(approach)), main= "No harvesting, constant resource")
approach <- vector(mode="numeric", length=runtime)
for (i in 1:runtime) {
  approach[i] <- sum(age_dist_1[1:longevity,i])
}
lines(approach[seq(round(add_genotypes/timescale),length(approach), by=4)], lty=2, pch=19, type="l", col="blue")
cols <- c(3,4,5,6)
for (j in 2:genotypes) {
  approach <- vector(mode="numeric", length=runtime)
  for (i in 1:runtime) {
    approach[i] <- sum(get(paste("age_dist_",j, sep=""))[1:longevity,i])
  }
  lines(approach[seq(round(add_genotypes/timescale),length(approach), by=4)], lty=cols[j-1], pch=19, col="red")
}
legend("topright",legend=c("all", "large maturer","small maturer"), lty=1:3, col=c("black","blue","red"))
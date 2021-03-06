# Dependent packages
library(pbapply)
library(truncnorm)

### COD MODEL

# run files containing functions for growth, reproduction and mortality processes

source("mortality_functions.R")

source("reproduction_functions.R")

source("growth_functions.R")

# Set longevity of the fish species, runtime of the model, and number of distinct genotypes
longevity_days = 9125
# set runtime to any arbitrary desired value, in units of days, that is a multiple of 365
runtime_days = 100*91.25 # 150*365
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
add_genotypes = 33440

# Set parameters
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
# Each timepoint consists of a matrix known as E of reversible masses, organized in rows for ages and columns for genotypes; a matrix of the same organization, l, for lengths; a matrix of the same organization for numbers of individuals, mu_exp (for "expressed"); and an array where each z-level is a genotype and the cells contain the proportion of each age-phenotype combo that has that genotype
timepoint1 <- list(matrix(data=0, nrow=longevity, ncol=genotypes), 
                   matrix(data=0, nrow=longevity, ncol=genotypes),
                   matrix(data=0, nrow=longevity, ncol=genotypes),
                   array(data=0, dim=c(longevity, genotypes, genotypes)))
# YOY individuals sizes and lengths - the model begins with a set number of YOY individuals
timepoint1[[1]][1,] <- 3.938
timepoint1[[2]][1,] <- 0.1102
# Setting initial number of individuals
timepoint1[[3]][1,] <- 100
# Genotypes initially correspond perfectly to phenotypes
timepoint1[[4]][1,,] <- diag(genotypes)
timepoints[[1]] <- timepoint1

# Iterate through time steps
pbmapply( function(x) {
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
    # Create vectors with length=number of new individuals, where each value is the individual's expressed phenotype (for "expressed") and their genotype (for "genotype")
    mu <- inheritance_object[[1]]
    expressed <- inheritance_object[[2]]
    genotype <- inheritance_object[[3]]
    if (length(expressed)!=length(genotype)) {
      print("bad")
    }
    # Filling in the genotypes of the new recruits
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







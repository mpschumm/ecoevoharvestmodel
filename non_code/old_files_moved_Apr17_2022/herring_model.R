results_inflection <- vector(mode="numeric",length=320)
results_LAM <- vector(mode="numeric",length=320)
results_AAM <- vector(mode="numeric",length=320)
results_L2 <- vector(mode="numeric",length=320)
results_Linf <- vector(mode="numeric",length=320)
results_popn <- vector(mode = "list", length = 320)
results_biomass <- vector(mode = "list", length = 320)

randos_vec <- rep(seq(0,42,by=6), each=40,times=1)*0.03
fishing_vec <- rep(2:9,each=5, times=8)*0.1

for (results_counter in 1:20) {
  
  # Dependent packages
  library(pbapply)
  library(truncnorm)
  
  ### HERRING MODEL
  
  # run files containing functions for growth, reproduction and mortality processes
  
  source("herring_mortality_functions.R")
  
  source("herring_reproduction_functions.R")
  
  source("herring_growth_functions.R")
  
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
  r = 3
  # The length at which the increasing function of lambda reaches its inflection point
  l_bar = 0.1
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
  e_E = 1
  # Efficiency with which energy is converted to structural mass
  e_S = 1
  # Multiplier and exponent for the cost of reproduction (in grams)
  # These have a different interpretation for herring
  r_0 = 0.7
  r_1 = 1
  # Effectively, reproductive investment (fraction of reversible mass devoted to reproduction)
  # This has a different interpretation for herring
  w= 0.5
  # Instantaneous base mortality rate, derived from daily instantaneous probability of survival for M=-0.2, multiplied to -1 to make positive
  m = log(0.999452)*(-1)
  # Instantaneous fishing mortality
  Fishing = 0 
  # For what fishing level to add when fishing is introduced partway through the simulation - put yearly instantaneous fishing mortality
  Fishing_add = fishing_vec[results_counter]
  Fishing_add = log((1-Fishing_add)^(1/365))*(-1)
  # Ricker
  # Based on Lynam et al. 2005 MEPS paper
  alpha = 1.16e-1
  b = 8.53e-7
  # Parameter setting size-selection ogive steepness
  q_f = 0.5
  # Minimum catch size for the species (in meters)
  mincatchsize = 0.25
  # Total population biomass
  b_t=0
  # Increased natural mortality rate at E/S = 0 or S = 0 
  # For roughly 0.5 year-based instantaneous mortality at ~1 yrs., Koster et al 2003
  m_p_max=log(0.99)*(-1)
  m_c_max=log(0.95)*(-1)
  # Egg packing constant
  packing = 3333333
  # Parameter for exponential decline in predation mortality with increased body length
  z_p=3.4
  # Amount of resource (initial)
  resource <- 3e+11
  # Resource carrying capacity
  resource_K <- 3e+11
  # Resource intrinsic growth rate
  resource_r <- 1.5
  # Resource variability
  r_SD =  randos_vec[results_counter]
  # Consumer functional response half-saturation constant
  K_half <- 3e+10
  # Additional herring parameters
  # Attack rate
  AR_1=150
  AR_2=1.7
  AR_3=0.5
  # Digestion
  H_1=4.8
  H_2=-0.74
  # Food conversion efficiency
  K_e=0.85
  # Optimal size for consuming part of resource
  M_opt= 327
  
  # Initializing the matrices for the first time point
  timepoints <- vector(mode = "list", length = runtime)
  # Each timepoint consists of a matrix known as E of reversible masses, organized in rows for ages and columns for genotypes; a matrix of the same organization, l, for lengths; a matrix of the same organization for numbers of individuals, mu_exp (for "expressed"); and an array where each z-level is a genotype and the cells contain the proportion of each age-phenotype combo that has that genotype
  timepoint1 <- list(matrix(data=0, nrow=longevity, ncol=genotypes), 
                     matrix(data=0, nrow=longevity, ncol=genotypes),
                     matrix(data=0, nrow=longevity, ncol=genotypes),
                     array(data=0, dim=c(longevity, genotypes, genotypes)))
  # YOY individuals sizes and lengths - the model begins with a set number of YOY individuals
  timepoint1[[1]][1,] <- 0.3712037
  timepoint1[[2]][1,] <- 0.057
  # Setting initial number of individuals
  timepoint1[[3]][1,] <- 10000
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
    mature <- (as.vector(timepoints[[runtime]][[1]][,tracker_counter])-as.vector((r_0*(structural_mass(timepoints[[runtime]][[2]][,tracker_counter]))^(r_1))))
    if (!(all(mature<=0))) {
    phenotype_averages[1,tracker_counter] <- try(which(mature==(mature[mature>0])[1]))
    phenotype_averages[2,tracker_counter] <- timepoints[[runtime]][[2]][phenotype_averages[1,tracker_counter],tracker_counter]
    phenotype_averages[3,tracker_counter] <- timepoints[[runtime]][[3]][phenotype_averages[1,tracker_counter],tracker_counter]
    }
  }
  
  phenotype_averages[is.na(phenotype_averages)] <- 0
  
  phenotype_averages[3,] <- phenotype_averages[3,]/sum(phenotype_averages[3,])
  
  results_inflection[results_counter] <- sum(c(1:genotypes)*(timepoints[[runtime-(365/timescale)]][[3]][1,]/sum(timepoints[[runtime-(365/timescale)]][[3]][1,])))
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
    r <<- seq(0, 5, by = (5-0)/(genotypes-1))
    r<<-matrix(rep(r,each=longevity),nrow=longevity)
  }
  new_E_l <<- list(list[[1]], list[[2]])
  new_mu_exp <<- list[[3]]
  # Generate vectors of random values, with some degree of temporal correlation phi
  times <- 1:round(timescale)
  mat <- as.matrix(dist(times, diag = T, upper= T))
  phi = 0.8
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
      new_E_l[[1]][1,] <- 0.3712037
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


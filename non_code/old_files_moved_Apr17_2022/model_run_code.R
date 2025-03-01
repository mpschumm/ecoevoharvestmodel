  results_inflection <- vector(mode="numeric",length=320)
  results_LAM <- vector(mode="numeric",length=320)
  results_AAM <- vector(mode="numeric",length=320)
  results_L2 <- vector(mode="numeric",length=320)
  results_Linf <- vector(mode="numeric",length=320)
  results_popn <- vector(mode = "list", length = 320)
  results_biomass <- vector(mode = "list", length = 320)
  
  randos_vec <- rep(seq(0,42,by=6), each=40,times=1)*0.03
  fishing_vec <- rep(2:9,each=5, times=8)*0.1
  
  for (results_counter in 1:320) {
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
    runtime_days =  100*365 # 150*365 # 20*365
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
    add_genotypes = 0.6*runtime_days
    
    # Set parameters
    # The steepness of the increase in the ratio lambda (reversible over structural mass)
    r = 6
    # The length at which the increasing function of lambda reaches its inflection point
    l_bar = 0.2
    # Scaling parameter for obtaining structural mass from standard length
    c_1 = 5787
    c_2 = 3
    # Multiplier for obtaining overall yearly energy intake as a function of structural mass
    p_0 = 0.1 
    # Exponent for obtaining overall yearly energy intake as a function of structural mass
    p_1 = 2/3
    # Maximum ratio lambda (reversible over structural mass) a fish can attain over its life
    lambda_max = 1.3
    # Minimum ratio lambda a fish can attain over its life
    lambda_min = 0
    # Multiplier for obtaining cost of maintenance of structural mass
    c_S = 0.003 
    # Multiplier for obtaining cost of maintenance of reversible mass
    c_E = 0.0003 
    # Efficiency with which energy is converted into reversible mass
    e_E = 0.9
    # Efficiency with which energy is converted to structural mass
    e_S = (1/3)
    # Multiplier and exponent for the cost of reproduction (in grams)
    r_0 = 6
    r_1 = 0.6
    # Effectively, reproductive investment (fraction of reversible mass devoted to reproduction)
    w= 0.6
    # Log of instantaneous base mortality rate, derived from daily instantaneous probability of survival for M=-0.2, multiplied to -1 to make positive
    m = log(0.999452)*(-1)
    # Instantaneous fishing mortality
    Fishing = 0
    # For what fishing level to add when fishing is introduced partway through the simulation - put yearly instantaneous fishing mortality
    Fishing_add = fishing_vec[results_counter]
    Fishing_add = log((1-Fishing_add)^(1/365))*(-1)
    # Ricker
    alpha = 8.44e-09
    b = 4.43e-15
    # Parameter setting size-selection ogive steepness
    q_f = 0.5
    # Minimum catch size for the species (in meters)
    mincatchsize = 0.5
    # Increased natural mortality rate at E/S = 0 or S = 0 
    # For roughly 0.5 year-based instantaneous mortality at ~1 yrs., Koster et al 2003
    m_p_max=log(0.9986)*(-1)
    m_c_max=log(0.95)*(-1)
    # Egg packing constant
    packing = 4.45e6
    # Parameter for exponential decline in condition-dependent mortality with improved body condition
    z_c=7
    # Parameter for exponential decline in predation mortality with increased body length
    z_p=8
    # Amount of resource (initial), set to reasonable starting values by calculating fish population resource consumption
    resource <- 3e+11
    # Resource carrying capacity
    resource_K <- 3e+11
    # Resource intrinsic growth rate
    resource_r <- 1.5
    # Resource variability
    r_SD = randos_vec[results_counter]
    # Consumer functional response half-saturation constant
    K_half <- 3e+10
    
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
    # Initialize resources vector
    resources <- vector(mode="numeric", length=runtime)
    # Initialize resource consumption vector
    resource_consumption <- vector(mode="numeric", length=runtime)
    
    # Iterate through time steps
    pbmapply( function(x) {
      timepoints[[x]] <<- model_run(timepoints[[x-1]], x)
    }, 2:runtime)
    
    phenotype_averages <- matrix(nrow=5, ncol=genotypes)
    
    for (tracker_counter in 1:genotypes) {
      mature<-(w*as.vector(timepoints[[runtime]][[1]][,tracker_counter])-as.vector((r_0*(structural_mass(timepoints[[runtime]][[2]][,tracker_counter]))^(r_1))))
      phenotype_averages[1,tracker_counter] <- which(mature==min(mature[mature>0]))
      phenotype_averages[2,tracker_counter] <- timepoints[[runtime]][[2]][phenotype_averages[1,tracker_counter],tracker_counter]
      phenotype_averages[3,tracker_counter] <- timepoints[[runtime]][[3]][phenotype_averages[1,tracker_counter],tracker_counter]
      phenotype_averages[4,tracker_counter] <- timepoints[[runtime]][[2]][9,tracker_counter]
      phenotype_averages[5,tracker_counter] <- timepoints[[runtime]][[2]][97,tracker_counter]
    }
    
    phenotype_averages[3,] <- phenotype_averages[3,]/sum(phenotype_averages[3,])
    
    results_inflection[results_counter] <- sum(c(1:genotypes)*(timepoints[[runtime-(365/timescale)]][[3]][1,]/sum(timepoints[[runtime-(365/timescale)]][[3]][1,])))
    results_AAM[results_counter] <- sum(phenotype_averages[3,]*phenotype_averages[1,])
    results_LAM[results_counter] <- sum(phenotype_averages[3,]*phenotype_averages[2,])
    results_L2[results_counter] <- sum(phenotype_averages[3,]*phenotype_averages[4,])
    results_Linf[results_counter] <- sum(phenotype_averages[3,]*phenotype_averages[5,])
    
    results_popn[[results_counter]]<-mapply(function(x) sum(timepoints[[x]][[3]]), seq((365/timescale)*10,runtime, by=365/timescale))
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
  phi = 0.8
  AR_mat <- as.matrix(phi^mat)
  L <- t(chol(AR_mat))
  ar_var <- L%*%matrix(rnorm(round(timescale)*round(timescale),0,1), ncol= round(timescale))*r_SD*sqrt(1-phi^2)
  random_values <- as.vector(ar_var[,1])
  mapply ( function(growing) {
    # New E and l values for each cohort-genotype combo
    new_E_l <<- get_pre_reproductive_size(new_E_l[[1]], new_E_l[[2]], resource, new_mu_exp)
    # Compute a random value for calculating r
    random_value <- random_values[growing]
    random_resource_K <- resource_K*(exp(random_value - (r_SD^2)/2))
    net_resource <- resource + resource*((random_resource_K-resource)/random_resource_K)*resource_r - new_E_l[[3]]
    resource <<- ifelse(net_resource>0, net_resource, 1)
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
      new_E_l[[1]][1,] <- 3.938
      new_E_l[[2]][1,] <- 0.1102
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





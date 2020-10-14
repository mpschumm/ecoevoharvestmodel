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
    # Converting from thousands recruits to recruits, dividing by 2 because Ricker makes both male and female fish
    N_0 = N_0 * 1000
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
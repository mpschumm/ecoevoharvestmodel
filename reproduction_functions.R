# REPRODUCTION FUNCTIONS

fecundity <- function(E, l, mu_exp) {
  # Calculate potential reproductive contribution in terms of biomass, with cost to reproduce substracted
  repro_biomass<-(w*E-r_0*structural_mass(l)^(r_1))*mu_exp
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
  # The matrix of fecundities for each cohort-phenotype combo is expanded into an array to multiply by the proportions in the mu array (where each z-level is a genotype, and the proportions are the proportion of individuals in that cohort-pheno combo that have that genotype)
  multiplier<-array(c(rep(repro_biomass_cohorts,ncol(repro_biomass_cohorts))), dim=c(dim(repro_biomass_cohorts),ncol(repro_biomass_cohorts)))
  genotype_contributions <- mu*multiplier
  # Produce a vector of total fecundity from each genotype
  genotype_contributions_vec <- apply(genotype_contributions, 3, sum)
  # Convert back to proportions
  genotype_contributions_vec <- genotype_contributions_vec/sum(genotype_contributions_vec)
  genotype_contributions_vec[is.na(genotype_contributions_vec)] <- 0
  # Multiply proportions by total individual offspring (or, half the number of total genotype copies being passed down)
  genotype_contributions_vec <- genotype_contributions_vec*N_0
  # Collapse into super-individuals when N_0 higher than 1000
  divider <- 1
  if (floor(log10(N_0)) > 2) {
    divider <- 10^(floor(log10(N_0)) - 2)
  }
  genotype_contributions_vec <- round(genotype_contributions_vec/divider)
  genotype_list <- mapply(function(z,y) rep(z,y), seq(1, ncol(repro_biomass_cohorts)), genotype_contributions_vec)
  genotype_list <- as.vector(unlist(genotype_list))
  SR_vector <- sample(genotype_list, length(genotype_list), replace=F)
  # Each genotype is paired with another, and averages are calculate to produce the offspring genotypes
  # Produce sequence of numbers corresponding to the genotypes of each of the super individuals
  genotype_list_avg <- round((genotype_list + SR_vector)/2)
  # Frequency distribution of genotypes
  genotype_contributions_vec_avg <- tabulate(genotype_list, nbins=genotypes)
  # Produce sequence of numbers corresponding to the phenotypes of each of the super individuals
  expression_list <- mapply(function(z,y) round_special(rtruncnorm(y,1,ncol(repro_biomass_cohorts),mean=z,sd=1)), seq(1, ncol(repro_biomass_cohorts)), genotype_contributions_vec_avg)
  expression_list <- as.vector(unlist(expression_list))
  expression_list <- expression_list[which(expression_list>0)]
  expression_vec <- tabulate(expression_list, nbins=genotypes)*divider
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


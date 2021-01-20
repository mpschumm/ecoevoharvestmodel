# REPRODUCTION FUNCTIONS

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
  # Could add sexual reproduction here
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
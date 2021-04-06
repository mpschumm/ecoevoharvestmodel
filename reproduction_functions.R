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
  # The 2D matrix (of fecundities of each cohort-phenotype combo; cohorts are rows, columns are phenotype) is expanded into an 3D array where each z-level is just a duplicate of the 2D matrix of fecundities, in order to then element-wise-multiply each cell by the corresponding proportions in the mu 3D array (in which each z-level is a genotype, and a cell is the proportion of individuals in that cohort-pheno combo that have that genotype)
  # The end result of the multiplication is a 3D array where each z level contains the number of new individuals spawned by each cohort-phenotype combo that have that z-level's corresponding genotype
  multiplier<-array(c(rep(repro_biomass_cohorts,ncol(repro_biomass_cohorts))), dim=c(dim(repro_biomass_cohorts),ncol(repro_biomass_cohorts)))
  genotype_contributions <- mu*multiplier
  # Produce a vector of total fecundity from each genotype, by summing on each z-level of the 3D array (i.e., each genotype)
  genotype_contributions_vec <- apply(genotype_contributions, 3, sum)
  # Convert back to proportions (proportions of total fecundity that are coming from parents of each genotype)
  genotype_contributions_vec <- genotype_contributions_vec/sum(genotype_contributions_vec)
  genotype_contributions_vec[is.na(genotype_contributions_vec)] <- 0
  # Multiply proportions by total individual offspring (or, half the number of total genotype copies being passed down)
  genotype_contributions_vec <- genotype_contributions_vec*N_0
  # Collapse into individuals into super-individuals when N_0 higher than 1000 (i.e., 200 offspring from genotype 1 are treated as 20 individuals when there ar between 1000-10,000 offspring total, or just 2 if there are 10,000-100,000, etc.)
  divider <- 1
  if (floor(log10(N_0)) > 2) {
    divider <- 10^(floor(log10(N_0)) - 2)
  }
  # Now that individuals have been collapsed, round any fractional fecundities to integers
  genotype_contributions_vec <- round(genotype_contributions_vec/divider)
  # Make a vector of the numerical values of each (super)individual's genotype, in numerical order (i.e., the number of elements is the same as the number of super individuals and each super individual's numerical genotype is an element)
  genotype_list <- mapply(function(z,y) rep(z,y), seq(1, ncol(repro_biomass_cohorts)), genotype_contributions_vec)
  # Turn from list into vector
  genotype_list <- as.vector(unlist(genotype_list))
  # Make a shuffled vector, but of the exact same length with the same elements
  shuffled_vector <- sample(genotype_list, length(genotype_list), replace=F)
  # Each genotype is now paired with another, and averages are calculate to produce a midpoint value from which the offspring phenotype will be determined (with some developmental noise added so the expressed phenotype does not exactly correspond to the expectation from the midpoint of the genotypes)
  # Although the phenotypes will be based on these genotype midpoints, the individual will only be tracked as having the genotype value in "genotype_list" and that will be the half of their diploid genome they will pass down if/when they reproduce in the future
  # Produce sequence of numbers corresponding to the midpoints of the genotypes of each of the super individuals (with each midpoint rounded to an integer)
  genotype_list_avg <- round((genotype_list + shuffled_vector)/2)
  # Frequency distribution of genotype midpoints
  genotype_contributions_vec_avg <- tabulate(genotype_list, nbins=genotypes)
  # Produce vector of numbers corresponding to the phenotypes of each of the super individuals:
  # Add developmental noise so the expressed phenotype does not exactly correspond to the expectation from the midpoint of the genotypes
  # For each possible value z of genotype midpoints, make y draws (with y = number of offspring with that genotype midpoint) from a Normal distribution centered at z (and truncated at either end of the possible range of values)
  expression_list <- mapply(function(z,y) round_special(rtruncnorm(y,1,ncol(repro_biomass_cohorts),mean=z,sd=1)), seq(1, ncol(repro_biomass_cohorts)), genotype_contributions_vec_avg)
  # Now each element is the value of the expressed phenotype of each offspring super-individual 
  # Go from list to vector
  expression_list <- as.vector(unlist(expression_list))
  # Get rid of zeroes (inserted in the phenotype vector when there were no offspring individuals having a possible genotype-midpoint value)
  expression_list <- expression_list[which(expression_list>0)]
  # Calculate the frequency distributions (numbers of offspring in each bin, with zeros for bins in the range but with no individuals) and multiply by "divider" to convert back to individuals from superindividuals
  expression_vec <- tabulate(expression_list, nbins=genotypes)*divider
  # Make this vector of frequencies the first row (youngest cohort) of the new "mu expressed" matrix 
  mu_exp[1,] <- expression_vec
  return(list(mu, expression_list, genotype_list, mu_exp))
}

scroll_through_phenotypes <- function(k, expressed, genotype) {
  # Create a row for each z-level in the new mu 3D array, where the row has the proportions of new recruits from each phenotype (column) that are in the genotype that this z-level of the array corresponds to
  row<-mapply(function(y) length(expressed[which(expressed==y & genotype==k)])/length(expressed[which(expressed==y)]), seq(1,genotypes))
  row[is.na(row)] <- 0
  return(row)
}

round_special <- function(value) {
  # A null value will be returned if there are no (super)individuals with that possible genotype midpoint value
  if (is.null(value)) {
    return(0)
  }
  if (!(is.null(value))) {
    return(round(value))
  }
}


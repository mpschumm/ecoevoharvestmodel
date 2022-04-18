#!/Users/mschumm/opt/anaconda3/bin/python

# Importing Python packages
import numpy as np

import math

import random

dists = {}

# Initialize matrix of columns (years) and rows (age classes)
for i in range(1,6):
    dists["age_dist_" + str(i)] = np.zeros((36, 401))

# Add (arbitrarily) 100 individuals to age=10 to start with
for i in range(1,6):
    (dists["age_dist_" + str(i)])[5,0] = 20

# Length at maturity
# Should be allowed to vary with evolution
L_mat = 25
# Von Bertalanffy growth rate
# Should be allowed to vary with evolution
k = 0.3
# Density dependence in growth (Eikeset et al. PNAS 2016)
b_k = -5.34e-14
# Asymptotic length
# Should be allowed to vary with evolution
L_inf = 40
# Allometric exponent for conversion of length to weight
c = 1
# Allometric multiplier for conversion of length to weight
b = 3
# Maturation ogive
q = 0.8
# Mortality coefficient
m = 0.4
# Cost of growth (experimental)
cost = 0.01
# Fishing coefficient
Fishing = 0
# Recruitment Beverton-Holt parameters
alpha = 84
beta = 1
# Number of different genotypes
genotypes = 6

# Obtain lengths at ages based on lengths of the same cohort at the past time point
def length_at_age (past_length_at_age, k, L_inf):
    return past_length_at_age + k*(L_inf - past_length_at_age)

# Density dependent version of length at age function
def dd_length_at_age (past_length_at_age, k, L_inf, t, B_t):
    return past_length_at_age + k*(L_inf - past_length_at_age)*(math.exp(B_t*b_k))*random.gauss(1,0.1)
       
# Convert length to weight
# Modified with FishBase values (!hard-coded!)
def weight_at_age (age, c, b):
    L = 0.005*(lengths_at_ages[age]**3.1)
    return(L)

# Calculate probability of maturation
def p_maturity_at_age (age, q):
    return 1 / (1 + math.exp(-q*(lengths_at_ages[age]-L_mat)))
    
# Calculate egg production for an individual at some age a from product of P(maturity) and weight
# (But, how does this end up in units of eggs? I don't think it does)
def egg_production (age):
    return weight_at_age(age, c, b)*p_maturity_at_age(age, q)

# Total egg production across all individuals, across all age classes
def egg_production_total (ages):
    total = 0
    for i in ages:
        total = total + age_dist[i,t] * (egg_production(i)/90718500000)
    return total

# Density dependent recruitment for this year
def recruitment (age, ages):
    egg_production_total = 0
    for g in range(1, genotypes):
        egg_production_total = egg_production_total + dists["egg_production_totals_"+str(g)]
    for g in range(1, genotypes):
        dists["relative_recruitment_" + str(g)] = dists["egg_production_totals_"+str(g)]/egg_production_total
    if age == 0:
        # currently Beverton-Holt recruitment – try this with other models in the future?
        N = ( (alpha*egg_production_total) / (1 + beta*egg_production_total) ) * 1000000000
    for g in range(1, genotypes):
        age_dist = dists["age_dist_"+str(g)]
        age_dist[0, t] = N*dists["relative_recruitment_" + str(g)]
        print(N*dists["relative_recruitment_" + str(g)])
        dists["age_dist_"+str(g)] = age_dist 
          
# Selective fishing mortality for a given length - rises from 0 to 1 with inflection at length = 10 for herring 
def selective_mortality (length_at_age):
    return 1 / (1 + math.exp(-0.8*(length_at_age - 30)))

# Calculate total number of individuals for an age class based on the cohort's size at the previous time point, mortality and fishing
def number_individuals (age, m, Fishing):
    if age > 0 :
        return age_dist[age-1, t-1] * math.exp(-1*m - Fishing*selective_mortality(lengths_at_ages[age])) * (1-(age_dist[age-1+18, t-1] - age_dist[age-2+18, t-2])*cost)
        
# Calculate total population biomass
def biomass_t (ages):
    b_t = 0
    for a in ages:
        b_t = b_t + age_dist[a, t-1]*weight_at_age(a, c, b)
    return(b_t)

### Run the model
# Initialize empty vector of lengths at different ages
for i in range(1,genotypes):
    dists["lengths_at_ages_" + str(i)] = np.zeros(18)

# Initial size is 4cm.
for i in range(1,genotypes):
    (dists["lengths_at_ages_" + str(i)])[0] = 4

for j in range(1,genotypes):
    for i in range(1,17):
        lengths_at_ages = (dists["lengths_at_ages_" + str(j)])
        lengths_at_ages[i] = length_at_age(lengths_at_ages[i-1],k,L_inf)
        (dists["lengths_at_ages_" + str(j)]) = lengths_at_ages

for i in range(1,genotypes):
    (dists["age_dist_" + str(i)])[18,0] = 4
    
lengths_at_ages = dists["lengths_at_ages_1"]

age_dist = dists["age_dist_1"]

# Initialize the mean total offspring to 0
mean_total_offspring = 0
# Loop through 390 years
for t in range(390):
    if t>200:
        Fishing = 15
    if t>0:
        for i in range(1,genotypes):
            dists["egg_production_totals_" + str(i)] = 0
        b_t=0
        for g in range(1, genotypes):
            age_dist = dists["age_dist_"+str(g)]
            lengths_at_ages = dists["lengths_at_ages_"+str(g)]
            b_t = b_t + biomass_t(range(1,17))
            dists["age_dist_"+str(g)]=age_dist
            dists["lengths_at_ages_"+str(g)] = lengths_at_ages
        for g in range(1, genotypes):
            if t>200:
                # L_mat = 21 + g*1
                k = 0.18 +g*0.04
            age_dist = dists["age_dist_"+str(g)]
            lengths_at_ages = dists["lengths_at_ages_"+str(g)]
            # For each year, loop through each of the non-zero ages
            for a in range(1,17):
                # Calculate the number of individuals in this age class
                age_dist[a,t] = number_individuals(a, m, Fishing)
                # Calculate the length of fish at this age
                value = age_dist[a-1+18, t-1]
                lengths_at_ages[a] = dd_length_at_age(value,k, L_inf, t, b_t)
            dists["egg_production_totals_" + str(g)] = egg_production_total(range(1,16))
            dists["age_dist_"+str(g)]=age_dist
            dists["lengths_at_ages_"+str(g)] = lengths_at_ages
        # For age 0, calculate recruitment 
        recruitment(0, range(1,16))
    # Record the sizes at ages for this year
    for j in range(1,genotypes):
        (dists["age_dist_" + str(j)])[range(18,36),t] = np.transpose(dists["lengths_at_ages_" + str(j)])

for i in range(1,genotypes):
    np.savetxt(("output_"+str(i)),dists["age_dist_"+str(i)])
    
### Cod model    

# Importing Python packages
import numpy as np

import math

dists = {}

# Initialize matrix of columns (years) and rows (age classes)
for i in range(1,6):
    dists["age_dist_" + str(i)] = np.zeros((44, 401))

# Add (arbitrarily) 100 individuals to age=10 to start with
for i in range(1,6):
    (dists["age_dist_" + str(i)])[5,0] = 20

# Length at maturity
# Should be allowed to vary with evolution
L_mat = 35
# Von Bertalanffy growth rate
# Should be allowed to vary with evolution
k = 0.12
# Density dependence in growth (Eikeset et al. PNAS 2016)
b_k = -5.34e-11
# Asymptotic length
# Should be allowed to vary with evolution
L_inf = 130
# Allometric exponent for conversion of length to weight
c = 7e-6
# Allometric multiplier for conversion of length to weight
b = 3
# Maturation ogive
q = 0.2
# Mortality coefficient
m = 0.2
# Cost of growth (experimental)
cost = 0.01
# Fishing coefficient
Fishing = 0
# Recruitment Beverton-Holt parameters
alpha = 0.5
beta = 1e-8
# Number of different genotypes
genotypes = 6

# Obtain lengths at ages based on lengths of the same cohort at the past time point
def length_at_age (past_length_at_age, k, L_inf):
    return past_length_at_age + k*(L_inf - past_length_at_age)

# Density dependent version of length at age function
def dd_length_at_age (past_length_at_age, k, L_inf, t, B_t):
    return past_length_at_age + k*(L_inf - past_length_at_age)*(math.exp(B_t*b_k))
       
# Convert length to weight
# Modified with FishBase values (!hard-coded!)
def weight_at_age (age, c, b):
    L = c*(lengths_at_ages[age]**b)
    return(L)

# Calculate probability of maturation
def p_maturity_at_age (age, q):
    return 1 / (1 + math.exp(-q*(lengths_at_ages[age]-L_mat)))
    
# Calculate egg production for an individual at some age a from product of P(maturity) and weight
# (But, how does this end up in units of eggs? I don't think it does)
def egg_production (age):
    return weight_at_age(age, c, b)*p_maturity_at_age(age, q)

# Total egg production across all individuals, across all age classes
def egg_production_total (ages):
    total = 0
    for i in ages:
        total = total + age_dist[i,t] * egg_production(i)
    return total

# Density dependent recruitment for this year
def recruitment (age, ages):
    egg_production_total = 0
    for g in range(1, genotypes):
        egg_production_total = egg_production_total + dists["egg_production_totals_"+str(g)]
    for g in range(1, genotypes):
        dists["relative_recruitment_" + str(g)] = dists["egg_production_totals_"+str(g)]/egg_production_total
    if age == 0:
        # currently Beverton-Holt recruitment – try this with other models in the future?
        N = ( (alpha*egg_production_total) / (1 + beta*egg_production_total) )
    for g in range(1, genotypes):
        age_dist = dists["age_dist_"+str(g)]
        age_dist[0, t] = N*dists["relative_recruitment_" + str(g)]
        print(N*dists["relative_recruitment_" + str(g)])
        dists["age_dist_"+str(g)] = age_dist 
          
# Selective fishing mortality for a given length - rises from 0 to 1 with inflection at length = 10 for herring 
def selective_mortality (length_at_age):
    return 1 / (1 + math.exp(-0.2*(length_at_age - 50)))

# Calculate total number of individuals for an age class based on the cohort's size at the previous time point, mortality and fishing
def number_individuals (age, m, Fishing):
    if age > 0 :
        return age_dist[age-1, t-1] * math.exp(-1*m - Fishing*selective_mortality(lengths_at_ages[age])) * (1-(age_dist[age-1+22, t-1] - age_dist[age-2+22, t-2])*cost)
        
# Calculate total population biomass
def biomass_t (ages):
    b_t = 0
    for a in ages:
        b_t = b_t + age_dist[a, t-1]*weight_at_age(a, c, b)
    return(b_t)

### Run the model
# Initialize empty vector of lengths at different ages
for i in range(1,genotypes):
    dists["lengths_at_ages_" + str(i)] = np.zeros(22)

# Initial size is 4cm.
for i in range(1,genotypes):
    (dists["lengths_at_ages_" + str(i)])[0] = 4

for j in range(1,genotypes):
    for i in range(1,21):
        lengths_at_ages = (dists["lengths_at_ages_" + str(j)])
        lengths_at_ages[i] = length_at_age(lengths_at_ages[i-1],k,L_inf)
        (dists["lengths_at_ages_" + str(j)]) = lengths_at_ages

for i in range(1,genotypes):
    (dists["age_dist_" + str(i)])[22,0] = 4
    
lengths_at_ages = dists["lengths_at_ages_1"]

age_dist = dists["age_dist_1"]

# Initialize the mean total offspring to 0
mean_total_offspring = 0
# Loop through 390 years
for t in range(390):
    if t>200:
        Fishing = 0.2
    if t>0:
        if b_t<0:
            break
        for i in range(1,genotypes):
            dists["egg_production_totals_" + str(i)] = 0
        b_t=0
        for g in range(1, genotypes):
            age_dist = dists["age_dist_"+str(g)]
            lengths_at_ages = dists["lengths_at_ages_"+str(g)]
            b_t = b_t + biomass_t(range(1,21))
            dists["age_dist_"+str(g)]=age_dist
            dists["lengths_at_ages_"+str(g)] = lengths_at_ages
        for g in range(1, genotypes):
            if t>200:
                # L_mat = 32 + g*1
                k = 0 + g*0.04
            age_dist = dists["age_dist_"+str(g)]
            lengths_at_ages = dists["lengths_at_ages_"+str(g)]
            # For each year, loop through each of the non-zero ages
            for a in range(1,21):
                # Calculate the number of individuals in this age class
                age_dist[a,t] = number_individuals(a, m, Fishing)
                # Calculate the length of fish at this age
                value = age_dist[a-1+22, t-1]
                lengths_at_ages[a] = dd_length_at_age(value,k, L_inf, t, b_t)
            dists["egg_production_totals_" + str(g)] = egg_production_total(range(1,20))
            dists["age_dist_"+str(g)]=age_dist
            dists["lengths_at_ages_"+str(g)] = lengths_at_ages
        # For age 0, calculate recruitment 
        recruitment(0, range(1,20))
    # Record the sizes at ages for this year
    for j in range(1,genotypes):
        (dists["age_dist_" + str(j)])[range(22,44),t] = np.transpose(dists["lengths_at_ages_" + str(j)])

for i in range(1,genotypes):
    np.savetxt(("cod_output_"+str(i)),dists["age_dist_"+str(i)])



    


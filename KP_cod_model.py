#!/Users/mschumm/opt/anaconda3/bin/python

# Importing Python packages
import numpy as np

import math

# Initialize matrix of columns (years) and rows (age classes)
age_dist = np.zeros((44, 401))

# Add (arbitrarily) 100 individuals to age=10 to start with
age_dist[10,0] = 100

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
# Allometric exponent for conversion of length to growth
c = 7e-6
# Allometric multiplier for conversion of length to growth
b = 3
# Maturation ogive
q = 0.2
# Mortality coefficient
m = 0.2
# Cost of growth (experimental)
cost = 0.001
# Fishing coefficient
Fishing = 0.2
# Recruitment Beverton-Holt parameters
alpha = 0.5
beta = 1e-8

# Obtain lengths at ages based on lengths of the same cohort at the past time point
def length_at_age (past_length_at_age, k, L_inf):
    return past_length_at_age + k*(L_inf - past_length_at_age)

# Density dependent version of length at age function
def dd_length_at_age (past_length_at_age, k, L_inf, t, B_t):
    return past_length_at_age + k*(L_inf - past_length_at_age)*(math.exp(B_t*b_k))
       
# Convert length to weight
def weight_at_age (age, c, b):
    L = lengths_at_ages[age]
    return(c*(L**b))

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
    if age == 0:
        # currently Beverton-Holt recruitment – try this with other models in the future?
        N = ( (alpha*egg_production_total(ages)) / (1 + beta*egg_production_total(ages)) )
    return N
   
# Selective fishing mortality for a given length - rises from 0 to 1 with inflection at length = 50 
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
lengths_at_ages = np.zeros(22)
# Initial size is 4cm.
lengths_at_ages[0] = 4
for i in range(1,21):
    lengths_at_ages[i] = length_at_age(lengths_at_ages[i-1],k,L_inf)
age_dist[22,0] = 4
# Initialize the mean total offspring to 0
mean_total_offspring = 0
# Loop through 390 years
for t in range(0,390):
    if t>0:
        # For each year, loop through each of the non-zero ages
        b_t = biomass_t(range(1,21))
        for a in range(1,21):
            # Calculate the number of individuals in this age class
            age_dist[a,t] = number_individuals(a, m, Fishing)
            # Calculate the length of fish at this age
            lengths_at_ages[a] = dd_length_at_age(age_dist[a-1+22, t-1],k, L_inf, t, b_t)
            # Lifetime fitness tracker: calculates the fitness for an individual in the model, at the final time step of the model, by summing (survivorship*offspring production) across ages
            if t == 388:
                # Calculate fraction surviving, l, at this age
                l = age_dist[a,t-1]/age_dist[0,t-1]
                print(egg_production(20))
                # Adding to the mean total offspring
                mean_total_offspring = mean_total_offspring + egg_production(a)*l
        # For age 0, calculate recruitment 
        age_dist[0, t] = recruitment(0, range(1,20))
    # Record the sizes at ages for this year
    age_dist[range(22,44),t] = np.transpose(lengths_at_ages)

print(age_dist.shape)

print(mean_total_offspring)

np.savetxt("output",age_dist)

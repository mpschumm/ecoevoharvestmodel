

r = 6
l_bar = 0
w = 0.674

reversible_masses <- vector(mode="numeric", length=100)
lengths <- vector(mode="numeric", length=100)
lambdas <- vector(mode="numeric", length=100)
Ps <- vector(mode="numeric", length=100)
repro_vec <- vector(mode="numeric", length=100)
reversible_masses[1] <- 3.938
lengths[1] <-0.1102
lambdas[1] <- reversible_masses[1]/(structural_mass(lengths[1]))
mortalities <- vector(mode="numeric", length=100)
for (age in 2:longevity) {
  new_E_l = get_pre_reproductive_size(reversible_masses[age-1],lengths[age-1], age-1)
  # lambdas[age] = intrinsic_lambda(lengths[age-1], l_bar, r, lambda_min, lambda_max)
  Ps[age] <- size_dependent_energy_intake(structural_mass(lengths[age-1]), p_0, p_1) - maintenance_cost(lengths[age-1], reversible_masses[age-1], c_S, c_E)
  reversible_masses[age] = new_E_l[1]
  lengths[age] = new_E_l[2]
  lambdas[age] = new_E_l[1]/(structural_mass(new_E_l[2]))
  repro_vec[age] <- max(0, (w*reversible_masses[age] - (r_0*structural_mass(lengths[age]))^(r_1)))
  mortalities[age] <- exp( -m + -Fishing*selective_mortality(lengths[age]) - m_c_max*exp(-z_c*lambdas[age]) - (m_p_max-m)*exp(-z_p*lengths[age]))
}

plot(lengths, type="l")
plot(reversible_masses, type="l")
plot(lambdas, type="l")
plot(Ps, type="l")
plot(repro_vec, type="l")


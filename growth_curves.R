

r = 10
l_bar = 0.2
w = 0.674

reversible_masses <- vector(mode="numeric", length=100)
lengths <- vector(mode="numeric", length=100)
reversible_masses[1] <- 0.3066
lengths[1] <- 0.0493
for (age in 2:longevity) {
  new_E_l = get_pre_reproductive_size(reversible_masses[age-1],lengths[age-1], age-1)
  reversible_masses[age] = new_E_l[1]
  lengths[age] = new_E_l[2]
}


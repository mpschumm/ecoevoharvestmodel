# Fishing mortality is size selective
selective_mortality <- function(length) {
  1 / (1 + exp(-q_f*(length - mincatchsize)*100))
}

# Mortality drops as E/S increases, as S increases
mortality <- function(age) {
  lambda <- ((age_dist[age+(longevity),t]))/structural_mass(age_dist[age+(2*longevity),t])
  length = age_dist[age+(2*longevity),t]
  if (age < longevity) {
    N_t = age_dist[age-1,t-1]*exp( -m + -Fishing*selective_mortality(length) - m_max*exp(-z_c*lambda) - (m_max-m)*exp(-z_p*length) )
  }
  if (age == longevity) {
    N_t = age_dist[age-1,t-1]*exp( -m + -Fishing*selective_mortality(length) - m_max*exp(-z_c*lambda) - (m_max-m)*exp(-z_p*length))
    N_t = age_dist[age,t-1]*exp( -m + -Fishing*selective_mortality(length) - m_max*exp(-z_c*lambda) - (m_max-m)*exp(-z_p*length))
  }
  return(N_t)
}

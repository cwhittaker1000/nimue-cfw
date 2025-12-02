################################################################################
### Vaccine model: deterministic ###############################################
################################################################################

### Initial setup ##############################################################
N_age <- user() # Number of age groups
N_vaccine <- user() # Number of vaccine groups
################################################################################

## In this version:
##  - i indexes age-groups
##  - j indexes vaccination status
##  - k indexes location

## Time output to line up with squire fitting infrastructure
time <- t
output(time) <- TRUE

### S: susceptible #############################################################

## Setting up initial conditions
S_0[,] <- user()
dim(S_0) <- c(N_age, N_vaccine)
initial(S[,]) <- S_0[i,j]
dim(S) <- c(N_age, N_vaccine)

## ODEs for S compartments
### Note: potentially need to replace vaccine_efficacy_infection_t with a constant vaccine_efficacy_infection
deriv(S[,1]) <- -(lambda[i] * vaccine_efficacy_infection_t[i,j] * S[i,j]) - (vr * vaccination_target[i] * S[i,j])
deriv(S[,2]) <- - (lambda[i] * vaccine_efficacy_infection_t[i,j] * S[i,j]) + (vr * vaccination_target[i] * S[i,j-1]) - (gamma_vaccine[j] * S[i,j])
deriv(S[,3:N_vaccine]) <- (gamma_vaccine[j-1] * S[i,j-1]) - (lambda[i] * vaccine_efficacy_infection_t[i,j] * S[i,j]) - (gamma_vaccine[j] * S[i,j])

## Outputs
output(S_overall[]) <- sum(S[i,])
dim(S_overall) <- N_age

################################################################################

### E : Latent ########################################################

## Setting up initial conditions
E_0[,] <- user()
dim(E_0) <- c(N_age, N_vaccine)
initial(E[,]) <- E_0[i,j]
dim(E) <- c(N_age, N_vaccine)

## Parameters related to E
gamma_E <- user() # rate of progression through latent infection

## ODEs for E compartments
### Note: potentially need to replace vaccine_efficacy_infection_t with a constant vaccine_efficacy_infection
deriv(E[,1]) <- (lambda[i] * vaccine_efficacy_infection_t[i,j] * S[i,j]) - (gamma_E * E[i,j]) - (vr * vaccination_target[i] * E[i,j])
deriv(E[,2]) <- (lambda[i] * vaccine_efficacy_infection_t[i,j] * S[i,j]) - (gamma_E * E[i,j]) + (vr * vaccination_target[i] * E[i,j-1]) - (gamma_vaccine[j] * E[i,j])
deriv(E[,3:N_vaccine]) <- (gamma_vaccine[j-1] * E[i,j-1]) + (lambda[i] * vaccine_efficacy_infection_t[i,j] * S[i,j]) - (gamma_E * E[i,j]) - (gamma_vaccine[j] * E[i,j])

## Outputs
output(E_overall[]) <- sum(E[i,])
dim(E_overall) <- N_age

################################################################################

### IMild: Unhospitalised infection ############################################

## Setting up initial conditions
IMild_0[,] <- user()
dim(IMild_0) <- c(N_age, N_vaccine)
initial(IMild[,]) <- IMild_0[i,j]
dim(IMild) <- c(N_age, N_vaccine)

## Parameters related to E
gamma_IMild <- user() # rate of progression from mild infection to recovery

## ODEs for E compartments
deriv(IMild[,1]) <- (gamma_E * E[i,j] * (1 - prob_hosp[i,j])) - (gamma_IMild * IMild[i,j])
deriv(IMild[,2]) <- (gamma_E * E[i,j] * (1 - prob_hosp[i,j])) - (gamma_IMild * IMild[i,j]) - (gamma_vaccine[j] * IMild[i,j])
deriv(IMild[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMild[i,j-1]) + (gamma_E * E[i,j] * (1 - prob_hosp[i,j])) - (gamma_IMild * IMild[i,j]) - (gamma_vaccine[j] * IMild[i,j])

## Outputs
output(IMild_overall[]) <- sum(IMild[i,])
dim(IMild_overall) <- N_age

################################################################################

### R: Recovered ####################################################

## Setting up initial conditions
R_0[,] <- user()
dim(R_0) <- c(N_age, N_vaccine)
initial(R[,]) <- R_0[i,j]
dim(R) <- c(N_age, N_vaccine)

## ODEs for R compartments
## Note - if at any point we make duration of hospitalisation diff for recovery vs dying people, need to have separate compartments, not just apportioning of flows
deriv(R[,1]) <- (gamma_IHosp * IHosp[i,j] * (1 - prob_death_hosp[i, j])) + (gamma_IMild * IMild[i,j]) - (vr * vaccination_target[i] * R[i,j])
deriv(R[,2]) <- (vr * vaccination_target[i] * R[i,j-1]) + (gamma_IHosp * IHosp[i,j] * (1 - prob_death_hosp[i, j])) + (gamma_IMild * IMild[i,j]) - (gamma_vaccine[j] * R[i,j])
deriv(R[,3:N_vaccine]) <- (gamma_vaccine[j-1] * R[i,j-1]) + (gamma_IHosp * IHosp[i,j] * (1 - prob_death_hosp[i, j])) + (gamma_IMild * IMild[i,j]) - (gamma_vaccine[j] * R[i,j])

## Outputs
output(R_overall[]) <- sum(R[i,])
dim(R_overall) <- N_age

################################################################################

### ICase (ICase1 & ICase2): To-be hospitalised infection ######################

## Setting up initial conditions
ICase_0[,] <- user()
dim(ICase_0) <- c(N_age, N_vaccine)
initial(ICase[,]) <- ICase_0[i,j]
dim(ICase) <- c(N_age, N_vaccine)

## Parameters related to ICase
gamma_ICase <- user() # rate of progression from symptom onset to requiring hospitalisation

## ODEs for ICase compartments
deriv(ICase[,1]) <- (gamma_E * E[i,j] * prob_hosp[i,j]) - (gamma_ICase * ICase[i,j])
deriv(ICase[,2]) <- (gamma_E * E[i,j] * prob_hosp[i,j]) - (gamma_ICase * ICase[i,j]) - (gamma_vaccine[j] * ICase[i,j])
deriv(ICase[,3:N_vaccine]) <- (gamma_vaccine[j-1] * ICase[i,j-1]) + (gamma_E * E[i,j] * prob_hosp[i,j]) - (gamma_ICase * ICase[i,j]) - (gamma_vaccine[j] * ICase[i,j])

## Outputs
output(ICase_overall[]) <- sum(ICase[i,])
dim(ICase_overall) <- N_age

################################################################################

### IHosp: Go to hospital #################################

## Setting up initial conditions
IHosp_0[,] <- user()
dim(IHosp_0) <- c(N_age, N_vaccine)
initial(IHosp[,]) <- IHosp_0[i,j]
dim(IHosp) <- c(N_age, N_vaccine)

## Parameters related to IHosp
gamma_IHosp <- user() # Note: rate of progression from hospitalisation to death ## Note that if we add different amounts of time in hospital for recovering vs dying people,
                      #       we'll need explicit compartments for each and a probability that shapes flow from ICase -> IHosp_Live & IHosp_Die

## ODEs for IHosp compartments
## Note: probs replace prob_severe_multi with prob_severe
deriv(IHosp[,1]) <- (gamma_ICase * ICase[i,j]) - (gamma_IHosp * IHosp[i,j]) # Note: if we make duration of stay different for dying vs recovering, need separate compartments, not just apportioning flow
deriv(IHosp[,2]) <- (gamma_ICase * ICase[i,j]) - (gamma_IHosp * IHosp[i,j]) - (gamma_vaccine[j] * IHosp[i,j])
deriv(IHosp[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IHosp[i,j-1]) + (gamma_ICase * ICase[i,j]) - (gamma_IHosp * IHosp[i,j]) - (gamma_vaccine[j] * IHosp[i,j])

################################################################################

### D: Dead ####################################################################

## Setting up initial conditions
D_0[,] <- user()
dim(D_0) <- c(N_age, N_vaccine)
initial(D[,]) <- D_0[i,j]
dim(D) <- c(N_age, N_vaccine)

## ODEs for D compartments
deriv(D[,1]) <- (gamma_IHosp * IHosp[i,j] * prob_death_hosp[i,j])
deriv(D[,2]) <- (gamma_IHosp * IHosp[i,j] * prob_death_hosp[i,j]) - (gamma_vaccine[j] * D[i,j])
deriv(D[,3:N_vaccine]) <- (gamma_vaccine[j-1] * D[i,j-1]) + (gamma_IHosp * IHosp[i,j] * prob_death_hosp[i,j]) - (gamma_vaccine[j] * D[i,j])

#################################################################################

################################################################################
### Vaccination capacity #######################################################
################################################################################
# Vaccination
# Vaccine prioritisation coverage matrix
N_prioritisation_steps <- user()
vaccine_coverage_mat[,] <- user()
dim(vaccine_coverage_mat) <- c(N_prioritisation_steps, N_age)

# Generating Vaccine Efficacy Over Time
vaccine_efficacy_infection_t[, ] <- interpolate(tt_vaccine_efficacy_infection, vaccine_efficacy_infection, "constant")
dim(vaccine_efficacy_infection_t) <- c(N_age, N_vaccine)
tt_vaccine_efficacy_infection[] <- user()
vaccine_efficacy_infection[, ,] <- user()
dim(tt_vaccine_efficacy_infection) <- user()
dim(vaccine_efficacy_infection) <- c(length(tt_vaccine_efficacy_infection), N_age, N_vaccine)

gamma_vaccine[] <- user() # Vector of vaccine progression parameters by vaccination status (First = 0 as handled separately as time-varying vaccination rate, Last = 0 as no progression from "previously vaccinated group)
dim(gamma_vaccine) <- N_vaccine

# Interpolation of vaccination rate over time
mv <- interpolate(tt_vaccine, max_vaccine, "constant")
tt_vaccine[] <- user()
max_vaccine[] <- user()
dim(tt_vaccine) <- user()
dim(max_vaccine) <- length(tt_vaccine)

# Track the proportion who have received vaccine in each age group

# Population size
pop_size[] <- sum(S[i,]) + sum(E[i,]) + sum(IMild[i,]) + sum(ICase[i,]) + sum(IHosp[i,]) + sum(R[i,])
dim(pop_size) <- N_age

# Proportion who have received vaccine
pr[] <- 1 - ((sum(S[i,1]) + sum(E[i,1]) + sum(IMild[i,1]) + sum(ICase[i,1]) + sum(IHosp[i,1]) + sum(R[i,1])) / (pop_size[i]))
dim(pr) <- N_age

# Isolate age groups below current target coverage which must be targeted
vaccination_target_mat[,] <- if(pr[j] < vaccine_coverage_mat[i,j]) 1 else 0
dim(vaccination_target_mat) <- c(N_prioritisation_steps, N_age)

vaccine_target_vec[] <- if(sum(vaccination_target_mat[i,]) == 0) 1 else 0
dim(vaccine_target_vec) <- N_prioritisation_steps
current_index <- min(sum(vaccine_target_vec) + 1, N_prioritisation_steps)

vaccination_target[] <- vaccination_target_mat[as.integer(current_index),i]
dim(vaccination_target) <- N_age

vr_temp[] <- S[i,1] * vaccination_target[i] + E[i,1] * vaccination_target[i] + + R[i,1] * vaccination_target[i]
dim(vr_temp) <- N_age
# Catch so vaccination rate does not exceed 1 if the number of people available for vaccination < number of vaccines
vr_den <- if(sum(vr_temp) <= mv) mv else sum(vr_temp)
vr <- if(mv == 0) 0 else mv / vr_den  # Vaccination rate to achieve capacity given number in vaccine-eligible population
################################################################################
################################################################################

################################################################################
### Hospitalisation and Transmission Dynamics Parameters
################################################################################

prob_hosp[,] <- user()
dim(prob_hosp) <- c(N_age, N_vaccine)

prob_death_hosp[,]  <- user()
dim(prob_death_hosp) <- c(N_age, N_vaccine)

rel_infectiousness[] <- user() # Relative infectiousness of age categories relative to maximum infectiousness age category
dim(rel_infectiousness) <- N_age

rel_infectiousness_vaccinated[,] <- user() # Relative infectiousness of vaccinated age categories relative to maximum infectiousness age category
dim(rel_infectiousness_vaccinated) <- c(N_age, N_vaccine)


################################################################################
### FOI and contact matrix #####################################################
################################################################################
# Generating Force of Infection
m[, ] <- interpolate(tt_matrix, mix_mat_set, "constant")
dim(m) <- c(N_age, N_age)
tt_matrix[] <- user()
mix_mat_set[, ,] <- user()
dim(tt_matrix) <- user()
dim(mix_mat_set) <- c(length(tt_matrix), N_age, N_age)

# Interpolation for beta
beta <- interpolate(tt_beta, beta_set, "constant")
tt_beta[] <- user()
beta_set[] <- user()
dim(tt_beta) <- user()
dim(beta_set) <- length(tt_beta)

# Generating Force of Infection
temp_rel[,] <- (IMild[i,j] * rel_infectiousness_vaccinated[i,j]) + (ICase[i,j] * rel_infectiousness_vaccinated[i,j])
temp[] <- sum(temp_rel[i,])
dim(temp_rel) <- c(N_age, N_vaccine)
dim(temp) <- c(N_age)

s_ij[,] <- m[i, j] * temp[j] * rel_infectiousness[j]
dim(s_ij) <- c(N_age, N_age)

lambda[] <- beta * sum(s_ij[i, ])
dim(lambda) <- N_age
################################################################################
################################################################################

################################################################################
### Output #####################################################################
################################################################################

# Unvaccinated
output(unvaccinated[]) <- sum(S[i,1]) + sum(E[i,1]) + sum(IMild[i,1]) + sum(ICase[i,1]) + sum(IHosp[i,1]) + sum(R[i,1]) + sum(D[i,1])
dim(unvaccinated) <- N_age
# Vaccinated
output(vaccinated[]) <- sum(S[i,2:N_vaccine]) + sum(E[i,2:N_vaccine]) + sum(IMild[i,2:N_vaccine]) + sum(ICase[i,2:N_vaccine]) + sum(IHosp[i,2:N_vaccine]) + sum(R[i,2:N_vaccine]) + sum(D[i,2:N_vaccine])
dim(vaccinated) <- N_age

output(N[]) <- pop_size[i] + sum(D[i,])
dim(N) <- N_age
################################################################################


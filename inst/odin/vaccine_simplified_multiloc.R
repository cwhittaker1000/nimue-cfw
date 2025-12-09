################################################################################
### Vaccine model: deterministic ###############################################
################################################################################

### Initial setup ##############################################################
N_age <- user() # Number of age groups
N_vaccine <- user() # Number of vaccine groups
N_locations <- user() # Number of locations to simulate
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
S_0[,,] <- user()
dim(S_0) <- c(N_age, N_vaccine,N_locations)
initial(S[,,]) <- S_0[i,j,k]
dim(S) <- c(N_age, N_vaccine,N_locations)

## ODEs for S compartments
### Note: potentially need to replace vaccine_efficacy_infection_t with a constant vaccine_efficacy_infection
deriv(S[,1,]) <- -(lambda[i,k] * vaccine_efficacy_infection[i,j] * S[i,j,k]) - (vr[k] * vaccination_target[i,k] * S[i,j,k])
deriv(S[,2,]) <- - (lambda[i,k] * vaccine_efficacy_infection[i,j] * S[i,j,k]) + (vr[k] * vaccination_target[i,k] * S[i,j-1,k]) - (gamma_vaccine[j] * S[i,j,k])
deriv(S[,3:N_vaccine,]) <- (gamma_vaccine[j-1] * S[i,j-1,k]) - (lambda[i,k] * vaccine_efficacy_infection[i,j] * S[i,j,k]) - (gamma_vaccine[j] * S[i,j,k])

## Outputs
output(S_overall[,]) <- sum(S[i,,j])
dim(S_overall) <- c(N_age, N_locations)

################################################################################

### E : Latent ########################################################

## Setting up initial conditions
E_0[,,] <- user()
dim(E_0) <- c(N_age, N_vaccine,N_locations)
initial(E[,,]) <- E_0[i,j,k]
dim(E) <- c(N_age, N_vaccine,N_locations)

## Parameters related to E
gamma_E <- user() # rate of progression through latent infection

## ODEs for E compartments
### Note: potentially need to replace vaccine_efficacy_infection_t with a constant vaccine_efficacy_infection
deriv(E[,1,]) <- (lambda[i,k] * vaccine_efficacy_infection[i,j] * S[i,j,k]) - (gamma_E * E[i,j,k]) - (vr[k] * vaccination_target[i,k] * E[i,j,k])
deriv(E[,2,]) <- (lambda[i,k] * vaccine_efficacy_infection[i,j] * S[i,j,k]) - (gamma_E * E[i,j,k]) + (vr[k] * vaccination_target[i,k] * E[i,j-1,k]) - (gamma_vaccine[j] * E[i,j,k])
deriv(E[,3:N_vaccine,]) <- (gamma_vaccine[j-1] * E[i,j-1,k]) + (lambda[i,k] * vaccine_efficacy_infection[i,j] * S[i,j,k]) - (gamma_E * E[i,j,k]) - (gamma_vaccine[j] * E[i,j,k])

## Outputs
output(E_overall[,]) <- sum(E[i,,j])
dim(E_overall) <- c(N_age, N_locations)

################################################################################

### IMild: Unhospitalised infection ############################################

## Setting up initial conditions
IMild_0[,,] <- user()
dim(IMild_0) <- c(N_age, N_vaccine, N_locations)
initial(IMild[,,]) <- IMild_0[i,j,k]
dim(IMild) <- c(N_age, N_vaccine, N_locations)

## Parameters related to E
gamma_IMild <- user() # rate of progression from mild infection to recovery

## ODEs for E compartments
deriv(IMild[,1,]) <- (gamma_E * E[i,j,k] * (1 - prob_hosp[i,j])) - (gamma_IMild * IMild[i,j,k])
deriv(IMild[,2,]) <- (gamma_E * E[i,j,k] * (1 - prob_hosp[i,j])) - (gamma_IMild * IMild[i,j,k]) - (gamma_vaccine[j] * IMild[i,j,k])
deriv(IMild[,3:N_vaccine,]) <- (gamma_vaccine[j-1] * IMild[i,j-1,k]) + (gamma_E * E[i,j,k] * (1 - prob_hosp[i,j])) - (gamma_IMild * IMild[i,j,k]) - (gamma_vaccine[j] * IMild[i,j,k])

## Outputs
output(IMild_overall[,]) <- sum(IMild[i,,j])
dim(IMild_overall) <-  c(N_age, N_locations)

################################################################################

### R: Recovered ####################################################

## Setting up initial conditions
R_0[,,] <- user()
dim(R_0) <- c(N_age, N_vaccine, N_locations)
initial(R[,,]) <- R_0[i,j,k]
dim(R) <- c(N_age, N_vaccine, N_locations)

## ODEs for R compartments
## Note - if at any point we make duration of hospitalisation diff for recovery vs dying people, need to have separate compartments, not just apportioning of flows
deriv(R[,1,]) <- (gamma_IHosp * IHosp[i,j,k] * (1 - prob_death_hosp[i, j])) + (gamma_IMild * IMild[i,j,k]) - (vr[k] * vaccination_target[i,k] * R[i,j,k])
deriv(R[,2,]) <- (vr[k] * vaccination_target[i,k] * R[i,j-1,k]) + (gamma_IHosp * IHosp[i,j,k] * (1 - prob_death_hosp[i, j])) + (gamma_IMild * IMild[i,j,k]) - (gamma_vaccine[j] * R[i,j,k])
deriv(R[,3:N_vaccine,]) <- (gamma_vaccine[j-1] * R[i,j-1,k]) + (gamma_IHosp * IHosp[i,j,k] * (1 - prob_death_hosp[i, j])) + (gamma_IMild * IMild[i,j,k]) - (gamma_vaccine[j] * R[i,j,k])

## Outputs
output(R_overall[,]) <- sum(R[i,,j])
dim(R_overall) <- c(N_age, N_locations)

################################################################################

### ICase (ICase1 & ICase2): To-be hospitalised infection ######################

## Setting up initial conditions
ICase_0[,,] <- user()
dim(ICase_0) <- c(N_age, N_vaccine, N_locations)
initial(ICase[,,]) <- ICase_0[i,j,k]
dim(ICase) <- c(N_age, N_vaccine, N_locations)

## Parameters related to ICase
gamma_ICase <- user() # rate of progression from symptom onset to requiring hospitalisation

## ODEs for ICase compartments
deriv(ICase[,1,]) <- (gamma_E * E[i,j,k] * prob_hosp[i,j]) - (gamma_ICase * ICase[i,j,k])
deriv(ICase[,2,]) <- (gamma_E * E[i,j,k] * prob_hosp[i,j]) - (gamma_ICase * ICase[i,j,k]) - (gamma_vaccine[j] * ICase[i,j,k])
deriv(ICase[,3:N_vaccine,]) <- (gamma_vaccine[j-1] * ICase[i,j-1,k]) + (gamma_E * E[i,j,k] * prob_hosp[i,j]) - (gamma_ICase * ICase[i,j,k]) - (gamma_vaccine[j] * ICase[i,j,k])

## Outputs
output(ICase_overall[,]) <- sum(ICase[i,,j])
dim(ICase_overall) <- c(N_age, N_locations)

################################################################################

### IHosp: Go to hospital #################################

## Setting up initial conditions
IHosp_0[,,] <- user()
dim(IHosp_0) <- c(N_age, N_vaccine, N_locations)
initial(IHosp[,,]) <- IHosp_0[i,j,k]
dim(IHosp) <- c(N_age, N_vaccine, N_locations)

## Parameters related to IHosp
gamma_IHosp <- user() # Note: rate of progression from hospitalisation to death ## Note that if we add different amounts of time in hospital for recovering vs dying people,
                      #       we'll need explicit compartments for each and a probability that shapes flow from ICase -> IHosp_Live & IHosp_Die

## ODEs for IHosp compartments
## Note: probs replace prob_severe_multi with prob_severe
deriv(IHosp[,1,]) <- (gamma_ICase * ICase[i,j,k]) - (gamma_IHosp * IHosp[i,j,k]) # Note: if we make duration of stay different for dying vs recovering, need separate compartments, not just apportioning flow
deriv(IHosp[,2,]) <- (gamma_ICase * ICase[i,j,k]) - (gamma_IHosp * IHosp[i,j,k]) - (gamma_vaccine[j] * IHosp[i,j,k])
deriv(IHosp[,3:N_vaccine,]) <- (gamma_vaccine[j-1] * IHosp[i,j-1,k]) + (gamma_ICase * ICase[i,j,k]) - (gamma_IHosp * IHosp[i,j,k]) - (gamma_vaccine[j] * IHosp[i,j,k])

## Outputs
output(IHosp_overall[,]) <- sum(IHosp[i,,j])
dim(IHosp_overall) <- c(N_age, N_locations)

################################################################################

### D: Dead ####################################################################

## Setting up initial conditions
D_0[,,] <- user()
dim(D_0) <- c(N_age, N_vaccine, N_locations)
initial(D[,,]) <- D_0[i,j,k]
dim(D) <- c(N_age, N_vaccine, N_locations)

## ODEs for D compartments
deriv(D[,1,]) <- (gamma_IHosp * IHosp[i,j,k] * prob_death_hosp[i,j])
deriv(D[,2,]) <- (gamma_IHosp * IHosp[i,j,k] * prob_death_hosp[i,j]) - (gamma_vaccine[j] * D[i,j,k])
deriv(D[,3:N_vaccine,]) <- (gamma_vaccine[j-1] * D[i,j-1,k]) + (gamma_IHosp * IHosp[i,j,k] * prob_death_hosp[i,j]) - (gamma_vaccine[j] * D[i,j,k])

## Outputs
output(D_overall[,]) <- sum(D[i,,j])
dim(D_overall) <- c(N_age, N_locations)

#################################################################################

################################################################################
### Vaccination capacity #######################################################
################################################################################

# Defining key vaccination inputs

## Vaccine prioritisation coverage matrix - assumed same across locations
N_prioritisation_steps <- user()
vaccine_coverage_mat[,] <- user()
dim(vaccine_coverage_mat) <- c(N_prioritisation_steps, N_age)

## Vaccine Efficacy against infection - assumed same across locations
vaccine_efficacy_infection[, ] <- user()
dim(vaccine_efficacy_infection) <- c(N_age, N_vaccine)
gamma_vaccine[] <- user() # Vector of vaccine progression parameters by vaccination status (First = 0 as handled separately as time-varying vaccination rate, Last = 0 as no progression from "previously vaccinated group)
dim(gamma_vaccine) <- N_vaccine

## Interpolation of vaccination rate over time
## Single loc version
# mv <- interpolate(tt_vaccine, max_vaccine, "constant")
# tt_vaccine[] <- user()
# max_vaccine[] <- user()
# dim(tt_vaccine) <- user()
# dim(max_vaccine) <- length(tt_vaccine)

tt_vaccine[] <- user()
dim(tt_vaccine) <- user()  # length T
max_vaccine_set[ , ] <- user()  # time x location
dim(max_vaccine_set) <- c(length(tt_vaccine), N_locations)
mv[] <- interpolate(tt_vaccine, max_vaccine_set, "constant")
dim(mv) <- N_locations

# Track the proportion who have received vaccine in each age group

## Population size
pop_size[,] <- sum(S[i,,j]) + sum(E[i,,j]) + sum(IMild[i,,j]) + sum(ICase[i,,j]) + sum(IHosp[i,,j]) + sum(R[i,,j])
dim(pop_size) <- c(N_age, N_locations)

## Proportion who have received vaccine
pr[,] <- 1 - ((sum(S[i,1,j]) + sum(E[i,1,j]) + sum(IMild[i,1,j]) + sum(ICase[i,1,j]) + sum(IHosp[i,1,j]) + sum(R[i,1,j])) / (pop_size[i,j]))
dim(pr) <- c(N_age, N_locations)

## Isolate age groups below current target coverage which must be targeted
vaccination_target_mat[,,] <- if(pr[j,k] < vaccine_coverage_mat[i,j]) 1 else 0
dim(vaccination_target_mat) <- c(N_prioritisation_steps, N_age, N_locations)

vaccine_target_vec[,] <- if(sum(vaccination_target_mat[i, , j]) == 0) 1 else 0 ## for each prioritisation and location, check whether all age-groups have been vaccinated to the level required
dim(vaccine_target_vec) <- c(N_prioritisation_steps, N_locations)
current_index[] <- min(sum(vaccine_target_vec[, i]) + 1, N_prioritisation_steps) ## picking the first prioritisation step which isn't 1 (i.e. vaccine target vec hasn't been reached across all age groups)
dim(current_index) <- N_locations

vaccination_target[,] <- vaccination_target_mat[as.integer(current_index[j]),i, j]
dim(vaccination_target) <- c(N_age, N_locations)
vr_temp[,] <- S[i,1,j] * vaccination_target[i, j] + E[i,1,j] * vaccination_target[i, j] + R[i,1,j] * vaccination_target[i, j]
dim(vr_temp) <- c(N_age, N_locations)
## vr_temp is total number of individuals in age group and location (across disease states) and disease available to be vaccinated,
## conditional on target coverage and prioritisation step we're at

# Catch so vaccination rate does not exceed 1 if the number of people available for vaccination < number of vaccines
vr_den[] <- if(sum(vr_temp[, i]) <= mv[i]) mv[i] else sum(vr_temp[, i])
dim(vr_den) <- N_locations
vr[] <- if(mv[i] == 0) 0 else mv[i] / vr_den[i]  # Vaccination rate to achieve capacity given number in vaccine-eligible population
                                           # this is really "what fraction of the vr_temp i.e. vaccinable pop, can be vaccinated in a timestep with mv
dim(vr) <- N_locations

######################################################################
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

## Defining the within-location mixing matrix (comprising a mixing matrix for each location hence N_age x N_age x N_locations)
mix_mat_set[, , ] <- user()
dim(mix_mat_set) <- c(N_age, N_age, N_locations)
within_loc_mixing[, ,] <- mix_mat_set[i,j,k]
dim(within_loc_mixing) <- c(N_age, N_age, N_locations)

## Defining the between location mixing matrix
##   q[i] = fraction of time residents of location i spend away from home
##   pi_travel[i, j] = probability (conditional on being away from i) that a trip from home i goes to destination j
##   Require: pi_travel[i, i] = 0, and sum_j pi_travel[i, j] = 1 for each i
q[] <- user()
dim(q) <- N_locations
pi_travel[,] <- user()
dim(pi_travel) <- c(N_locations, N_locations)  # home (i), dest (j)

## Defining between-location mixing matrix (comprising a mixing matrix for between locations, initially assuming no age-specific mobility)

# Interpolation for beta ## Double check I HAVE THE MATRIX BETA_SET THE RIGHT WAY ROUND HERE I THINK BASED ON https://mrc-ide.github.io/odin/articles/odin.html#interpolating-functions
tt_beta[] <- user()
beta_set[,] <- user()
beta[] <- interpolate(tt_beta, beta_set, "constant")
dim(tt_beta) <- user()
dim(beta_set) <- c(length(tt_beta), N_locations)
dim(beta) <- N_locations

# Generating Force of Infection

## 1) Infectiousness by age & location
temp_rel[,,] <- (IMild[i,j,k] * rel_infectiousness_vaccinated[i,j]) + (ICase[i,j,k] * rel_infectiousness_vaccinated[i,j])
dim(temp_rel) <- c(N_age, N_vaccine, N_locations)
temp[,] <- sum(temp_rel[i,,j])
dim(temp) <- c(N_age, N_locations)

## 2) Age-structured mixing within each *physical* location
s_ij[,,] <- within_loc_mixing[i, j, k] * temp[j,k] * rel_infectiousness[j]
dim(s_ij) <- c(N_age, N_age, N_locations)

## 3) Local FOI: FOI experienced when physically in location i
##    (susceptible age = i, location = j)
lambda_local[,] <- beta[j] * sum(s_ij[i, ,j])
dim(lambda_local) <- c(N_age, N_locations)

## 4) Travel FOI: define travel lambda here
## For each age i and *home* location j:
##   lambda_travel[i, j] = sum_k pi_travel[j, k] * lambda_local[i, k]
## i.e. average FOI over destinations, conditional on being away
travel_weighted[,,] <- lambda_local[i, k] * pi_travel[j, k]
dim(travel_weighted) <- c(N_age, N_locations, N_locations)
lambda_travel[,] <- sum(travel_weighted[i, j, ])
dim(lambda_travel) <- c(N_age, N_locations)

## 5 Total FOI : Add local and travel lambda here
lambda[,] <- (1 - q[j]) * lambda_local[i, j] + q[j] * lambda_travel[i, j]
dim(lambda) <- c(N_age, N_locations)

################################################################################
### Output #####################################################################
################################################################################

# Unvaccinated
output(unvaccinated[,]) <- sum(S[i,1,j]) + sum(E[i,1,j]) + sum(IMild[i,1,j]) + sum(ICase[i,1,j]) + sum(IHosp[i,1,j]) + sum(R[i,1,j]) + sum(D[i,1,j])
dim(unvaccinated) <- c(N_age, N_locations)
# Vaccinated
output(vaccinated[,]) <- sum(S[i,2:N_vaccine,j]) + sum(E[i,2:N_vaccine,j]) + sum(IMild[i,2:N_vaccine,j]) + sum(ICase[i,2:N_vaccine,j]) + sum(IHosp[i,2:N_vaccine,j]) + sum(R[i,2:N_vaccine,j]) + sum(D[i,2:N_vaccine,j])
dim(vaccinated) <- c(N_age, N_locations)

output(N[,]) <- pop_size[i,j] + sum(D[i,,j])
dim(N) <- c(N_age, N_locations)
################################################################################


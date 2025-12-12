################################################################################
### Global metapopulation model: stochastic version ############################
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

################################################################################
### Time step (for odin.dust discrete-time model) ##############################
################################################################################
dt <- user(1)              # length of one time step (in days)
initial(time) <- 0
update(time) <- time + dt

################################################################################
################################################################################
### S: susceptible #############################################################
################################################################################
################################################################################

### Note: need to check what to do given competing risks of being vaccinated vs transition to E etc

## Setting up initial conditions
S_0[,,] <- user()
dim(S_0) <- c(N_age, N_vaccine,N_locations)
initial(S[,,]) <- S_0[i,j,k]
dim(S) <- c(N_age, N_vaccine,N_locations)

################################################################################
## 1) Probabilities of transitions out of S
################################################################################
## Infection probability for each (age, vaccine, location) over dt
p_S_inf[,,] <- 1 - exp(-lambda[i, k] * vaccine_efficacy_infection[i, j] * dt)
dim(p_S_inf) <- c(N_age, N_vaccine, N_locations)

## Vaccination probability for unvaccinated susceptibles (j = 1) - i.e. note this only applies vaccination stratum 1
p_S_vacc[ , ] <- 1 - exp(-vr[j] * vaccination_target[i, j] * dt) ## consider making this not stochastic
dim(p_S_vacc) <- c(N_age, N_locations)

## Progression probability between vaccine strata (j -> j+1)
p_S_prog[,,] <- 1 - exp(-gamma_vaccine[j] * dt)
dim(p_S_prog) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 2) Draw transition counts (binomial approximations)
################################################################################
## Infection from all S strata
n_S_inf[,,] <- rbinom(S[i, j, k], p_S_inf[i, j, k])
dim(n_S_inf) <- c(N_age, N_vaccine, N_locations)

## Vaccination from S[,1,] into S[,2,] - use remaining S[,1,] after infection to avoid double-counting
### have to make sure here that we don't inadvertently make this a negative number - perhaps a check for S[i, 1, k] - n_S_inf[i, 1, k] > 0
n_S_vacc[, ] <- rbinom(S[i, 1, j] - n_S_inf[i, 1, j], p_S_vacc[i, j])
dim(n_S_vacc) <- c(N_age, N_locations)

## Progression between vaccine strata j -> j+1 (j >= 2)
## Use remaining S after infection; gamma_vaccine[1] should be 0 so no S1->S2 progression
n_S_prog[,,] <- rbinom(S[i, j, k] - n_S_inf[i, j, k], p_S_prog[i, j, k])
dim(n_S_prog) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 3) Updates for S compartments
################################################################################
## j = 1: unvaccinated stratum (unvaccinated)
## Outflows: infection, vaccination
update(S[, 1, ]) <- S[i, 1, k] - n_S_inf[i, 1, k] - n_S_vacc[i, k]

## j = 2: vaccinated stratum 1 (vaccinated but not protected)
## Inflows: vaccination from S[,1,]
## Outflows: infection, progression to S[,3,]
update(S[, 2, ]) <- S[i, 2, k] - n_S_inf[i, 2, k] - n_S_prog[i, 2, k] + n_S_vacc[i, k]

## j = 3 : vaccinated stratum 2 (vaccinated and protected)
## Inflows: progression from previous stratum
## Outflows: infection, progression to next stratum
update(S[, 3, ]) <- S[i, 3, k] - n_S_inf[i, 3, k] - n_S_prog[i, 3, k] + n_S_prog[i, 2, k]

## j = 4 : vaccinated stratum 3 (vaccinated and protection has waned)
## Inflows: progression from previous stratum
## Outflows: infection only (no further vaccine progression)
update(S[, 4, ]) <- S[i, 4, k] - n_S_inf[i, 4, k] + n_S_prog[i, 3, k]

################################################################################
## 4) Outputs e.g. S_overall
################################################################################
initial(S_overall[, ]) <- 0
update(S_overall[, ]) <- sum(S[i,,j])
dim(S_overall) <- c(N_age, N_locations)

################################################################################
################################################################################
### E : Latent #################################################################
################################################################################
################################################################################

## Setting up initial conditions
E_0[,,] <- user()
dim(E_0) <- c(N_age, N_vaccine, N_locations)
initial(E[,,]) <- E_0[i, j, k]
dim(E) <- c(N_age, N_vaccine, N_locations)

## Parameters related to E
gamma_E <- user() # rate of progression through latent infection

################################################################################
## 1) Probabilities of transitions out of E
################################################################################

## Latent progression probability (E -> infectious states) over dt
p_E_latent[,,] <- 1 - exp(-gamma_E * dt)
dim(p_E_latent) <- c(N_age, N_vaccine, N_locations)

## Vaccination probability for latent E in stratum 1 (E[,1,] -> E[,2,])
## Same vaccination hazard as S[,1,]
p_E_vacc[, ] <- 1 - exp(-vr[j] * vaccination_target[i, j] * dt)
dim(p_E_vacc) <- c(N_age, N_locations)

## Vaccine-status progression probability between latent strata (j -> j+1)
## Uses the same gamma_vaccine vector as S
p_E_prog[,,] <- 1 - exp(-gamma_vaccine[j] * dt)
dim(p_E_prog) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 2) Draw transition counts (binomial approximations)
################################################################################

## Latent progression out of E (to IMild/ICase etc - to be consumed downstream)
n_E_latent[,,] <- rbinom(E[i, j, k], p_E_latent[i, j, k])
dim(n_E_latent) <- c(N_age, N_vaccine, N_locations)

## Vaccination from E[,1,] into E[,2,]
## Use remaining E[,1,] after latent progression to avoid double-counting
n_E_vacc[, ] <- rbinom(E[i, 1, j] - n_E_latent[i, 1, j], p_E_vacc[i, j])
dim(n_E_vacc) <- c(N_age, N_locations)

## Vaccine-status progression between latent strata j -> j+1 (j >= 2)
## Use remaining E after latent progression; gamma_vaccine[1] should be 0, so no E1->E2 via gamma_vaccine
n_E_prog[,,] <- rbinom(E[i, j, k] - n_E_latent[i, j, k], p_E_prog[i, j, k])
dim(n_E_prog) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 3) Updates for E compartments
################################################################################
## Remember: inflows from S are n_S_inf[i, j, k] from the S block

## j = 1: unvaccinated latent stratum
## Inflows: new infections from S[,1,]
## Outflows: latent progression, vaccination to E[,2,]
update(E[, 1, ]) <- E[i, 1, k] +
  n_S_inf[i, 1, k] -    # inflow from S1
  n_E_latent[i, 1, k] - # outflow to infectious states
  n_E_vacc[i, k]       # outflow to E2

## j = 2: latent, vaccinated but not yet protected
## Inflows: new infections from S[,2,], vaccination from E[,1,]
## Outflows: latent progression, vaccine-status progression to E[,3,]
update(E[, 2, ]) <- E[i, 2, k] +
  n_S_inf[i, 2, k] +    # inflow from S2
  n_E_vacc[i, k] -      # inflow from E1
  n_E_latent[i, 2, k] - # outflow to infectious states
  n_E_prog[i, 2, k]     # outflow to E3

## j = 3: latent, vaccinated and protected
## Inflows: new infections from S[,3,], vaccine-status progression from E[,2,]
## Outflows: latent progression, vaccine-status progression to E[,4,]
update(E[, 3, ]) <- E[i, 3, k] +
  n_S_inf[i, 3, k] +     # inflow from S3
  n_E_prog[i, 2, k] -    # inflow from E2
  n_E_latent[i, 3, k] -  # outflow to infectious states
  n_E_prog[i, 3, k]      # outflow to E4

## j = 4: latent, vaccinated with waned protection (last stratum)
## Inflows: new infections from S[,4,], vaccine-status progression from E[,3,]
## Outflows: latent progression only (no further vaccine progression)
update(E[, 4, ]) <- E[i, 4, k] +
  n_S_inf[i, 4, k] +     # inflow from S4
  n_E_prog[i, 3, k] -    # inflow from E3
  n_E_latent[i, 4, k]    # outflow to infectious states

################################################################################
## 4) Outputs e.g. E_overall
################################################################################
initial(E_overall[, ]) <- 0
update(E_overall[, ]) <- sum(E[i,,j])
dim(E_overall) <- c(N_age, N_locations)

################################################################################
################################################################################
### IMild: Unhospitalised infection ############################################
################################################################################
################################################################################

## Setting up initial conditions
IMild_0[,,] <- user()
dim(IMild_0) <- c(N_age, N_vaccine, N_locations)
initial(IMild[,,]) <- IMild_0[i, j, k]
dim(IMild) <- c(N_age, N_vaccine, N_locations)

## Parameters related to IMild
gamma_IMild <- user() # rate of progression from mild infection to recovery

################################################################################
## 1) Inflows from E: split latent exits into IMild vs ICase
################################################################################

## n_E_latent[i,j,k] was drawn in the E block:
##   n_E_toIMild  = n_E_latent - n_E_toICase
n_E_toIMild[,,] <- n_E_latent[i, j, k] - n_E_toICase[i, j, k]
dim(n_E_toIMild) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 2) Probabilities of transitions out of IMild
################################################################################

## Recovery probability from IMild over dt
p_IMild_rec[,,] <- 1 - exp(-gamma_IMild * dt)
dim(p_IMild_rec) <- c(N_age, N_vaccine, N_locations)

## Vaccine-status progression probability between IMild strata (j -> j+1)
## Uses the same gamma_vaccine as S and E
p_IMild_prog[,,] <- 1 - exp(-gamma_vaccine[j] * dt)
dim(p_IMild_prog) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 3) Draw IMild transition counts
################################################################################

## Recovery counts from IMild
n_IMild_rec[,,] <- rbinom(IMild[i, j, k], p_IMild_rec[i, j, k])
dim(n_IMild_rec) <- c(N_age, N_vaccine, N_locations)

## Vaccine-status progression counts j -> j+1 (j >= 2)
## Use remaining IMild after recovery; gamma_vaccine[1] should be 0, so no IMild1->IMild2 via gamma_vaccine
n_IMild_prog[,,] <- rbinom(IMild[i, j, k] - n_IMild_rec[i, j, k], p_IMild_prog[i, j, k])
dim(n_IMild_prog) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 4) Updates for IMild compartments
################################################################################
## j = 1: unvaccinated mild infections
## Inflows: n_E_toIMild from E[,1,]
## Outflows: recovery only
update(IMild[, 1, ]) <- IMild[i, 1, k] +
  n_E_toIMild[i, 1, k] -
  n_IMild_rec[i, 1, k]

## j = 2: mild, vaccinated but not protected
## Inflows: n_E_toIMild from E[,2,], progression from IMild[,1,]
## Outflows: recovery, progression to IMild[,3,]
update(IMild[, 2, ]) <- IMild[i, 2, k] +
  n_E_toIMild[i, 2, k] +
  n_IMild_prog[i, 1, k] -
  n_IMild_rec[i, 2, k] -
  n_IMild_prog[i, 2, k]

## j = 3: mild, vaccinated and protected
## Inflows: n_E_toIMild from E[,3,], progression from IMild[,2,]
## Outflows: recovery, progression to IMild[,4,]
update(IMild[, 3, ]) <- IMild[i, 3, k] +
  n_E_toIMild[i, 3, k] +
  n_IMild_prog[i, 2, k] -
  n_IMild_rec[i, 3, k] -
  n_IMild_prog[i, 3, k]

## j = 4: mild, vaccinated with waned protection (last stratum)
## Inflows: n_E_toIMild from E[,4,], progression from IMild[,3,]
## Outflows: recovery only (no further vaccine progression)
update(IMild[, 4, ]) <- IMild[i, 4, k] +
  n_E_toIMild[i, 4, k] +
  n_IMild_prog[i, 3, k] -
  n_IMild_rec[i, 4, k]

################################################################################
## 5) Outputs e.g. IMild_overall
################################################################################

initial(IMild_overall[, ]) <- 0
update(IMild_overall[, ]) <- sum(IMild[i,,j])
dim(IMild_overall) <- c(N_age, N_locations)

################################################################################
################################################################################
### R: Recovered ###############################################################
################################################################################
################################################################################

## Setting up initial conditions
R_0[,,] <- user()
dim(R_0) <- c(N_age, N_vaccine, N_locations)
initial(R[,,]) <- R_0[i, j, k]
dim(R) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 1) IHosp exits: split into recovered vs dead ################################
################################################################################
## NOTE: these will be used by both R (recovered) and D (dead)

## Total exits from IHosp
n_IHosp_exit[,,] <- rbinom(IHosp[i, j, k], p_IHosp_exit[i, j, k])
dim(n_IHosp_exit) <- c(N_age, N_vaccine, N_locations)

## Among exits, split into death vs recovery
n_IHosp_toD[,,] <- rbinom(n_IHosp_exit[i, j, k], prob_death_hosp[i, j])
dim(n_IHosp_toD) <- c(N_age, N_vaccine, N_locations)

n_IHosp_toR[,,] <- n_IHosp_exit[i, j, k] - n_IHosp_toD[i, j, k]
dim(n_IHosp_toR) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 2) Probabilities of transitions out of R ####################################
################################################################################

## Vaccination probability for R[,1,] -> R[,2,]  (same hazard as S[,1,])
p_R_vacc[, ] <- 1 - exp(-vr[j] * vaccination_target[i, j] * dt)
dim(p_R_vacc) <- c(N_age, N_locations)

## Vaccine-status progression probability between R strata (j -> j+1)
## Uses same gamma_vaccine as S, E, IMild
p_R_prog[,,] <- 1 - exp(-gamma_vaccine[j] * dt)
dim(p_R_prog) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 3) Draw R transition counts (binomial approximations) #######################
################################################################################

## Vaccination from R[,1,] into R[,2,]
n_R_vacc[, ] <- rbinom(R[i, 1, j], p_R_vacc[i, j])
dim(n_R_vacc) <- c(N_age, N_locations)

## Vaccine-status progression R[,j] -> R[,j+1] (j >= 2)
n_R_prog[,,] <- rbinom(R[i, j, k], p_R_prog[i, j, k])
dim(n_R_prog) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 4) Updates for R compartments ###############################################
################################################################################
## Inflows:
##  - from IMild: n_IMild_rec[i, j, k]  (defined in IMild block)
##  - from IHosp: n_IHosp_toR[i, j, k]  (defined above)

## j = 1: unvaccinated recovered
## Inflows: from IMild[,1,], from IHosp[,1,]
## Outflows: vaccination to R[,2,]
update(R[, 1, ]) <- R[i, 1, k] +
  n_IMild_rec[i, 1, k] +
  n_IHosp_toR[i, 1, k] -
  n_R_vacc[i, k]

## j = 2: recovered, vaccinated but not protected
## Inflows: from IMild[,2,], from IHosp[,2,], vaccination from R[,1,]
## Outflows: progression to R[,3,]
update(R[, 2, ]) <- R[i, 2, k] +
  n_IMild_rec[i, 2, k] +
  n_IHosp_toR[i, 2, k] +
  n_R_vacc[i, k] -
  n_R_prog[i, 2, k]

## j = 3: recovered, vaccinated and protected
## Inflows: from IMild[,3,], from IHosp[,3,], progression from R[,2,]
## Outflows: progression to R[,4,]
update(R[, 3, ]) <- R[i, 3, k] +
  n_IMild_rec[i, 3, k] +
  n_IHosp_toR[i, 3, k] +
  n_R_prog[i, 2, k] -
  n_R_prog[i, 3, k]

## j = 4: recovered, vaccinated with waned protection (last stratum)
## Inflows: from IMild[,4,], from IHosp[,4,], progression from R[,3,]
## Outflows: none (no further vaccine progression; gamma_vaccine[4] = 0)
update(R[, 4, ]) <- R[i, 4, k] +
  n_IMild_rec[i, 4, k] +
  n_IHosp_toR[i, 4, k] +
  n_R_prog[i, 3, k]

################################################################################
## 5) Outputs e.g. R_overall ###################################################
################################################################################
initial(R_overall[, ]) <- 0
update(R_overall[, ]) <- sum(R[i,,j])
dim(R_overall) <- c(N_age, N_locations)

################################################################################
################################################################################
### ICase: To-be hospitalised infection ########################################
################################################################################
################################################################################

## Setting up initial conditions
ICase_0[,,] <- user()
dim(ICase_0) <- c(N_age, N_vaccine, N_locations)
initial(ICase[,,]) <- ICase_0[i, j, k]
dim(ICase) <- c(N_age, N_vaccine, N_locations)

## Parameters related to ICase
gamma_ICase <- user() # rate of progression from symptom onset to requiring hospitalisation

################################################################################
## 1) Inflows from E: n_E_toICase
################################################################################

n_E_toICase[,,] <- rbinom(n_E_latent[i, j, k], prob_hosp[i, j])
dim(n_E_toICase) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 2) Probabilities of transitions out of ICase
################################################################################

## Progression from ICase -> IHosp over dt
p_ICase_toIHosp[,,] <- 1 - exp(-gamma_ICase * dt)
dim(p_ICase_toIHosp) <- c(N_age, N_vaccine, N_locations)

## Vaccine-status progression within ICase strata (j -> j+1)
p_ICase_prog[,,] <- 1 - exp(-gamma_vaccine[j] * dt)
dim(p_ICase_prog) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 3) Draw ICase transition counts
################################################################################

## Progression into IHosp
n_ICase_toIHosp[,,] <- rbinom(ICase[i, j, k], p_ICase_toIHosp[i, j, k])
dim(n_ICase_toIHosp) <- c(N_age, N_vaccine, N_locations)

## Vaccine-status progression j -> j+1 within ICase
## Use remaining ICase after progression to IHosp
n_ICase_prog[,,] <- rbinom(ICase[i, j, k] - n_ICase_toIHosp[i, j, k], p_ICase_prog[i, j, k])
dim(n_ICase_prog) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 4) Updates for ICase compartments
################################################################################
## Inflows: from E via n_E_toICase[i,j,k]
## Outflows: to IHosp (n_ICase_toIHosp), to next vaccine stratum (n_ICase_prog)

## j = 1: unvaccinated to-be-hospitalised infections
## Inflows: from E[,1,]
## Outflows: to IHosp; (gamma_vaccine[1] should be 0 so no ICase1->ICase2 via gamma_vaccine)
update(ICase[, 1, ]) <- ICase[i, 1, k] +
  n_E_toICase[i, 1, k] -
  n_ICase_toIHosp[i, 1, k]

## j = 2: to-be-hospitalised, vaccinated but not protected
## Inflows: from E[,2,], progression from ICase[,1,]
## Outflows: to IHosp, progression to ICase[,3,]
update(ICase[, 2, ]) <- ICase[i, 2, k] +
  n_E_toICase[i, 2, k] +
  n_ICase_prog[i, 1, k] -
  n_ICase_toIHosp[i, 2, k] -
  n_ICase_prog[i, 2, k]

## j = 3: to-be-hospitalised, vaccinated and protected
## Inflows: from E[,3,], progression from ICase[,2,]
## Outflows: to IHosp, progression to ICase[,4,]
update(ICase[, 3, ]) <- ICase[i, 3, k] +
  n_E_toICase[i, 3, k] +
  n_ICase_prog[i, 2, k] -
  n_ICase_toIHosp[i, 3, k] -
  n_ICase_prog[i, 3, k]

## j = 4: to-be-hospitalised, vaccinated with waned protection (last stratum)
## Inflows: from E[,4,], progression from ICase[,3,]
## Outflows: to IHosp only (no further vaccine progression; gamma_vaccine[4] = 0)
update(ICase[, 4, ]) <- ICase[i, 4, k] +
  n_E_toICase[i, 4, k] +
  n_ICase_prog[i, 3, k] -
  n_ICase_toIHosp[i, 4, k]

################################################################################
## 5) Outputs e.g. ICase_overall
################################################################################

initial(ICase_overall[, ]) <- 0
update(ICase_overall[, ]) <- sum(ICase[i,,j])
dim(ICase_overall) <- c(N_age, N_locations)

################################################################################
################################################################################
### IHosp: Go to hospital (stochastic, odin.dust) ##############################
################################################################################
################################################################################

## Setting up initial conditions
IHosp_0[,,] <- user()
dim(IHosp_0) <- c(N_age, N_vaccine, N_locations)
initial(IHosp[,,]) <- IHosp_0[i, j, k]
dim(IHosp) <- c(N_age, N_vaccine, N_locations)

## Parameters related to IHosp
gamma_IHosp <- user()  # rate of leaving IHosp (to R or D)

################################################################################
## 1) Probabilities of transitions out of IHosp
################################################################################

## Probability of exiting hospital (to R or D) over dt
p_IHosp_exit[,,] <- 1 - exp(-gamma_IHosp * dt)
dim(p_IHosp_exit) <- c(N_age, N_vaccine, N_locations)

## Vaccine-status progression probability within IHosp (j -> j+1)
p_IHosp_prog[,,] <- 1 - exp(-gamma_vaccine[j] * dt)
dim(p_IHosp_prog) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 2) Draw IHosp transition counts
################################################################################

## Vaccine-status progression within IHosp (j -> j+1)
## Use remaining IHosp after exit
n_IHosp_prog[,,] <- rbinom(IHosp[i, j, k] - n_IHosp_exit[i, j, k], p_IHosp_prog[i, j, k])
dim(n_IHosp_prog) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 3) Updates for IHosp compartments
################################################################################
## Inflows: from ICase via n_ICase_toIHosp[i,j,k]
## Outflows: exit to R/D (n_IHosp_exit), progression to next vaccine stratum (n_IHosp_prog)

## j = 1: unvaccinated hospitalised
## Inflows: n_ICase_toIHosp[,1,]
## Outflows: exit to R/D
update(IHosp[, 1, ]) <- IHosp[i, 1, k] +
  n_ICase_toIHosp[i, 1, k] -
  n_IHosp_exit[i, 1, k]

## j = 2: hospitalised, vaccinated but not protected
## Inflows: n_ICase_toIHosp[,2,], progression from IHosp[,1,]
## Outflows: exit to R/D, progression to IHosp[,3,]
update(IHosp[, 2, ]) <- IHosp[i, 2, k] +
  n_ICase_toIHosp[i, 2, k] +
  n_IHosp_prog[i, 1, k] -
  n_IHosp_exit[i, 2, k] -
  n_IHosp_prog[i, 2, k]

## j = 3: hospitalised, vaccinated and protected
## Inflows: n_ICase_toIHosp[,3,], progression from IHosp[,2,]
## Outflows: exit to R/D, progression to IHosp[,4,]
update(IHosp[, 3, ]) <- IHosp[i, 3, k] +
  n_ICase_toIHosp[i, 3, k] +
  n_IHosp_prog[i, 2, k] -
  n_IHosp_exit[i, 3, k] -
  n_IHosp_prog[i, 3, k]

## j = 4: hospitalised, vaccinated with waned protection (last stratum)
## Inflows: n_ICase_toIHosp[,4,], progression from IHosp[,3,]
## Outflows: exit to R/D only (no further vaccine progression; gamma_vaccine[4] = 0)
update(IHosp[, 4, ]) <- IHosp[i, 4, k] +
  n_ICase_toIHosp[i, 4, k] +
  n_IHosp_prog[i, 3, k] -
  n_IHosp_exit[i, 4, k]

################################################################################
## 4) Outputs e.g. IHosp_overall
################################################################################

initial(IHosp_overall[, ]) <- 0
update(IHosp_overall[, ]) <- sum(IHosp[i,,j])
dim(IHosp_overall) <- c(N_age, N_locations)

################################################################################
### D: Dead (stochastic, odin.dust) ###########################################
################################################################################

## Setting up initial conditions
D_0[,,] <- user()
dim(D_0) <- c(N_age, N_vaccine, N_locations)
initial(D[,,]) <- D_0[i, j, k]
dim(D) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 1) Vaccine-status progression probabilities within D
################################################################################

## Vaccine-status progression probability in D strata (j -> j+1)
## Uses same gamma_vaccine vector as S/E/IMild/R/IHosp
p_D_prog[,,] <- 1 - exp(-gamma_vaccine[j] * dt)
dim(p_D_prog) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 2) Draw D progression counts
################################################################################

## Vaccine-status progression D[,j] -> D[,j+1]
n_D_prog[,,] <- rbinom(D[i, j, k], p_D_prog[i, j, k])
dim(n_D_prog) <- c(N_age, N_vaccine, N_locations)

################################################################################
## 3) Updates for D compartments
################################################################################
## Inflows: from IHosp via n_IHosp_toD[i,j,k] (defined in IHosp block)
## Outflows: vaccine-status progression within D (conceptually possible, but
##           you may keep gamma_vaccine[last] = 0 so last stratum is absorbing)

## j = 1: unvaccinated deaths
## Inflows: deaths from IHosp[,1,]
## Outflows: progression to D[,2,] (if gamma_vaccine[1] > 0; normally 0)
update(D[, 1, ]) <- D[i, 1, k] +
  n_IHosp_toD[i, 1, k] -
  n_D_prog[i, 1, k]

## j = 2: deaths, vaccinated but not protected
## Inflows: deaths from IHosp[,2,], progression from D[,1,]
## Outflows: progression to D[,3,]
update(D[, 2, ]) <- D[i, 2, k] +
  n_IHosp_toD[i, 2, k] +
  n_D_prog[i, 1, k] -
  n_D_prog[i, 2, k]

## j = 3: deaths, vaccinated and protected
## Inflows: deaths from IHosp[,3,], progression from D[,2,]
## Outflows: progression to D[,4,]
update(D[, 3, ]) <- D[i, 3, k] +
  n_IHosp_toD[i, 3, k] +
  n_D_prog[i, 2, k] -
  n_D_prog[i, 3, k]

## j = 4: deaths, vaccinated with waned protection (last stratum)
## Inflows: deaths from IHosp[,4,], progression from D[,3,]
## Outflows: none if gamma_vaccine[4] = 0 (recommended), so this is absorbing
update(D[, 4, ]) <- D[i, 4, k] +
  n_IHosp_toD[i, 4, k] +
  n_D_prog[i, 3, k]

################################################################################
## 4) Outputs e.g. D_overall
################################################################################
initial(D_overall[, ]) <- 0
update(D_overall[, ]) <- sum(D[i,,j])
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

## Vaccination rate over time

## Per-day max vaccine capacity: rows = day, cols = location
max_vaccine_day[, ] <- user()
N_days <- user()
dim(max_vaccine_day) <- c(N_days, N_locations)

## Map continuous time -> day index (1-based), capped at N_days
day_index <- if (as.integer(time) >= N_days) N_days else as.integer(time) + 1L

## Current per-location capacity at this time step
mv[] <- max_vaccine_day[day_index, i]
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

## Defining the between-location mixing matrices
##   q_flight[i]     = fraction of time residents of location i spend away via flights
##   q_non_flight[i] = fraction of time residents of location i spend away via non-flight (gravity) travel
##   pi_travel_flight[i, j]      = probability (conditional on being away via flights from i)
##                                 that a trip from home i goes to destination j
##   pi_travel_non_flight[i, j]  = probability (conditional on being away via non-flight travel from i)
##                                 that movement from home i goes to destination j
## Require:
##   pi_travel_flight[i, i]     = 0, sum_j pi_travel_flight[i, j]     = 1 for each i with q_flight[i] > 0
##   pi_travel_non_flight[i, i] = 0, sum_j pi_travel_non_flight[i, j] = 1 for each i with q_non_flight[i] > 0
q_flight[] <- user()
dim(q_flight) <- N_locations
q_non_flight[] <- user()
dim(q_non_flight) <- N_locations
pi_travel_flight[,] <- user()
dim(pi_travel_flight) <- c(N_locations, N_locations)   # home (i), dest (j)
pi_travel_non_flight[,] <- user()
dim(pi_travel_non_flight) <- c(N_locations, N_locations)  # home (i), dest (j)

## Defining between-location mixing matrix (comprising a mixing matrix for between locations, initially assuming no age-specific mobility)

# Interpolation for beta ## Double check I HAVE THE MATRIX BETA_SET THE RIGHT WAY ROUND HERE I THINK BASED ON https://mrc-ide.github.io/odin/articles/odin.html#interpolating-functions
# beta schedule: rows = day (1..N_days), cols = location (1..N_locations)
beta_day[, ] <- user()
dim(beta_day) <- c(N_days, N_locations)
beta[] <- beta_day[day_index, i]
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

## 4a)  Travel FOI via flight-based mixing
## For each age i and *home* location j:
##   lambda_travel_flight[i, j] = sum_k pi_travel_flight[j, k] * lambda_local[i, k]
## i.e. average FOI over destinations, conditional on being away
travel_weighted_flight[,,] <- lambda_local[i, k] * pi_travel_flight[j, k]
dim(travel_weighted_flight) <- c(N_age, N_locations, N_locations)
lambda_travel_flight[,] <- sum(travel_weighted_flight[i, j, ])
dim(lambda_travel_flight) <- c(N_age, N_locations)

## 4b) Travel FOI via non-flight (gravity-based) mixing
## For each age i and *home* location j:
##   lambda_travel_non_flight[i, j] =
##       sum_k pi_travel_non_flight[j, k] * lambda_local[i, k]
## where pi_travel_non_flight[j, k] is a gravity-derived kernel (home j -> dest k)
travel_weighted_non_flight[,,] <- lambda_local[i, k] * pi_travel_non_flight[j, k]
dim(travel_weighted_non_flight) <- c(N_age, N_locations, N_locations)
lambda_travel_non_flight[,] <- sum(travel_weighted_non_flight[i, j, ])
dim(lambda_travel_non_flight) <- c(N_age, N_locations)

## 5) Total FOI: local + flight + non-flight components  ## NEW
## q_flight[j]     = fraction of time residents of j spend away via flights
## q_non_flight[j] = fraction of time residents of j spend away via non-flight travel
## Total away fraction for j: q_flight[j] + q_non_flight[j]
lambda[,] <- (1 - (q_flight[j] + q_non_flight[j])) * lambda_local[i, j] +
  q_flight[j] * lambda_travel_flight[i, j] +
  q_non_flight[j] * lambda_travel_non_flight[i, j]
dim(lambda) <- c(N_age, N_locations)

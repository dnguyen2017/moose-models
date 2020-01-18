# run regressions on simulated data, make vital rate functions, define functions to create matrices
source("structured_moose_calculations.R") 

# params for moose survival
cow.params <- c("surv.int" = unname(coef(surv_cow)[1]),
                  "surv.slope" = unname(coef(surv_cow)[2])
)
calf.params <- c("surv.int" = unname(coef(surv_calf)[1]),
                 "surv.slope" = unname(coef(surv_calf)[2])
)
# note, calves.multi is the multinomial reg model for the calving function

tick.params <- c("pEF" = 0.5,      # proportion of tick load that becomes engorged female
                 "t.surv" = 0.5,   # tick survival 
                 "l.surv" = 0.3,   # larval survival
                 "no.eggs" = 5000, # eggs per EF
                 "maxMu" = 50,     # maximum mean of the burden dist'n. In thousands of ticks.
                 "saturationRate" = 10^9 # how quickly mu(larvae) -> maxMu
)

### simulation model code

structured_moose_model <- function (initCalf, # init no. female calf moose
                                    initCow, # init no. female moose
                                    q0, # init no. questing larvae
                                    tfinal, # max number of years
                                    calf_surv,
                                    cow_surv,
                                    calving_model,
                                    tick_params,
                                    harvest, # annual proportion of moose harvested
                                    xmin = 0, # minimum burden
                                    xmax = 150, # maximum burden
                                    xsize = 150, # size of iteration matrix
                                    xincrement = 1000) { # increments of burden classes      
  
  # create meshpoints for mid-point rule approx
  # burden can take values within [L,U]
  # NOTE: bins are in 1000 tick increments
  L <- xmin   # lower limit of burden
  U <- xmax # upper limit of burden
  m <- xsize # size of the iteration matrix
  h <- (U-L)/m
  
  # meshpts are the midpoints of the 1000 increment bins used for for numerical integration
  # scaledpts are the true no. of ticks for each bin midpoint and is used for calc. of vital rates
  meshpts <- L + (1:m)*h - h/2
  scaledpts <- meshpts*xincrement
  
  # tick specific adjustments
  meshpts_ticks <- floor(meshpts) # since negative binomial only takes integers
  scaledpts_ticks <- scaledpts - (0.5 * xincrement)  # so that lowest burden class doesn't recruit larvae
  
  # init vectors to store data
  calfPop <- c(initCalf, rep(0, tfinal - 1))
  cowPop <- c(initCow, rep(0, tfinal - 1))
  larvaePop <- c(q0, rep(0, tfinal - 1))
  larvae_from_calf <- c(NA, rep(0, tfinal - 1)) 
  
  # save burden structure seperatly for calf and cow
  calfBurden <- matrix(0, nrow = tfinal, ncol = m) # after acquisition
  calfBurdenWK <- calfBurden                       # after winter kill
  
  cowBurden <- calfBurden
  cowBurdenWK <- calfBurden
  
  # create proj matrices (except for attachment matrix, make that in the simulation loop)
  survMat <- survival_matrix_str(scaledpts, 
                                 calf_parameters = calf_surv,
                                 cow_parameters = cow_surv)
  recrMat <- calving_matrix_str(scaledpts, calving_model)  
  
  poolMat <- pool_matrix_str(meshpts)
  
  larvaeMat <- larvae_matrix_str(scaledpts_ticks, tick_params) 
  
  # simulations
  t <- 1  
  while(t <= tfinal) {
    # mean burden (argument for attachment matrix)
    meanBurden <- getMu(larvaePop[t], tick_params["maxMu"], tick_params["saturationRate"])
    
    # initialize host popvec and attachment matrix
    pop_0 <- c(c(calfPop[t], rep(0, m - 1)), # calf
               c(cowPop[t], rep(0, m - 1)))  # cow
    aMat <- attach_matrix_str(meshpts_ticks, 
                                 size_calf = 3, mu_calf = meanBurden,
                                 size_cow = 3, mu_cow = meanBurden) # what if I didn't hardcode size = 3?
    
    # tick attachment
    pop_t1 <- aMat %*% pop_0            # note, lose a bit of the population since tickMat doesn't sum to 1
    # calfBurden[t, ] <- as.vector(pop_t1)[1:m]
    # cowBurden[t, ] <- as.vector(pop_t1)[m+1:(2*m)]
    
    # harvest
    pop_harvested <- (1 - harvest) * pop_t1 # harvest rate is independant of burden and all ticks on them die
    
    # burden dependent moose survival
    pop_t2 <- survMat %*% pop_harvested                    
    
    # calfBurdenWK[t, ] <- as.vector(pop_t2)[1:m]
    # cowBurdenWK[t, ] <- as.vector(pop_t2)[m+1:(2*m)]
    
    # questing larvae
    q_new <- larvaeMat %*% pop_t2  # pools into first element
    larvaePop[t + 1] <- sum(as.vector(q_new))     # save total larval recruitment from all moose stages
    larvae_from_calf[t + 1] <- as.vector(q_new)[1] # save calf contribution seperately too
    
    # burden dependent calf recruitment and maturation of last years calves into cow
    pop_t3 <- recrMat %*% pop_t2
    
    # pool moose classes
    pop_new <- as.vector(poolMat %*% pop_t3)  # pools into first element
    calfPop[t+1] <- pop_new[1]            # saves pop value of calves
    cowPop[t+1] <- pop_new[m + 1]            # saves pop value of cows
    
    t <- t + 1
    if (t %% 10 == 0) print(paste("time = ", t, "; aMat = ", sum(aMat), sep = ""))
  }
  
  # merge and return data matrices
  data_out <- tibble("t" = 0:tfinal,
                     "cows" = cowPop,
                     "calves" = calfPop,
                     "larvae" = larvaePop,
                     "larvae_calf" = larvae_from_calf,
                     "harvest" = harvest)
}

out <- structured_moose_model(initCalf = 50, # init no. female calf moose
                              initCow = 100, # init no. female moose
                              q0 = 10^6, # init no. questing larvae
                              tfinal = 100, # max number of years
                              calf_surv = calf.params,
                              cow_surv = cow.params,
                              calving_model = calves.multi,
                              tick_params = tick.params,
                              harvest = 0, # annual proportion of moose harvested
                              xmin = 0, # minimum burden
                              xmax = 150, # maximum burden
                              xsize = 150, # size of iteration matrix
                              xincrement = 1000)

ggplot(out) + 
  geom_line(aes(t, cows), col = "blue") +
  geom_line(aes(t, calves), col = "red")

ggplot(out) + 
  geom_line(aes(t, larvae), col = "blue") +
  geom_line(aes(t, larvae_calf), col = "red")
  

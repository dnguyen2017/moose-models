library(tidyverse)
library(Matrix)

### -----------------------------------------------------
### Regressions to estimate burden dependnent vital rates
### -----------------------------------------------------
# moose_data <- read_csv("data/simulated_moose_data.csv") # old cow only data
moose_data <- read_csv("data/simulated_structured_moose_data.csv") # calf-cow structured data

### survival ~ burden using logistic model 
surv_cow <- glm(survival ~ burden, family = binomial, data = filter(moose_data, stage == "cow"))
surv_calf <- glm(survival ~ burden, family = binomial, data = filter(moose_data, stage == "calf"))

# summary(surv_calf)

### number of calves ~ burden using multinomial regression
calves.multi <- nnet::multinom(calves ~ burden, 
                               data = filter(moose_data, stage == "cow")) # note, the filter is not strictly necessary, since all fecundity data for calves are NA

# summary(calves.multi)

## -------------------------------------------------
## define vital rate functions
## -------------------------------------------------

# note, vector of parameters is defined in moose_simulations.R

# survival rate function
surv.x = function(x, params){
  u = exp(params["surv.int"] + params["surv.slope"] * x)
  return(u/(1+u))
}

# calving function
# get probabilities using predict() function
# becuase of this the second argument is model, not params
recr.x <- function(x, model) {
  preds <- predict(model, tibble(burden = x), type = "prob") # return named vec of probabilities for "0", "1", or "2" calves
  return((1/2) * (2 * unname(preds["2"]) + unname(preds["1"]) )) # return expected number of female calves
}

# For plotting only
# same as above but takes vector argument for x 
recr.x.vector <- function(x, model) {
  preds <- predict(model, tibble(burden = x), type = "prob") # return named vec of probabilities for "0", "1", or "2" calves
  return((1/2) *(2 * unname(preds[,"2"]) + unname(preds[,"1"]) )) # return expected number of female calves
}


# tick recruitment function
tick.recr.x <- function(x, params) {
  # p(is engorged female) * p(survive) * p(larval survival) * E(number eggs per EF) * burden
  return(params["pEF"] * params["t.surv"] * params["l.surv"] * params["no.eggs"] * x) # could also estimate no. eggs from poisson distn
}

# calculate mu (the mean of the tick distribution) based on the abundance of questing larvae
# plug the output of getMu into attach_matrix during simulations
getMu <- function(qlarvae, max, satrate) {
  mu <- max * (qlarvae/(qlarvae + satrate))
  return(mu)
}


### -----------------------------------------------------
### calf-cow structured matrices
### -----------------------------------------------------

survival_matrix_str <- function(values,
                                calf_parameters,
                                cow_parameters){
  # dimension of final matrix
  len <- 2 * length(values)
  
  # calculate surv probability for midpoint of each bin
  surv_calf <- unname(sapply(values, function(x) surv.x(x, calf_parameters)))  
  surv_cow <- unname(sapply(values, function(x) surv.x(x, cow_parameters)))  
  
  # survival matrix has survival probabilities on the diagonal, zero elsewhere
  surv_mat <- sparseMatrix(i = 1:len, j = 1:len, 
                           x = c(surv_calf, surv_cow),
                           dims = list(len, len))  
  return(surv_mat)
}

# survival_matrix_str <- function(values,
#                                 calf_parameters,
#                                 cow_parameters){
#   # dims of matrix per stage
#   len1 <- length(values)
#   # dimension of final matrix
#   len2 <- len1 *2 
#   
#   # indices to fill with survivals
#   row_index <- c((len1+1):len2, # for calf survival into adult
#                  (len1+1):len2) # for cow survival
#   col_index <- 1:len2
#   
#   # calculate surv probability for midpoint of each bin
#   surv_calf <- unname(sapply(values, function(x) surv.x(x, calf_parameters)))  
#   surv_cow <- unname(sapply(values, function(x) surv.x(x, cow_parameters)))  
#   
#   # survival matrix has survival probabilities on the diagonal, zero elsewhere
#   surv_mat <- sparseMatrix(i = row_index, j = col_index, 
#                            x = c(surv_calf, surv_cow),
#                            dims = list(len2, len2))  
#   return(surv_mat)
# }


## calving matrix
calving_matrix_str <- function(values,      # calc. expected number of calves for these values
                               model) {     # regression model to calculate probabilities
  # dims of matrix per stage
  len1 <- length(values)
  # dimension of final matrix
  len2 <- len1 *2 
  
  # calculate expected no. calves for midpoint of each bin
  calving_x <- sapply(values, function(x) recr.x(x, model)) 
  calf_reproduction <- rep(0, length(values))
  
  # first row gets expected calves ~ burden since all calves are born into tick free class
  calving_mat <- sparseMatrix(i = rep(1, len2), j = 1:len2,   
                              x = c(calf_reproduction, calving_x), dims = list(len2,len2))
  
  
  # create maturation matrix: calf -> cow and cow -> cow
  # {{0, 0},{D1, D1}} <- a 2 X 2 block matrix where D1 denotes a diagonal matrix of 1s with dimension len x len
  # indices to fill with survivals
  row_index <- c((len1+1):len2, # for calf survival into adult
                 (len1+1):len2) # for cow survival
  col_index <- 1:len2
  
  maturation_mat <- sparseMatrix(i = row_index, j = col_index,
                           x = rep(1, len2), dims = list(len2,len2))
  return(calving_mat + maturation_mat)
}

# output looks right
calving_matrix_str(1:10,
                   calves.multi)

## matrix to pool all moose into lowest burden class (e.g., all ticks detached)
pool_matrix_str <- function(values){
  
  # dims of matrix per stage
  len1 <- length(values)
  # dimension of final matrix
  len2 <- len1 *2 
  
  # indices to fill with 1s
  row_index <- c(rep(1, len1), # pools calves
                 rep(len1 + 1, len1)) # pools cows
  col_index <- 1:len2
  
  # creates block diagonal matrix from two matrices of size len1 x len1
  # with 1s on the first row and zeros elsewhere
  pool_mat <- sparseMatrix(i = row_index, j = col_index, 
                           x = rep(1, len2), 
                           dims = list(len2, len2))
  return(pool_mat)
}

# works correctly
# pool_matrix_str(1:5) %*% c(1:5,1:5)

### tick matrices

# tick attachment matrix
# assume it all happens in one step and has negative binomial distn
attach_matrix_str <- function(values, 
                              size_calf, mu_calf, # burden dist pars for calf
                              size_cow, mu_cow) {# burden dist pars for cow
  
  # dims of matrix per stage
  len1 <- length(values)
  # dimension of final matrix
  len2 <- len1 * 2 
  
  # indices to fill with 1s
  row_index <- 1:len2
  col_index <- c(rep(1, len1), # calf burden prob
                 rep(len1 + 1, len1)) # cow burden prob
  
  
  # calculate probabilities of getting x burden
  calf_x <- dnbinom(values, size = size_calf, mu = mu_calf)
  cow_x <- dnbinom(values, size = size_cow, mu = mu_cow)
  
  # make attachment
  aMat <- sparseMatrix(i = row_index, j = col_index,
                       x = c(calf_x, cow_x),
                       dims = list(len2, len2))
  return(aMat)
}

# # works right
# attach_matrix_str(1:5,
#                   3, 5,
#                   3, 5) %*% 
#   c(c(5,0,0,0,0), c(5,0,0,0,0))

# larval recruitment matrix
# assumes that ticks on calves and cows have same recruitment

larvae_matrix_str <- function(values,    # use vector of true burdens (scaledpts) - 500 (so that lowest burden class doesn't produce larvae)
                              parameters) {
  # dims of matrix per stage
  len1 <- length(values)
  # dimension of final matrix
  len2 <- len1 *2 
  
  # indices to fill with 1s
  row_index <- c(rep(1, len1), # pools calves
                 rep(len1 + 1, len1)) # pools cows
  col_index <- 1:len2
  
  # calculate tick recruitment
  larvae_from_calf <- unname(sapply(values, function(x) tick.recr.x(x, parameters)))
  larvae_from_cow <- larvae_from_calf
  
  larvae_mat <- sparseMatrix(i = row_index, j = col_index, 
                             x = c(larvae_from_calf, larvae_from_cow),
                             dims = list(len2, len2))
  return(larvae_mat)
}

# works correctly
# larvae_matrix_str(1:5,
#                   tick.params)

# alternatively, can also just represent as a vector. remember to get use %*% not regular multiplication
larvae_vector_str <- function(values,    # use vector of true burdens (scaledpts) - 500 (so that lowest burden class doesn't produce larvae)
                              parameters) {
  larvae_x <- unname(sapply(values, function(x) tick.recr.x(x, parameters)))
  return(c(larvae_x, larvae_x))
}

### -----------------------------------------------------
### simulation model code
### -----------------------------------------------------

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
  meanCalfBurden <- rep(NA, tfinal)
  meanCowBurden <- rep(NA, tfinal)
  
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
    # print(as.vector(pop_t1)[1:m])
    calfBurden[t, ] <- as.vector(pop_t1)[1:m]
    cowBurden[t, ] <- as.vector(pop_t1)[(m+1):(2*m)]
    
    # calculate mean burdens pre-harvest (since this is what would be sampled in check station data)
    meanCalfBurden[t+1] <- (calfBurden[t, ] * scaledpts)/sum(calfBurden[t,])
    meanCowBurden[t+1] <- (cowBurden[t, ] * scaledpts)/sum(cowBurden[t,])
    
    # harvest
    pop_harvested <- (1 - harvest) * pop_t1 # harvest rate is independant of burden and all ticks on them die
    
    # burden dependent moose survival
    pop_t2 <- survMat %*% pop_harvested                    
    
    calfBurdenWK[t, ] <- as.vector(pop_t2)[1:m]
    cowBurdenWK[t, ] <- as.vector(pop_t2)[(m+1):(2*m)]
    
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
                     "moose_cows" = cowPop,
                     "moose_calves" = calfPop,
                     "larvae_cows" = larvaePop,
                     "larvae_calves" = larvae_from_calf,
                     "harvest" = harvest,
                     "moose_calvesBurden" = meanCalfBurden,
                     "moose_cowsBurden" = meanCowBurden)
}

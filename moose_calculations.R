library(tidyverse)

### -----------------------------------------------------
### Regressions to estimate burden dependnent vital rates
### -----------------------------------------------------
moose_data <- read_csv("data/simulated_moose_data.csv")

### survival ~ burden using logistic model 
surv.glm <- glm(survival ~ burden, family = binomial, data = moose_data)
# summary(surv.glm)

### number of calves ~ burden using multinomial regression
calves.multi <- nnet::multinom(calves ~ burden, data = moose_data)
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


## -------------------------------------------------
## define functions to make proj. matrices
## -------------------------------------------------

### moose matrices

## survival matrix

survival_matrix <- function(values,      # calculate survival probs for these values of the structure variable
                            parameters){ # parameter estimates from regressions
  
  # calculate surv probability for midpoint of each bin
  surv_x <- unname(sapply(values, function(x) surv.x(x, parameters)))  
  # survival matrix has survival probabilities on the diagonal, zero elsewhere
  surv_mat <- diag(surv_x)  
  return(surv_mat)
} 

## calving matrix
calving_matrix <- function(values,      # calc. expected number of calves for these values
                           model) {     # regression model to calculate probabilities
  len <- length(values)
  
  calving_x <- sapply(values, function(x) recr.x(x, model)) # calculate expected no. calves for midpoint of each bin
  calving_mat <- matrix(0, nrow = len, ncol = len)  # create empty recruitment matrix
  calving_mat[1, ] <- calving_x # first column gets expected calves ~ burden since all calves are born into tick free class
  calving_mat <- calving_mat + diag(rep(1, len)) # add diagonal of ones since all mortality accounted for in surv_mat
  
  return(calving_mat)
}


## matrix to pool all moose into lowest burden class (e.g., all ticks detached)
pool_matrix <- function(values){
  pool_mat <- matrix(0, nrow = length(values), ncol = length(values))
  pool_mat[1,] <- rep(1, length(values))
  return(pool_mat)
}


### tick matrices

# tick attachment matrix
# assume it all happens in one step and has negative binomial distn
attach_matrix <- function(values, size, mu) {
  
  aMat <- matrix(0, nrow = length(values), ncol = length(values))
  # since all moose start at 0, only first col matters
  aMat[,1] <- dnbinom(values, size = size, mu = mu) 
  return(aMat)
}

# calculate mu (the mean of the tick distribution) based on the abundance of questing larvae
# plug the output of getMu into attach_matrix during simulations
getMu <- function(qlarvae, max, satrate) {
  mu <- max * (qlarvae/(qlarvae + satrate))
  return(mu)
}

# larval recruitment matrix

larvae_matrix <- function(values,    # use vector of true burdens (scaledpts) - 500 (so that lowest burden class doesn't produce larvae)
                          parameters) {
  larvae_mat <- matrix(0, nrow = length(values), ncol = length(values))
  larvae_mat[1,] <- unname(sapply(values, function(x) tick.recr.x(x, parameters)))
  return(larvae_mat)
}

### simulation model code

moose_model <- function (n0, # init no. female moose
                         q0, # init no. questing larvae
                         tfinal, # max number of years
                         moose_surv,
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
  meshpts_ticks <- floor(meshpts) # since negadtive binomial only takes integers
  scaledpts_ticks <- scaledpts - (0.5 * xincrement)  # so that lowest burden class doesn't recruit larvae
  
  # init vectors to store data
  moosePop <- c(n0, rep(0, tfinal - 1))
  larvaePop <- c(q0, rep(0, tfinal - 1))
  burdenDist <- matrix(0, nrow = tfinal, ncol = m)
  burdenPreDemo <- burdenDist
  
  # create proj matrices (except for attachment matrix, make that in the simulation loop)
  survMat <- survival_matrix(scaledpts, moose_surv)
  recrMat <- calving_matrix(scaledpts, calving_model)  
  
  demoMat <- recrMat %*% survMat # moose survival and calving
  
  poolMat <- pool_matrix(meshpts)
  
  larvaeMat <- larvae_matrix(scaledpts_ticks, tick_params) 
  
  # simulations
  t <- 1  
  while(t <= tfinal) {
    # mean burden (argument for attachment matrix)
    meanBurden <- getMu(larvaePop[t], tick_params["maxMu"], tick_params["saturationRate"])
    
    # initialize host popvec and attachment matrix
    pop_0 <- c(moosePop[t], rep(0, m - 1))
    aMat <- attach_matrix(meshpts_ticks, size = 3, mu = meanBurden) # what if I didn't hardcode size = 3?
    
    # tick attachment
    pop_t1 <- aMat %*% pop_0            # note, lose a bit of the population since tickMat doesn't sum to 1
    burdenPreDemo[t, ] <- pop_t1
    
    # harvest
    pop_harvested <- (1 - harvest) * pop_t1 # harvest rate is independant of burden and all ticks on them die
    
    # burden dependent moose demography
    pop_t2 <- demoMat %*% pop_harvested                    
    burdenDist[t , ] <- pop_t2       # save burden distn
    
    # questing larvae
    q_new <- (larvaeMat %*% pop_t2)[1]  # pools into first element
    larvaePop[t + 1] <- q_new           # save larval value 
    
    # pool moose classes
    pop_new <- (poolMat %*% pop_t2)[1]  # pools into first element
    moosePop[t+1] <- pop_new            # saves pop value
    
    t <- t + 1
    if (t %% 10 == 0) print(paste("time = ", t, "; aMat = ", sum(aMat), sep = ""))
  }
  
  # merge and return data matrices
  data_out <- tibble("t" = 0:tfinal,
                     "moose" = moosePop,
                     "larvae" = larvaePop,
                     "harvest" = harvest)
}

### trash heap

# here, consider burden classes of 1000 ticks
## note! I've floored meshpts since nbinom needs integer values
#tick_mat <- attach_matrix(floor(meshpts), 3, 33)

# tick recruitment matrix
#larvae_mat <- matrix(0, nrow = length(meshpts), ncol = length(meshpts))
# use scaledpts - 500, so that no larvae are produced from lowest burden class
#larvae_mat[1,] <- unname(sapply(scaledpts - 500, function(x) tickRec(x, params.est)))


# ## -------------------------------------------------
# ## Create projection matrices
# ## -------------------------------------------------
# 
# # create meshpoints for mid-point rule approx
# # burden can take values within [L,U]
# # NOTE: bins are in 1000 tick increments
# L <- 0   # lower limit of burden
# U <- 150 # upper limit of burden
# m <- 150 # size of the iteration matrix
# h <- (U-L)/m
# 
# # meshpts are the midpoints of the 1000 increment bins used for for numerical integration
# # scaledpts are the true no. of ticks for each bin midpoint and is used for calc. of vital rates
# meshpts <- L + (1:m)*h - h/2 
# scaledpts <- meshpts*1000 
# 
# ### moose matrices
# 
# ## survival matrix
# 
# surv_x <- unname(sapply(scaledpts, function(x) surv.x(x, params.est)))  # calculate surv probability for midpoint of each bin
# surv_mat <- diag(surv_x)  # proj. matrix for survival has survival probabilities on the diagonal
# 
# ## calving matrix
# 
# recruit_x <- sapply(scaledpts, function(x) recr.x(x,calves.multi)) # calculate expected no. calves for midpoint of each bin
# recr_mat <- matrix(0, nrow = length(meshpts), ncol = length(meshpts))  # create empty recruitment matrix
# recr_mat[,1] <- recruit_x # first column gets expected calves ~ burden since all calves are born into tick free class
# recr_mat <- recr_mat + diag[rep(1, length(meshpts) )] # add diagonal of ones since all mortality accounted for in surv_mat
# 
# # matrix to pool all moose into lowest burden class (e.g., all ticks detached)
# pool_mat <- matrix(0, nrow = length(meshpts), ncol = length(meshpts))
# pool_mat[1,] <- rep(1, len)
# 
# 
# ### tick matrices
# 
# # here, consider burden classes of 1000 ticks
# ## note! I've floored meshpts since nbinom needs integer values
# tick_mat <- attachMat(floor(meshpts), 3, 33)
# 
# # tick recruitment matrix
# larvae_mat <- matrix(0, nrow = length(meshpts), ncol = length(meshpts))
# # use scaledpts - 500, so that no larvae are produced from lowest burden class
# larvae_mat[1,] <- unname(sapply(scaledpts - 500, function(x) tick.recr.x(x, params.est)))

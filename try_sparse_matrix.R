library(Matrix)
library(tidyverse)

# denseMat <- diag(c(1,1,1))
# sparseMat <- Matrix(diag(c(1,1,1)), sparse = TRUE)
# 
# denseMat %*% sparseMat
# 
# sparseMat + sparseMat
# 
# # can construct sparse matrices by specifying three vectors - i, j, and x - where (i,j) specifies the entry with non-zero value x.
# set.seed(123)
# size <- 10^3
# x <- rnorm(0, 1, n=size)
# denseMat <- diag(x)
# sparseMat <- sparseMatrix(i = 1:size, j = 1:size, x = x, dims = list(size, size))
# 
# object.size(denseMat)
# object.size(sparseMat)
# 
# 17504/8000216 * 100 # percent memory usage of sparse matrix compared to dense matrix


## -------------------------------------------------
## define functions to make proj. matrices
## -------------------------------------------------

### moose matrices

## survival matrix

survival_matrix_sparse <- function(values,      # calculate survival probs for these values of the structure variable
                            parameters){ # parameter estimates from regressions
  len <- length(values)
  
  # calculate surv probability for midpoint of each bin
  surv_x <- unname(sapply(values, function(x) surv.x(x, parameters)))  
  
  # survival matrix has survival probabilities on the diagonal, zero elsewhere
  surv_mat <- sparseMatrix(i = 1:len, j = 1:len, x = surv_x,
                           dims = list(len, len))  
  return(surv_mat)
} 

## calving matrix
calving_matrix_sparse <- function(values,      # calc. expected number of calves for these values
                           model) {     # regression model to calculate probabilities
  len <- length(values)
  # calculate expected no. calves for midpoint of each bin
  calving_x <- sapply(values, function(x) recr.x(x, model)) 
  # first column gets expected calves ~ burden since all calves are born into tick free class
  calving_mat <- sparseMatrix(i = rep(1, len), j = 1:len,   
                              x = calving_x, dims = list(len,len))
  # add diagonal of ones since all mortality accounted for in surv_mat
  diag_mat <- sparseMatrix(i = 1:len, j = 1:len,
                              x = rep(1, len), dims = list(len,len))
  return(calving_mat + diag_mat)
}

## matrix to pool all moose into lowest burden class (e.g., all ticks detached)
pool_matrix_sparse <- function(values){
  len <- length(values)
  pool_mat <- sparseMatrix(i = rep(1, len), j = 1:len, x = rep(1, len), dims = list(len, len))
  return(pool_mat)
}


### tick matrices

# tick attachment matrix
# assume it all happens in one step and has negative binomial distn
attach_matrix_sparse <- function(values, size, mu) {
  len <- length(values)
  aMat <- sparseMatrix(i = 1:len, j = rep(1, len),
                       x = dnbinom(values, size = size, mu = mu),
                       dims = list(len, len))
  return(aMat)
}

# larval recruitment matrix

larvae_matrix_sparse <- function(values,    # use vector of true burdens (scaledpts) - 500 (so that lowest burden class doesn't produce larvae)
                          parameters) {
  len <- length(values)
  larvae_x <- unname(sapply(values, function(x) tick.recr.x(x, parameters)))
  larvae_mat <- sparseMatrix(i = rep(1, len), j = 1:len, 
               x = larvae_x, dims = list(len, len))
  
  return(larvae_mat)
}
# alternatively, can also just represent as a vector. remember to get use %*% not regular multiplication
larvae_vector <- function(values,    # use vector of true burdens (scaledpts) - 500 (so that lowest burden class doesn't produce larvae)
                          parameters) {
  larvae_x <- unname(sapply(values, function(x) tick.recr.x(x, parameters)))
  return(larvae_x)
}

### simulation model code

moose_model_sparse <- function (n0, # init no. female moose
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
  survMat <- survival_matrix_sparse(scaledpts, moose_surv)
  recrMat <- calving_matrix_sparse(scaledpts, calving_model)  
  
  demoMat <- recrMat %*% survMat # moose survival and calving
  
  poolMat <- pool_matrix_sparse(meshpts)
  
  larvaeMat <- larvae_matrix_sparse(scaledpts_ticks, tick_params) 
  
  # simulations
  t <- 1  
  while(t <= tfinal) {
    # mean burden (argument for attachment matrix)
    meanBurden <- getMu(larvaePop[t], tick_params["maxMu"], tick_params["saturationRate"])
    
    # initialize host popvec and attachment matrix
    pop_0 <- c(moosePop[t], rep(0, m - 1))
    aMat <- attach_matrix_sparse(meshpts_ticks, size = 3, mu = meanBurden) # what if I didn't hardcode size = 3?
    
    # tick attachment
    pop_t1 <- aMat %*% pop_0            # note, lose a bit of the population since tickMat doesn't sum to 1
    burdenPreDemo[t, ] <- as.vector(pop_t1)
    
    # harvest
    pop_harvested <- (1 - harvest) * pop_t1 # harvest rate is independant of burden and all ticks on them die
    
    # burden dependent moose demography
    pop_t2 <- demoMat %*% pop_harvested                    
    burdenDist[t , ] <- as.vector(pop_t2)    # save burden distn
    
    # questing larvae
    q_new <- (larvaeMat %*% pop_t2)[1]  # pools into first element
    larvaePop[t + 1] <- as.vector(q_new)           # save larval value 
    
    # pool moose classes
    pop_new <- (poolMat %*% pop_t2)[1]  # pools into first element
    moosePop[t+1] <- as.vector(pop_new)            # saves pop value
    
    t <- t + 1
    if (t %% 10 == 0) print(paste("time = ", t, "; aMat = ", sum(aMat), sep = ""))
  }
  
  # merge and return data matrices
  data_out <- tibble("t" = 0:tfinal,
                     "moose" = moosePop,
                     "larvae" = larvaePop,
                     "harvest" = harvest)
}



# ### -----------------------------------------------------
# ### Regressions to estimate burden dependnent vital rates
# ### -----------------------------------------------------
# moose_data <- read_csv("data/simulated_moose_data.csv")
# 
# ### survival ~ burden using logistic model 
# surv.glm <- glm(survival ~ burden, family = binomial, data = moose_data)
# # summary(surv.glm)
# 
# ### number of calves ~ burden using multinomial regression
# calves.multi <- nnet::multinom(calves ~ burden, data = moose_data)
# # summary(calves.multi)
# 
# ## -------------------------------------------------
# ## define vital rate functions
# ## -------------------------------------------------
# 
# # note, vector of parameters is defined in moose_simulations.R
# 
# # survival rate function
# surv.x = function(x, params){
#   u = exp(params["surv.int"] + params["surv.slope"] * x)
#   return(u/(1+u))
# }
# 
# # calving function
# # get probabilities using predict() function
# # becuase of this the second argument is model, not params
# recr.x <- function(x, model) {
#   preds <- predict(model, tibble(burden = x), type = "prob") # return named vec of probabilities for "0", "1", or "2" calves
#   return((1/2) * (2 * unname(preds["2"]) + unname(preds["1"]) )) # return expected number of female calves
# }
# 
# # For plotting only
# # same as above but takes vector argument for x 
# recr.x.vector <- function(x, model) {
#   preds <- predict(model, tibble(burden = x), type = "prob") # return named vec of probabilities for "0", "1", or "2" calves
#   return((1/2) *(2 * unname(preds[,"2"]) + unname(preds[,"1"]) )) # return expected number of female calves
# }
# 
# 
# # tick recruitment function
# tick.recr.x <- function(x, params) {
#   # p(is engorged female) * p(survive) * p(larval survival) * E(number eggs per EF) * burden
#   return(params["pEF"] * params["t.surv"] * params["l.surv"] * params["no.eggs"] * x) # could also estimate no. eggs from poisson distn
# }
# 
# 
# # calculate mu (the mean of the tick distribution) based on the abundance of questing larvae
# # plug the output of getMu into attach_matrix during simulations
# getMu <- function(qlarvae, max, satrate) {
#   mu <- max * (qlarvae/(qlarvae + satrate))
#   return(mu)
# }
# 
# 
# # params for moose survival
# moose.params <- c("surv.int" = unname(coef(surv.glm)[1]),
#                   "surv.slope" = unname(coef(surv.glm)[2])
# )
# # note, calves.multi is the multinomial reg model for the calving function
# 
# tick.params <- c("pEF" = 0.5,      # proportion of tick load that becomes engorged female
#                  "t.surv" = 0.5,   # tick survival 
#                  "l.surv" = 0.3,   # larval survival
#                  "no.eggs" = 5000, # eggs per EF
#                  "maxMu" = 50,     # maximum mean of the burden dist'n. In thousands of ticks.
#                  "saturationRate" = 10^9 # how quickly mu(larvae) -> maxMu
# )
# 
# # check if output makes
# inputs <- seq(0, 10, 1)
# 
# survival_matrix_sparse(inputs, moose.params) # survival along diagonal
# 
# calving_matrix_sparse(inputs, calves.multi) # recruitment on first row + diag of 1s
# 
# pool_matrix_sparse(inputs) # first row of 1s
# 
# attach_matrix_sparse(inputs, 3, 5) # negbinom probs on first column
# 
# larvae_matrix_sparse(inputs, tick.params)
# larvae_vector(inputs, tick.params)

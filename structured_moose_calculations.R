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

# output looks right
# survival_matrix_str(1:10,
#                     calf_parameters = calf.params,
#                     cow_parameters = moose.params)


## calving matrix
calving_matrix_str <- function(values,      # calc. expected number of calves for these values
                               model) {     # regression model to calculate probabilities
  # dimension of final matrix
  len <- 2 * length(values)
  
  # calculate expected no. calves for midpoint of each bin
  calving_x <- sapply(values, function(x) recr.x(x, model)) 
  calf_reproduction <- rep(0, length(values))
  
  # first row gets expected calves ~ burden since all calves are born into tick free class
  calving_mat <- sparseMatrix(i = rep(1, len), j = 1:len,   
                              x = c(calf_reproduction, calving_x), dims = list(len,len))
  # add diagonal of ones since all mortality accounted for in surv_mat
  diag_mat <- sparseMatrix(i = 1:len, j = 1:len,
                           x = rep(1, len), dims = list(len,len))
  return(calving_mat + diag_mat)
}

# output looks right
# calving_matrix_str(1:10,
#                    calves.multi)

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

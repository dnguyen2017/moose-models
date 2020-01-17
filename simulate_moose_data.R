## -------------------------------------
## create hypothetical data set for moose-tick ipm 
## burden structure only (no moose stage structure)
## -------------------------------------
library(tidyverse)
set.seed(1)
# samuel (2004, pg 34): average of 33,000 ticks per hide of 214 moose from western Canada. 
burdens <- rnbinom(n = 200, size = 3, mu = 33000)

# Hist presented in Samuel 2004 looks similar to mine.
hist(burdens, breaks = seq(0, 140000, 10000))
quantile(burdens)

logEqn <- function(x, p0, beta1){
  # p0 == max probability, i.e., probability == p0 when x = 0
  # to get beta0, solve p0 == exp(beta0)/(1+exp(beta0)) -> beta0 == log(1/(1-p0)) 
  beta0 <- log(1/(1-p0)) 
  exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))
}

# survival data: P(survive|burden = x)
p0_s = 0.9 # P(survive|burden = 0)
# reasonable results using beta1 == 0.00005 (when burdens <- rbinom(n = 200, size = 3, mu = 33000))
beta1_s <- -0.00005

prob_surv <- logEqn(burdens, p0_s, beta1_s)
range(prob_surv)

# generate survival outcomes
survival <- rbinom(n=length(burdens), size=1, prob=prob_surv)
data <- data.frame(burdens, prob_surv, survival)
names(data) <- c("burden", "psurvival", "survival")

# fecundity data: P(calves = y|burden = x), y[0, 1.4]
# note: should really be P(calf = 0, 1, or 2| survive and burden = x), but that is too complicated for now
# P(calf = x | burden = 0)
p0_f2 <- 0.2 # calf = 2
p0_f1 <- 0.79 # calf = 1
p0_f0 <- 1 - p0_f2 - p0_f1 # calf = 0

# effect of burden on calf production
beta1_f2 <- -0.0001  # calf = 2
beta1_f1 <- -0.00005  # calf = 1

# calculate the probabilities of 0, 1, or 2 calves using log link
prob_twin <- logEqn(burdens, p0_f2, beta1_f2)
prob_single <- logEqn(burdens, p0_f1, beta1_f1)
prob_none <- 1 - (prob_twin + prob_single)

# should sum to total number of moose = 200
sum(prob_twin, prob_single, prob_none)

data$ptwin <- prob_twin
data$psingle <- prob_single
data$pnone <- prob_none

# problem with p(calf = 0| burden = 0) < 0
# ad hoc fix: convert all negative probs to zero
data[prob_none < 0,]
data[prob_none < 0,]$pnone <- 0

# generate calving outcomes: calves == 2, 1, 0 for twins, singletons, and none respectivelly
# Procedure:
# apply the burden specific probabilities for twinning, singleton, and no calvs to rmultinom
# rmultinom returns a vector with length == number of probabilities supplied (i.e., 3)
# there is a 1 in the returned vector indicating which outcome happened.
# so, return the index of the event using which.max which returns values 1:3, then subtract 1 to get the number of calves 0:2
data$calves <- sapply(1:200, function(x) which.max(rmultinom(n = 1, size = 1, prob = c(data$pnone[x],data$psingle[x],data$ptwin[x]))) - 1)

# save simulated data
# write_csv(data,"data/simulated_moose_data.csv")

### add calf winter-survival data
calf_burdens <- rnbinom(n = 200, size = 3, mu = 33000)

s_calf_0 <- 0.8 # w/i ballpark of High winter survival reported in Table 3-5 of Henry Jones' thesis
beta1_s_calf <- -0.000075

prob_surv_calf <- logEqn(calf_burdens, s_calf_0, beta1_s_calf)
range(prob_surv)

# compare survival probs
# plot(prob_surv ~ burdens, ylim = c(0,1), col = "red")
# points(prob_surv_calf ~ calf_burdens, ylim = c(0,1), col = "blue")


# make df of calf data
calf_data <- tibble(burden = calf_burdens,
                    psurvival = prob_surv_calf,
                    survival = rbinom(n=length(calf_burdens), size=1, prob=prob_surv_calf), # generate survival outcomes
                    ptwin = 0,
                    psingle = 0,
                    pnone = 0,
                    calves = 0,
                    stage = "calf"
                    )

# add col to data for stage
data$stage <- "cow"
# merge data
structured_data <- rbind(data, calf_data)

# save simulated stage-structured data
write_csv(structured_data,"data/simulated_structured_moose_data.csv")

#ggplot(structured_data, aes(x = burden, y = psurvival)) + geom_point() + facet_wrap(~stage, nrow = 2)

## save true vital rates
## used for evaluating regression model estimates
inputs <- 1:140000

# true survival probabilities
prob_surv_cow <- logEqn(inputs, p0_s, beta1_s)
prob_surv_calf <- logEqn(inputs, s_calf_0, beta1_s_calf)

# true calving probabilities
prob_twin <- logEqn(inputs, p0_f2, beta1_f2)
prob_single <- logEqn(inputs, p0_f1, beta1_f1)
prob_none <- 1 - (prob_twin + prob_single)

true_probs <- 
  tibble(burden = inputs,
         cow_survival = prob_surv_cow,
         calf_survival = prob_surv_calf,
         prob_twin = prob_twin,
         prob_single = prob_single,
         prob_none = prob_none) %>%
  # fix negative calving probabilities
  mutate(prob_none = ifelse(prob_none < 0, 0, prob_none))
  
true_probs <-
  true_probs %>%
  # normalize calving probabilities
  mutate(total_prob = prob_twin + prob_single + prob_none,
         prob_twin = prob_twin/total_prob,
         prob_single = prob_single/total_prob,
         prob_none = prob_none/total_prob,
         new_total_prob = prob_twin + prob_single + prob_none)
# # confirm that calving probs are correctly normalized. important for making comparisons to estimated probs later on
# sum(near(true_probs$new_total_prob, 1)) 

# save df of true_probs
true_probs %>%
  select(-c(total_prob, new_total_prob)) %>%
write_csv("data/true_moose_probabilities.csv")

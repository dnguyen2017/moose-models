### harvest sensitivity analysis

# run regressions on simulated data, make vital rate functions, define functions to create matrices
source("moose_calculations.R") 

# params for moose survival
moose.params <- c("surv.int" = unname(coef(surv.glm)[1]),
                  "surv.slope" = unname(coef(surv.glm)[2])
)
# note, calves.multi is the multinomial reg model for the calving function

tick.params <- c("pEF" = 0.5,      # proportion of tick load that becomes engorged female
                 "t.surv" = 0.5,   # tick survival 
                 "l.surv" = 0.3,   # larval survival
                 "no.eggs" = 5000, # eggs per EF
                 "maxMu" = 50,     # maximum mean of the burden dist'n. In thousands of ticks.
                 "saturationRate" = 10^9 # how quickly mu(larvae) -> maxMu
)

# analysis of harvesting rates
harvest_rates <- seq(0, 1, 0.1)

# simulate for different harvest rates
h_sensitivity <- lapply(harvest_rates, function(x) moose_model(n0 = 100,
                                                               q0 = 10^6,
                                                               tfinal = 100,
                                                               moose_surv = moose.params,
                                                               calving_model = calves.multi,
                                                               tick_params = tick.params,
                                                               harvest = x))
# extract output of harvest rate analysis
h_out <- bind_rows(h_sensitivity)

# make new df with larvae and moose abundances gathered so that I can make faceted plots
h_out <- 
  h_out %>% 
  pivot_longer(c("moose", "larvae"),
               names_to = "spp",
               values_to = "abundance")

# save harvest analysis
write_csv(h_out,"data/harvest_sensitivity_output.csv")

# # facet plots
# # when ticks are wiped out
# h_controlled <- c(0.4, 0.5)
# h_long %>% 
#   filter(round(harvest,1) %in% h_controlled) %>% # need to round harvest variable due to numerical differences (not sure why they popped up)
#   ggplot() + geom_point(aes(x = t, y = abundance)) +
#   facet_grid(spp ~ harvest, scales = "free")
# 
# # when moose are regulated by ticks
# h_regulated <- 0.0 #c(0.2, 0.3)
# h_long %>% 
#   filter(round(harvest,1) %in% h_regulated) %>% # need to round harvest variable due to numerical differences (not sure why they popped up)
#   ggplot() + geom_point(aes(x = t, y = abundance)) +
#   facet_grid(spp ~ harvest, scales = "free")

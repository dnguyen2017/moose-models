### harvest sensitivity analysis

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

# analysis of harvesting rates
harvest_rates <- seq(0, 1, 0.1)

# simulate for different harvest rates
h_sensitivity <- lapply(harvest_rates, function(x) structured_moose_model(initCalf = 50, # init no. female calf moose
                                                                          initCow = 100, # init no. female moose
                                                                          q0 = 10^6, # init no. questing larvae
                                                                          tfinal = 100, # max number of years
                                                                          calf_surv = calf.params,
                                                                          cow_surv = cow.params,
                                                                          calving_model = calves.multi,
                                                                          tick_params = tick.params,
                                                                          harvest = x, # annual proportion of moose harvested
                                                                          xmin = 0, # minimum burden
                                                                          xmax = 150, # maximum burden
                                                                          xsize = 150, # size of iteration matrix
                                                                          xincrement = 1000))
# extract output of harvest rate analysis
h_out <- bind_rows(h_sensitivity)


# save harvest analysis
write_csv(h_out,"data/structured_harvest_sensitivity_output.csv")

# try plotting some stuff

# moose faceted by harvest
h_out %>%
  ggplot() +
  geom_line(aes(t, cows), col = "blue") +
  geom_line(aes(t, calves), col = "red") +
  facet_wrap(~harvest, scales = "free")
# larvae faceted by harvest
h_out %>%
  ggplot() +
  geom_line(aes(t, larvae), col = "blue") +
  geom_line(aes(t, larvae_calf), col = "red") +
  facet_wrap(~harvest, scales = "free")

# make new df with larvae and moose abundances gathered so that I can make faceted plots
# gather moose stages
h_new <-
  h_out %>%
  pivot_longer(c("cows", "calves"),
               names_to = "moose_stage",
               values_to = "abundance")
# line
h_new %>%
  ggplot() +
  geom_line(aes(t, abundance, col = moose_stage)) +
  facet_wrap(~harvest, scales = "free")
# area
h_new %>%
  ggplot() +
  geom_area(aes(t, abundance, fill = moose_stage)) +
  facet_wrap(~harvest, scales = "free")

# now gather tick sources
h_new <-
  h_out %>%
  pivot_longer(c("larvae", "larvae_calf"),
               names_to = "tick_source",
               values_to = "abundance")
# line
h_new %>%
  ggplot() +
  geom_line(aes(t, abundance, col = tick_source)) +
  facet_wrap(~harvest, scales = "free")

# area
h_new %>%
  ggplot() +
  geom_area(aes(t, abundance, fill = tick_source)) +
  facet_wrap(~harvest, scales = "free")

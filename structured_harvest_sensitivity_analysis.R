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

# reformat and plot data

# get tick contribution from each moose stage. 
# Note, I lose the init value of larvae, since the first value of larvae_calf is NA, so larvae[1] - larvae_calf[1] == NA
df <-
  h_out %>% 
  mutate(larvae_cows = larvae_cows - larvae_calves)

# make two new cols: spp [moose, larvae] & stage [cows, calves]
# this allows me to facet by species and then get stage distribution of moose population
# and the absolute contribution of each moose stage to questing tick recruitment
df <- 
  df %>%
  pivot_longer(cols = c("moose_cows","moose_calves","larvae_cows","larvae_calves"),
               names_sep = "_",
               names_to = c("spp", "stage"),
               values_to = "abundance")

# generate a plot for each harvesting rate
library(viridis)

p_list <-
  lapply(harvest_rates, function(x) 
  filter(df, near(harvest, x) & t < 50 ) %>%
         ggplot(mapping = aes(x = t, y = abundance, fill = stage)) +
         geom_area(colour = "white") +
         labs(title = paste(x*100, "% annual harvest")) +
         theme_bw() + scale_color_viridis(discrete = T, alpha = 0.8) +
         facet_wrap(~spp, scales = "free", nrow = 2)
)

# can now plot each individual harvest rate plot
p_list[[4]] + 
  scale_fill_viridis(discrete = T, alpha = c(1,0.5) ) # can specify level specific alpha values

filter(df, near(harvest, rate) & t < 20 ) %>%
  ggplot(mapping = aes(x = t, y = abundance, fill = stage)) +
  geom_area(colour = "white") +
  scale_fill_viridis(discrete = T, alpha = 0.5) +
  facet_grid(spp ~ harvest, scales = "free")

# try out methods fpr
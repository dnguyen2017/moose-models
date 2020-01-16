# compare sparse matrix implementation and dense matrix implementation
# sparse matrices implemented using Matrix package

source("moose_calculations.R") # dense matrix
source("try_sparse_matrix.R") # sparse matrix

## set model parameters

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

# simulation with dense matrices
out_dense <- moose_model(n0 = 100,
                   q0 = 10^6,
                   tfinal = 100,
                   moose_surv = moose.params,
                   calving_model = calves.multi,
                   tick_params = tick.params,
                   harvest = 0.2)

# simulation with sparse matrices
out_sparse <- moose_model_sparse(n0 = 100,
                         q0 = 10^6,
                         tfinal = 100,
                         moose_surv = moose.params,
                         calving_model = calves.multi,
                         tick_params = tick.params,
                         harvest = 0.2)

# outputs are the same
identical(out_dense, out_sparse)

# compare size of dense vs. sparse matrix
inputs <- seq(0, 1500)
dense_mat <- calving_matrix(inputs, calves.multi)
sparse_mat <- calving_matrix_sparse(inputs, calves.multi)

object.size(sparse_mat)[1]/object.size(dense_mat)[1]
object.size(sparse_mat)[1]

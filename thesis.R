
# Install packages:
install.packages("devtools")
devtools::install_github("hesim-dev/hesim")
devtools::install_github('chjackson/msm')
install.packages("msm")
install.packages("hesim")
install.packages("data.table") 
library("msm")
library("hesim")
library("data.table")
library(dplyr)

# import data
# getwd()
library(readr)
dataall_f <- read_csv("Desktop/Thesis code/model data/dataall_f.csv", 
                      col_types = cols(...1 = col_skip()))
dataall <- dataall_f
head(dataall)
dataall <- dataall %>% 
  mutate_at(c('female', 'age', 'patient_id', 'time', 'strategy_id', 'state_id'), as.numeric)
summary(dataall)

# define global variables, target patients, states, and strategies
n_patients <- 1000
patients <- data.table(
  patient_id = 1:n_patients,
  age = rnorm(n_patients, mean = 45, sd = 7),
  female = rbinom(n_patients, size = 1, prob = .51)
)

states <- data.table(
  state_id = c(1,2,3),
  state_name = c("Hyper","Acute","Stable")
)

strategies <- data.table(
  strategy_id = 1:3,
  strategy_name = c("SOC", "Box", "Box1"),
  soc = c(1,0,0),
  box = c(0,1,0),
  box1 = c(0,0,1)
)

hesim_dat <- hesim_data(
  strategies = strategies,
  patients = patients,
  states = states
  )

print(hesim_dat)

# parameterization
# initialize transition probability matrix

q <- rbind ( c(0.9762, 0.0113, 0,  0.0125),
             c(0, 0, 0.8944, 0.0156), 
             c(0, 0.375, 0.6125, 0.0125),
              c(0, 0, 0, 0)
             )

# fit multi-state markov model

# diffs <- aggregate(cbind(minDiff=time)~patient_id, FUN=function(x) min(diff(x)),data=dataall)
# dataall1 <- merge(dataall,diffs,by='patient_id',all.x=T)

msm_fit <- msm(state_id ~ time, subject = patient_id,
               data = subset(dataall, strategy_name == "SOC"),
               exacttimes = TRUE,
               covariates = ~ age + female, 
               qmatrix = q, gen.inits = FALSE,
               control = list(fnscale=4000),
               fixedpars = c(1,2,3,4,5,6,7,8)
               )

# get transition probability matrix
pmatrix.msm(msm_fit, t=1)   # 1 year

# set relative risk (RR) for interventions
params_rr <- list( 
  lrr_12_est = c(soc = log(1), box = log(.61), box1 = log(.8)),
  lrr_12_se = c(soc = 0, box = .02, box1 = .02),
  lrr_13_est = c(soc = log(1), box = log(.61), box1 = log(.8)),
  lrr_13_se = c(soc = 0, box = .02, box1 = .02),
  lrr_14_est = c(soc = log(1), box = log(.61), box1 = log(.8)),
  lrr_14_se = c(soc = 0, box = .02, box1 = .02)
  )

n_samples <- 1000 # set PSA sample number

params_rr_def <- define_rng({
  list(
    rr_12 = lognormal_rng(lrr_12_est, lrr_12_se),
    rr_13 = lognormal_rng(lrr_13_est, lrr_13_se),
    rr_14 = lognormal_rng(lrr_14_est, lrr_14_se)
    )}, n = n_samples)
params_rr_rng <- eval_rng(params_rr_def, params_rr)

# predicting transition intensity matrices for the PSA sample cohort
survmods_data <- expand(hesim_dat, by = c("strategies", "patients"))
transmod_data <- survmods_data

transmod_data <- transmod_data %>% 
  mutate_at(c('female', 'age', 'patient_id', 'strategy_id'), as.numeric)

summary(transmod_data)
qmat_soc <- qmatrix(msm_fit,
                    newdata = subset(transmod_data, strategy_name == "SOC"),
                      # transmod_data[strategy_name == "SOC"],
                    uncertainty = "normal", n = n_samples)
dim(qmat_soc)
qmat_soc[,, 1]

# transition probability matrices are derived from the transition intensity matrices with the matrix exponential
cycle_len <- 1
pmat_soc <- expmat(qmat_soc, t = cycle_len)
pmat_soc[,, 1]

# create identifiers for the array of transition matrices
tpmat_id <- tpmatrix_id(transmod_data,  n_samples)
head(tpmat_id)

# predict RR for each row of the input data and PSA samples
xbeta <- function(x, beta) c(x %*% t(beta))
x_rr <- as.matrix(transmod_data[, .(soc, box, box1)])
rr <- cbind(xbeta(x_rr, params_rr_rng$rr_12),
            xbeta(x_rr, params_rr_rng$rr_13),
            xbeta(x_rr, params_rr_rng$rr_14)
            )
head(rr)

pmat <- apply_rr(pmat_soc, rr = rr,
               index = list(c(1, 2), c(1,3), c(1,4)
                            ))
tprobs <- tparams_transprobs(pmat, tpmat_id)

# final disease model
transmod <- CohortDtstmTrans$new(params = tprobs,
                                 cycle_length = cycle_len)

# construct utility and cost model
utility_tbl <- stateval_tbl(
  data.table(state_id = states$state_id,
             mean = c(0.86, -0.133, 0.72),
             se = c(0.02, 0, 0.01)
  ),
  dist = "beta")
print(utility_tbl)

medcost_tbl <- stateval_tbl(
  data.table(state_id = states$state_id,
             mean = c(2150.3, 19430.65, 2961.06),
             se = c(6.32, 15.81, 6.32)
  ),
  dist = "gamma")
print(medcost_tbl)

drugcost_dt <- data.table(strategy_id = c(1,2,3,1,2,3,1,2,3),
                          state_id = c(1,1,1,2,2,2,3,3,3),
                          time_id = c(1,1,1,1,1,1,1,1,1),
                          time_start = c(0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00),
                          time_stop = c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf),
                          est = c(0,135,0,0,0,0,0,135,0)
                          )
drugcost_tbl <- stateval_tbl(
  drugcost_dt,
  dist = "fixed")
print(drugcost_tbl)

# final utility model
utilitymod <- create_StateVals(utility_tbl, n = n_samples,
                               hesim_data = hesim_dat)
# final cost model
drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples,
                                time_reset = TRUE, hesim_data = hesim_dat)
medcostmod <- create_StateVals(medcost_tbl, n = n_samples,
                               hesim_data = hesim_dat)
costmods <- list(Drug = drugcostmod,
                 Medical = medcostmod)

# health economic model integration
cdtstm <- CohortDtstm$new(trans_model = transmod,
                          utility_model = utilitymod,
                          cost_models = costmods)

# simulate state probability
cdtstm$sim_stateprobs(n_cycles = 30/cycle_len)

# simulate QALYs
cdtstm$sim_qalys(dr = .04)

# simulate costs
cdtstm$sim_costs(dr = .04)

# simulation results
ce_sim_cdtstm <- cdtstm$summarize()
labs_cohort <- get_labels(hesim_dat)
summary(ce_sim_cdtstm, labels = labs_cohort) %>%
  format()

# cost-effectiveness analysis
ce_sim_cdtstm

wtp <- seq(0, 60000, 500)
cea_cdtstm <- cea(ce_sim_cdtstm, dr_costs = .04, dr_qalys = .04, k = wtp)
cea_pw_cdtstm <- cea_pw(ce_sim_cdtstm, comparator = 1,
                        dr_qalys = .04, dr_costs = .04,
                        k = wtp)

# incremental cost-effectiveness ratio
icer(cea_pw_cdtstm, k = 53397.50, labels = labs_cohort) %>%
  format()

# uncertainty visualization
plot_ceplane(cea_pw_cdtstm, labels = labs_cohort)
plot_ceac(cea_cdtstm, labels = labs_cohort)
plot_ceac(cea_pw_cdtstm, labels = labs_cohort)
plot_ceaf(cea_cdtstm, labels = labs_cohort)
plot_evpi(cea_cdtstm, labels = labs_cohort)

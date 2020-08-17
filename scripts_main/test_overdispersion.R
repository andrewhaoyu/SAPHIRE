#test the potential overdispersion on Wuhan's data
code_root= "/n/holystore01/LABS/xlin/Lab/hzhang/SAPHIRE/"

setwd(paste0(code_root, "scripts_main"))
source(paste0(code_root, "R/init_cond.R"))
init_sets_list=get_init_sets_list(r0 = 0.23)
source(paste0(code_root, "R/fun_SEIRsimu.R"))
mcmc_pars_estimate=read.table("../output/pars_est_run_main_analysis.txt", header=T)

library(cairoDevice)
pars_estimate = mcmc_pars_estimate 
file_name = "main_analysis"
init_settings = init_sets_list
panel_B_R_ylim = 4
estN_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 5)[, "Onset_expect"])
estN <- round(apply(estN_mat, 1, mean), 0)
onset = init_sets_list$daily_new_case
plot(onset,estN)
outcome <- ((onset-estN)^2-onset)/estN
model <- lm(outcome~estN-1)
summary(model)

#possion model Deviance calculation
  logL <- sum(log(dpois(onset, estN)))
  logL_full <- sum(log(dpois(onset, onset)))
D <- 2*(logL_full-logL)

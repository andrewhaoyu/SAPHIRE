rm(list = ls())
## IMPORTANT: Please set code_root variable properly. 
## code_root should be set to the directory where the repository README file is located. 
## For more information, please read the repository README file
code_root="/n/holystore01/LABS/xlin/Lab/hzhang/SAPHIRE/"

setwd(paste0(code_root, "scripts_main"))
run_id = "main_analysis_nb"

library(readr)
library(coda)
library(cairoDevice)

# if (!file.exists("../output/pars_est_run_main_analysis.txt") | 
#     !file.exists("../output/pars_est_run_main_analysis_rep1.txt") |
#     !file.exists("../output/pars_est_run_main_analysis_rep2.txt")) {
#   stop("Outputs from main analysis cannot be found.\n
#       Probably scripts_main/Run_SEIR_main_analysis.R has not yet been run or code_root is not set correctly.")
# }

pars_estimate_main=read.table(paste0("../output/pars_est_run_",run_id,".txt"), header=T)
pars_estimate_main_rep1=read.table(paste0("../output/pars_est_run_",run_id,"_rep1.txt"), header=T)
pars_estimate_main_rep2=read.table(paste0("../output/pars_est_run_",run_id,"_rep2.txt"), header=T)

mcmc_main=mcmc(data=pars_estimate_main)
mcmc_rep1=mcmc(data=pars_estimate_main_rep1)
mcmc_rep2=mcmc(data=pars_estimate_main_rep2)

mcmc_3traj=mcmc.list(mcmc_main,mcmc_rep1,mcmc_rep2)

gelman.diag(mcmc_3traj)
# our results: multivariate psrf = 1
p_vec <- c(4001:14000)
plot_par_3traj = function(par_name, plotmath_name) {
  # red-like: #BC3C29, blue-like: #0072B5, orange-like: #E18727
  plot(p_vec, pars_estimate_main[p_vec, par_name], type="l", col="#0072B5", main=plotmath_name, xlab="", ylab="")
  points(p_vec, pars_estimate_main_rep1[p_vec, par_name], type="l", col="#BC3C29")
  points(p_vec, pars_estimate_main_rep2[p_vec, par_name], type="l", col="#E18727")
}

png(paste0("../output/mcmc_convergence_",run_id,".png"), width=10, height=5,res = 300, units = "in")
par(mfrow=c(2, 4), mar=c(2, 2, 3.2, 0.9))
plot_par_3traj("b12", expression(b[12]))
plot_par_3traj("b3", expression(b[3]))
plot_par_3traj("b4", expression(b[4]))
plot_par_3traj("b5", expression(b[5]))
plot_par_3traj("r12", expression(r[12]))
plot_par_3traj("delta3", expression(delta[3]))
plot_par_3traj("delta4", expression(delta[4]))
plot_par_3traj("delta5", expression(delta[5]))
dev.off()



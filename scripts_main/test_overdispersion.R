#test the potential overdispersion on Wuhan's data
#generate ascertainment rate for Wuhan's data
code_root= "/n/holystore01/LABS/xlin/Lab/hzhang/SAPHIRE/"

setwd(paste0(code_root, "scripts_main"))
source(paste0(code_root, "R/init_cond.R"))
init_sets_list=get_init_sets_list(r0 = 0.23)
source(paste0(code_root, "R/fun_SEIRsimu.R"))
mcmc_pars_estimate=read.table("../output/pars_est_run_main_analysis.txt", header=T)
init_settings = init_sets_list
init_settings$days_to_fit <- 1:68
library(vioplot)
transform_var_main_5stage=function(pars) {
  b_vec <- pars[1:4]
  b_vec <- c(b_vec[1], b_vec[1], b_vec[2:4])
  ##
  r12 <- pars[5]
  r3 <- 1 / (1 + (1 - r12) / (r12 * exp(pars[6])))
  r4 <- 1 / (1 + (1 - r3) / (r3 * exp(pars[7])))
  r5 <- 1 / (1 + (1 - r4) / (r4 * exp(pars[8])))
  r_vec <- c(r12,r12,r3,r4,r5)
  
  return(list(b_vec, r_vec))
}
r_mat <- mcmc_pars_estimate[,5:8]
colnames(r_mat) <- c("r12","r3","r4","r5")
for(l in 1:nrow(mcmc_pars_estimate)){
  r_mat[l,] <- transform_var_main_5stage(mcmc_pars_estimate[l,])[[2]][c(1,3,4,5)] 
}

r_est = colMeans(r_mat)
r_est_low = apply(r_mat,2,function(x){quantile(x,0.025)})
r_est_high = apply(r_mat,2,function(x){quantile(x,0.975)})

mean(as.matrix(r_mat))
quantile(as.matrix(r_mat),0.025)
quantile(as.matrix(r_mat),0.975)
analysis_result_table <- r_est
pla = 2
for(l in 1:length(r_est)){
  analysis_result_table[l] =  paste0(round(r_est[l],pla)," (",
         round(r_est_low[l],pla),
         "-",
         round(r_est_high[l],pla)
         ,")")
}
analysis_result_table_temp = analysis_result_table
#write.csv(analysis_result_table,file = "../output/analysis_result_table.csv")
# 
# 
# ##
# onset_obs_all <- init_settings$daily_new_case_all
# ptime <- 1:length(onset_obs_all)
# estAIP_mat <- apply(mcmc_pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 5)[, c("I", "A", "P")])
# estI_mat <- estAIP_mat[ptime, ]
# estA_mat <- estAIP_mat[ptime + length(ptime), ]
# estP_mat <- estAIP_mat[ptime + length(ptime) * 2, ]
# estI_mean <- apply(estI_mat, 1, mean)
# estA_mean <- apply(estA_mat, 1, mean)
# estP_mean <- apply(estP_mat, 1, mean)
# estAIP_dat <- rbind(estI_mean, estA_mean, estP_mean)
# library(lubridate)
# ascertainment_rate = estAIP_dat[1,]/(estAIP_dat[1,]+estAIP_dat[2,])
# mydate <- as.Date(c(paste("Jan", 1:31), paste("Feb", 1:29), paste("Mar", 1:8)),"%b%d")
# 
# data <- data.frame(date= mydate,ascertainment_rate)
# library(ggplot2)
# library(ggthemes)
# ggplot(data) + geom_point(aes(mydate,ascertainment_rate))+
#   theme_Publication()+
#   
# 
# plot(ascertainment_rate)
# barpos <- barplot(estAIP_dat, col = c("#BC3C29FF", "#FFDC91FF", "#0072B5FF"), xlab = "", ylab = "", border = "NA")
# mtext("Date (2020)", side = 1, line  = 3, cex = 1.01)
# mtext("No. of active infectious cases", side = 2, line = 3, cex = 1.01)
# axis(1, at = barpos[seq(1, 68, 11)], labels = mydate[seq(1, 68, 11)])
# legend("topleft", legend = c("Presymptomatic (P)", "Unascertained (A)", "Ascertained (I)"), fill = c("#0072B5FF", "#FFDC91FF", "#BC3C29FF"), bty = "n")
# 
# library(cairoDevice)
# pars_estimate = mcmc_pars_estimate 
# file_name = "main_analysis"
# init_settings = init_sets_list
# panel_B_R_ylim = 4
# estN_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 5)[, "Onset_expect"])
# estN <- round(apply(estN_mat, 1, mean), 0)
# onset = init_sets_list$daily_new_case
# plot(onset,estN)
# outcome <- ((onset-estN)^2-onset)/estN
# model <- lm(outcome~estN-1)
# summary(model)
# 
# #possion model Deviance calculation
#   logL <- sum(log(dpois(onset, estN)))
#   logL_full <- sum(log(dpois(onset, onset)))
# D <- 2*(logL_full-logL)








setwd(paste0(code_root, "scripts_main"))
source(paste0(code_root, "R/init_cond_new.R"))
init_sets_list=get_init_sets_list(r0 = 0.23)
source(paste0(code_root, "R/fun_SEIRsimu.R"))
mcmc_pars_estimate=read.table("../output/pars_est_run_main_analysis_nb.txt", header=T)
colMeans(mcmc_pars_estimate)
apply(mcmc_pars_estimate,2,function(x){quantile(x,0.025)})
apply(mcmc_pars_estimate,2,function(x){quantile(x,0.975)})
r_mat <- mcmc_pars_estimate[,5:8]
colnames(r_mat) <- c("r12","r3","r4","r5")
for(l in 1:nrow(mcmc_pars_estimate)){
  r_mat[l,] <- transform_var_main_5stage(mcmc_pars_estimate[l,])[[2]][c(1,3,4,5)] 
}

r_est = colMeans(r_mat)
r_est_low = apply(r_mat,2,function(x){quantile(x,0.025)})
r_est_high = apply(r_mat,2,function(x){quantile(x,0.975)})

mean(as.matrix(r_mat))
quantile(as.matrix(r_mat),0.025)
quantile(as.matrix(r_mat),0.975)
analysis_result_table <- r_est
pla = 2
for(l in 1:length(r_est)){
  analysis_result_table[l] =  paste0(round(r_est[l],pla)," (",
                                     round(r_est_low[l],pla),
                                     "-",
                                     round(r_est_high[l],pla)
                                     ,")")
}
analysis_result_table = rbind(analysis_result_table_temp,analysis_result_table)
write.csv(analysis_result_table,file = "../output/analysis_result_table_all.csv")




library(cairoDevice)
pars_estimate = mcmc_pars_estimate 
init_settings = init_sets_list
estN_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 5)[, "Onset_expect"])
estN <- round(apply(estN_mat, 1, mean), 0)
onset = init_sets_list$daily_new_case
plot(onset,estN)
outcome <- ((onset-estN)^2-onset)/estN
model <- lm(outcome~estN-1)
summary(model)
pars = colMeans(mcmc_pars_estimate)
phi = pars[length(pars)]
#p = phi/(phi+as.numeric(ypred))
# meant to suppress warnings when ypred is negative


#possion model Deviance calculation
logL <- sum(dnbinom(x = as.numeric(onset), 
                     size = phi,
                     mu = estN,log=T))
logL_full <- sum(dnbinom(x = as.numeric(onset), 
                         size = phi,
                         mu = onset,log=T))
D <- 2*(logL_full-logL)

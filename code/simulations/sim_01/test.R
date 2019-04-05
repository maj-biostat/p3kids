

suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rbenchmark))
suppressPackageStartupMessages(library(ggplot2))
library(lme4)

ggplot2::theme_set(ggplot2::theme_bw())
ggplot2::theme_update(text = element_text(size = 10))
ggplot2::theme_update(legend.position = "top")
# ggplot2::theme_update(legend.title = element_blank())
ggplot2::theme_update(axis.text.x = element_text(size = 10))
ggplot2::theme_update(axis.text.y = element_text(size = 10))

prob_to_odd <- cmpfun(function(x){
  return(x/(1-x))
})

inv_logit <- cmpfun(function(x){
  return(exp(x)/(1+exp(x)))
})



create_base <- cmpfun(function(hospts, n_clinic, t){
  # Set up the data 
  # repeat observations of participants within clinics within hospt
  id <- unlist(lapply(n_clinic, function(x) 1:x))
  
  dt <- expand.grid(t, 1:max(n_clinic))
  names(dt) <- c("time", "id")
  
  ds <- data.frame(hosp_id = 1:hospts, start_tx = rep(1:4, each = 2))
  
  d0 <- data.frame(id = rep(id, length = hospts * sum(n_clinic)))
  d0$hosp_id <- rep(1:hospts, each = sum(n_clinic))
  d0$clinic_id <- rep(unlist(lapply(1:length(n_clinic), 
                                    function(x) rep(x, n_clinic[x]))), 
                      len = nrow(d0))
  
  d0 <- merge(d0, dt, by = "id")
  d0 <- merge(d0, ds, by = "hosp_id")
  d0$tx_active <- as.numeric(d0$time >= d0$start_tx)
  
  # calendar trend
  d0$k1 <- as.numeric(d0$time == 1)
  d0$k2 <- as.numeric(d0$time == 2)
  d0$k3 <- as.numeric(d0$time == 3)
  d0$k4 <- as.numeric(d0$time == 4)
  
  # exposure trend
  d0$tx_durn <- ave(d0$tx_active, d0$hosp_id, d0$clinic_id, d0$id, FUN=cumsum)
  d0$tx1 <- as.numeric(d0$tx_durn == 1)
  d0$tx2 <- as.numeric(d0$tx_durn == 2)
  d0$tx3 <- as.numeric(d0$tx_durn == 3)
  d0$tx4 <- as.numeric(d0$tx_durn == 4)
  
  d0 <- d0[order(d0$hosp_id, d0$clinic_id, d0$id, d0$time),
           c("hosp_id", "clinic_id", "id", "time", "tx_active", "tx_durn",
             paste0("k", 1:4), paste0("tx", 1:4))]

  rownames(d0) <- NULL
  d0
})





# # total number of hospitals
# hospts <- 8
# n_clinic <- c(3, 2, 4)
# # Baseline probs of vac by clinic
# p0 <- c(0.7, 0.3, 0.05)
# # Timepoints
# t <- 0:4

gendat <- cmpfun(function(){
  
  # total number of hospitals
  hospts <- 8
  # Number of patients within each clinic
  # length(n) = number of clinics
  # Cardiology, Respitarory, Immunosuppressed, Neuro, Other
  n_clinic <- c(50, 400, 200, 200, 100)
  # Baseline probs of vac by clinic
  p0 <- c(0.45, 0.3, 0.05, 0.35, 0.2)
  # Timepoints
  t <- 0:4
  
  sig_hosp <- 0.05
  sig_clinic <- 0.2
  sig_person <- 0.5
  
  d0 <- create_base(hospts, n_clinic, t)
  
  gensub <- function(x, hosp_int, clinic_int, p_baseline, clinic_odds0){
    
    # each subject has a random intercept
    subj_int <- rnorm(1, clinic_int, sig_person)
    
    # eta <- log(rep(clinic_odds0, length(t))) + rep(hosp_int, length(t)) + 
    #   rep(clinic_int, length(t)) + rep(subj_int, length(t)) 
    
    data.table(hosp_int = hosp_int,
               clinic_int = clinic_int,
               subj_int = rep(subj_int, length(t)),
               p_baseline = p_baseline)
  }
  
  
  
  l <- list()
  for(i in 1:hospts){
    hosp_int <- rnorm(1, 0, sig_hosp)
    
    # For each clinic, generate data sized by clinic with each clinic
    # having a random intercept
    for(j in 1:length(n_clinic)){
      
      clinic_int <- rnorm(1, hosp_int, sig_clinic)
      clinic_odds0 <- prob_to_odd(p0[j])
      
      l <- c(l, lapply(1:n_clinic[j], 
                       gensub,
                       hosp_int, clinic_int, p0[j], clinic_odds0))
      
      
    }
  }
  
  d <- as.data.table(cbind(d0, rbindlist(l)))
  
  d$eta <- log(prob_to_odd(p0[d0$clinic_id])) + 
    d$hosp_int + 
    d$clinic_int + 
    d$subj_int +
    #0.05*d$k1*d$time + 0.05*d$k2*d$time + 0.05*d$k3*d$time + 0.05*d$k4*d$time +
    #0.25*d$tx1*d$tx_durn + 0.25*d$tx2*d$tx_durn + 0.25*d$tx3*d$tx_durn + 0.25*d$tx4*d$tx_durn 
    1.05*d$tx1 + 1.05*d$tx2 + 1.05*d$tx3 + 1.05*d$tx4
  
  #head(inv_logit(d$eta))
  d$p <- inv_logit(d$eta)
  d$y <- rbinom(nrow(d), 1, prob = d$p)
  d
})

# sig_hosp <- 0.05
# sig_clinic <- 0.2
# sig_person <- 0.5

# d <- gendat()
# 
# names(d)
# 
# dfig <- d[, .(prop = mean(y)), keyby = .(hosp_id, clinic_id, time, tx_active)] 
# 
# dfig <- merge(dfig, data.frame(hosp_id = 1:hospts, 
#                                start_tx = rep(1:4, each = 2)), by = "hosp_id")
# 
# ds1 <- merge(expand.grid(hosp_id = 1:hospts, 
#                  clinic_id = 1:length(n_clinic)),
#              data.frame(hosp_id = 1:hospts, 
#                         start_tx = rep(1:4, each = 2)), by = "hosp_id")
# 
# ggplot(dfig, aes(x = time, y = prop, group = paste0(hosp_id, ":", clinic_id)))+
#   geom_point(size = 0.4, alpha = 0.2) + 
#   geom_line() +
#   scale_x_continuous(lim = c(0, 5))+
#   scale_y_continuous(lim = c(0, 1), breaks = c(0, 0.5, 1)) +
#   facet_grid(paste0("Start tx ", start_tx) ~ paste0("Clinic ", clinic_id))+
#   geom_vline(data = ds1, aes(xintercept = start_tx), linetype = 2)
# 
# 
# 
# names(d)


library(doParallel)
library(foreach)




cl <- makeCluster(parallel::detectCores() - 2, outfile="")
registerDoParallel(cl)

results <- foreach(i = 1:1000,
                   .errorhandling = 'pass',
                   .packages = c("data.table", "lme4")
                   #.export=c("gendat", "create_base", "prob_to_odd", "inv_logit")
                   #.options.snow=opts,
) %dopar%{
  
  # i = 1
  #flog.info("Starting trial: i = %s", i)
  set.seed(10000000 + i)
 
cat(paste0("Iteration ", i, "\n")) 

  d <- gendat()
  l1 <- glmer(y ~ k1 + k2 + k3 + k4 + tx_active + (1|hosp_id) + (1|clinic_int) + (1|subj_int), 
    data = d, family = binomial)
  # format(object.size(d), units = "MB")
  # format(object.size(l1), units = "MB")
  # format(object.size(s), units = "b")
  
  s <- list(coef = summary(l1)$coef, varcorr = VarCorr(l1))
  d <- NULL
  l1 <- NULL
  
  return(s)
}


stopCluster(cl)

beepr::beep() 
rdsfilename <- paste0("out/res-",format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".RDS")

saveRDS(list(results = results,
  w = warnings()), 
  rdsfilename)



# res <- readRDS("out/res-2019-04-04-23-20-13.RDS")
# m <- res$results
# win <- unlist(lapply(1:length(m), function(x) summary(results[[x]])$coef["tx_active", "Pr(>|z|)"] < 0.05 ))
# mean(win)

# l1 <- glmer(y ~ k1 + k2 + k3 + k4 + tx_active + (1|hosp_id), data = d, family = binomial)
# summary(l1)
# 
# julia_command("using MixedModels"); 
# dj <- JuliaObject(d)
# julia_assign("dj", d)
# julia_command("f = @formula(y ~ k1 + k2 + k3 + k4 + tx_active + (1|hosp_id));")
# l2 <- julia_command("gm1 = fit(GeneralizedLinearMixedModel, f, dj, Bernoulli())"); 
# l2 <- julia_eval("gm1")


## Use optimization to estimate single-drug fitnesses
## August 26, 2019
## V. Kim

argv = commandArgs(trailingOnly = TRUE)
# argv[1] = strain 1 data directory
# argv[2] = strain 2 data directory
# argv[3] = strain 3 data directory
# argv = c("~/Documents/embl/screen_K", "~/Documents/embl/screen_M",
# "~/Documents/embl/screen_D")
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tools)

datadir = "~/Documents/embl/gitlab/species_conservation/data/"
if(!dir.exists(datadir)) {
  datadir = "data/"
}

figdir = "~/Documents/embl/gitlab/species_conservation/figures/qc/"
if(!dir.exists(figdir)) {
  figdir = "figures/qc/"
}

cons_scriptdir = "~/Documents/embl/gitlab/species_conservation/conservation/"
if(!dir.exists(cons_scriptdir)) cons_scriptdir = "conservation/"
source(paste0(cons_scriptdir, "fun_cons.R"))
# source the script with core functions
core_scriptdir = "~/Documents/embl/gitlab/drug_comb_screen/dataprep/"
if(!dir.exists(core_scriptdir)) core_scriptdir = "/g/huber/users/vkim/gitlab/parallel/staph-main/dataprep/"
source(paste0(core_scriptdir, "code_screen_core_fun.R"))
# downstream processing of screen results
source(paste0(core_scriptdir, "downstream_fun.R"))

strains = unlist(lapply(argv, function(x) gsub("(.+)_([A-Z])", "\\2", x)))

get.fitness.data <- function(home_dir, trimtp=6, maxfit=1) {
  
  batches = grep("batch", list.files(home_dir), value = T)
  
  alldata = lapply(batches, function(b) {
    # load annotated data
    batch <- list.files(paste0(home_dir, b))
    batch_annot = lapply(paste0(home_dir, b, "/", batch), read.csv, header=T,
                         check.names = FALSE)
    names(batch_annot) <- file_path_sans_ext(batch)
    
    params = rjson::fromJSON(file = paste0("params/",b, ".json" ))
    batch_annot = compute_plate_stats(annotated = batch_annot,
                                      params = params)
    batch_annot = plyr::ldply(batch_annot, data.frame, .id = "plate")
    batch_annot = filter_data(batch_annot, params = params)
    batch_annot
  })
  names(alldata) = batches
  
  
  alldata_df = plyr::ldply(alldata, data.frame, .id=NULL)
  alldata_df = dplyr::distinct(alldata_df, plate, Well, Drug, Concentration,
                               Replicate, run, donor, donor_conc, biol_rep, batch,
                               flag_fixed, AUCnorm_bug_rob, correl)
  
  alldata_df
}

get.top.perc2 <- function(data, p=0.9) {
  top.wells = group_by(data, batch, plate) %>%
    filter(fit_comb > quantile(fit_comb, p, na.rm=T)) %>%
    group_by(Well) %>%
    dplyr::summarise(n=n()) %>% filter(n > quantile(n, p))
  
  annot = distinct(data, Well, Drug, flag_fixed)
  
  top.drugs = inner_join(top.wells, annot) %>%
    arrange(desc(n))
  
  don.estim = inner_join(data, top.wells) %>%
    group_by(batch, biol_rep, plate, donor, donor_real_conc) %>%
    dplyr::summarise(fit_don = median(fit_comb))
  
  don.estim
  
}


get.top.wells <- function(data, p=0.9) {
  top.wells = group_by(data, batch, plate) %>%
    filter(fit_comb > quantile(fit_comb, p, na.rm=T)) %>%
    group_by(Well) %>%
    dplyr::summarise(n=n()) %>% filter(n > quantile(n, p))
  
  top.wells
  
}


fitdata_all = lapply(1:3, function(i) {
  AUCdf =  get.fitness.data(home_dir = paste0(argv[i], 
                                              "/screen",
                                              strains[i],
                                              "_annotated_reads/"))
  AUCdf
})
names(fitdata_all) = strains
fitdata_df = plyr::ldply(fitdata_all, data.frame, .id = "strain")

donconc = lapply(1:3, function(i) {
  read.table(file.path(argv[i], paste0("donconc", strains[i], ".txt")),
             sep = "\t", header = T, stringsAsFactors = F)
  
})



optimize_epsilon <- function(batchdata, norm="L2") {
  # objective function
  Fcn <- function(x) {
    donors = as.matrix(x[1:ncol(f_rd)])
    recips = as.matrix(x[(ncol(f_rd)+1):length(x)])
    
    f_exp = recips %*% t(donors)
    
    sum((f_rd - f_exp)^2)
  }
  # gradient of the objectie function
  gradFcn <- function(x) {
    donors = as.matrix(x[1:ncol(f_rd)])
    recips = as.matrix(x[(ncol(f_rd)+1):length(x)])
    
    f_exp = recips %*% t(donors)
    f_d = matrix(rep(donors, length(recips)), nrow = length(recips), byrow = T)
    f_r = matrix(rep(recips, length(donors)), ncol = length(donors))
    c(colSums(-2*f_r * (f_rd - f_exp)), rowSums(-2*f_d * (f_rd - f_exp)))
    
  }
  
  Fcn_L1 <- function(x) {
    donors = as.matrix(x[1:ncol(f_rd)])
    recips = as.matrix(x[(ncol(f_rd)+1):length(x)])
    
    f_exp = recips %*% t(donors)
    
    sum(abs(f_rd - f_exp))
  }
  
  grad_L1 <- function(x) {
    donors = as.matrix(x[1:ncol(f_rd)])
    recips = as.matrix(x[(ncol(f_rd)+1):length(x)])
    
    f_exp = recips %*% t(donors)
    f_d = matrix(rep(donors, length(recips)), nrow = length(recips), byrow = T)
    f_r = matrix(rep(recips, length(donors)), ncol = length(donors))
    c(colSums(-f_r * sign(f_rd - f_exp)), rowSums(-f_d * sign(f_rd - f_exp)))
    
  }
  
  # filter out TSB plates and non-drug wells
  batchdata = filter(batchdata, !Drug %in% c('BUG', 'ONE', 'NOT'))
  batchdata = filter(batchdata, donor != "TSB")
  
  # construct f_rd matrix
  combfit = distinct(batchdata, Drug, Concentration, Replicate, donor, donor_conc, biol_rep, AUCnorm_bug_rob)
  # f_rd matrix
  f_rd = mutate(combfit, rown = paste(Drug, Concentration, Replicate, sep = "_"),
                coln = paste(donor, donor_conc, biol_rep, sep = "_")) %>%
    select(-c(Drug, Concentration, Replicate, donor, donor_conc, biol_rep)) %>%
    distinct(rown, coln, .keep_all = T) %>%
    tidyr::spread(key = coln, value = AUCnorm_bug_rob)
  
  rownames(f_rd) = f_rd$rown
  f_rd = as.matrix(f_rd[,-1])
  
  
  set.seed(100); xinit = runif(n= ncol(f_rd)+nrow(f_rd))
  
  if(norm == "L2") {
    res = optim(xinit, Fcn, gradFcn, method = "L-BFGS-B",
                lower = rep(0, length(xinit)), upper = rep(1, length(xinit)),
                control = list(maxit = 1e5))
  } else {
    res = optim(xinit, Fcn_L1, grad_L1, method = "L-BFGS-B",
                lower = rep(0, length(xinit)), upper = rep(1, length(xinit)),
                control = list(maxit = 1e5))
  }
  
  # check if the result converged
  stopifnot(res$convergence == 0)
  
  # now check the Bliss scores for L2
  donfit = tibble(fit_don = res$par[1:ncol(f_rd)],
                  donconc = colnames(f_rd))
  recfit = tibble(fit_rec = res$par[(ncol(f_rd)+1):length(res$par)],
                  recconcrep = rownames(f_rd))
  
  reccombs = inner_join(mutate(combfit, recconcrep = paste(Drug, Concentration, Replicate, sep="_"),
                               donconc = paste(donor, donor_conc, biol_rep, sep="_")), recfit, by = "recconcrep")
  donreccombs = inner_join(reccombs, donfit, by = "donconc")
  donreccombs = mutate(donreccombs, eps = AUCnorm_bug_rob - fit_rec * fit_don)
  
  donreccombs
}


newman = filter(fitdata_df, strain == "K")
epsnewman = lapply(unique(newman$batch), function(b) {
  batchdata = filter(newman, batch == b)
  epsbatch = optimize_epsilon(batchdata, norm = "L2")
  epsbatch
})
names(epsnewman) = unique(newman$batch)
epsnewman = plyr::ldply(epsnewman, .id = "batch")
epsnewman$comb = purrr::map2_chr(as.character(epsnewman$Drug),
                                 as.character(epsnewman$donor),
                                 ~ paste(sort(c(.x, .y)), collapse = "_"))

pval_newman = epsnewman %>% group_by(comb) %>%
  mutate(bernoulli = eps > median(epsnewman$eps)) %>%
  dplyr::summarise(sumbern = sum(bernoulli),
                   n=n()) %>%
  group_by(comb) %>%
  dplyr::mutate(pval = binom.test(x=sumbern,
                                  n=n)$p.value)
pval_newman = ungroup(pval_newman) %>%
  mutate(padj = p.adjust(pval, method = "BH"))

compute.covariates <- function(eps_df) {
  eps_df = mutate(eps_df, 
                  recconc_num = gsub("([A-Z]+_)(.+)", "\\2", recconc),
                  donconc_num = gsub("([A-Z]+_)(.+)", "\\2", donconc),
                  combid = paste(donor, donor_real_conc, Drug, Concentration, sep = "_"))
  
  eps_df = mutate(eps_df,
                  unique = paste(batch, biol_rep, Replicate, donor, sep = "_"))
  
  eps_df = distinct(eps_df, unique, combid, .keep_all = T)
  
  
  combs = unique(eps_df$comb)
  medpbliss = numeric()
  cov.fit = numeric()
  
  for(i in 1:length(combs)) {
    d1 = sub("([A-Z]+)_(.+)", "\\1", combs[i])
    comb_df = filter(eps_df, comb == combs[i])
    
    comb_df = comb_df %>% mutate(d1_conc = ifelse(donor == d1, 
                                                  as.character(donconc),
                                                  as.character(recconc)),
                                 d2_conc = ifelse(donor == d1, as.character(recconc),
                                                  as.character(donconc)))
    # estimate effect size as overall effect of median polished
    # matrix with (9 combs x n replicates)
    comb_eps_mat = comb_df %>%
      reshape2::dcast(d1_conc + d2_conc  ~ unique, mean, value.var = "eps")
    
    eps_heatmat = as.matrix(comb_eps_mat[,3:ncol(comb_eps_mat)])
    rownames(eps_heatmat) = paste(comb_eps_mat$d1_conc, comb_eps_mat$d2_conc, sep = "_")
    
    noise = sd(as.numeric(eps_heatmat), na.rm = T)
    cov.fit = c(cov.fit, noise)
    mdp = medpolish(eps_heatmat, maxiter = 100, na.rm = T)
    eps.medp = mdp$overall
    medpbliss = c(medpbliss, eps.medp)
    
  }
  
  data.frame(comb = combs,
             noise = cov.fit,
             medpbliss = medpbliss)
  
  
}

donreal = inner_join(distinct(newman, batch, biol_rep, run, donor, donor_conc),
           donconc[[1]], by = c("donor", "donor_conc", "run"))
eps_df = inner_join(mutate(epsnewman, batch = as.integer(batch)), donreal, by = c("donor", "donor_conc", "batch", "biol_rep")) %>%
  select(-c(donconc, recconcrep)) %>%
  mutate(recconc = paste(Drug, Concentration, sep = "_"),
         donconc = paste(donor, donor_real_conc, sep = "_"))
medp = compute.covariates(eps_df)

pval_newman = inner_join(pval_newman, medp, by = "comb")

#View(filter(pval_newman, padj <= 0.05 & abs(medpbliss) >= 0.05))


dsm = filter(fitdata_df, strain == "M")
epsdsm = lapply(unique(dsm$batch), function(b) {
  batchdata = filter(dsm, batch == b)
  epsbatch = optimize_epsilon(batchdata, norm = "L2")
  epsbatch
})
names(epsdsm) = unique(dsm$batch)
epsdsm = plyr::ldply(epsdsm, .id = "batch")
epsdsm$comb = purrr::map2_chr(as.character(epsdsm$Drug),
                                 as.character(epsdsm$donor),
                                 ~ paste(sort(c(.x, .y)), collapse = "_"))

pval_dsm = epsdsm %>% group_by(comb) %>%
  mutate(bernoulli = eps > median(epsdsm$eps)) %>%
  dplyr::summarise(sumbern = sum(bernoulli),
                   n=n()) %>%
  group_by(comb) %>%
  dplyr::mutate(pval = binom.test(x=sumbern,
                                  n=n)$p.value)
pval_dsm = ungroup(pval_dsm) %>%
  mutate(padj = p.adjust(pval, method = "BH"))

donreal = inner_join(distinct(dsm, batch, biol_rep, run, donor, donor_conc),
                     donconc[[2]], by = c("donor", "donor_conc", "run"))
eps_df = inner_join(mutate(epsdsm, batch = as.numeric(as.character(batch)),
                           donor = as.character(donor)),
                    donreal, by = c("donor", "donor_conc", "batch", "biol_rep")) %>%
  select(-c(donconc, recconcrep)) %>%
  mutate(recconc = paste(Drug, Concentration, sep = "_"),
         donconc = paste(donor, donor_real_conc, sep = "_"))
medp = compute.covariates(eps_df)
pval_dsm = inner_join(pval_dsm, medp, by = "comb")
#View(filter(pval_dsm, padj <= 0.05 & abs(medpbliss) >= 0.05))


bsub = filter(fitdata_df, strain == "D")
epsbsub = lapply(unique(bsub$batch), function(b) {
  batchdata = filter(bsub, batch == b)
  epsbatch = optimize_epsilon(batchdata, norm = "L2")
  epsbatch
})
names(epsbsub) = unique(bsub$batch)
epsbsub = plyr::ldply(epsbsub, .id = "batch")

epsbsub$comb = purrr::map2_chr(as.character(epsbsub$Drug),
                              as.character(epsbsub$donor),
                              ~ paste(sort(c(.x, .y)), collapse = "_"))

pval_bsub = epsbsub %>% group_by(comb) %>%
  mutate(bernoulli = eps > median(epsbsub$eps)) %>%
  dplyr::summarise(sumbern = sum(bernoulli),
                   n=n()) %>%
  group_by(comb) %>%
  dplyr::mutate(pval = binom.test(x=sumbern,
                                  n=n)$p.value)
pval_bsub = ungroup(pval_bsub) %>%
  mutate(padj = p.adjust(pval, method = "BH"))

donreal = inner_join(distinct(bsub, batch, biol_rep, run, donor, donor_conc),
                     donconc[[3]], by = c("donor", "donor_conc", "run"))
eps_df = inner_join(mutate(epsbsub, batch = as.numeric(as.character(batch)),
                           donor = as.character(donor)),
                    donreal, by = c("donor", "donor_conc", "batch", "biol_rep")) %>%
  select(-c(donconc, recconcrep)) %>%
  mutate(recconc = paste(Drug, Concentration, sep = "_"),
         donconc = paste(donor, donor_real_conc, sep = "_"))
medp = compute.covariates(eps_df)
pval_bsub = inner_join(pval_bsub, medp, by = "comb")
#View(filter(pval_bsub, padj <= 0.05 & abs(medpbliss) >= 0.05))



write.table(pval_newman, file = paste0(datadir, "Newman-L2.txt"),
            sep = "\t", quote = F, row.names = F,
            col.names = T)

write.table(pval_dsm, file = paste0(datadir, "DSM20231-L2.txt"),
            sep = "\t", quote = F, row.names = F,
            col.names = T)

write.table(pval_bsub, file = paste0(datadir, "Bsubtilis-L2.txt"),
            sep = "\t", quote = F, row.names = F,
            col.names = T)

# set.seed(100); xinit = runif(n= ncol(f_rd)+nrow(f_rd))
# res = optim(xinit, Fcn, control = list(maxit = 1e6))
# 
# set.seed(100); xinit = runif(n= ncol(f_rd)+nrow(f_rd))
# res = optim(xinit, Fcn, gradFcn, method = "BFGS",
#             control = list(maxit = 1e4))

# L-BFGS is the only optimzation method with box constraint
set.seed(100); xinit = runif(n= ncol(f_rd)+nrow(f_rd))
res_L2 = optim(xinit, Fcn, gradFcn, method = "L-BFGS-B",
            lower = rep(0, length(xinit)), upper = rep(1, length(xinit)),
            control = list(maxit = 1e4))

res_L1 = optim(xinit, Fcn_L1, grad_L1, method = "L-BFGS-B",
            lower = rep(0, length(xinit)), upper = rep(1, length(xinit)),
            control = list(maxit = 1e4))

# compare with experimentally determined initial values
# recipient fitness
recfit = lapply(fitdata_all, function(x) {
  tsb = filter(x, donor == "TSB")
  tsb = dplyr::filter(tsb, !run %in% c(4, 15, 18, 52))

  recfit = dplyr::group_by(tsb, batch, Drug, Concentration, Replicate, biol_rep) %>%
    dplyr::summarise(fit_rec = median(AUCnorm_bug_rob))

  recfit
})
recfit_df = plyr::ldply(recfit, .id = "strain")


topwells = c(0.75, 0.75, 0.9)
# donor fitneess
donfit = lapply(1:3, function(i) {
  combs = dplyr::filter(fitdata_all[[i]], donor != "TSB" & !Drug %in% c("ONE", "BUG", "NOT"))
  combs = dplyr::inner_join(combs, donconc[[i]])
  combs = dplyr::rename(combs, fit_comb = AUCnorm_bug_rob)

  topfit = get.top.perc2(combs, p=topwells[i])
  topfit
})
names(donfit) = strains
donfit_df = plyr::ldply(donfit, data.frame, .id = "strain")



estimate_single_fit <- function(batchdata, norm="L2", mode = "donor") {
  # objective function
  Fcn <- function(x) {
    donors = as.matrix(x[1:ncol(f_rd)])
    recips = as.matrix(x[(ncol(f_rd)+1):length(x)])
    
    f_exp = recips %*% t(donors)
    
    sum((f_rd - f_exp)^2)
  }
  # gradient of the objectie function
  gradFcn <- function(x) {
    donors = as.matrix(x[1:ncol(f_rd)])
    recips = as.matrix(x[(ncol(f_rd)+1):length(x)])
    
    f_exp = recips %*% t(donors)
    f_d = matrix(rep(donors, length(recips)), nrow = length(recips), byrow = T)
    f_r = matrix(rep(recips, length(donors)), ncol = length(donors))
    c(colSums(-2*f_r * (f_rd - f_exp)), rowSums(-2*f_d * (f_rd - f_exp)))
    
  }
  
  Fcn_L1 <- function(x) {
    donors = as.matrix(x[1:ncol(f_rd)])
    recips = as.matrix(x[(ncol(f_rd)+1):length(x)])
    
    f_exp = recips %*% t(donors)
    
    sum(abs(f_rd - f_exp))
  }
  
  grad_L1 <- function(x) {
    donors = as.matrix(x[1:ncol(f_rd)])
    recips = as.matrix(x[(ncol(f_rd)+1):length(x)])
    
    f_exp = recips %*% t(donors)
    f_d = matrix(rep(donors, length(recips)), nrow = length(recips), byrow = T)
    f_r = matrix(rep(recips, length(donors)), ncol = length(donors))
    c(colSums(-f_r * sign(f_rd - f_exp)), rowSums(-f_d * sign(f_rd - f_exp)))
    
  }
  
  # filter out TSB plates and non-drug wells
  batchdata = filter(batchdata, !Drug %in% c('BUG', 'ONE', 'NOT'))
  batchdata = filter(batchdata, donor != "TSB")
  
  # construct f_rd matrix
  combfit = distinct(batchdata, Drug, Concentration, Replicate, donor, donor_conc, biol_rep, AUCnorm_bug_rob)
  # f_rd matrix
  f_rd = mutate(combfit, rown = paste(Drug, Concentration, Replicate, sep = "_"),
                coln = paste(donor, donor_conc, biol_rep, sep = "_")) %>%
    select(-c(Drug, Concentration, Replicate, donor, donor_conc, biol_rep)) %>%
    distinct(rown, coln, .keep_all = T) %>%
    tidyr::spread(key = coln, value = AUCnorm_bug_rob)
  
  rownames(f_rd) = f_rd$rown
  f_rd = as.matrix(f_rd[,-1])

  set.seed(100); xinit = runif(n= ncol(f_rd)+nrow(f_rd))
  
  if(norm == "L2") {
    res = optim(xinit, Fcn, gradFcn, method = "L-BFGS-B",
                lower = rep(0, length(xinit)), upper = rep(1, length(xinit)),
                control = list(maxit = 1e5))
  } else {
    res = optim(xinit, Fcn_L1, grad_L1, method = "L-BFGS-B",
                lower = rep(0, length(xinit)), upper = rep(1, length(xinit)),
                control = list(maxit = 1e5))
  }
  
  # check if the result converged
  stopifnot(res$convergence == 0)
  
  if(mode == "donor") {
    tibble(fit_don = res$par[1:ncol(f_rd)],
           donconc = colnames(f_rd))
  } else {
    tibble(fit_rec = res$par[(ncol(f_rd)+1):length(res$par)],
           recconcrep = rownames(f_rd))
  }

}

newman = filter(fitdata_df, strain == "K")
donfitnewman = lapply(unique(newman$batch), function(b) {
  batchdata = filter(newman, batch == b)
  epsbatch = estimate_single_fit(batchdata, norm = "L2", mode = "donor")
  epsbatch
})
names(donfitnewman) = unique(newman$batch)
donfitnewman_L2 = plyr::ldply(donfitnewman, .id = "batch")

donfitnewman = lapply(unique(newman$batch), function(b) {
  batchdata = filter(newman, batch == b)
  epsbatch = estimate_single_fit(batchdata, norm = "L1", mode = "donor")
  epsbatch
})
names(donfitnewman) = unique(newman$batch)
donfitnewman_L1 = plyr::ldply(donfitnewman, .id = "batch")

donfitnewman_L2$norm = "L2"
donfitnewman_L1$norm = "L1"

donfitnewman = bind_rows(donfitnewman_L1, donfitnewman_L2)


donfitnewman = tidyr::separate(donfitnewman, donconc,
                               c("donor", "donor_conc", "biol_rep"), sep="_")

donreal = inner_join(distinct(newman, batch, biol_rep, run, donor, donor_conc),
                     donconc[[1]], by = c("donor", "donor_conc", "run"))
donfitnewman = inner_join(mutate(donfitnewman, donor_conc = as.numeric(donor_conc),
                  batch = as.numeric(as.character(batch)),
                  biol_rep = as.numeric(biol_rep)), 
           donreal,
           by = c("donor", "donor_conc", "batch", "biol_rep")) 

donfitnewman = rename(donfitnewman, fit_don_est = fit_don)
donfitnewman$strain = "K"

p1 = ggplot(inner_join(donfitnewman, donfit_df),
       aes(x = fit_don, y = fit_don_est,
           color = norm)) +
  geom_point() +
  facet_wrap(~ batch) +
  scale_color_manual(values = c("#a0cc78", "#207561")) + 
  geom_abline(slope = 1, color = "grey", linetype="dashed") +
  theme_bw() +
  labs(x = "Donor fitness based on top 25% growing wells",
       y = "Estimated donor fitness",
       color = "Norm")
ggsave(p1, filename = paste0(figdir, "donfit-estim-Newman.pdf"))



donfitdsm = lapply(unique(dsm$batch), function(b) {
  batchdata = filter(dsm, batch == b)
  epsbatch = estimate_single_fit(batchdata, norm = "L2", mode = "donor")
  epsbatch
})
names(donfitdsm) = unique(dsm$batch)
donfitdsm_L2 = plyr::ldply(donfitdsm, .id = "batch")

donfitdsm = lapply(unique(dsm$batch), function(b) {
  batchdata = filter(dsm, batch == b)
  epsbatch = estimate_single_fit(batchdata, norm = "L1", mode = "donor")
  epsbatch
})
names(donfitdsm) = unique(dsm$batch)
donfitdsm_L1 = plyr::ldply(donfitdsm, .id = "batch")

donfitdsm_L2$norm = "L2"
donfitdsm_L1$norm = "L1"

donfitdsm = bind_rows(donfitdsm_L1, donfitdsm_L2)


donfitdsm = tidyr::separate(donfitdsm, donconc,
                               c("donor", "donor_conc", "biol_rep"), sep="_")

donreal = inner_join(distinct(dsm, batch, biol_rep, run, donor, donor_conc),
                     donconc[[2]], by = c("donor", "donor_conc", "run"))
donfitdsm = inner_join(mutate(donfitdsm, donor_conc = as.numeric(donor_conc),
                                 batch = as.numeric(as.character(batch)),
                                 biol_rep = as.numeric(biol_rep)), 
                          donreal,
                          by = c("donor", "donor_conc", "batch", "biol_rep")) 

donfitdsm = rename(donfitdsm, fit_don_est = fit_don)
donfitdsm$strain = "M"

p2 = ggplot(inner_join(donfitdsm, donfit_df),
            aes(x = fit_don, y = fit_don_est,
                color = norm)) +
  geom_point() +
  facet_wrap(~ batch) +
  scale_color_manual(values = c("#a0cc78", "#207561")) + 
  geom_abline(slope = 1, color = "grey", linetype="dashed") +
  theme_bw() +
  labs(x = "Donor fitness based on top 25% growing wells",
       y = "Estimated donor fitness",
       color = "Norm")
ggsave(p2, filename = paste0(figdir, "donfit-estim-DSM20231.pdf"))

donfitbsub = lapply(unique(bsub$batch), function(b) {
  batchdata = filter(bsub, batch == b)
  epsbatch = estimate_single_fit(batchdata, norm = "L2", mode = "donor")
  epsbatch
})
names(donfitbsub) = unique(bsub$batch)
donfitbsub_L2 = plyr::ldply(donfitbsub, .id = "batch")

donfitbsub = lapply(unique(bsub$batch), function(b) {
  batchdata = filter(bsub, batch == b)
  epsbatch = estimate_single_fit(batchdata, norm = "L1", mode = "donor")
  epsbatch
})
names(donfitbsub) = unique(bsub$batch)
donfitbsub_L1 = plyr::ldply(donfitbsub, .id = "batch")

donfitbsub_L2$norm = "L2"
donfitbsub_L1$norm = "L1"

donfitbsub = bind_rows(donfitbsub_L1, donfitbsub_L2)


donfitbsub = tidyr::separate(donfitbsub, donconc,
                               c("donor", "donor_conc", "biol_rep"), sep="_")

donreal = inner_join(distinct(bsub, batch, biol_rep, run, donor, donor_conc),
                     donconc[[3]], by = c("donor", "donor_conc", "run"))
donfitbsub = inner_join(mutate(donfitbsub, donor_conc = as.numeric(donor_conc),
                                 batch = as.numeric(as.character(batch)),
                                 biol_rep = as.numeric(biol_rep)), 
                          donreal,
                          by = c("donor", "donor_conc", "batch", "biol_rep")) 

donfitbsub = rename(donfitbsub, fit_don_est = fit_don)
donfitbsub$strain = "D"

p3 = ggplot(inner_join(donfitbsub, donfit_df),
            aes(x = fit_don, y = fit_don_est,
                color = norm)) +
  geom_point() +
  facet_wrap(~ batch) +
  scale_color_manual(values = c("#a0cc78", "#207561")) + 
  geom_abline(slope = 1, color = "grey", linetype="dashed") +
  theme_bw() +
  labs(x = "Donor fitness based on top 25% growing wells",
       y = "Estimated donor fitness",
       color = "Norm") +
  theme(aspect.ratio = 1)
ggsave(p3, filename = paste0(figdir, "donfit-estim-Bsubtilis.pdf"))



recfitnewman = lapply(unique(newman$batch), function(b) {
  batchdata = filter(newman, batch == b)
  epsbatch = estimate_single_fit(batchdata, norm = "L2", mode = "recipient")
  epsbatch
})
names(recfitnewman) = unique(newman$batch)
recfitnewman_L2 = plyr::ldply(recfitnewman, .id = "batch")

recfitnewman = lapply(unique(newman$batch), function(b) {
  batchdata = filter(newman, batch == b)
  epsbatch = estimate_single_fit(batchdata, norm = "L1", mode = "recipient")
  epsbatch
})
names(recfitnewman) = unique(newman$batch)
recfitnewman_L1 = plyr::ldply(recfitnewman, .id = "batch")

recfitnewman_L2$norm = "L2"
recfitnewman_L1$norm = "L1"

recfitnewman = bind_rows(recfitnewman_L1, recfitnewman_L2)

recfitnewman = tidyr::separate(recfitnewman, recconcrep,
                               c("Drug", "Concentration", "Replicate"), sep="_")

recfitnewman = rename(recfitnewman, fit_rec_est = fit_rec)
recfitnewman$strain = "K"

p4 = ggplot(inner_join(mutate(recfitnewman,
                              batch = as.numeric(as.character(batch)),
                              Concentration = as.numeric(Concentration),
                              Replicate = as.numeric(Replicate)),
                       recfit_df),
            aes(x = fit_rec, y = fit_rec_est,
                color = norm)) +
  geom_point() +
  facet_wrap(~ batch) +
  scale_color_manual(values = c("#7ec4e7", "#2c6ea0")) + 
  geom_abline(slope = 1, color = "grey", linetype="dashed") +
  theme_bw() +
  labs(x = "TSB recipient fitness",
       y = "Estimated recipient fitness",
       color = "Norm")
ggsave(p4, filename = paste0(figdir, "recfit-estim-Newman.pdf"))

p5 = ggplot(inner_join(mutate(recfitnewman,
                              batch = as.numeric(as.character(batch)),
                              Concentration = as.numeric(Concentration),
                              Replicate = as.numeric(Replicate)),
                       recfit_df) %>% filter(norm == "L2"),
            aes(x = fit_rec, y = fit_rec_est)) +
  geom_point() +
  facet_wrap(~ batch) +
  #scale_color_manual(values = c("#7ec4e7", "#2c6ea0")) + 
  geom_abline(slope = 1, color = "grey", linetype="dashed") +
  ggrepel::geom_text_repel(aes(label = ifelse(abs(fit_rec_est - fit_rec) > 0.3,
                                              paste(Drug, Concentration, sep = "_"), "")),
                           segment.colour = NA) + 
  theme_bw() +
  labs(x = "TSB recipient fitness",
       y = "Estimated recipient fitness",
       color = "Norm")

ggsave(p5, filename = paste0(figdir, "recfit-estim-Newman-outliers.pdf"))

#---
recfitdsm = lapply(unique(dsm$batch), function(b) {
  batchdata = filter(dsm, batch == b)
  epsbatch = estimate_single_fit(batchdata, norm = "L2", mode = "recipient")
  epsbatch
})
names(recfitdsm) = unique(dsm$batch)
recfitdsm_L2 = plyr::ldply(recfitdsm, .id = "batch")

recfitdsm = lapply(unique(dsm$batch), function(b) {
  batchdata = filter(dsm, batch == b)
  epsbatch = estimate_single_fit(batchdata, norm = "L1", mode = "recipient")
  epsbatch
})
names(recfitdsm) = unique(dsm$batch)
recfitdsm_L1 = plyr::ldply(recfitdsm, .id = "batch")

recfitdsm_L2$norm = "L2"
recfitdsm_L1$norm = "L1"

recfitdsm = bind_rows(recfitdsm_L1, recfitdsm_L2)

recfitdsm = tidyr::separate(recfitdsm, recconcrep,
                               c("Drug", "Concentration", "Replicate"), sep="_")

recfitdsm = rename(recfitdsm, fit_rec_est = fit_rec)
recfitdsm$strain = "M"

p6 = ggplot(inner_join(mutate(recfitdsm,
                              batch = as.numeric(as.character(batch)),
                              Concentration = as.numeric(Concentration),
                              Replicate = as.numeric(Replicate)),
                       recfit_df),
            aes(x = fit_rec, y = fit_rec_est,
                color = norm)) +
  geom_point() +
  facet_wrap(~ batch) +
  scale_color_manual(values = c("#7ec4e7", "#2c6ea0")) + 
  geom_abline(slope = 1, color = "grey", linetype="dashed") +
  theme_bw() +
  labs(x = "TSB recipient fitness",
       y = "Estimated recipient fitness",
       color = "Norm")
ggsave(p6, filename = paste0(figdir, "recfit-estim-DSM20231.pdf"))

p7 = ggplot(inner_join(mutate(recfitdsm,
                              batch = as.numeric(as.character(batch)),
                              Concentration = as.numeric(Concentration),
                              Replicate = as.numeric(Replicate)),
                       recfit_df) %>% filter(norm == "L2"),
            aes(x = fit_rec, y = fit_rec_est)) +
  geom_point() +
  facet_wrap(~ batch) +
  #scale_color_manual(values = c("#7ec4e7", "#2c6ea0")) + 
  geom_abline(slope = 1, color = "grey", linetype="dashed") +
  ggrepel::geom_text_repel(aes(label = ifelse(abs(fit_rec_est - fit_rec) > 0.3,
                                              paste(Drug, Concentration, sep = "_"), "")),
                           segment.colour = NA) + 
  theme_bw() +
  labs(x = "TSB recipient fitness",
       y = "Estimated recipient fitness",
       color = "Norm")

ggsave(p7, filename = paste0(figdir, "recfit-estim-DSM20231-outliers.pdf"))

#----
recfitbsub = lapply(unique(bsub$batch), function(b) {
  batchdata = filter(bsub, batch == b)
  epsbatch = estimate_single_fit(batchdata, norm = "L2", mode = "recipient")
  epsbatch
})
names(recfitbsub) = unique(bsub$batch)
recfitbsub_L2 = plyr::ldply(recfitbsub, .id = "batch")

recfitbsub = lapply(unique(bsub$batch), function(b) {
  batchdata = filter(bsub, batch == b)
  epsbatch = estimate_single_fit(batchdata, norm = "L1", mode = "recipient")
  epsbatch
})
names(recfitbsub) = unique(bsub$batch)
recfitbsub_L1 = plyr::ldply(recfitbsub, .id = "batch")

recfitbsub_L2$norm = "L2"
recfitbsub_L1$norm = "L1"

recfitbsub = bind_rows(recfitbsub_L1, recfitbsub_L2)

recfitbsub = tidyr::separate(recfitbsub, recconcrep,
                            c("Drug", "Concentration", "Replicate"), sep="_")

recfitbsub = rename(recfitbsub, fit_rec_est = fit_rec)
recfitbsub$strain = "D"

p8 = ggplot(inner_join(mutate(recfitbsub,
                              batch = as.numeric(as.character(batch)),
                              Concentration = as.numeric(Concentration),
                              Replicate = as.numeric(Replicate)),
                       recfit_df),
            aes(x = fit_rec, y = fit_rec_est,
                color = norm)) +
  geom_point() +
  facet_wrap(~ batch) +
  scale_color_manual(values = c("#7ec4e7", "#2c6ea0")) + 
  geom_abline(slope = 1, color = "grey", linetype="dashed") +
  theme_bw() +
  labs(x = "TSB recipient fitness",
       y = "Estimated recipient fitness",
       color = "Norm") +
  theme(aspect.ratio = 1)
ggsave(p8, filename = paste0(figdir, "recfit-estim-Bsubtilis.pdf"))

p9 = ggplot(inner_join(mutate(recfitbsub,
                              batch = as.numeric(as.character(batch)),
                              Concentration = as.numeric(Concentration),
                              Replicate = as.numeric(Replicate)),
                       recfit_df) %>% filter(norm == "L2"),
            aes(x = fit_rec, y = fit_rec_est)) +
  geom_point() +
  facet_wrap(~ batch) +
  #scale_color_manual(values = c("#7ec4e7", "#2c6ea0")) + 
  geom_abline(slope = 1, color = "grey", linetype="dashed") +
  ggrepel::geom_text_repel(aes(label = ifelse(abs(fit_rec_est - fit_rec) > 0.25,
                                              paste(Drug, Concentration, sep = "_"), "")),
                           segment.colour = NA) + 
  theme_bw() +
  labs(x = "TSB recipient fitness",
       y = "Estimated recipient fitness",
       color = "Norm") +
  theme(aspect.ratio = 1)

ggsave(p9, filename = paste0(figdir, "recfit-estim-Bsubtilis-outliers.pdf"))



donfit_all =bind_rows(donfitnewman, donfitdsm, donfitbsub)
p10 = ggplot(mutate(donfit_all, strain = get_strain_name(strain)),
       aes(x = fit_don_est, fill = norm, group = norm))+
  geom_histogram(position = "identity",bins = 50, alpha = 0.7) +
  scale_fill_manual(values = c("#a0cc78", "#207561")) + 
  facet_wrap(~ strain) +
  theme_classic() +
  labs(x = "Estimated donor fitness") +
  theme(aspect.ratio = 1)
ggsave(p10, filename = paste0(figdir, "estim-donor-distrib.pdf"))


recfit_all = bind_rows(recfitnewman, recfitdsm, recfitbsub)
p11 = ggplot(mutate(recfit_all, strain = get_strain_name(strain)),
             aes(x = fit_rec_est, fill = norm, group = norm))+
  geom_histogram(position = "identity",bins = 50, alpha = 0.7) +
  scale_fill_manual(values = c("#7ec4e7", "#2c6ea0")) + 
  facet_wrap(~ strain) +
  theme_classic() +
  labs(x = "Estimated recipient fitness") +
  theme(aspect.ratio = 1)
ggsave(p11, filename = paste0(figdir, "estim-recip-distrib.pdf"))




estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

blissK = ggplot(epsnewman, aes(x = eps)) +
  geom_density(fill = "blue", colour=NA, alpha=0.3) + 
  theme_classic() + 
  geom_vline(xintercept = median(epsnewman$eps),
             colour = "red") + 
  geom_vline(xintercept = estimate_mode(epsnewman$eps),
             colour = "black",
             linetype = "dotted") + 
  labs(x = paste("Bliss scores in Newman"),
       y = "Count",
       title = paste("Distribution of Bliss scores in Newman. Median:",
                     format(median(epsnewman$eps), digits = 2),
                     "Mode:", format(estimate_mode(epsnewman$eps),digits = 2))) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(blissK, filename = paste0(figdir, "bliss-estim-fitness-Newman.pdf"))


blissM = ggplot(epsdsm, aes(x = eps)) +
  geom_density(fill = "blue", colour=NA, alpha=0.3) + 
  theme_classic() + 
  geom_vline(xintercept = median(epsdsm$eps),
             colour = "red") + 
  geom_vline(xintercept = estimate_mode(epsdsm$eps),
             colour = "black",
             linetype = "dotted") + 
  labs(x = paste("Bliss scores in DSM20231"),
       y = "Count",
       title = paste("Distribution of Bliss scores in DSM20231. Median:",
                     format(median(epsdsm$eps), digits = 2),
                     "Mode:", format(estimate_mode(epsdsm$eps),digits = 2))) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(blissM, filename = paste0(figdir, "bliss-estim-fitness-DSM20231.pdf"))


blissD = ggplot(epsbsub, aes(x = eps)) +
  geom_density(fill = "blue", colour=NA, alpha=0.3) + 
  theme_classic() + 
  geom_vline(xintercept = median(epsbsub$eps),
             colour = "red") + 
  geom_vline(xintercept = estimate_mode(epsbsub$eps),
             colour = "black",
             linetype = "dotted") + 
  labs(x = paste("Bliss scores in B. subtilis"),
       y = "Count",
       title = paste("Distribution of Bliss scores in B. subtilis. Median:",
                     format(median(epsbsub$eps), digits = 2),
                     "Mode:", format(estimate_mode(epsbsub$eps),digits = 2))) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(blissD, filename = paste0(figdir, "bliss-estim-fitness-Bsubtilis.pdf"))


# use Medpolish overall effect to summarise the Bliss scores of a combination
## plot the summarized distribution

pval_newman$strain = "Newman"
pval_dsm$strain = "DSM20231"
pval_bsub$strain = "B. subtilis"


pval_all = bind_rows(pval_newman, pval_dsm, pval_bsub)


qqnorm(scale(epsbsub$eps))
qqline(scale(epsbsub$eps), col = "steelblue", lwd =2)

qqnorm(scale(pval_all$medpbliss))
qqline(scale(pval_all$medpbliss), col = "steelblue", lwd = 2)

ggplot(mutate(pval_all, medpbliss = scale(medpbliss)), aes(x = medpbliss)) +
  geom_density(fill = "blue", colour=NA, alpha=0.3) + 
  theme_classic() +
  xlim(c(-7,7))

fdr_out = fdrtool::fdrtool(pval_all$medpbliss)


pval_all$qval = fdr_out$qval


#View(filter(pval_all, qval <= 0.05 & abs(medpbliss) > 0.03))







## New statistical test for the combinatorial screen
## September 3, 2019
## V. Kim

argv = commandArgs(trailingOnly = TRUE)
# argv[1] = strain 1 data directory
# argv[2] = strain 2 data directory
# argv[3] = strain 3 data directory
argv = c("~/Documents/embl/screen_K", "~/Documents/embl/screen_M",
"~/Documents/embl/screen_D")
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tools)
library(locfdr)

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

# UNCOMMENT to speed up computations
saveRDS(fitdata_df, file = paste0(datadir, "fitdata-grampos.rds"))
fitdata_df = readRDS(paste0(datadir, "fitdata-grampos.rds"))


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
                   n=n(),
                   eps.med = median(eps)) %>%
  group_by(comb) %>%
  dplyr::mutate(pval = binom.test(x=sumbern,
                                  n=n)$p.value)
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
                   n=n(),
                   eps.med = median(eps)) %>%
  group_by(comb) %>%
  dplyr::mutate(pval = binom.test(x=sumbern,
                                  n=n)$p.value)
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
                   n=n(),
                   eps.med = median(eps)) %>%
  group_by(comb) %>%
  dplyr::mutate(pval = binom.test(x=sumbern,
                                  n=n)$p.value)

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

# use Medpolish overall effect to summarise the Bliss scores of a combination
## plot the summarized distribution
pval_newman$strain = "Newman"
pval_dsm$strain = "DSM20231"
pval_bsub$strain = "B. subtilis"
pval_all = bind_rows(pval_newman, pval_dsm, pval_bsub)

# qqnorm(scale(epsbsub$eps))
# qqline(scale(epsbsub$eps), col = "steelblue", lwd =2)
# 
# ggplot(pval_all, aes(x = medpbliss)) +
#   geom_density(fill = "blue", colour=NA, alpha=0.3) + 
#   theme_classic() 
# 
# ggplot(pval_all, aes(x = medpbliss, y = eps.med)) +
#   geom_point() + 
#   theme_classic() +
#   geom_abline(slope = 1, color = "red")
# 
# fdr_out = fdrtool::fdrtool(pval_all$medpbliss)
# pval_all$qval = fdr_out$qval

#View(filter(pval_all, qval <= 0.05 & abs(medpbliss) > 0.03))
# w = locfdr(zz = pval_all$medpbliss, bre = 300, df = 8)
# # add local FDR
# pval_all$lfdr = w$fdr
# xleft = w$mat[which.min(abs(w$mat[,'Fdrleft'] - 0.1)), 'x']
# xright = w$mat[which.min(abs(w$mat[,'Fdrright'] - 0.1)), 'x']
# fdrleft = w$mat[which.min(abs(w$mat[,'Fdrleft'] - 0.1)), 'fdr']
# fdrright = w$mat[which.min(abs(w$mat[,'Fdrright'] - 0.1)), 'fdr']
# 
# p = ggplot(as.data.frame(w$mat),
#        aes(x = x,
#            y = Fdrleft)) +
#   geom_point(color = "#ffbb33") +
#   geom_line(data = as.data.frame(w$mat),
#             aes(x = x, y = fdr)) +
#   theme_bw() +
#   geom_hline(yintercept = 0.1, color = "#0066cc", linetype = "dashed") +
#   annotate("text", x = -0.25, y = 0.12, label = "FDR = 0.1") + 
#   geom_vline(xintercept = w$mat[which.min(abs(w$mat[,'Fdrleft'] - 0.1)), 'x'],
#              linetype = "dotted") +
#   annotate("text", x = -0.13, y = 0.34, label = paste("fdr =", 
#                                                      format(w$mat[which.min(abs(w$mat[,'Fdrleft'] - 0.1)), 'fdr'], digits = 2)))+
#   geom_point(data = as.data.frame(w$mat),
#              aes(x = x,
#                  y = Fdrright),
#              color = "#595959") +
#   geom_vline(xintercept = w$mat[which.min(abs(w$mat[,'Fdrright'] - 0.1)), 'x'],
#              linetype = "dotted") +
#   labs(x = "Effect size",
#        y = "Local (solid) and tail-area (points) FDR",
#        title = paste(nrow(filter(pval_all, (lfdr < fdrleft & medpbliss < 0) | (lfdr < fdrright & medpbliss > 0))), "interactions in Gram-positive")) + 
#   annotate("text", x = 0.14, y = 0.26, label = paste("fdr =", 
#                                                       format(w$mat[which.min(abs(w$mat[,'Fdrright'] - 0.1)), 'fdr'], digits = 2)))+
#   theme(plot.title = element_text(hjust = 0.5))
# ggsave(p, filename = "~/Documents/embl/presentations/2019/20190904_statistical_test_drugcombsreen/figures/locfdr-FDR.pdf",
#        width = 9, height = 7)


# test per strain as each combination is repeated 3 times if 
# the summarised Bliss scores are considered all together
w_newman = locfdr(zz = pval_newman$medpbliss, bre = 200)
pval_newman$lfdr = w_newman$fdr

xleft = w_newman$mat[which.min(abs(w_newman$mat[,'Fdrleft'] - 0.1)), 'x']
xright = w_newman$mat[which.min(abs(w_newman$mat[,'Fdrright'] - 0.1)), 'x']
fdrleft = w_newman$mat[which.min(abs(w_newman$mat[,'Fdrleft'] - 0.1)), 'fdr']
fdrright = w_newman$mat[which.min(abs(w_newman$mat[,'Fdrright'] - 0.1)), 'fdr']

pval_newman$fdrleft = fdrleft
pval_newman$fdrright = fdrright

write.table(select(pval_newman, -c(sumbern, n, strain)),
            file = paste0(datadir, "strains/Newman-lfdr.txt"),
            sep = "\t", quote = F, row.names = F)

p = ggplot(as.data.frame(w_newman$mat),
         aes(x = x,
             y = Fdrleft)) +
geom_point(color = "#ffbb33") +
geom_line(data = as.data.frame(w_newman$mat),
          aes(x = x, y = fdr)) +
theme_bw() +
geom_hline(yintercept = 0.1, color = "#0066cc", linetype = "dashed") +
annotate("text", x = -0.2, y = 0.12, label = "FDR = 0.1") + 
geom_vline(xintercept = xleft,
           linetype = "dotted") +
annotate("text", x = -0.12, y = 0.24, label = paste("fdr =", 
                                                    format(fdrleft, digits = 2)))+
geom_point(data = as.data.frame(w_newman$mat),
           aes(x = x,
               y = Fdrright),
           color = "#595959") +
geom_vline(xintercept = xright,
           linetype = "dotted") +
labs(x = "Effect size",
     y = "Local (solid) and tail-area (points) FDR",
     title = paste(nrow(filter(pval_newman, (lfdr < fdrleft & medpbliss < 0) | (lfdr < fdrright & medpbliss > 0))), "interactions in Newman")) + 
annotate("text", x = 0.12, y = 0.33, label = paste("fdr =", 
                                                   format(fdrright, digits = 2)))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(p, filename = "~/Documents/embl/presentations/2019/20190904_statistical_test_drugcombsreen/figures/locfdr-FDR-newman.pdf",
     width = 9, height = 7)


w_dsm = locfdr(zz = pval_dsm$medpbliss, bre = 200)
pval_dsm$lfdr = w_dsm$fdr
xleft = w_dsm$mat[which.min(abs(w_dsm$mat[,'Fdrleft'] - 0.1)), 'x']
xright = w_dsm$mat[which.min(abs(w_dsm$mat[,'Fdrright'] - 0.1)), 'x']
fdrleft = w_dsm$mat[which.min(abs(w_dsm$mat[,'Fdrleft'] - 0.1)), 'fdr']
fdrright = w_dsm$mat[which.min(abs(w_dsm$mat[,'Fdrright'] - 0.1)), 'fdr']

pval_dsm$fdrleft = fdrleft
pval_dsm$fdrright = fdrright

write.table(select(pval_dsm, -c(sumbern, n, strain)),
            file = paste0(datadir, "strains/DSM20231-lfdr.txt"),
            sep = "\t", quote = F, row.names = F)

p = ggplot(as.data.frame(w_dsm$mat),
           aes(x = x,
               y = Fdrleft)) +
  geom_point(color = "#ffbb33") +
  geom_line(data = as.data.frame(w_dsm$mat),
            aes(x = x, y = fdr)) +
  theme_bw() +
  geom_hline(yintercept = 0.1, color = "#0066cc", linetype = "dashed") +
  annotate("text", x = -0.2, y = 0.12, label = "FDR = 0.1") + 
  geom_vline(xintercept = xleft,
             linetype = "dotted") +
  annotate("text", x = -0.12, y = 0.25, label = paste("fdr =", 
                                                      format(fdrleft, digits = 2)))+
  geom_point(data = as.data.frame(w_dsm$mat),
             aes(x = x,
                 y = Fdrright),
             color = "#595959") +
  geom_vline(xintercept = xright,
             linetype = "dotted") +
  labs(x = "Effect size",
       y = "Local (solid) and tail-area (points) FDR",
       title = paste(nrow(filter(pval_dsm, (lfdr < fdrleft & medpbliss < 0) | (lfdr < fdrright & medpbliss > 0))), "interactions in DSM20231")) + 
  annotate("text", x = 0.14, y = 0.27, label = paste("fdr =", 
                                                     format(fdrright, digits = 2)))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(p, filename = "~/Documents/embl/presentations/2019/20190904_statistical_test_drugcombsreen/figures/locfdr-FDR-DSM20231-with-human.pdf",
       width = 9, height = 7)

combs_common = intersect(pval_newman$comb, pval_dsm$comb)
pval_dsm_nohum = filter(pval_dsm, comb %in% combs_common)

w_dsm = locfdr(zz = pval_dsm_nohum$medpbliss, bre = 200)
pval_dsm_nohum$lfdr = w_dsm$fdr
xleft = w_dsm$mat[which.min(abs(w_dsm$mat[,'Fdrleft'] - 0.1)), 'x']
xright = w_dsm$mat[which.min(abs(w_dsm$mat[,'Fdrright'] - 0.1)), 'x']
fdrleft = w_dsm$mat[which.min(abs(w_dsm$mat[,'Fdrleft'] - 0.1)), 'fdr']
fdrright = w_dsm$mat[which.min(abs(w_dsm$mat[,'Fdrright'] - 0.1)), 'fdr']

pval_dsm_nohum$fdrleft = fdrleft
pval_dsm_nohum$fdrright = fdrright

p = ggplot(as.data.frame(w_dsm$mat),
           aes(x = x,
               y = Fdrleft)) +
  geom_point(color = "#ffbb33") +
  geom_line(data = as.data.frame(w_dsm$mat),
            aes(x = x, y = fdr)) +
  theme_bw() +
  geom_hline(yintercept = 0.1, color = "#0066cc", linetype = "dashed") +
  annotate("text", x = -0.2, y = 0.12, label = "FDR = 0.1") + 
  geom_vline(xintercept = xleft,
             linetype = "dotted") +
  annotate("text", x = -0.12, y = 0.25, label = paste("fdr =", 
                                                      format(fdrleft, digits = 2)))+
  geom_point(data = as.data.frame(w_dsm$mat),
             aes(x = x,
                 y = Fdrright),
             color = "#595959") +
  geom_vline(xintercept = xright,
             linetype = "dotted") +
  labs(x = "Effect size",
       y = "Local (solid) and tail-area (points) FDR",
       title = paste(nrow(filter(pval_dsm_nohum, (lfdr < fdrleft & medpbliss < 0) | (lfdr < fdrright & medpbliss > 0))), "interactions in DSM20231")) + 
  annotate("text", x = 0.14, y = 0.27, label = paste("fdr =", 
                                                     format(fdrright, digits = 2)))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(p, filename = "~/Documents/embl/presentations/2019/20190904_statistical_test_drugcombsreen/figures/locfdr-FDR-DSM20231.pdf",
       width = 9, height = 7)


w_bsub = locfdr(zz = pval_bsub$medpbliss, bre = 200)
pval_bsub$lfdr = w_bsub$fdr
xleft = w_bsub$mat[which.min(abs(w_bsub$mat[,'Fdrleft'] - 0.1)), 'x']
xright = w_bsub$mat[which.min(abs(w_bsub$mat[,'Fdrright'] - 0.1)), 'x']
fdrleft = w_bsub$mat[which.min(abs(w_bsub$mat[,'Fdrleft'] - 0.1)), 'fdr']
fdrright = w_bsub$mat[which.min(abs(w_bsub$mat[,'Fdrright'] - 0.1)), 'fdr']

pval_bsub$fdrleft = fdrleft
pval_bsub$fdrright = fdrright

write.table(select(pval_bsub, -c(sumbern, n, strain)),
            file = paste0(datadir, "strains/B.subtilis-lfdr.txt"),
            sep = "\t", quote = F, row.names = F)

p = ggplot(as.data.frame(w_bsub$mat),
           aes(x = x,
               y = Fdrleft)) +
  geom_point(color = "#ffbb33") +
  geom_line(data = as.data.frame(w_bsub$mat),
            aes(x = x, y = fdr)) +
  theme_bw() +
  geom_hline(yintercept = 0.1, color = "#0066cc", linetype = "dashed") +
  annotate("text", x = -0.2, y = 0.12, label = "FDR = 0.1") + 
  geom_vline(xintercept = xleft,
             linetype = "dotted") +
  annotate("text", x = -0.12, y = fdrright, label = paste("fdr =", 
                                                      format(fdrleft, digits = 2)))+
  geom_point(data = as.data.frame(w_bsub$mat),
             aes(x = x,
                 y = Fdrright),
             color = "#595959") +
  geom_vline(xintercept = xright,
             linetype = "dotted") +
  labs(x = "Effect size",
       y = "Local (solid) and tail-area (points) FDR",
       title = paste(nrow(filter(pval_bsub, (lfdr < fdrleft & medpbliss < 0) | (lfdr < fdrright & medpbliss > 0))), "interactions in B. subtilis")) + 
  annotate("text", x = 0.14, y = fdrright, label = paste("fdr =", 
                                                     format(fdrright, digits = 2)))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(p, filename = "~/Documents/embl/presentations/2019/20190904_statistical_test_drugcombsreen/figures/locfdr-FDR-Bsubtilis.pdf",
       width = 9, height = 7)

# pval_all = mutate(pval_all, binmedpbliss = ifelse(medpbliss > 0, 
#                                                   medpbliss *sumbern / n,
#                                                   medpbliss * (n-sumbern)/n))
# 
# w = locfdr(zz = pval_all$binmedpbliss, bre = 300, df = 12, pct0 = 0.1)

# pval_all = mutate(pval_all,
#                   pi0 = ifelse(medpbliss > 0,
#                                1 - (sumbern + 200)/(n+2000),
#                                1 - (n - sumbern + 200)/(n+2000)))
# pval_all = mutate(pval_all, lfdr_update = lfdr * pi0 / w$fp0["mlest", "p0"])
# View(filter(pval_all, lfdr < .2))

pval_all_common = bind_rows(pval_newman,
                            pval_dsm_nohum,
                            pval_bsub)

pval_gold = filter(pval_all_common, (lfdr < fdrleft & medpbliss < 0) | (lfdr < fdrright & medpbliss > 0))
pval_gold = filter(pval_gold, gsub("(.+)_(.+)", "\\1", comb) != gsub("(.+)_(.+)", "\\2", comb))

saveRDS(pval_gold, file = paste0(datadir, "finalset-nohuman-drugs.rds"))

pval_all = bind_rows(pval_newman,
                            pval_dsm,
                            pval_bsub)
pval_gold = filter(pval_all, (lfdr < fdrleft & medpbliss < 0) | (lfdr < fdrright & medpbliss > 0))
pval_gold = filter(pval_gold, gsub("(.+)_(.+)", "\\1", comb) != gsub("(.+)_(.+)", "\\2", comb))
saveRDS(pval_gold, file = paste0(datadir, "finalset-with-human.rds"))


# normalized rate
pval_all_common = filter(pval_all_common, gsub("(.+)_(.+)", "\\1", comb) != gsub("(.+)_(.+)", "\\2", comb))


syn.obs = nrow(filter(pval_all_common, lfdr < fdrleft & medpbliss < 0)) 
syn.det = nrow(filter(pval_all_common, medpbliss < 0))

ant.obs = nrow(filter(pval_all_common, lfdr < fdrright & medpbliss > 0))
ant.det = nrow(filter(pval_all_common, medpbliss > 0))

bpl = ggplot(data.frame(rate = c(syn.obs/syn.det, ant.obs/ant.det),
                        type = c("Synergy", "Antagonism")),
             aes(x = type, y = rate, fill = type)) +
  geom_bar(stat = 'identity', width = 0.7) +
  scale_fill_manual(values = c("#FFCC33", "#009999")) + 
  labs(x = "",
       y = "Normalized rate (observed / detectable)",
       title = "Normalized rate in Gram-positive strains",
       fill = "") + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16))

ggsave(bpl, filename ="~/Documents/embl/presentations/2019/20190904_statistical_test_drugcombsreen/figures/normalizedrate-grampos.pdf",
       width = 5, height = 6)

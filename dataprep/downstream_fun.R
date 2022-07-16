## Functions for downstream processing of the
## normalized well AUCs and aggregation of the
## screen data
PLATE_EXCL = c("ASA_1_D_1_52_1_06", "ASA_2_D_1_52_2_06",
               "ASA_3_D_1_52_3_06", "PEN_2_D_1_58_11_06",
               "TSB_1_K_2_25_16_03", "TSB_2_K_2_25_17_03",
               "TSB_3_K_2_25_18_03", "CER_2_K_1_6_2_01", "GEN_1_K_1_1_7_01", 
               "LIN_1_K_1_2_13_01", "DOX_1_K_1_7_7_01",
               "AZM_2_K_2_19_2_03", "TEL_1_K_2_24_10_03",
               "MUP_1_K_1_64_7_05", "FUS_1_K_2_64_13_05",
               "FUS_2_K_2_64_14_05",  "ADEP_2_K_2_65_11_05",
               "ORI_1_M_2_42_4_05",  "TEC_3_D_1_79_6_07",
               "TSB_1_D_1_51_16_06", "TSB_2_D_1_51_17_06", 
               "TSB_3_D_1_51_18_06", "CIL_3_M_1_87_15_09",
               "MON_3_M_1_87_12_09", "GEF_2_M_1_82_14_09",
               "CLO_2_M_1_82_8_09", "TRO_2_M_1_86_14_09",
               "TAM_2_M_2_83_8_09", "TSB_1_M_1_88_16_09",
               "TSB_2_M_1_88_17_09", "TSB_3_M_1_88_18_09",
               "HEX_1_M_1_78_7_08", "HEX_2_M_1_78_8_08",
               "HEX_3_M_1_78_9_08", "OLE_1_M_1_87_7_09",
               "OLE_2_M_1_87_8_09", "OLE_3_M_1_87_9_09",
               "PHE_1_M_1_88_1_09", "PHE_2_M_1_88_2_09",
               "PHE_3_M_1_88_3_09")

# exclude BUG wells that grow worse than the first quartile of bug wells
# in the same run
exclude_bug_wells <- function(annotated) {
  # remove non-growing BUG wells
  bugs = lapply(annotated, function(x) filter(x, Drug == "BUG"))
  bugs <- lapply (bugs, AUC_calc_ODno1st)
  bugs_df = plyr::ldply(bugs, data.frame, .id = NULL) %>%
    dplyr::select(-c(Time_h, Time_points, OD, OD_no1st)) %>%
    distinct(Well, plate_name, AUC_well_ODno1st, .keep_all = T) 
  
  bug_quant = dplyr::group_by(bugs_df, run) %>% 
    dplyr::summarise(AUCq1 = quantile(AUC_well_ODno1st, 0.25, na.rm = T))
  
  dead_bugs = dplyr::inner_join(bugs_df, bug_quant) %>%
    dplyr::filter(AUC_well_ODno1st < AUCq1) %>%
    dplyr::select(Well, plate_name) %>%
    dplyr::filter(plate_name %in% c("CTX_1_K_1_2_7_01",
                                    "TEL_3_M_2_27_12_0",
                                    "TMP_3_M_2_27_15_0",
                                    "CLR_3_M_1_28_6_03",
                                    "MMC_1_M_2_39_4_04",
                                    "VAN_1_M_2_47_16_05",
                                    "LOS_1_M_1_78_4_0",
                                    "LOS_2_M_1_78_5_0",
                                    "CLA_2_D_2_84_7_07",
                                    "OLE_1_M_1_87_7_09",
                                    "OLE_2_M_1_87_8_09",
                                    "OLE_3_M_1_87_9_09"))
  
  annotated = lapply(annotated, function(df) {
    pl = as.character(unique(df$plate_name))
    if(nrow(dead_bugs)) {
      excl_well = as.character(dplyr::filter(dead_bugs, plate_name == pl)$Well)
      if(length(excl_well)) df = dplyr::filter(df, ! Well %in% excl_well)
    }
    df
  })
  
  annotated
}

compute_plate_stats <- function(annotated, params) {
  # trim AUC curves at this time point
  trimtp = params$trimtp
  # exclude tp
  if ("excludetp" %in% names(params)) excludetp = params$excludetp

  max_hour = ifelse(strain %in% c("K", "M"), 8, ifelse(strain == 'P', 5, 7.8))
  annotated = lapply(annotated, function(x) subset(x, x$Time_h < max_hour))
  
  # exclude time points if needed
  if(exists("excludetp")) {
    annotated = lapply(annotated, function(x) {
      if(!unique(x$run) %in% c(82,83)) subset(x, !x$Time_points %in% excludetp)
      else x
    })
  }
  annotated <- exclude_bug_wells(annotated)
  annotated <- lapply (annotated, fix_nan)
  annotated <- lapply (annotated, AUC_calc_ODno1st)
  annotated = lapply(annotated, function(x) {
    dplyr::group_by(x, Well) %>%
      dplyr::mutate(time_diff = abs(Time_h - max_hour) ) %>%
      dplyr::filter(time_diff <= min(time_diff)) %>%
      dplyr::select(-time_diff)
  })
  annotated = lapply(annotated, POC_bugwells_AUC)
  annotated
  
}

compute_plate_stats_OD <- function(annotated, params, strain) {
  # trim AUC curves at this time point
  trimtp = params$trimtp
  # exclude tp
  if ("excludetp" %in% names(params)) excludetp = params$excludetp
 
  max_hour = ifelse(strain %in% c("K", "M"), 8, ifelse(strain == 'P', 5, 7.8))
  # this has to be added just to make sure that run 70 gets trimmed at the right time point
  annotated = lapply(annotated, function(x) subset(x, x$Time_h < max_hour))
  
  # exclude time points if needed
  if(exists("excludetp")) {
    annotated = lapply(annotated, function(x) {
      if(!unique(x$run) %in% c(82,83)) subset(x, !x$Time_points %in% excludetp)
      else x
    })
  }
  annotated <- exclude_bug_wells(annotated)
  annotated <- lapply (annotated, fix_nan)
  annotated = lapply(annotated, function(x) {
    dplyr::group_by(x, Well) %>%
      dplyr::mutate(time_diff = abs(Time_h - max_hour) ) %>%
      dplyr::filter(time_diff <= min(time_diff)) %>%
      dplyr::select(-time_diff) %>%
      ungroup()
  })
  annotated = lapply(annotated, POC_bugwells_OD)
  annotated
  
}

# remove plates with poor correlation
get_cor <- function(data) {
  rep_cor_drugs = inner_join(filter(data, Replicate== 1 & donor != "TSB"),
                             filter(data, Replicate == 2 & donor != "TSB"), 
                             by = c("plate", "Drug", "Concentration", "biol_rep", "batch"))
  
  rep_cor_drugs = group_by(rep_cor_drugs, plate) %>%
    dplyr::summarise(correl = cor(ODnorm_bug_rob.x, ODnorm_bug_rob.y))
  
  
  rep_cor_tsb = inner_join(filter(data, Replicate== 1 & donor == "TSB"),
                           filter(data, Replicate == 2 & donor == "TSB"),
                           by = c("plate", "Drug", "Concentration", "biol_rep", "batch"))
  rep_cor_tsb = group_by(rep_cor_tsb, plate) %>%
    dplyr::summarise(correl = cor(ODnorm_bug_rob.x, ODnorm_bug_rob.y))
  
  rep_cor = bind_rows(rep_cor_drugs, rep_cor_tsb)
  
  rep_cor
}

filter_data <- function(ODdf, params) {
  # trim OD curves at this time point
  trimtp = params$trimtp
  # maximum normalized fitness (BUG control: fitness = 1)
  maxfit = params$maxfit
  # exclude tp
  if ("excludetp" %in% names(params)) excludetp = params$excludetp
  # minimum medium of the plate (to remove dead plates)
  minmedplate = params$minmedplate
  # maximum medium of the plate (to remove overgrowing plates)
  maxmedplate = params$maxmedplate
  # minimum IQR of a plate
  iqrplate = params$iqrplate
  # minimum correlation
  mincor = params$mincor
  topwells = params$topwells
  
  # exclude some plates
  write(PLATE_EXCL, file = "manually.excluded")
  ODdf = dplyr::filter(ODdf, !plate %in% PLATE_EXCL)
  # exclude runs 15 and 18 (with dead control wells)
  ODdf = subset(ODdf, ! run %in% c(15, 18))
  
  plate_mean <- dplyr::group_by(ODdf, plate) %>%
    dplyr::summarize(med_plate = median(ODnorm_bug_rob),
                     iqr_plate = EnvStats::iqr(ODnorm_bug_rob))
  
  
  rmed = as.character(dplyr::filter(plate_mean,
                                    med_plate < minmedplate & iqr_plate < iqrplate)$plate)
  
  # filter dead plates
  plate_mean = dplyr::filter(plate_mean, med_plate >= minmedplate & iqr_plate >= iqrplate)
  
  message("Filtered out the following dead plates:")
  cat(rmed, sep="\n")
  
  ODdf = inner_join(ODdf, plate_mean, by = "plate")
  
  # subset overgrowing wells
  ODdf.overgr = dplyr::filter(ODdf, med_plate > maxmedplate) %>%
    dplyr::filter(donor != "TSB" & Drug != "BUG") %>%
    dplyr::filter(ODnorm_bug_rob > maxmedplate)
  
  rep_cor = get_cor(ODdf)
  rep_cor = inner_join(rep_cor, plate_mean)
  
  message("Filtered out the plates with low replicate correlation:")
  cat(as.character(dplyr::filter(rep_cor, correl < mincor & med_plate > 0.25)$plate), sep="\n")
  
  rep_cor = dplyr::filter(rep_cor, correl >= mincor | med_plate <= 0.25)
  ODdf = dplyr::inner_join(ODdf, rep_cor)
  
  # set an upper bound on fitness (ODnorm_bug_rob)
  ODdf = dplyr::mutate(ODdf,
                        ODnorm_bug_rob = ifelse(ODnorm_bug_rob <= maxfit, ODnorm_bug_rob, maxfit))
  
  ODdf
  
}

# the first for code_inter_scores_new.R
get.top.perc <- function(data, p) {
  top.wells = group_by(data, plate) %>%
    filter(fit_comb > quantile(fit_comb, p, na.rm=T)) %>%
    group_by(Well) %>%
    dplyr::summarise(n=n()) %>% filter(n > quantile(n, p))
  
  annot = distinct(data, Well, Drug, flag_fixed)
  
  top.drugs = inner_join(top.wells, annot) %>%
    arrange(desc(n))
  
  don.estim = inner_join(data, top.wells) %>%
    group_by(plate, donor, donor_conc) %>%
    dplyr::summarise(fit_don = median(fit_comb))
  
  don.estim
  
}


#-----------------------------------------------------
# Only for checkerboards
# the second version is for checkerboard donor fitness
get.top.perc2 <- function(data, p) {
  top.wells = group_by(data, batch, plate) %>%
    filter(AUCnorm_bug_rob > quantile(AUCnorm_bug_rob, p, na.rm=T)) %>%
    group_by(batch, Well) %>%
    dplyr::summarise(n=n()) %>% filter(n > quantile(n, p))
  
  
  don.estim = inner_join(data, top.wells) %>%
    group_by(batch, plate, donor, donor_conc) %>%
    dplyr::summarise(fit.don.med = median(AUCnorm_bug_rob))
  
  don.estim
  
}

# have to equalize mask_excluded_wells with filter_data
mask_excluded_wells <- function(data, params) {
  # trim OD curves at this time point
  trimtp = params$trimtp
  # maximum normalized fitness (BUG control: fitness = 1)
  maxfit = params$maxfit
  # exclude tp
  if ("excludetp" %in% names(params)) excludetp = params$excludetp
  # minimum medium of the plate (to remove dead plates)
  minmedplate = params$minmedplate
  # maximum medium of the plate (to remove overgrowing plates)
  maxmedplate = params$maxmedplate
  # minimum IQR of a plate
  iqrplate = params$iqrplate
  # minimum correlation
  mincor = params$mincor
  topwells = params$topwells
  
  # set an upper bound on fitness (ODnorm_bug_rob)
  data = dplyr::mutate(data,
                             ODnorm_bug_rob = ifelse(ODnorm_bug_rob < maxfit,
                                                      ODnorm_bug_rob, maxfit))
  data = filter(data, ! run %in% c(15, 18))
  
  dead_plates = dplyr::group_by(data, batch, plate) %>%
    dplyr::filter(median(ODnorm_bug_rob) <= minmedplate,
                  EnvStats::iqr(ODnorm_bug_rob) <= iqrplate) %>%
    distinct(batch, plate)
  
  overgr_wells = dplyr::filter(data, donor != "TSB" & Drug != "BUG") %>% 
    dplyr::group_by(batch, plate) %>%
    dplyr::filter(median(ODnorm_bug_rob) >= maxmedplate) %>%
    dplyr::filter(ODnorm_bug_rob >= maxmedplate) %>%
    distinct(batch, plate, Well)
  
  rep_cor = get_cor(data)
  low_cor = dplyr::filter(rep_cor, correl < mincor)
  
  data = dplyr::mutate(data, ODnorm_bug_rob = ifelse(plate %in% dead_plates$plate,
                                                      NA, ODnorm_bug_rob))
  data = dplyr::mutate(data, 
                       ODnorm_bug_rob = ifelse(plate %in% overgr_wells$plate & 
                                                  Well %in% overgr_wells$Well, NA, ODnorm_bug_rob))
  data = dplyr::mutate(data,
                       ODnorm_bug_rob = ifelse(plate %in% low_cor$plate,
                                                NA, ODnorm_bug_rob))
  data = dplyr::mutate(data,
                       ODnorm_bug_rob = ifelse(plate %in% PLATE_EXCL,
                                                NA, ODnorm_bug_rob))
  data
  
}

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
  combfit = distinct(batchdata, Drug, Concentration, Replicate, donor, donor_conc, biol_rep, fit_comb)
  # f_rd matrix
  f_rd = mutate(combfit, rown = paste(Drug, Concentration, Replicate, sep = "_"),
                coln = paste(donor, donor_conc, biol_rep, sep = "_")) %>%
    dplyr::select(-c(Drug, Concentration, Replicate, donor, donor_conc, biol_rep)) %>%
    distinct(rown, coln, .keep_all = T) %>%
    tidyr::spread(key = coln, value = fit_comb)
  
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
  donreccombs = mutate(donreccombs, eps = fit_comb - fit_rec * fit_don)
  
  donreccombs
}

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
  combfit = distinct(batchdata, Drug, Concentration, Replicate, donor, donor_conc, biol_rep, ODnorm_bug_rob)
  # f_rd matrix
  f_rd = mutate(combfit, rown = paste(Drug, Concentration, Replicate, sep = "_"),
                coln = paste(donor, donor_conc, biol_rep, sep = "_")) %>%
    dplyr::select(-c(Drug, Concentration, Replicate, donor, donor_conc, biol_rep)) %>%
    distinct(rown, coln, .keep_all = T) %>%
    tidyr::spread(key = coln, value = ODnorm_bug_rob)
  
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
    tibble(ODnorm_bug_rob = res$par[1:ncol(f_rd)],
           donconc = colnames(f_rd))
  } else {
    tibble(ODnorm_bug_rob = res$par[(ncol(f_rd)+1):length(res$par)],
           recconcrep = rownames(f_rd))
  }
  
}

#------------------------------------------
# Defunct methods (not in use anymore)
get_plate_median <- function(time6) {
  plate_mean <- dplyr::group_by(time6, plate) %>%
    dplyr::summarize(med_plate = median(ODnorm_bug_rob),
                     iqr_plate = EnvStats::iqr(ODnorm_bug_rob))
  # filter dead plates
  plate_mean = dplyr::filter(plate_mean, med_plate > minmedplate & iqr_plate > iqrplate)
  plate_mean
}
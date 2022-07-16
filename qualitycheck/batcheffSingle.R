# Single Drug Fitness of the two strains
argv = commandArgs(trailingOnly = TRUE)
# argv[1] = strain data directory

strain = gsub("(.+)_([A-Z])", "\\2", argv[1])

params = rjson::fromJSON(file = "params.json")

pacman::p_load(dplyr, reshape2, plyr, stringr, gtools, Rmisc,
               data.table, RColorBrewer, zoo, hexbin, tools,
               directlabels, grid, ggpubr, gridExtra,
               readr, ggplot2, stringi, ggrepel, viridis,
               scales, tidyr, sigr, lemon, here, smoothmest,EnvStats)


# trim AUC curves at this time point
trimtp = params$trimtp
# maximum normalized fitness (BUG control: fitness = 1)
maxfit = params$maxfit
# minimum medium of the plate (to remove dead plates)
minmedplate = params$minmedplate
# maximum medium of the plate (to remove overgrowing plates)
maxmedplate = params$maxmedplate
# minimum IQR of a plate
iqrplate = params$iqrplate
# minimum plate-level replicate correlation
mincor = params$mincor

figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/integrated/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/integrated/")
  dir.create(figdir)
}

outdir = "~/Documents/embl/gitlab/drug_comb_screen/data/"


script_dir = "~/Documents/embl/gitlab/drug_comb_screen/dataprep/"
if(!dir.exists(script_dir)) script_dir = "dataprep/"

# source the script with core functions
source(paste0(script_dir, "code_screen_core_fun.R"))

compute_plate_stats <- function(annotated) {
  #trim at 6th time point
  annotated = lapply(annotated, function(x) subset(x, x$Time_points <= trimtp))
  
  # remove non-growing BUG wells in batch 5 from plate VAN_1
  ## THIS IS BAD CODING. REWRITE THIS IN A MORE USEABLE/SCALABLE MANNER
  if("VAN_1_M_2_47_16_05" %in% names(annotated)) {
    rm.bug = dplyr::filter(annotated$VAN_1_M_2_47_16_05, Drug == "BUG") %>% 
      dplyr::group_by(Well) %>%
      dplyr::summarise(meanbug = mean(OD_no1st)) %>%
      dplyr::filter(meanbug < 0.1)
    
    annotated$VAN_1_M_2_47_16_05 = dplyr::filter(annotated$VAN_1_M_2_47_16_05,
                                                 !Well %in% rm.bug$Well)
  }
  ## END bad code
  
  # DO basic stats: mean median etc and normalisation
  annotated <- lapply (annotated, POC_bugwells)
  annotated <- lapply (annotated, fix_nan)
  annotated <- lapply (annotated, AUC_calc_ODno1st)
  #proper ranking function (120 is 15*8 time points) --> wells ranked according to AUCs
  annotated <- lapply (annotated, POC_top4)
  annotated <- lapply (annotated, fix_nan)
  annotated <- lapply (annotated, AUC_calc_poc_top)
  annotated <- lapply (annotated, AUC_calc_poc_bug)
  #this are the normalisation functions, per plate
  #by median of BUG controls
  annotated <- lapply (annotated, AUC_med_bug)
  #or by median of top 4% wells
  annotated <- lapply (annotated, AUC_med_top)
  
  #calculate robust mean of BUG wells and normalise for that
  robustMean <- function(x){
    if(length(x) == 1){return(x)}
    else if (unique(x) == 0) {return (mean(x))}
    else{
      return(smhuber(x)$mu)
    } 
  }
  
  POC_bugwells <- function (df) {
    plate_mean_bug <- dplyr::filter(df, Drug == "BUG") %>%
      dplyr::group_by(Time_points) %>%
      dplyr::summarise(mean_AUC_bug = mean(AUC_well_ODno1st), med_AUC_bug = median(AUC_well_ODno1st),
                       mean_AUC_bug_rob = robustMean(AUC_well_ODno1st))
    
    df = inner_join(df, plate_mean_bug, by = "Time_points") %>% mutate(AUCnorm_bug_mean = AUC_well_ODno1st/mean_AUC_bug,
                                                                       AUCnorm_bug_rob = AUC_well_ODno1st/mean_AUC_bug_rob,
                                                                       AUCnorm_bug_med = AUC_well_ODno1st/med_AUC_bug)
    return(df)
  }
  
  annotated = lapply(annotated, POC_bugwells)
  
  annotated
}

get.fitness.data <- function(home_dir) {
  
  batches = grep("batch", list.files(home_dir), value = T)
  
  alldata = lapply(batches, function(b) {
    # load annotated data
    batch <- list.files(paste0(home_dir, b))
    # process batch 1 data
    batch_annot = lapply(paste0(home_dir, b, "/", batch), read.csv, header=T,
                         check.names = FALSE)
    names(batch_annot) <- file_path_sans_ext(batch)
    batch_annot = compute_plate_stats(annotated = batch_annot)
    
    batch_annot = plyr::ldply(batch_annot, data.frame, .id = "plate")
    batch_annot
  })
  
  names(alldata) = batches
  
  alldata_df = plyr::ldply(alldata, data.frame, .id=NULL)
  
  # set an upper bound on fitness (AUCnorm_bug_rob)
  alldata_df = dplyr::mutate(alldata_df,
                             AUCnorm_bug_rob = ifelse(AUCnorm_bug_rob < maxfit,
                                                      AUCnorm_bug_rob, maxfit))
  
  alldata_df = dplyr::distinct(alldata_df, plate, plate_name, Well, Drug, Concentration,
                               Replicate, run, donor, donor_conc, biol_rep, batch,
                               flag_fixed, AUCnorm_bug_rob)
  
  alldata_df
}

get_cor <- function(data) {
  rep_cor_drugs = inner_join(filter(data, Replicate== 1 & donor != "TSB"),
                             filter(data, Replicate == 2 & donor != "TSB"), 
                             by = c("plate_name", "Drug", "Concentration", "biol_rep", "batch"))
  
  rep_cor_drugs = group_by(rep_cor_drugs, batch, plate_name) %>%
    dplyr::summarise(correl = cor(AUCnorm_bug_rob.x, AUCnorm_bug_rob.y))
  
  
  rep_cor_tsb = inner_join(filter(data, Replicate== 1 & donor == "TSB"),
                           filter(data, Replicate == 2 & donor == "TSB"),
                           by = c("plate_name", "Drug", "Concentration", "biol_rep", "batch"))
  rep_cor_tsb = group_by(rep_cor_tsb, batch, plate_name) %>%
    dplyr::summarise(correl = cor(AUCnorm_bug_rob.x, AUCnorm_bug_rob.y))
  
  rep_cor = bind_rows(rep_cor_drugs, rep_cor_tsb)
  
  rep_cor
}


get_plate_median <- function(time6) {
  plate_mean <- dplyr::group_by(time6, plate_name) %>%
    dplyr::summarize(med_plate = median(AUCnorm_bug_rob),
                     iqr_plate = iqr(AUCnorm_bug_rob))
  # filter dead plates
  plate_mean = dplyr::filter(plate_mean, med_plate > minmedplate & iqr_plate > iqrplate)
  plate_mean
}


get.top.perc <- function(data, p) {
  top.wells = group_by(data, batch, plate_name) %>%
    filter(AUCnorm_bug_rob > quantile(AUCnorm_bug_rob, p, na.rm=T)) %>%
    group_by(batch, Well) %>%
    dplyr::summarise(n=n()) %>% filter(n > quantile(n, p))
  
  
  don.estim = inner_join(data, top.wells) %>%
    group_by(batch, plate_name, donor, donor_conc) %>%
    dplyr::summarise(fit.don.med = median(AUCnorm_bug_rob))
  
  don.estim
  
}


filter_fitness_data <- function(data) {
  nondead_plates = dplyr::group_by(data, batch, plate_name) %>%
    dplyr::filter(median(AUCnorm_bug_rob) > minmedplate,
                  iqr(AUCnorm_bug_rob) > iqrplate) %>%
    distinct(batch, plate_name)
  
  
  data_filt = dplyr::inner_join(data, nondead_plates)
  
  overgr_wells = dplyr::filter(data, donor != "TSB" & Drug != "BUG") %>% 
    dplyr::group_by(batch, plate_name) %>%
    dplyr::filter(median(AUCnorm_bug_rob) >= maxmedplate) %>%
    dplyr::filter(AUCnorm_bug_rob >= maxmedplate) %>%
    distinct(batch, plate_name, Well)
  
  data_filt = dplyr::anti_join(data_filt, overgr_wells)
  
  rep_cor = get_cor(data_filt)
  high_cor = dplyr::filter(rep_cor, correl > 0.5)
  data_filt = dplyr::inner_join(data_filt, high_cor)
  

  data_filt
  
}


get_single_fitness <- function(data_filt, p=0.8) {
  combs = dplyr::filter(data_filt, donor != "TSB" & !Drug %in% c("ONE", "BUG", "NOT"))
  donfit = get.top.perc(combs, p=p)

  combs = dplyr::select(combs, Well, Drug, Concentration, Replicate, 
                        #run, # only for M
                        donor, donor_conc, plate_name, AUCnorm_bug_rob, biol_rep)

  donrec = inner_join(combs, donfit )

  # filter "noisy" wells
  noisy.wells = dplyr::group_by(donrec, batch, Well) %>%
    dplyr::summarise(combdon.cor = cor(AUCnorm_bug_rob, fit.don.med)) %>%
    dplyr::filter(combdon.cor > mincor)

  donrec = dplyr::inner_join(donrec, noisy.wells)
  
  donrec = dplyr::mutate(donrec, recconc = paste(Drug, Concentration, sep="_"))
  
  # cor.text = distinct(donrec, batch, recconc, combdon.cor, Replicate)
  # 
  # 
  # pdf(file = paste0(figdir, "slope-fitting", strain, ".pdf"),
  #     width = 12, height = 12)
  # pages = ceiling(nrow(distinct(donrec, recconc, Drug, batch)) / 36)
  # for (i in seq_len(pages)) {
  #   page = ggplot(donrec, aes(x = fit.don.med,
  #                             y= AUCnorm_bug_rob, 
  #                             colour = factor(Replicate))) +
  #     geom_point() +
  #     xlim(c(0,1)) + ylim(c(0,1)) + 
  #     geom_text(data = filter(cor.text, Replicate == 1), aes(label = round(combdon.cor, digits = 2)),
  #               x = 0.1, y=1) + 
  #     geom_text(data = filter(cor.text, Replicate == 2), aes(label = round(combdon.cor, digits = 2)),
  #               x = 0.1, y=0.8) + 
  #     facet_wrap_paginate(~ batch + recconc, ncol = 6, nrow = 6,
  #                         page = i) +
  #     theme_bw() +
  #     theme(aspect.ratio = 1)
  #   print(page)
  # }
  # 
  # dev.off()

  # continue here
  donrec = dplyr::group_by(donrec, batch, recconc, Replicate) %>%
    do(data.frame(fit.rec.est = coef(lm(AUCnorm_bug_rob ~ fit.don.med,data = .))[2])) %>%
    dplyr::left_join(donrec,.) %>%
    dplyr::distinct(batch, plate_name, Drug, Concentration, donor, donor_conc, Replicate, .keep_all=T)
  
  
  # here need to restrict the estimated slope fitness of recipients to (0,1)
  donrec = dplyr::mutate(donrec, 
                             fit.rec.est = ifelse(fit.rec.est < maxfit, fit.rec.est, maxfit))
  donrec
  
}

data.fit = get.fitness.data(home_dir = paste0(argv[1], 
                                            "/screen",
                                            strain,
                                            "_annotated_reads/"))

data_filt = filter_fitness_data(data = data.fit)


# load donor concentrations
donconc = read.table(file.path(argv[1], paste0("donconc", strain, ".txt")),
           sep = "\t", header = T, stringsAsFactors = F)


donconc = inner_join(donconc, distinct(data_filt, run, batch) )

donrec = get_single_fitness(data_filt = data_filt)

donrec = inner_join(donrec, donconc)

library(ggforce)
pdf(file = paste0(figdir, "dose-repsone-don-estim", strain, ".pdf"),
    width = 12, height = 12)
pages = ceiling(length(unique(donrec$donor)) / 16)
                
for (i in seq_len(pages)){
  page = ggplot(donrec, aes(x = donor_real_conc, 
                                           y = fit.don.med,
                                           colour = factor(batch),
                                           shape = factor(biol_rep))) +
    geom_point() +
    labs(x = 'Donor concentration [uM]',
         y = "Estimated donor fitness",
         colour = "Batch",
         shape = "Biological replicate") + 
    facet_wrap_paginate(~ donor, scales = "free_x", nrow = 4, ncol = 4, page = i) +
    theme_bw() +
    theme(aspect.ratio = 1)
  
  print(page)

}

dev.off()

donfit = dplyr::distinct(donrec, donor, donor_conc,
                         donor_real_conc, biol_rep, fit.don.med, plate_name)

cross_biol_rep <- function(df, x, y, by =c("donor", "donor_real_conc")) {
  rep1 = dplyr::filter(df, biol_rep == x)
  rep2 = dplyr::filter(df, biol_rep == y)
  
  repcor =inner_join(rep1, rep2, by = by)
  repcor
}

# matrix of all combinations
cmat = utils::combn(1:length(unique(donfit$biol_rep)),2)
repcor = lapply(1:ncol(cmat), function(i) cross_biol_rep(donfit, x = cmat[1,i], y = cmat[2,i]))
repcor = plyr::ldply(repcor, 'data.frame')


pears.cor = cor(repcor$fit.don.med.x, repcor$fit.don.med.y)

don.est.cor = ggplot(repcor, aes(x = fit.don.med.x,
                   y = fit.don.med.y)) +
  geom_point() +
  geom_label(aes(label = paste(donor, round(donor_real_conc,digits = 2),sep="_"))) + 
  geom_abline(slope = 1, linetype = "dotted") + 
  labs(x = "Estimated donor fitness (biol_rep = x)",
       y = "Estimated donor fitness (biol_rep = y)",
       title = paste("Donor estimate correlation across biol rep: r =", format(pears.cor, digits = 2), ",strain", strain)) + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

ggsave(don.est.cor, filename = paste0(figdir, "donor-estimate-correlation-", strain, ".pdf"))


# do the same for ONE fitness

donors = subset (data_filt, Drug == "ONE")

donors = inner_join(donconc, donors)

onefit = dplyr::group_by(donors,batch, biol_rep,donor, donor_real_conc) %>%
  dplyr::summarise(fit_one = median(AUCnorm_bug_rob))


repcor = lapply(1:ncol(cmat), function(i) cross_biol_rep(onefit, x = cmat[1,i], y = cmat[2,i]))
repcor = plyr::ldply(repcor, 'data.frame')

pears.cor = cor(repcor$fit_one.x, repcor$fit_one.y)

don.one.cor = ggplot(repcor, aes(x = fit_one.x,
                                 y = fit_one.y)) +
  geom_point() +
  geom_label(aes(label = paste(donor, round(donor_real_conc,digits = 2),sep="_"))) + 
  geom_abline(slope = 1, linetype = "dotted") + 
  labs(x = "ONE donor fitness (biol_rep = x)",
       y = "ONE donor fitness (biol_rep = y)",
       title = paste("Donor estimate correlation across biol rep: r =", format(pears.cor, digits = 2), ",strain", strain)) + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

ggsave(don.one.cor, filename = paste0(figdir, "donor-ONE-correlation-", strain, ".pdf"))


pdf(file = paste0(figdir, "dose-repsone-ONE-", strain, ".pdf"),
    width = 12, height = 12)
pages = ceiling(length(unique(onefit$donor)) / 16)

for (i in seq_len(pages)){
  page = ggplot(onefit, aes(x = donor_real_conc, 
                            y = fit_one,
                            colour = factor(batch),
                            shape = factor(biol_rep))) +
    geom_point() +
    labs(x = 'Donor concentration [uM]',
         y = "ONE donor fitness",
         colour = "Batch",
         shape = "Biological replicate") + 
    facet_wrap_paginate(~ donor, scales = "free_x", nrow = 4, ncol = 4, page = i) +
    theme_bw() +
    theme(aspect.ratio = 1)
  
  print(page)
  
}

dev.off()


# do the same for estimated recipient fitness
pdf(file = paste0(figdir, "dose-repsone-rec-estim", strain, ".pdf"),
    width = 12, height = 12)
pages = ceiling(length(unique(donrec$Drug)) / 16)

for (i in seq_len(pages)){
  page = ggplot(donrec, aes(x = Concentration, 
                            y = fit.rec.est,
                            colour = factor(batch),
                            shape = factor(biol_rep))) +
    geom_point() +
    labs(x = 'Recipient concentration [uM]',
         y = "Estimated recipient fitness",
         colour = "Batch",
         shape = "Biological replicate") + 
    facet_wrap_paginate(~ Drug, scales = "free_x", nrow = 4, ncol = 4, page = i) +
    theme_bw() +
    theme(aspect.ratio = 1)
  
  print(page)
  
}

dev.off()

recfit = distinct(donrec, Drug, Concentration,  Replicate, batch, biol_rep, fit.rec.est, plate_name)
recfit = dplyr::group_by(recfit, Drug, Concentration, batch, biol_rep) %>%
  dplyr::summarise(fit.rec.est.av = mean(fit.rec.est))

cross_batches <- function(df, x, y, by =c("Drug", "Concentration")) {
  rep1 = dplyr::filter(df, batch == x)
  rep2 = dplyr::filter(df, batch == y)
  
  repcor =inner_join(rep1, rep2, by = by)
  repcor
}

# matrix of all combinations
cmat = utils::combn(unique(recfit$batch),2)


repcor = lapply(1:ncol(cmat), function(i) cross_batches(recfit, 
                                                         x = cmat[1,i],
                                                         y = cmat[2,i]))
repcor = plyr::ldply(repcor, 'data.frame')

pears.cor = cor(repcor$fit.rec.est.av.x, repcor$fit.rec.est.av.y)



rec.est.cor = ggplot(repcor, aes(x = fit.rec.est.av.x,
                                 y = fit.rec.est.av.y)) +
  geom_point() +
  # ggrepel::geom_text_repel(aes(label = ifelse(abs(fit.rec.est.av.x - fit.rec.est.av.y) > 0.25 & fit.rec.est.av.x < 0.9 & fit.rec.est.av.y < 0.9,
  #                               paste(Drug, round(Concentration,digits = 2),sep="_"),
  #                               ""))) + 
  geom_abline(slope = 1, linetype = "dotted") + 
  labs(x = "Estimated recipient fitness (batch = x)",
       y = "Estimated recipient fitness (batch = y)",
       title = paste("Recipient estimate correlation across batches: r =", format(pears.cor, digits = 2), ",strain", strain)) + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

ggsave(rec.est.cor, filename = paste0(figdir, "recipient-estimate-correlation-", strain, ".pdf"))


# do the same for TSB (batch fitness)



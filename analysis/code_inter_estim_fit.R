#get new rob mean data and from this extract
#rec fitness as mean rob norm AUC of a well across all tsb plates
#don fitness as mean rob norm AUC of ONE wells
#comb fitness as rob norm AUC of the 2 replicates in each combination plate

argv = commandArgs(trailingOnly = TRUE)

# argv[1] = data/ directory
# argv[2] = batch num
# argv[3] = strain

params = rjson::fromJSON(file = "params.json")

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

mincor = params$mincor



pacman::p_load(dplyr, reshape2, plyr, stringr, gtools, Rmisc, data.table, RColorBrewer, zoo, hexbin, tools, directlabels, grid, ggpubr, gridExtra,
               readr, ggplot2, stringi, ggrepel, viridis, scales, tidyr, sigr, lemon, here, smoothmest,EnvStats)

strain = paste0("screen", argv[3])

script_dir = "~/Documents/embl/gitlab/drug_comb_screen/dataprep/"
if(!dir.exists(script_dir)) script_dir = "dataprep/"

figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/interactions/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/", strain, "/batch", argv[2], "/interactions/")
  dir.create(figdir)
}

home_dir = "~/Documents/embl/screen_M/screenM_annotated_reads/batch5/"
if(!dir.exists(home_dir)) {
  home_dir = file.path(argv[1], 
                       paste0(strain, "_annotated_reads"),
                       paste0("batch", argv[2], "/"))
}


# directory for writing epsilon's (interaction scores)
eps_dir = "~/Documents/embl/screen_M/interactions/batch5/"
if(!dir.exists(eps_dir)) {
  eps_dir = file.path(argv[1], "interactions",
                      paste0("batch", argv[2], "/"))
  dir.create(eps_dir)
}


# source the script with core functions
source(paste0(script_dir, "code_screen_core_fun.R"))

#list them
# load annotated data
list_all_plates <- list.files(home_dir)

#read them
annotated <- lapply(paste0(home_dir, list_all_plates), read.csv, header=T,
                    check.names = FALSE)
names(annotated) <- file_path_sans_ext(list_all_plates)

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

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

fix_nan <- function (x) {
  x[is.nan(x)] <- 0
  return(x)
}
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

time6 = plyr::ldply(lapply(annotated, function(x) {
  x[x$Time_points == 1,]
}), 'data.frame')


# exclude runs 15 and 18 (with dead control wells)
time6 = subset(time6, run !=15 & run!=18)

plate_mean <- dplyr::group_by(time6, plate_name) %>%
  dplyr::summarize(med_plate = median(AUCnorm_bug_rob),
                   iqr_plate = iqr(AUCnorm_bug_rob))

rmed = as.character(dplyr::filter(plate_mean, med_plate < minmedplate & iqr_plate < iqrplate)$plate_name)

# filter dead plates
plate_mean = dplyr::filter(plate_mean, med_plate >= minmedplate & iqr_plate >= iqrplate)

message("Filtered out the following dead plates:")
cat(rmed, sep="\n")

time6 = inner_join(time6, plate_mean, by = "plate_name")

# subset overgrowing wells
time6.overgr = dplyr::filter(time6, med_plate > maxmedplate) %>%
               dplyr::filter(donor != "TSB" & Drug != "BUG") %>%
               dplyr::filter(AUCnorm_bug_rob > maxmedplate)


if(nrow(time6.overgr)) {
  drugs.overgr = ggplot(time6.overgr, aes(x = factor(Drug),
                                          fill = Drug)) + geom_bar() +
    labs(x = "Drug", title = paste("Overgrowing recipients in batch", argv[2])) +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle=90))
  ggsave(drugs.overgr, filename = paste0(figdir, "overgrowing-recipients.pdf"))
  
  # remove overgrowing wells
  time6 = dplyr::filter(time6, !(plate_name %in% time6.overgr$plate_name & Well %in% time6.overgr$Well))
}

# remove plates with poor correlation
get_cor <- function(data) {
  rep_cor_drugs = inner_join(filter(data, Replicate== 1 & donor != "TSB"),
                             filter(data, Replicate == 2 & donor != "TSB"), 
                             by = c("plate_name", "Drug", "Concentration", "biol_rep", "batch"))
  
  rep_cor_drugs = group_by(rep_cor_drugs, plate_name) %>%
    dplyr::summarise(correl = cor(AUCnorm_bug_rob.x, AUCnorm_bug_rob.y))
  
  
  rep_cor_tsb = inner_join(filter(data, Replicate== 1 & donor == "TSB"),
                           filter(data, Replicate == 2 & donor == "TSB"),
                           by = c("plate_name", "Drug", "Concentration", "biol_rep", "batch"))
  rep_cor_tsb = group_by(rep_cor_tsb, plate_name) %>%
    dplyr::summarise(correl = cor(AUCnorm_bug_rob.x, AUCnorm_bug_rob.y))
  
  rep_cor = bind_rows(rep_cor_drugs, rep_cor_tsb)
  
  rep_cor
}

rep_cor = get_cor(time6)

message("Filtered out the plates with low replicate correlation:")
cat(as.character(dplyr::filter(rep_cor, correl < mincor)$plate_name), sep="\n")

rep_cor = dplyr::filter(rep_cor, correl >= mincor)
time6 = dplyr::inner_join(time6, rep_cor)


#new columns
time6$recid = interaction(time6$Drug, time6$flag_fixed, sep = "_")
time6$recconc = interaction(time6$Drug, time6$Concentration, sep = "_")

# set an upper bound on fitness (AUCnorm_bug_rob)
time6 = dplyr::mutate(time6,
                      AUCnorm_bug_rob = ifelse(AUCnorm_bug_rob < maxfit, AUCnorm_bug_rob, maxfit))

#-------------------------get donor fitness from top 20% growing wells---------------------------
get.top.perc <- function(data, p) {
  top.wells = group_by(data, plate_name) %>%
    filter(AUCnorm_bug_rob > quantile(AUCnorm_bug_rob, p, na.rm=T)) %>%
    group_by(Well) %>%
    dplyr::summarise(n=n()) %>% filter(n > quantile(n, p))
  
  annot = distinct(data, Well, Drug, flag_fixed)
  
  top.drugs = inner_join(top.wells, annot) %>%
    arrange(desc(n))
  
  don.estim = inner_join(data, top.wells) %>%
    group_by(plate_name, donor, donor_conc) %>%
    dplyr::summarise(fit.don.med = median(AUCnorm_bug_rob))
  
  don.estim
  
}

combs = dplyr::filter(time6, donor != "TSB" & !Drug %in% c("ONE", "BUG", "NOT"))
donfit = get.top.perc(combs, p=0.75)

combs = dplyr::select(combs, Drug, Concentration, Replicate,
                      recconc, recid, donor, donor_conc, plate_name, AUCnorm_bug_rob)

donrec = inner_join(combs, donfit )

rmed_wells = dplyr::group_by(donrec, recconc) %>%
  dplyr::summarise(combdon.cor = cor(AUCnorm_bug_rob, fit.don.med)) %>%
  dplyr::filter(combdon.cor < mincor)

message("Filtered out the following noisy wells:")
cat(as.character(rmed_wells$recconc), sep="\n")

# filter "noisy" wells
well.cor = dplyr::group_by(donrec, recconc) %>%
  dplyr::summarise(combdon.cor = cor(AUCnorm_bug_rob, fit.don.med)) %>%
  dplyr::filter(combdon.cor >= mincor)

donrec = inner_join(donrec, well.cor)

# continue here
donreccomb = dplyr::group_by(donrec, Drug, recconc, Replicate) %>%
  do(data.frame(fit.rec.est = coef(lm(AUCnorm_bug_rob ~ fit.don.med,data = .))[2])) %>%
  dplyr::left_join(donrec,.) %>%
  dplyr::distinct(plate_name, Drug, Concentration, donor, donor_conc, Replicate, .keep_all=T)


# here need to restrict the estimated slope fitness of recipients to (0,1)
donreccomb = dplyr::mutate(donreccomb, 
                           fit.rec.est = ifelse(fit.rec.est < maxfit, fit.rec.est, maxfit))

#epsilon
donreccomb = dplyr::mutate(donreccomb, multipl = fit.don.med * fit.rec.est)
donreccomb = dplyr::mutate(donreccomb, eps = AUCnorm_bug_rob - multipl)

#get a format suitable for Wilcoxon
donreccomb$combid_broad = paste0(donreccomb$Drug, "_", donreccomb$donor)

donreccomb$donconc = paste0(donreccomb$donor, "_", donreccomb$donor_conc)
donreccomb$donid = paste0(donreccomb$donor, "_", donreccomb$donor_conc)
# fix the donor flag issue
donreccomb$combid = paste0(donreccomb$donid, "_", donreccomb$recconc)

donreccomb = dplyr::rename(donreccomb, fit_comb = AUCnorm_bug_rob)


#write it down
write.table (donreccomb,
             paste0(eps_dir, "interactions.tsv"),
             sep = "\t", row.names = FALSE)



# check estimated vs measured single drug fitness
tsb = subset(time6, donor == "TSB")
# exclude TSB plates from run 4
tsb = dplyr::filter(tsb, !run %in% c(4, 15, 18, 52))

recs = cbind.data.frame(tsb$Drug, tsb$Concentration, tsb$flag_fixed, 
                        tsb$Drug_Flag_Conc, tsb$recid,
                        tsb$AUC_well_ODno1st, tsb$AUCnorm_bug_rob, tsb$AUCnorm_bug_mean, tsb$AUCnorm_bug_med,
                        tsb$plate_name, tsb$run, tsb$batch, tsb$position)


colnames(recs) <- c("Drug", "recconc", "recflag", "Drug_Flag_Conc", "recid",
                    "AUC_well_ODno1st", "AUCnorm_bug_rob", "AUCnorm_bug", "AUCnorm_bug_med",
                    "plate_name", "run", "batch", "position")

recs$recconc = interaction(recs$Drug, recs$recconc, sep = "_")


tsbfit = dplyr::group_by(recs, Drug, recconc) %>%
  dplyr::summarise(fit_tsb = median(AUCnorm_bug_rob))


tsb.vs.estim = inner_join(donreccomb, tsbfit)


est.check = ggplot(distinct(tsb.vs.estim, Drug, recconc, .keep_all = T), aes(x = fit_tsb, y = fit.rec.est)) + 
  geom_point() +
  labs(x = "Median TSB fitness",
       y = "Estimated recipient fitness") + 
  geom_label(aes(label = recid)) + 
  xlim(c(0,1)) + ylim(c(0,1)) + 
  geom_abline(slope = 1, linetype = "dotted") +
  theme_bw()

ggsave(est.check, filename = paste0(figdir, "Estimated-vs-TSB.pdf"))


# donor fitness vs top 20% growing wells
donors = subset (time6, Drug == "ONE")

onefit = dplyr::group_by(donors, donor, donor_conc) %>%
  dplyr::summarise(fit_one = median(AUCnorm_bug_rob))

one.vs.estim = inner_join(donreccomb, onefit)

est.don.check = ggplot(distinct(one.vs.estim, donor, donor_conc, .keep_all = T), aes(x = fit_one,
                                                                             y = fit.don.med)) + 
  geom_point() +
  labs(x = "Median ONE fitness",
       y = "Estimated donor fitness") + 
  xlim(c(0,1)) + ylim(c(0,1)) + 
  geom_abline(slope = 1, linetype = "dotted") +
  geom_label(aes(label = donid)) + 
  theme_bw()

ggsave(est.don.check, filename = paste0(figdir, "Estimated-vs-ONE.pdf"))

#---alternative estimate of recipient fitness from combo plates (like Ana)
# combs = subset(time6, donor != "TSB")
# 
# plate = combs %>% dplyr::group_by(plate_name) %>%
#   filter(Drug == "ONE") %>%
#   dplyr::summarize(mean_donor = median(AUC_well_ODno1st)) %>%
#   dplyr::ungroup()
# 
# combs = inner_join(plate,combs, by = 'plate_name')
# 
# combs$donid = interaction(combs$donor, combs$donor_conc, sep = "_")
# 
# comb_df = data.table(combs)
# 
# 
# #get this estimate of recipient fitness as Ana did: the slope of the best fit between the fitness of the combination vs the fitness of the single donor
# rec_fit = comb_df[,list(intercept=coef(lm(AUCnorm_bug_rob~mean_donor))[1], coef=coef(lm(AUCnorm_bug_rob~mean_donor))[2]),by=Drug_Flag_Conc]
# 
# rec_fit = inner_join (rec_fit, recs, by = "Drug_Flag_Conc")
# rec_fit = filter(rec_fit, !Drug %in% c("BUG", "NOT"))
# 
# pdf(paste0(figdir, "recTSBvscombs.pdf"), width = 10, height = 10)
# 
# ggplot(distinct(rec_fit, Drug, recconc, .keep_all = T), aes(x= AUCr_norm_bug_rob, y = coef)) +
#   geom_point() +
#   coord_cartesian() +
#   labs(x = "rec from TSB", y = "rec from best fit combs", 
#        title="Fitness from tsb or combs") +
#   geom_abline(slope = 1) +
#   #geom_text_repel(data = data, aes(label=recid), size = 3) +
#   #geom_text_repel(data=subset(distinct(final, Drug, recconc, .keep_all = T), coef < 0.75) , aes(label=Drug), size = 3) +
#   geom_smooth(method = "lm", color="red",
#               formula = y ~ x, linetype="dashed", size = 0.5) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5))+
#   #use this to play with lines and fetch outliers
#   #geom_abline(slope = 0.998, intercept = 0.1) +
#   stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE)
# 
# dev.off()
# 
# 
# #decided to make this plot for all of them (good quality control) --> this has to be included in QC
# library(ggforce)
# comb_df = filter(comb_df, !Drug %in% c("BUG", "NOT", "ONE"))
# 
# pdf(paste0(figdir, "fitcombvs_fitdon.pdf"), width = 20, height = 20)
# pages = ceiling(length(unique(comb_df$Drug_Flag_Conc)) / 36)
# 
# for (i in seq_len(pages)){
#   page = ggplot(comb_df,aes(x = mean_donor, y = AUCnorm_bug_rob)) +
#     facet_wrap_paginate(~ Drug_Flag_Conc, ncol =  6, nrow = 6, page = i) +
#     geom_point(aes(colour = factor(Replicate)), alpha = 0.7) +
#     labs(x = "fit donor", y = "fit comb",
#          color = "Replicate",
#          title="Fit Comb vs Fit donor") +
#     coord_cartesian() +
#     geom_smooth(method = 'lm', colour = 'red', formula = y~x, linetype = 'dashed', size = 0.5) +
#     stat_smooth_func(geom="text",method='lm',hjust = 0, parse = T) +
#     theme_bw() + 
#     theme(aspect.ratio = 1,
#           axis.text.x = element_text(size=11, angle=60, hjust=1),
#           axis.text.y = element_text(size=11),
#           strip.text.x = element_text(size=15),
#           panel.spacing = unit(1, "lines"))
#   print(page)
# }
# 
# dev.off()
# 
# 



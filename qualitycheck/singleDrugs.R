# Single Drug Fitness of the two strains
argv = commandArgs(trailingOnly = TRUE)
# argv[1] = strain 1 data directory
# argv[2] = strain 2 data directory
# argv[3] = strain 3 data directory

strain_1 = gsub("(.+)_([A-Z])", "\\2", argv[1])
strain_2 = gsub("(.+)_([A-Z])", "\\2", argv[2])
strain_3 = gsub("(.+)_([A-Z])", "\\2", argv[3])


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
if(!dir.exists(outdir)) {
  outdir = paste0("data/", strain_1, "-vs-", strain_2, "/")
  if(!dir.exists(outdir)) dir.create(outdir)
}

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
  
  alldata_df = dplyr::distinct(alldata_df, plate, Well, Drug, Concentration,
                        Replicate, run, donor, donor_conc, biol_rep, batch,
                        flag_fixed, AUCnorm_bug_rob)
  
  alldata_df
}

get_plate_median <- function(time6) {
  plate_mean <- dplyr::group_by(time6, plate) %>%
    dplyr::summarize(med_plate = median(AUCnorm_bug_rob),
                     iqr_plate = iqr(AUCnorm_bug_rob))
  
  # filter dead plates
  plate_mean = dplyr::filter(plate_mean, med_plate > minmedplate & iqr_plate > iqrplate)
  
  plate_mean
}

get_cor <- function(data, strain) {
  rep_cor_drugs = inner_join(filter(data, Replicate== 1 & donor != "TSB"),
                             filter(data, Replicate == 2 & donor != "TSB"), by = c("plate", "Drug", "Concentration", "biol_rep", "batch"))
  
  rep_cor_drugs = group_by(rep_cor_drugs, plate) %>%
    dplyr::summarise(correl = cor(AUCnorm_bug_rob.x, AUCnorm_bug_rob.y))
  
  rep_cor_drugs$strain = paste("Double-drug plates", strain)
  
  
  rep_cor_tsb = inner_join(filter(data, Replicate== 1 & donor == "TSB"),
                           filter(data, Replicate == 2 & donor == "TSB"), by = c("plate", "Drug", "Concentration", "biol_rep", "batch"))
  rep_cor_tsb = group_by(rep_cor_tsb, plate) %>%
    dplyr::summarise(correl = cor(AUCnorm_bug_rob.x, AUCnorm_bug_rob.y))
  
  rep_cor_tsb$strain = paste("TSB", strain)
  
  rep_cor = bind_rows(rep_cor_drugs, rep_cor_tsb)
  
  rep_cor
}


data.1 = get.fitness.data(home_dir = paste0(argv[1], 
                                            "/screen",
                                            strain_1,
                                            "_annotated_reads/"))
# dead plates already filtered out
plate_med = get_plate_median(data.1)

data.1 = inner_join(data.1, plate_med, by = "plate")

# subset overgrowing wells
data.1.overgr = dplyr::filter(data.1, med_plate >= maxmedplate) %>%
  dplyr::filter(donor != "TSB" & Drug != "BUG") %>%
  dplyr::filter(AUCnorm_bug_rob >= maxmedplate)

drugs.overgr = ggplot(data.1.overgr, aes(x = factor(Drug),
                                        fill = Drug)) + geom_bar() +
  labs(x = "Drug", title = paste("Overgrowing recipients in strain", strain_1)) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle=90))
ggsave(drugs.overgr, filename = paste0(figdir, "overgrowing-recipients-", strain_1, ".pdf"))



# remove overgrowing wells
data.1 = dplyr::filter(data.1, !(plate %in% data.1.overgr$plate & Well %in% data.1.overgr$Well))



rec.tsb1 = filter(data.1, donor == "TSB" & !(Drug %in% c("BUG", "NOT", "ONE")))

# make nice plots for these
rec.overall = ggplot(mutate(rec.tsb1,
              biol_rep = paste("biol_rep", biol_rep)),
       aes(x = AUCnorm_bug_rob, fill = factor(Replicate))) + 
  geom_histogram() +
  labs(fill = "Replicate",
       title = paste("Recipient (TSB) fitness in strain", strain_1)) + 
  facet_wrap(~ biol_rep) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

ggsave(rec.overall, 
       filename = paste0(figdir, 'recipient-distrib-', strain_1, ".pdf"),
       width = 10, height = 6)




pdf(paste0(figdir, "rec-tsb-distrib-by-drug-biolrep-", strain_1, ".pdf"),
    width = 10, height = 10)
pages = ceiling(length(unique(rec.tsb1$Drug)) / 16)

for (i in seq_len(pages)){
  page = ggplot(rec.tsb1,
                aes(x = AUCnorm_bug_rob, fill = factor(biol_rep))) + 
    geom_histogram() +
    labs(title = "TSB fitness distribution by drug and biol_rep",
         fill = "biol_rep") + 
    ggforce::facet_wrap_paginate(~ Drug, ncol=4, nrow=4, page=i) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))
  print(page)
}
dev.off()

pdf(paste0(figdir, "rec-tsb-distrib-by-drug-rep-", strain_1, ".pdf"),
    width = 10, height = 10)
pages = ceiling(length(unique(rec.tsb1$Drug)) / 16)

for (i in seq_len(pages)){
  page = ggplot(rec.tsb1,
                aes(x = AUCnorm_bug_rob, fill = factor(Replicate))) + 
    geom_histogram() +
    labs(title = "TSB fitness distribution by drug and Replicate",
         fill = "Replicate") + 
    ggforce::facet_wrap_paginate(~ Drug, ncol=4, nrow=4, page=i) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))
  print(page)
}
dev.off()


# estimate recipient fitness from the slope

calc_recfit_from_slope <- function(alldata, don) {
  alldata = dplyr::mutate(alldata, donor_flag = paste(donor, donor_conc, sep="_"),
                          recid = paste(Drug, flag_fixed, sep="_"))
  combs = subset(alldata, donor %in% don & Drug != "ONE")
  # combs = dplyr::select(combs, Drug, recconc, recid, donor_flag, AUCnorm_bug_rob)
  # donor fitnesses
  donfit = dplyr::filter(alldata, Drug=="ONE") %>%
    dplyr::group_by(donor_flag, batch, biol_rep, Replicate) %>%
    dplyr::summarise(fit_don_av = mean(AUCnorm_bug_rob)) %>%
    dplyr::rename(recid = donor_flag)
  
  donrec = inner_join(combs, donfit, by = c("recid", "batch", "biol_rep", "Replicate"))
  
  rec_fit_slope = dplyr::group_by(donrec, donor_flag, batch, biol_rep, Replicate) %>%
    do(data.frame(fit_rec_av = coef(lm(AUCnorm_bug_rob ~ fit_don_av,data = .))[2])) %>%
    left_join(donrec,.) %>%
    distinct(donor_flag, batch, biol_rep, Replicate, fit_rec_av) %>%
    dplyr::rename(recid = donor_flag) %>%
    dplyr::mutate(Drug = gsub("(^.+)(_.+)", "\\1", recid)) %>%
    dplyr::mutate(recconc = recid)
  
  rec_fit_slope
  
}


dons = unique(as.character(data.1$donor))
dons = dons[dons!= "TSB"]
rec_slope1 = calc_recfit_from_slope(alldata = data.1,
                                    don = dons)

rec_slope1 = dplyr::group_by(rec_slope1, biol_rep, recid, batch) %>%
  dplyr::mutate(fit_rec_av =mean(fit_rec_av)) %>%
  dplyr::distinct(recid, batch, biol_rep, fit_rec_av)

tsb.av.1 = dplyr::group_by(rec.tsb1, batch, Drug, flag_fixed) %>%
  dplyr::filter(AUCnorm_bug_rob < maxmedplate) %>%
  dplyr::summarise(fit_tsb_av = median(AUCnorm_bug_rob)) %>%
  dplyr::mutate(recid = paste(Drug, flag_fixed, sep="_"))



rec_comp1 = inner_join(rec_slope1, tsb.av.1)

tsb.vs.slope = ggplot(rec_comp1,
       aes(x = fit_rec_av,
           y = fit_tsb_av,
           colour = factor(batch))) + 
  geom_point() +
  labs(title = "Estimated recipient fitness vs TSB fitness: overgrowing wells removed from TSB",
       x = "Estimated recipient fitness",
       y = "Median TSB recipient fitness",
       colour = "batch") + 
  facet_wrap(~ Drug) + 
  theme_bw() + 
  geom_abline(slope=1, linetype = "dotted") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(tsb.vs.slope, filename = paste0(figdir, "Estimated-TSB-strain", strain_1, ".pdf"),
       width = 10, height = 9)


# single donor fitness
one.1 = filter(data.1, Drug=="ONE" & donor!= "TSB")

one_fit = ggplot(one.1, aes(x = AUCnorm_bug_rob)) +
  geom_histogram() +
  labs(title = paste("Single donor (ONE) fitness distributuon in", strain_1) ) + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(one_fit, filename = paste0(figdir, "ONE-fitness-distrib-", strain_1, ".pdf"))

perc = c(0.8, 0.85, 0.9, 0.95, 0.98)
pdf(paste0(figdir, "one-vs-topwells-", strain_1, ".pdf"),
    width = 24, height = 24, onefile = T)
for(p in perc) {
  top.percent = group_by(data.1, plate) %>%
    mutate(top.perc = quantile(AUCnorm_bug_rob, p, na.rm=T)) %>%
    mutate(top.p = ifelse(AUCnorm_bug_rob >= top.perc, T, F))
  
  top.perc = filter(top.percent, top.p) %>%
    dplyr::select(plate, AUCnorm_bug_rob) %>%
    dplyr::summarise(medAUCtop = median(AUCnorm_bug_rob))
  
  one.vs.perc = inner_join(group_by(one.1, plate) %>%
               mutate(med.donor = median(AUCnorm_bug_rob)), 
             top.perc)
  
  
  g = ggplot(one.vs.perc,
         aes(x = med.donor, y = medAUCtop,
             colour = factor(donor_conc))) +
    geom_point(size = 3) +
    labs(x = "Median ONE fitness",
         y = paste("Median AUCnorm of top", (1-p)*100, "%"),
         colour = "Donor concentration",
         title = paste0("Top ", (1-p)*100, "% growing well fitness vs ONE well fitness")) + 
    geom_abline(slope = 1, linetype="dashed") + 
    facet_wrap(~ donor) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(g)
  
  
}
dev.off()


# in one plot

get.top.perc <- function(p) {
  top.percent = group_by(data.1, plate) %>%
    mutate(top.perc = quantile(AUCnorm_bug_rob, p, na.rm=T)) %>%
    mutate(top.p = ifelse(AUCnorm_bug_rob >= top.perc, T, F))
  
  top.perc = filter(top.percent, top.p) %>%
    dplyr::select(plate, AUCnorm_bug_rob) %>%
    dplyr::summarise(medAUCtop = median(AUCnorm_bug_rob))
  
  one.vs.perc = inner_join(group_by(one.1, plate) %>%
                             mutate(med.donor = median(AUCnorm_bug_rob)), 
                           top.perc)
  
  one.vs.perc = mutate(one.vs.perc, which.perc = paste0("top ", (1-p)*100, "%"))
  one.vs.perc
  
}

perc = c(0.8, 0.85, 0.9, 0.95)


top.perc.df = lapply(perc, get.top.perc)
top.perc.df = plyr::ldply(top.perc.df, data.frame)

library(wesanderson)
g = ggplot(top.perc.df,
           aes(x = med.donor, y = medAUCtop,
               colour = factor(which.perc, levels = c("top 5%", "top 10%", "top 15%", "top 20%")))) +
  geom_point(size = 2) +
  scale_color_manual(values = brewer.pal(n=9, "Blues")[c(4,6,8,9)]) + 
  labs(x = "Median ONE fitness",
       y = paste("Median AUCnorm of top growing wells"),
       colour = "",
       title = paste0("Top growing well fitness vs ONE well fitness")) + 
  geom_abline(slope = 1, linetype="dashed") + 
  facet_wrap(~ donor) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(g, filename = paste0(figdir, "one-vs-topwells-", strain_1, ".pdf"),
       width = 16, height = 16)


pdf(paste0(figdir, "one-vs-topwells-", strain_1, ".pdf"),
    width = 24, height = 24, onefile = T)
for(p in perc) {
  top.percent = group_by(data.1, plate) %>%
    mutate(top.perc = quantile(AUCnorm_bug_rob, p, na.rm=T)) %>%
    mutate(top.p = ifelse(AUCnorm_bug_rob >= top.perc, T, F))
  
  top.perc = filter(top.percent, top.p) %>%
    dplyr::select(plate, AUCnorm_bug_rob) %>%
    dplyr::summarise(medAUCtop = median(AUCnorm_bug_rob))
  
  one.vs.perc = inner_join(group_by(one.1, plate) %>%
                             mutate(med.donor = median(AUCnorm_bug_rob)), 
                           top.perc)
  
}
dev.off()

# direct donor measurements: compare with combination fitness
one.med1 = group_by(one.1, plate) %>%
  mutate(med.donor = median(AUCnorm_bug_rob))

donor_fit1 = dplyr::select(one.med1, plate, donor, donor_conc, biol_rep, batch, med.donor)
donor_fit1 = dplyr::distinct(donor_fit1)
donor_fit1 = dplyr::rename(donor_fit1, don_fit = med.donor)
drug.vs.donor1 = inner_join(data.1, donor_fit1)
drug.vs.donor1 = dplyr::filter(drug.vs.donor1, !Drug %in% c("BUG", "NOT", "ONE"))


comb.vs.donor1 = ggplot(drug.vs.donor1, aes(x = don_fit, y = AUCnorm_bug_rob, colour = Well)) +
  scale_color_viridis_d() + 
  geom_point(alpha=0.3) + 
  facet_wrap(~ Drug) +
  labs(x = "Donor fitness",
       y = "Combination fitness",
       title = paste("Combination vs donor fitness by recipient, colored by well in", strain_1)) + 
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

ggsave(comb.vs.donor1, filename = paste0(figdir, "comb-vs-donor-", strain_1, ".pdf"),
       width = 12, height = 12)


# combination vs top 20% growing wells
top.perc = get.top.perc(p=0.8)

top_fit1 = dplyr::select(top.perc, plate, donor, donor_conc, biol_rep, batch, medAUCtop)
top_fit1 = dplyr::rename(top_fit1, don_fit = medAUCtop)
drug.vs.top1 = inner_join(data.1, top_fit1)
drug.vs.top1 = dplyr::filter(drug.vs.top1, !Drug %in% c("BUG", "NOT", "ONE"))


comb.vs.top1 = ggplot(drug.vs.top1, aes(x = don_fit, y = AUCnorm_bug_rob, colour = Well)) +
  scale_color_viridis_d() + 
  geom_point(alpha=0.3) + 
  facet_wrap(~ Drug) +
  labs(x = "Median fitness of top 20% growing wells of a plate",
       y = "Combination fitness",
       title = paste("Combination vs top 20% by recipient, colored by well in", strain_1)) + 
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

ggsave(comb.vs.top1, filename = paste0(figdir, "comb-vs-top-", strain_1, ".pdf"),
       width = 12, height = 12)




drugs = unique(data.1$donor)
drugs = as.character(drugs[drugs!= "TSB"])

pdf(paste0(figdir, "plate-plots-", strain_1, ".pdf"),
    onefile = T)
for(d in drugs) {
  # make a plot in a loop for each drug
  don_data = filter(data.1, donor == d)
  g = ggplot(don_data,
         aes(x = factor(biol_rep),
             y = AUCnorm_bug_rob,
             fill = factor(donor_conc))) +
    labs(fill = "Donor concentration",
         x = "biol_rep",
         title = paste(d, "plate effects in strain", strain_1)) +
  geom_boxplot() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  print(g)
}
dev.off()



# AND DO THE same for strain M

data.2 = get.fitness.data(home_dir = paste0(argv[2], 
                                            "/screen",
                                            strain_2,
                                            "_annotated_reads/"))

# dead plates already filtered out
plate_med = get_plate_median(data.2)

data.2 = inner_join(data.2, plate_med, by = "plate")

# subset overgrowing wells
data.2.overgr = dplyr::filter(data.2, med_plate >= maxmedplate) %>%
  dplyr::filter(donor != "TSB" & Drug != "BUG") %>%
  dplyr::filter(AUCnorm_bug_rob >= maxmedplate)

drugs.overgr = ggplot(data.2.overgr, aes(x = factor(Drug),
                                         fill = Drug)) + geom_bar() +
  labs(x = "Drug", title = paste("Overgrowing recipients in strain", strain_2)) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle=90))
ggsave(drugs.overgr, filename = paste0(figdir, "overgrowing-recipients-", strain_2, ".pdf"))



# remove overgrowing wells
data.2 = dplyr::filter(data.2, !(plate %in% data.2.overgr$plate & Well %in% data.2.overgr$Well))


rec.tsb2 = filter(data.2, donor == "TSB" & !(Drug %in% c("BUG", "NOT", "ONE")))


dons = unique(as.character(data.2$donor))
dons = dons[dons!= "TSB"]
rec_slope2 = calc_recfit_from_slope(alldata = data.2,
                                    don = dons)

tsb.av.2 = dplyr::group_by(rec.tsb2, Replicate, biol_rep, batch, Drug, flag_fixed) %>%
  dplyr::filter(AUCnorm_bug_rob < maxmedplate) %>%
  dplyr::summarise(fit_tsb_av = mean(AUCnorm_bug_rob)) %>%
  dplyr::mutate(recid = paste(Drug, flag_fixed, sep="_"))



rec_comp2 = inner_join(rec_slope2, tsb.av.2)

tsb.vs.slope = ggplot(rec_comp2,
                      aes(x = fit_rec_av,
                          y = fit_tsb_av,
                          colour = factor(biol_rep))) + 
  geom_point() +
  labs(title = "Estimated recipient fitness vs TSB fitness: overgrowing wells removed from TSB",
       x = "Estimated recipient fitness",
       y = "Mean TSB recipient fitness",
       colour = "biol_rep") + 
  facet_wrap(~ Drug) + 
  theme_bw() + 
  geom_abline(slope=1, linetype = "dotted") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(tsb.vs.slope, filename = paste0(figdir, "Estimated-TSB-strain", strain_2, ".pdf"),
       width = 10, height = 9)



# single donor fitness
one.2 = filter(data.2, Drug=="ONE" & donor!= "TSB")


perc = c(0.8, 0.85, 0.9, 0.95)

get.top.perc <- function(p) {
  top.percent = group_by(data.2, plate) %>%
    mutate(top.perc = quantile(AUCnorm_bug_rob, p, na.rm=T)) %>%
    mutate(top.p = ifelse(AUCnorm_bug_rob >= top.perc, T, F))
  
  top.perc = filter(top.percent, top.p) %>%
    dplyr::select(plate, AUCnorm_bug_rob) %>%
    dplyr::summarise(medAUCtop = median(AUCnorm_bug_rob))
  
  one.vs.perc = inner_join(group_by(one.2, plate) %>%
                             mutate(med.donor = median(AUCnorm_bug_rob)), 
                           top.perc)
  
  one.vs.perc = mutate(one.vs.perc, which.perc = paste0("top ", (1-p)*100, "%"))
  one.vs.perc
  
}


top.perc.df = lapply(perc, get.top.perc)
top.perc.df = plyr::ldply(top.perc.df, data.frame)

library(wesanderson)
g = ggplot(top.perc.df,
           aes(x = med.donor, y = medAUCtop,
               colour = factor(which.perc, levels = c("top 5%", "top 10%", "top 15%", "top 20%")))) +
  geom_point(size = 2) +
  scale_color_manual(values = brewer.pal(n=9, "Blues")[c(4,6,8,9)]) + 
  labs(x = "Median ONE fitness",
       y = paste("Median AUCnorm of top growing wells"),
       colour = "",
       title = paste0("Top growing well fitness vs ONE well fitness in strain ", strain_2)) + 
  geom_abline(slope = 1, linetype="dashed") + 
  facet_wrap(~ donor) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(g, filename = paste0(figdir, "one-vs-topwells-", strain_2, ".pdf"),
       width = 16, height = 16)


# direct donor measurements: compare with combination fitness
one.med2 = group_by(one.2, plate) %>%
  mutate(med.donor = median(AUCnorm_bug_rob))

donor_fit2 = dplyr::select(one.med2, plate, donor, donor_conc, biol_rep, batch, med.donor)
donor_fit2 = dplyr::distinct(donor_fit2)
donor_fit2 = dplyr::rename(donor_fit2, don_fit = med.donor)
drug.vs.donor2 = inner_join(data.2, donor_fit2)
drug.vs.donor2 = dplyr::filter(drug.vs.donor2, !Drug %in% c("BUG", "NOT", "ONE"))


comb.vs.donor2 = ggplot(drug.vs.donor2, aes(x = don_fit, y = AUCnorm_bug_rob, colour = Well)) +
  scale_color_viridis_d() + 
  geom_point(alpha=0.3) + 
  facet_wrap(~ Drug) +
  labs(x = "Donor fitness",
       y = "Combination fitness",
       title = paste("Combination vs donor fitness by recipient, colored by well in", strain_2)) + 
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

ggsave(comb.vs.donor2, filename = paste0(figdir, "comb-vs-donor-", strain_2, ".pdf"),
       width = 12, height = 12)


# combination vs top 20% growing wells
top.perc = get.top.perc(p=0.8)

top_fit2 = dplyr::select(top.perc, plate, donor, donor_conc, biol_rep, batch, medAUCtop)
top_fit2 = dplyr::rename(top_fit2, don_fit = medAUCtop)
drug.vs.top2 = inner_join(data.2, top_fit2)
drug.vs.top2 = dplyr::filter(drug.vs.top2, !Drug %in% c("BUG", "NOT", "ONE"))


comb.vs.top2 = ggplot(drug.vs.top2, aes(x = don_fit, y = AUCnorm_bug_rob, colour = Well)) +
  scale_color_viridis_d() + 
  geom_point(alpha=0.3) + 
  facet_wrap(~ Drug) +
  labs(x = "Median fitness of top 20% growing wells of a plate",
       y = "Combination fitness",
       title = paste("Combination vs top 20% by recipient, colored by well in", strain_2)) + 
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

ggsave(comb.vs.top2, filename = paste0(figdir, "comb-vs-top-", strain_2, ".pdf"),
       width = 12, height = 12)



trimtp = 10
data.3 = get.fitness.data(home_dir = paste0(argv[3], 
                                            "/screen",
                                            strain_3,
                                            "_annotated_reads/"))

# dead plates already filtered out
plate_med = get_plate_median(data.3)

data.3 = inner_join(data.3, plate_med, by = "plate")

# subset overgrowing wells
data.3.overgr = dplyr::filter(data.3, med_plate >= maxmedplate) %>%
  dplyr::filter(donor != "TSB" & Drug != "BUG") %>%
  dplyr::filter(AUCnorm_bug_rob >= maxmedplate)

drugs.overgr = ggplot(data.3.overgr, aes(x = factor(Drug),
                                         fill = Drug)) + geom_bar() +
  labs(x = "Drug", title = paste("Overgrowing recipients in strain", strain_3)) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle=90))
ggsave(drugs.overgr, filename = paste0(figdir, "overgrowing-recipients-", strain_3, ".pdf"))



# remove overgrowing wells
data.3 = dplyr::filter(data.3, !(plate %in% data.3.overgr$plate & Well %in% data.3.overgr$Well))

rec.tsb3 = filter(data.3, donor == "TSB" & !(Drug %in% c("BUG", "NOT", "ONE")))

dons = unique(as.character(data.3$donor))
dons = dons[dons!= "TSB"]
rec_slope3 = calc_recfit_from_slope(alldata = data.3,
                                    don = dons)

tsb.av.3 = dplyr::group_by(rec.tsb3, Replicate, biol_rep, batch, Drug, flag_fixed) %>%
  dplyr::filter(AUCnorm_bug_rob < maxmedplate) %>%
  dplyr::summarise(fit_tsb_av = mean(AUCnorm_bug_rob)) %>%
  dplyr::mutate(recid = paste(Drug, flag_fixed, sep="_"))



rec_comp3 = inner_join(rec_slope3, tsb.av.3)

tsb.vs.slope = ggplot(rec_comp3,
                      aes(x = fit_rec_av,
                          y = fit_tsb_av,
                          colour = factor(biol_rep))) + 
  geom_point() +
  labs(title = "Estimated recipient fitness vs TSB fitness: overgrowing wells removed from TSB",
       x = "Estimated recipient fitness",
       y = "Mean TSB recipient fitness",
       colour = "biol_rep") + 
  facet_wrap(~ Drug) + 
  theme_bw() + 
  geom_abline(slope=1, linetype = "dotted") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(tsb.vs.slope, filename = paste0(figdir, "Estimated-TSB-strain", strain_3, ".pdf"),
       width = 10, height = 9)



# single donor fitness
one.3 = filter(data.3, Drug=="ONE" & donor!= "TSB")


perc = c(0.8, 0.85, 0.9, 0.95)

get.top.perc <- function(p) {
  top.percent = group_by(data.3, plate) %>%
    mutate(top.perc = quantile(AUCnorm_bug_rob, p, na.rm=T)) %>%
    mutate(top.p = ifelse(AUCnorm_bug_rob >= top.perc, T, F))
  
  top.perc = filter(top.percent, top.p) %>%
    dplyr::select(plate, AUCnorm_bug_rob) %>%
    dplyr::summarise(medAUCtop = median(AUCnorm_bug_rob))
  
  one.vs.perc = inner_join(group_by(one.3, plate) %>%
                             mutate(med.donor = median(AUCnorm_bug_rob)), 
                           top.perc)
  
  one.vs.perc = mutate(one.vs.perc, which.perc = paste0("top ", (1-p)*100, "%"))
  one.vs.perc
  
}


top.perc.df = lapply(perc, get.top.perc)
top.perc.df = plyr::ldply(top.perc.df, data.frame)

library(wesanderson)
g = ggplot(top.perc.df,
           aes(x = med.donor, y = medAUCtop,
               colour = factor(which.perc, levels = c("top 5%", "top 10%", "top 15%", "top 20%")))) +
  geom_point(size = 2) +
  scale_color_manual(values = brewer.pal(n=9, "Blues")[c(4,6,8,9)]) + 
  labs(x = "Median ONE fitness",
       y = paste("Median AUCnorm of top growing wells"),
       colour = "",
       title = paste0("Top growing well fitness vs ONE well fitness in strain ", strain_3)) + 
  geom_abline(slope = 1, linetype="dashed") + 
  facet_wrap(~ donor) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(g, filename = paste0(figdir, "one-vs-topwells-", strain_3, ".pdf"),
       width = 16, height = 16)


# direct donor measurements: compare with combination fitness
one.med3 = group_by(one.3, plate) %>%
  mutate(med.donor = median(AUCnorm_bug_rob))

donor_fit3 = dplyr::select(one.med3, plate, donor, donor_conc, biol_rep, batch, med.donor)
donor_fit3 = dplyr::distinct(donor_fit3)
donor_fit3 = dplyr::rename(donor_fit3, don_fit = med.donor)
drug.vs.donor3 = inner_join(data.3, donor_fit3)
drug.vs.donor3 = dplyr::filter(drug.vs.donor3, !Drug %in% c("BUG", "NOT", "ONE"))


comb.vs.donor3 = ggplot(drug.vs.donor3, aes(x = don_fit, y = AUCnorm_bug_rob, colour = Well)) +
  scale_color_viridis_d() + 
  geom_point(alpha=0.3) + 
  facet_wrap(~ Drug) +
  labs(x = "Donor fitness",
       y = "Combination fitness",
       title = paste("Combination vs donor fitness by recipient, colored by well in", strain_3)) + 
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

ggsave(comb.vs.donor3, filename = paste0(figdir, "comb-vs-donor-", strain_3, ".pdf"),
       width = 12, height = 12)


# combination vs top 20% growing wells
top.perc = get.top.perc(p=0.8)

top_fit3 = dplyr::select(top.perc, plate, donor, donor_conc, biol_rep, batch, medAUCtop)
top_fit3 = dplyr::rename(top_fit3, don_fit = medAUCtop)
drug.vs.top3 = inner_join(data.3, top_fit3)
drug.vs.top3 = dplyr::filter(drug.vs.top3, !Drug %in% c("BUG", "NOT", "ONE"))


comb.vs.top3 = ggplot(drug.vs.top3, aes(x = don_fit, y = AUCnorm_bug_rob, colour = Well)) +
  scale_color_viridis_d() + 
  geom_point(alpha=0.3) + 
  facet_wrap(~ Drug) +
  labs(x = "Median fitness of top 20% growing wells of a plate",
       y = "Combination fitness",
       title = paste("Combination vs top 20% by recipient, colored by well in", strain_3)) + 
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

ggsave(comb.vs.top3, filename = paste0(figdir, "comb-vs-top-", strain_3, ".pdf"),
       width = 12, height = 12)


rep_cor1 = get_cor(data = data.1, strain = strain_1)
rep_cor2 = get_cor(data = data.2, strain = strain_2)
rep_cor3 = get_cor(data = data.3, strain = strain_3)

rep_cor = bind_rows(rep_cor1, rep_cor2, 
                    filter(rep_cor3, correl > 0.3))

low_cor_plates = filter(rep_cor, correl < 0.5)

write.table(low_cor_plates, file = "low-correlation-plates.txt",
            row.names = F, col.names = T, quote = F)

cor_plot = ggplot(rep_cor,
       aes(y=correl, x=factor(strain,
                              levels = c("Double-drug plates K","TSB K","Double-drug plates M", "TSB M",
                                         "Double-drug plates D", "TSB D" )),
           fill = factor(strain,
                         levels = c("Double-drug plates K","TSB K","Double-drug plates M", "TSB M",
                                    "Double-drug plates D", "TSB D" )))) + 
  scale_fill_manual(values = brewer.pal(6, "Paired")) + 
  geom_boxplot() + 
  labs(x = "", 
       y = "Pearson correlation",
       title = "Plate-level replicate well correlation") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

ggsave(cor_plot,
       filename = paste0(figdir, "correlation-boxplot-all-strains.pdf"),
       width = 8, height = 6)


rep_cor = bind_rows(rep_cor1, rep_cor2)

cor_plot = ggplot(rep_cor,
                  aes(y=correl, x=factor(strain,
                                         levels = c("Double-drug plates K","TSB K","Double-drug plates M", "TSB M")),
                      fill = factor(strain,
                                    levels = c("Double-drug plates K","TSB K","Double-drug plates M", "TSB M")))) + 
  scale_fill_manual(values = brewer.pal(6, "Paired")) + 
  geom_boxplot() + 
  labs(x = "", 
       y = "Pearson correlation",
       title = "Plate-level replicate well correlation") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

ggsave(cor_plot,
       filename = paste0(figdir, "correlation-boxplot-all-K-and-M.pdf"),
       width = 8, height = 6)




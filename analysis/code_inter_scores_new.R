#get new rob mean data and from this extract
#rec fitness as mean rob norm AUC of a well across all tsb plates
#don fitness as mean rob norm AUC of ONE wells
#comb fitness as rob norm AUC of the 2 replicates in each combination plate
argv = commandArgs(trailingOnly = TRUE)
# argv[1] = data/ directory
# argv[2] = batch num
# argv[3] = strain
# argv = c("/Volumes/gitlab/parallel/pneumo-main/data/screen_P",
#          10, "P")
params = rjson::fromJSON(file = paste0("params/batch", argv[2], ".json" ))

pacman::p_load(dplyr, reshape2, plyr, stringr, gtools, Rmisc,
               data.table, RColorBrewer, zoo, hexbin, tools, 
               directlabels, grid, ggpubr, gridExtra,
               readr, ggplot2, stringi, ggrepel, viridis,
               scales, tidyr, sigr, lemon, here, smoothmest,EnvStats)

strain = paste0("screen", argv[3])

script_dir = "~/Documents/embl/gitlab/drug_comb_screen/dataprep/"
if(!dir.exists(script_dir)) script_dir = "dataprep/"

figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/interactions/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/", strain, "/batch", argv[2], "/interactions/")
  dir.create(figdir)
}

home_dir = file.path(argv[1], 
                     paste0(strain, "_annotated_reads"),
                     paste0("batch", argv[2], "/"))

# directory for writing epsilon's (interaction scores)
eps_dir = file.path(argv[1], "interactions",
                    paste0("batch", argv[2], "/"))
if(!dir.exists(eps_dir)) dir.create(eps_dir)

# source the script with core functions
source(paste0(script_dir, "code_screen_core_fun.R"))
# downstream processing of screen results
source(paste0(script_dir, "downstream_fun.R"))

# load annotated data
list_all_plates <- list.files(home_dir)

#read them
annotated <- lapply(paste0(home_dir, list_all_plates), read.csv, header=T,
                    check.names = FALSE)
names(annotated) <- file_path_sans_ext(list_all_plates)
annotated = compute_plate_stats_OD(annotated, params, strain=argv[3])


ODdf = plyr::ldply(lapply(annotated, function(x) {
  dplyr::distinct(x, plate_name, Well, Drug, Concentration, Replicate, flag_fixed,
                  run, donor, donor_conc, biol_rep, batch, ODnorm_bug_rob)
}), 'data.frame', .id = "plate")
ODdf = filter_data(ODdf, params)

# load real donor concentrations
donconc = read.table(file.path(argv[1], paste0("donconc", argv[3], ".txt")),
                     sep = "\t", header = T, stringsAsFactors = F)

combs = dplyr::filter(ODdf, donor != "TSB" & !Drug %in% c("ONE", "BUG", "NOT"))
if(nrow(dplyr::inner_join(combs, donconc)) < length(unique(combs$plate))) {
  warning("Missing real donor concentrations")
  message("Using donor_conc as donor_real_conc")
  donconc = dplyr::distinct(combs, donor, run, donor_conc) %>%
    dplyr::mutate(donor_real_conc = donor_conc)
}
combs = dplyr::inner_join(combs, donconc)

combs = dplyr::rename(combs, fit_comb = ODnorm_bug_rob)
donreccomb = optimize_epsilon(combs, norm = "L2")

#-------recipient fitness------

#extract tsb plates and get recipient behaviour from here
tsb = subset(ODdf, donor == "TSB")
# exclude TSB plates from run 4
tsb = dplyr::filter(tsb, !run %in% c(4, 15, 18, 52))

recfit = dplyr::group_by(tsb, Well, Drug, Concentration, Replicate, flag_fixed) %>%
  dplyr::summarise(fit_rec = median(ODnorm_bug_rob))

#-------------------------get donor fitness alone---------------------------
donors = subset(ODdf, Drug == "ONE")
donors = dplyr::inner_join(donors, donconc)

onefit = dplyr::group_by(donors, donor, donor_real_conc, biol_rep, plate) %>%
  dplyr::summarise(fit_don = median(ODnorm_bug_rob))


donconc_map = dplyr::select(donors, donor, donor_conc, donor_real_conc, biol_rep) %>%
  dplyr::distinct()
topfit = dplyr::inner_join(donreccomb,  donconc_map) %>%
  dplyr::distinct(donor, donor_real_conc, biol_rep, fit_don)


one.vs.top = ggplot(inner_join(onefit, topfit, by = c("donor", "donor_real_conc", "biol_rep")),
       aes(x=fit_don.x, y=fit_don.y)) +
  labs(x = "ONE donor fitness",
       y = "Estimated donor fitness (L2-norm)",
       title = paste("Estimated donor fitness vs ONE well fitneess in strain",
                     argv[3], "batch", argv[2])) + 
  geom_point() +
  coord_equal() +
  xlim(c(0,1)) + ylim(c(0,1)) +
  ggrepel::geom_text_repel(aes(label = ifelse(abs(fit_don.x - fit_don.y) > 0.1 & fit_don.y < 0.8,
                                              paste(donor, donor_real_conc, sep="_"), "")),
                           segment.colour = NA) +
  geom_abline(slope = 1, linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(one.vs.top, filename = paste0(figdir, "ONE-vs-top.pdf"))

# join combination, recipient and donor fitness
# reccomb = dplyr::inner_join(combs, recfit)
# donreccomb = inner_join(reccomb, topfit)

# filter noisy wells
# wellcor = dplyr::group_by(donreccomb, Well, Drug, Concentration) %>%
#   dplyr::mutate(combdon.cor = cor(fit_comb, fit_don)) %>%
#   dplyr::distinct(Well, Drug, Concentration, Replicate, fit_rec, combdon.cor) %>%
#   dplyr::mutate(row = gsub("^([A-Z])([0-9]+)", "\\1", Well),
#                 col = gsub("^([A-Z])([0-9]+)", "\\2", Well))
# 
# wellcor.plot = ggplot(wellcor,
#        aes(x = factor(col, levels = 1:24),
#            y = factor(row, levels = rev(LETTERS[1:16])), fill = combdon.cor)) +
#   geom_tile(colour = "white", size = 0.25) +
#   scale_fill_gradient(low="blue", high = "white") +
#   labs(x = "", y = "",
#        title = paste("Noisy wells: combination-donor correlation for", argv[3], "batch", argv[2]),
#        fill = "fitcomb-fitdon correlation") + 
#   coord_fixed() + 
#   theme_grey(base_size=12)
# ggsave(wellcor.plot, filename = paste0(figdir, "noisy-wells-mincor05.pdf"))

# rep1 = dplyr::filter(donreccomb, Replicate == 1)
# rep2 = dplyr::filter(donreccomb, Replicate == 2)
# 
# combrepcor = dplyr::inner_join(rep1, rep2,
#                                by = c("Drug", "Concentration", "donor", "donor_conc")) %>%
#   dplyr::group_by(Drug, Concentration) %>%
#   dplyr::summarise(repcor = cor(fit_comb.x, fit_comb.y))
# 
# noisywells = dplyr::filter(combrepcor, repcor < 0.5)
# # noisywells = inner_join(noisywells, group_by(recfit, Drug, Concentration) %>%
# #   dplyr::summarise(fit_rec = mean(fit_rec))) %>%
# #   filter(fit_rec > 0.1)
# 
# write.table(noisywells,
#           file = paste0(eps_dir, "noisy-wells-removed.tsv"),
#           sep = "\t", quote = F, row.names = F)

#epsilon
donreccomb$multipl <- donreccomb$fit_rec * donreccomb$fit_don
#donreccomb$eps <- donreccomb$fit_comb - donreccomb$multipl

#get a format suitable for Wilcoxon
donreccomb$combid_broad = paste0(donreccomb$Drug, "_", donreccomb$donor)

donreccomb$recconc = interaction(donreccomb$Drug, donreccomb$Concentration, sep = "_")

donreccomb = dplyr::inner_join(donreccomb,  donconc_map)

donreccomb$donconc = interaction(donreccomb$donor, donreccomb$donor_real_conc, sep = "_")
donreccomb$combid = paste0(donreccomb$donconc, "_", donreccomb$recconc)

donreccomb = dplyr::distinct(donreccomb, combid, combid_broad,
                             donconc, recconc, Replicate, biol_rep,
                             donor, Drug, fit_comb, fit_don, fit_rec,
                             multipl, eps)
#write it
write.table (donreccomb,
             paste0(eps_dir, "interactions.tsv"),
             sep = "\t", row.names = FALSE)

estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

# estimate recipient fitness from the slope of fit_comb vs fit_don
estim.recfit = dplyr::distinct(donreccomb, Drug, recconc, Replicate, fit_rec) %>%
  dplyr::rename(fit.rec.est = fit_rec) %>%
  mutate(Concentration = as.numeric(gsub("(.+)_([0-9]+)", "\\2", recconc)) )

tsb.vs.estim = inner_join(estim.recfit, recfit)

est.check = ggplot(distinct(tsb.vs.estim, Drug, recconc, .keep_all = T), aes(x = fit_rec, y = fit.rec.est)) + 
  geom_point() +
  labs(x = "TSB fitness",
       y = "Estimated recipient fitness",
       title = "Estimated vs TSB recipient fitness, batch", argv[2]) + 
  geom_text(aes(label = ifelse(abs(fit.rec.est - fit_rec) > 0.25, paste(Drug, Concentration, sep="_"), ""))) + 
  xlim(c(0,1)) + ylim(c(0,1)) + 
  geom_abline(slope = 1, linetype = "dotted") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(est.check, filename = paste0(figdir, "batch", argv[2], "-estimated-vs-TSB.pdf"))

#tsb.vs.estim = dplyr::mutate(tsb.vs.estim, diff = fit.rec.est - fit_rec)

null.dist = ggplot(donreccomb, aes(x = eps)) +
  geom_density(fill = "blue", colour=NA, alpha=0.3) + 
  theme_classic() + 
  geom_vline(xintercept = median(donreccomb$eps),
             colour = "red") + 
  geom_vline(xintercept = estimate_mode(donreccomb$eps),
             colour = "black",
             linetype = "dotted") + 
  labs(x = paste("Bliss scores in strain",
                 strain),
       y = "Count",
       title = paste("Strain", argv[3], "batch", argv[2],
                     "Bliss distrib. Median:",
                     format(median(donreccomb$eps), digits = 2),
                     "Mode:", format(estimate_mode(donreccomb$eps),digits = 2))) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(null.dist, filename = paste0(figdir, "batch", argv[2], "-eps-distribution.pdf"))

# reorder the levels of the 'recconc' factor variable
newlevels = sort(levels(donreccomb$recconc))
donreccomb = dplyr::mutate(donreccomb, recconc = factor(recconc, levels = newlevels))

# #decided to make this plot for all of them (good quality control) --> this has to be included in QC
library(ggforce)
pdf(paste0(figdir, "fitcombvs_fitdon.pdf"), width = 20, height = 20)
pages = ceiling(length(unique(donreccomb$recconc)) / 36)
for (i in seq_len(pages)){
  page = ggplot(donreccomb, aes(x = fit_don, y = fit_comb)) +
    facet_wrap_paginate(~ recconc, ncol =  6, nrow = 6, page = i) +
    geom_point(aes(colour = factor(Replicate)), alpha = 0.7) +
    labs(x = "Donor fitness", y = "Combination fitness",
         color = "Replicate",
         title="Combination vs donor fitness") +
    coord_cartesian() +
    geom_smooth(method = 'lm', colour = 'red', formula = y~x, linetype = 'dashed', size = 0.5) +
    stat_smooth_func(geom="text",method='lm',hjust = 0, parse = T) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(size=11, angle=60, hjust=1),
          axis.text.y = element_text(size=11),
          strip.text.x = element_text(size=15),
          panel.spacing = unit(1, "lines"))
  print(page)
}

dev.off()







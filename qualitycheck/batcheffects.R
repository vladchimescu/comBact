# Assess batch, run, stacking, shaking effects
# on fitness and reproducibility
argv = commandArgs(trailingOnly = TRUE)
# argv[1] = data/ directory
# argv[2] = strain
pacman::p_load(dplyr, reshape2, plyr, stringr, gtools, Rmisc,
               data.table, RColorBrewer, zoo, hexbin, tools,
               directlabels, grid, ggpubr, gridExtra,
               readr, ggplot2, stringi, ggrepel, viridis,
               scales, tidyr, sigr, lemon, here, smoothmest,EnvStats)

strain = paste0("screen", argv[2])

script_dir = "~/Documents/embl/gitlab/drug_comb_screen/dataprep/"
if(!dir.exists(script_dir)) script_dir = "dataprep/"

figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/batcheffects/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/", strain, "/batcheffects/")
  dir.create(figdir, recursive = T)
}

# this directory has all the pre-processed batches from the same strain
home_dir = "~/Documents/embl/screen_M/screenM_annotated_reads/"
if(!dir.exists(home_dir)) {
  home_dir = file.path(argv[1], 
                       paste0(strain, "_annotated_reads/"))
}

# source the script with core functions
source(paste0(script_dir, "code_screen_core_fun.R"))
# downstream processing of screen results
source(paste0(script_dir, "downstream_fun.R"))

# now we have the complete annotated data set from batch 1 and 2
# we can first load them separately and merge them into one object
# to compare fitness of perfect replicates

#========LOAD BATCH DATA=============================

batches = grep("batch", list.files(home_dir), value = T)

alldata = lapply(batches, function(b) {
  # load annotated data
  batch <- list.files(paste0(home_dir, b))
  # process batch 1 data
  batch_annot = lapply(paste0(home_dir, b, "/", batch), read.csv, header=T,
                       check.names = FALSE)
  names(batch_annot) <- file_path_sans_ext(batch)
  
  params = rjson::fromJSON(file = paste0("params/",b, ".json" ))
  batch_annot = compute_plate_stats(annotated = batch_annot,
                                    params = params)
  batch_annot = plyr::ldply(batch_annot, data.frame, .id = "plate")
  batch_annot
})

names(alldata) = batches


bugs = lapply(alldata, function(x) filter(x, Drug == "BUG"))
bugs_df = plyr::ldply(bugs, data.frame, .id = NULL)
bugs_df = distinct(bugs_df, plate, Drug, Concentration, Replicate, donor, donor_conc, .keep_all = T)

# let's check now the bug wells in batches 1 and 2
bug_by_run = ggplot(bugs_df, aes(x = factor(run),
                    y = AUC_well_ODno1st,
                    fill = factor(run))) + 
  geom_boxplot() + 
  scale_fill_viridis(discrete = T, option='D') + 
  theme_bw() + 
  labs(x = "Experimental run",
       title = "BUG well fitness by run") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(bug_by_run, filename = paste0(figdir, "bug-by-run.pdf"))


# remove the dead runs
bug_by_pos = ggplot(bugs_df, aes(x = factor(run),
                                 y = AUC_well_ODno1st,
                                 fill = factor(position))) + 
  geom_boxplot(position = position_dodge(preserve = "single"),
               lwd=0.05, outlier.shape = NA) + 
  scale_fill_viridis(discrete = T, option='B') + 
  #scale_x_discrete(breaks = c(1:14,16:17)) +
  theme_bw() + 
  labs(x = "Experimental run",
       fill = "Position",
       title = "BUG well fitness by run (colored by position)") + 
  coord_polar(start=-pi/8) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") 
ggsave(bug_by_pos, 
       filename = paste0(figdir, "bug-by-position.pdf"),
       width = 7, height = 7)


# compare TSB plates 
tsb = lapply(alldata, function(x) filter(x, donor == "TSB"))
tsb_all = plyr::ldply(tsb, data.frame, .id=NULL)

## exclude bad runs 15 and 18 (as we already know that these
## would produce bad results)
tsb_all = filter(tsb_all, ! run %in% c(15,18))

run.ef.all = ggplot(tsb_all, aes(x = factor(run), 
                             y = AUCnorm_bug_rob,
                             fill = factor(position))) + 
  scale_fill_viridis(discrete = T, option='A') + 
  geom_boxplot() + 
  scale_y_continuous(limits = quantile(tsb_all$AUCnorm_bug_rob, c(0.01, 0.99))) +
  labs(x = "Experimental run",
       fill = "Position",
       title = "TSB Plates split by run") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(run.ef.all, 
       filename = paste0(figdir, "TSB-by-run.pdf"),
       width = 9)


# within batch concordance for TSB plates
#tsb1 = filter(tsb_all, batch == 1)
rep1 = filter(tsb_all, Replicate == 1)
rep2 = filter(tsb_all, Replicate == 2)

tsb_rep = inner_join(rep1, rep2,
                     by = c("plate","Drug", "Concentration",
                            "Time_points", "run", "batch"))
tsb_rep = mutate(tsb_rep, OD_diff = OD.x - OD.y,
                 AUC_diff = AUCnorm_bug_rob.x -AUCnorm_bug_rob.y)

tsb_rep = filter(tsb_rep, Drug != "BUG")



batch_num = sort(gsub("batch", "", batches))

pdf(paste0(figdir, "TSB-repcor.pdf"),
    width = 10, height = 6,
    onefile = T)
for(i in batch_num) {
  p = ggplot(filter(tsb_rep, batch == i),
             aes(x = Drug, y = OD_diff,
                 colour = Drug)) +
    geom_boxplot(width=0.5) + 
    theme_classic() + 
    labs(x = 'Recipient',
         y = "Replicate OD difference",
         title = paste("TSB plate replicate concordance in batch", i)) + 
    geom_hline(yintercept = 0,
               linetype = "dotted") + 
    scale_y_continuous(limits = quantile(filter(tsb_rep, batch==i)$OD_diff, c(0.01, 0.99))) +
    theme(legend.position = 'none',
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle=90),
          plot.title = element_text(hjust = 0.5))
  print(p)
}

dev.off()

pdf(paste0(figdir,  "TSB-repTp.pdf"),
    onefile = T)

for (i in batch_num) {
  p = ggplot(filter(tsb_rep, batch == i),
             aes(x = factor(Time_points),
                 y = OD_diff,
                 color = factor(Time_points))) + 
    geom_boxplot() +
    scale_color_viridis(discrete = T) +
    labs(x = "Time points",
         y = "Replicate OD difference",
         title = paste("TSB replicate concordance by run in batch", i)) + 
    facet_wrap(~ run) + 
    theme_bw() + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1)
  print(p)
}

dev.off()


pdf(paste0(figdir,  "TSB-AUC_replicates.pdf"),
    width = 10, height = 6,
    onefile = T)

for (i in batch_num) {
  p = ggplot(filter(tsb_rep, batch == i & Time_points == 6),
             aes(x = Drug,
                 y = AUC_diff,
                 colour = Drug)) +
    geom_boxplot() + 
    labs(x = '',
         y = "Replicate AUCnorm_bug_rob difference",
         title = paste("TSB replicate AUC_norm_bug diff. by recipient in batch", i)) + 
    theme_classic() + 
    geom_hline(yintercept = 0,
               linetype = "dotted") + 
    scale_y_continuous(limits = quantile(filter(tsb_rep, batch==i & Time_points == 6)$AUC_diff, c(0.01, 0.99))) +
    theme(legend.position = 'none',
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle=90),
          plot.title = element_text(hjust = 0.5))
  print(p)
}

dev.off()



# check the concordance using scatterplots
tsb_sct = ggplot(tsb_rep,
       aes(x = OD.x,
           y = OD.y,
           colour = factor(run))) + 
  labs(x = "OD replicate 1",
       y = "OD replciate 2",
       title= "TSB OD replicate scatter by run") + 
  geom_point() + 
  geom_abline(slope = 1, colour = "black") + 
  theme_bw() + 
  facet_wrap(~ run) + 
  coord_equal() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

ggsave(tsb_sct, filename = paste0(figdir,
                                   "TSB-scatter-OD.pdf"),
       width = 10, height = 10)


tsb_sct = ggplot(filter(tsb_rep, Time_points == 6),
                 aes(x = AUCnorm_bug_rob.x,
                     y = AUCnorm_bug_rob.y,
                     colour = factor(run))) + 
  labs(x = "AUCnorm_bug_rob replicate 1",
       y = "AUCnrom_bug_rob replciate 2",
       title= "TSB normalized AUC replicate scatter by run") + 
  geom_point() + 
  geom_abline(slope = 1, colour = "black") + 
  theme_bw() + 
  facet_wrap(~ run) + 
  coord_equal() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

ggsave(tsb_sct, filename = paste0(figdir,
                                  "TSB-scatter-AUC.pdf"),
       width = 10, height = 10)


#-----------ACROSS BATCH COMPARISON-----------------------
stopifnot(length(batches) > 1)
# we can also compare the measurements across batches
# for TSB plates (as these are always replicates of each other)
cm = combinations(n = length(batches), v = batch_num, r=2)

pdf(paste0(figdir, "TSB-across-batches.pdf"),
    width = 12, height = 9,
    onefile = T)
for (i in 1:nrow(cm)) {
  tsb_both = inner_join(filter(tsb_all, batch == cm[i,1]),
                        filter(tsb_all, batch == cm[i,2]),
                        by = c("Drug", "Concentration"))
  
  p = ggplot(distinct(tsb_both,
                                  plate.x, plate.y,
                                  Drug, Concentration, .keep_all = T),
                         aes(x = AUCnorm_bug_rob.x,
                             y = AUCnorm_bug_rob.y)) + 
    geom_point(alpha = 0.5) + 
    geom_abline(slope = 1, colour = "red") + 
    facet_wrap(~ run.x + run.y, scales = "free") + 
    theme_bw() + 
    labs(x = paste("AUCnorm_bug_rob, batch", cm[i,1]),
         y = paste("AUCnorm_bug_rob, batch", cm[i,2]),
         title = paste("Across batch TSB plate comparison (run by run)",
                       "batch", cm[i,1], "vs", cm[i,2])) + 
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
}
dev.off()


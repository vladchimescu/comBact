## Compare perfect replicates

# Assess batch, run, stacking, shaking effects
# on fitness and reproducibility

pacman::p_load(dplyr, reshape2, plyr, stringr, gtools, Rmisc,
               data.table, RColorBrewer, zoo, hexbin, tools,
               directlabels, grid, ggpubr, gridExtra,
               readr, ggplot2, stringi, ggrepel, viridis,
               scales, tidyr, sigr, lemon, here, smoothmest,EnvStats)

script_dir = "~/Documents/embl/gitlab/drug_comb_screen/dataprep/"
if(!dir.exists(script_dir)) script_dir = "dataprep/"



figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/batcheffects/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/batch", argv[2], "/batcheffects/")
  dir.create(figdir, recursive = T)
}

home_dir = "~/Documents/embl/screen_K/screenK1_annotated_reads/"
if(!dir.exists(home_dir)) {
  home_dir = file.path(argv[1], "screenK1_annotated_reads",
                       paste0("batch", argv[2], "/"))
}


# source the script with core functions
source(paste0(script_dir, "code_screen_core_fun.R"))

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

fix_nan <- function (x) {
  x[is.nan(x)] <- 0
  return(x)
}

compute_plate_stats <- function(annotated) {
  #trim at 6th time point
  annotated = lapply(annotated, function(x) subset(x, x$Time_points <= 6))
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


# now we have the complete annotated data set from batch 1 and 2
# we can first load them separately and merge them into one object
# to compare fitness of perfect replicates

#========LOAD ANNOTATED DATA=============================

# load annotated data
batch1 <- list.files(paste0(home_dir, "batch1"))
batch2 <- list.files(paste0(home_dir, "batch2"))


# now we have the complete annotated data set from batch 1 and 2
# we can first load them separately and merge them into one object
# to compare fitness of perfect replicates
#========LOAD ANNOTATED DATA=============================
# load annotated data
batch1 <- list.files(paste0(home_dir, "batch1"))
batch2 <- list.files(paste0(home_dir, "batch2"))

# extract the donor names
don1 = gsub("([A-Z]+)(_.+)", "\\1", batch1)
don2 = gsub("([A-Z]+)(_.+)", "\\1", batch2)

# check plate names for common donors
common = intersect(don1, don2)


# subset the files in batches 1 and 2 only to those that have common donors
batch1 = batch1[don1 %in% common[common != "TSB"]]
batch2 = batch2[don2 %in% common[common != "TSB"]]


# process batch 1 data
batch1_annot = lapply(paste0(home_dir, "batch1/", batch1), read.csv, header=T,
                      check.names = FALSE)
names(batch1_annot) <- file_path_sans_ext(batch1)
batch1_annot = compute_plate_stats(annotated = batch1_annot)

# process batch 2 data
batch2_annot = lapply(paste0(home_dir, "batch2/", batch2), read.csv, header=T,
                      check.names = FALSE)
names(batch2_annot) <- file_path_sans_ext(batch2)
batch2_annot = compute_plate_stats(annotated = batch2_annot)


select.cols <- function(df) {
  dplyr::select(df, plate, Time_h, OD, Drug,
                Concentration, Replicate, donor,
                donor_conc, position,
                batch, Time_points,
                run, OD_no1st,
                AUCnorm_bug_rob,
                AUC_well_ODno1st,
                mean_AUC_bug)
}


batch1_df = plyr::ldply(batch1_annot, data.frame, .id = "plate")
batch2_df = plyr::ldply(batch2_annot, data.frame, .id = "plate")

batch2_df = dplyr::filter(batch2_df, !run %in% c(15, 18))
batch1_df = dplyr::filter(batch1_df, run != 2)

# remove BUG and ONE
batch1_df = dplyr::filter(batch1_df, !Drug %in% c("BUG", "ONE", "NOT"))
batch2_df = dplyr::filter(batch2_df, !Drug %in% c("BUG", "ONE", "NOT"))

b1_data = select.cols(batch1_df) %>% group_by(donor, donor_conc, Drug, Concentration, Replicate) %>%
  dplyr::mutate(AUCnorm_bug_rob = mean(AUCnorm_bug_rob),
                AUC_well_ODno1st = mean(AUC_well_ODno1st)) %>%
  dplyr::select(-c(Time_h, Time_points, OD, OD_no1st)) %>% 
  distinct(.keep_all = T)

b2_data = select.cols(batch2_df) %>% group_by(donor, donor_conc, Drug, Concentration, Replicate) %>%
  dplyr::mutate(AUCnorm_bug_rob = mean(AUCnorm_bug_rob),
                AUC_well_ODno1st = mean(AUC_well_ODno1st)) %>%
  dplyr::select(-c(Time_h, Time_points, OD, OD_no1st)) %>% 
  distinct(.keep_all = T)

all_data = inner_join(b1_data, b2_data,
           by = c("Drug", "donor", "donor_conc", "Concentration"))

# first check if the measurements correlate in general
# for all Donor-Recipient_Donor_conc-Recip_conc combinations

rep_plates = ggplot(all_data,
       aes(x = AUCnorm_bug_rob.x,
           y = AUCnorm_bug_rob.y)) +
  geom_point() + 
  theme_bw() + 
  geom_abline(slope = 1, colour = "red") + 
  labs(title = "Replicate plate comparison",
       x = "AUCnorm_bug_rob (batch 1)",
       y = "AUCnorm_bug_rob (batch 2)") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ run.x + run.y,
             ncol = 2, scales="free") 

ggsave(rep_plates, 
       filename = paste0(figdir, "replicatePlatesAUCbugnorm.pdf"))



all_data = dplyr::mutate(all_data,
                        log.AUC.ratio = log2(AUCnorm_bug_rob.x / AUCnorm_bug_rob.y))


rep.ratio = ggplot(all_data, aes(x = log.AUC.ratio)) + 
  geom_histogram(bins = 100) +
  labs(x = "Log2 (AUC_batch1 / AUC_batch2)",
        y = "Count") + 
  theme_bw()

ggsave(rep.ratio,
       filename = paste0(figdir, "repratio-acrossbatch.pdf"))


# data with large discrepancy
bad = filter(all_data, log.AUC.ratio > 2 | log.AUC.ratio < - 2)

setdiff(common, unique(bad$donor))

bad_don = ggplot(bad, aes(x = donor,
                fill = donor)) + 
  geom_bar() + theme_classic() + 
  labs(x = "Donor",
       y = "Count in data with poor reproducibility") + 
  theme(legend.position = "none")
ggsave(bad_don,
       filename = paste0(figdir, "bad-donors.pdf"))

# what about recipients

bad_rec = ggplot(bad, aes(x = Drug,
                fill = Drug)) + 
  geom_bar() + theme_classic() + 
  labs(x = "Recipient",
       y = "Count in data with poor reproducibility") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90))
ggsave(bad_rec,
       filename = paste0(figdir, "bad-recipients.pdf"))




bad_run1 = ggplot(bad, aes(x = factor(run.x),
                          fill = factor(run.x))) +
  geom_bar() + theme_classic() + 
  labs(x = "Experimental run in batch 1 replicates") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90))


bad_run2 = ggplot(bad, aes(x = factor(run.y),
                           fill = factor(run.y))) +
  geom_bar() + theme_classic() + 
  labs(x = "Experimental run in batch 2 replicates") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90))



bad_pos1 = ggplot(bad, aes(x = factor(position.x),
                           fill = factor(position.x))) +
  geom_bar() + theme_classic() + 
  labs(x = "Position in batch 1 replicates") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90)) + 
  facet_wrap(~ run.x)


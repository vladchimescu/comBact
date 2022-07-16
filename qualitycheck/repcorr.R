# replicate correlation

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

compute_plate_stats <- function(annotated, trim = F) {
  # DO basic stats: mean median etc and normalisation
  
  POC_bugwells <- function (df) {
    plate_mean_bug <- dplyr::group_by(df, Time_points) %>%
      dplyr::filter(Drug == "BUG") %>%
      dplyr::summarize(mean_plate_bug = mean(OD), med_plate_bug = median(OD), sd_plate_bug = sd(OD), se_bug = sd(OD)/sqrt(length(OD)),
                       min_plate_bug= min(OD), max_plate_bug = max (OD))
    
    df = inner_join(df, plate_mean_bug, by = "Time_points") %>% mutate(POC_bug = ((OD/mean_plate_bug)*100))
    return(df)
  }
  
  annotated <- lapply (annotated, POC_bugwells)
  annotated <- lapply (annotated, fix_nan)
  
  if(trim) {
    # trim the time point 5 onwards
    annotated <- lapply(annotated, function(x) {
      subset(x, Time_points %in% 1:5)
    })
  }
  
  #basic function to calculate AUC for one well
  AUC_calc_ODno1st <- function (x) {
    AUCs <- dplyr::group_by(x,wellID) %>%
      dplyr::summarise(AUC_well_ODno1st = DescTools::AUC(Time_h, OD,
                                                         method = "spline"))
    x= inner_join(x, AUCs, by = "wellID")
    return(x)
  }
  
  annotated <- lapply (annotated, AUC_calc_ODno1st)
  #proper ranking function (120 is 15*8 time points) --> wells ranked according to AUCs
  
  #proper ranking function (120 is 15*8 time points)
  POC_top4 <- function (df) {
    ranked_AUC <- df %>% mutate(AUC_rank = BiocGenerics::rank(AUC_well_ODno1st, ties.method = "min")) %>%
      arrange(desc(AUC_rank))
    top15 <- ranked_AUC[1:120, ]
    
    plate_mean_top15 <- dplyr::group_by(top15, Time_points) %>%
      dplyr::summarize(mean_plate_top = mean(OD), med_plate_top = median(OD), sd_plate_top = sd(OD), se_top = sd(OD)/sqrt(length(OD)),
                       min_plate_top= min(OD), max_plate_top = max (OD))
    
    df = inner_join(df, plate_mean_top15, by = "Time_points") %>% mutate(POC_top = ((OD/mean_plate_top)*100))
    return(df)
    
  }
  
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
batch1_annot = compute_plate_stats(annotated = batch1_annot, trim = F)

# process batch 2 data
batch2_annot = lapply(paste0(home_dir, "batch2/", batch2), read.csv, header=T,
                      check.names = FALSE)
names(batch2_annot) <- file_path_sans_ext(batch2)
batch2_annot = compute_plate_stats(annotated = batch2_annot)

rep_corr = function (x) {
  
  x = x[(x$Drug != "ONE" | x$Drug != "NOT"), ]
  rep1 = subset(x, Replicate == 1)
  rep2 = subset(x, Replicate == 2)
  rep1 = dplyr::rename (rep1, OD_r1 = OD,
                        OD_no1st_r1 = OD_no1st,
                        AUCnorm_r1 = AUCnorm_bug_rob,
                        AUC_well_ODno1st_r1 = AUC_well_ODno1st)
  rep2 = dplyr::rename (rep2, OD_r2 = OD,
                        OD_no1st_r2 = OD_no1st,
                        AUCnorm_r2 = AUCnorm_bug_rob,
                        AUC_well_ODno1st_r2 = AUC_well_ODno1st)
  rep1$Replicate = NULL
  rep2$Replicate = NULL
  rep1$wellID = NULL
  rep2$wellID = NULL
  allrep = merge(rep1, rep2, by = c("Drug_Flag_Conc",
                                    "Time_points",
                                    "position", "flag_fixed",
                                    "batch", "run", "Drug",
                                    "plate_name", "donor",
                                    "donor_conc"))
  return(allrep)
}

allrep = lapply(batch1_annot, rep_corr)
allrep = lapply(allrep, function (x) transform(x, slope = x$OD_r1/x$OD_r2))
allrep = lapply(allrep, function (x) transform(x, ratio_AUC = x$AUCnorm_r1/x$AUCnorm_r2))
cor_b1 = plyr::ldply(allrep, data.frame, .id = "plate")


label_b1 = cor_b1 %>% group_by(donor) %>%
  dplyr::summarise(pos = quantile(slope, 0.6))

OD_rep_b1 = ggplot(cor_b1, aes(x = donor,
                   y = slope,
                   colour = donor)) +
  geom_boxplot() + 
  geom_text(data = label_b1,
            aes(label = donor,
                y = pos)) + 
  theme_classic() + 
  scale_x_discrete(position = "top") + 
  labs(x = '',
       y = "Replicate OD ratio",
       title = "Replicate concordance by donor in batch 1") + 
  theme(legend.position = 'none',
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

ggsave(OD_rep_b1, 
       filename = paste0(figdir, "ODrep-batch1.pdf"),
       width = 12, height = 8)




## plot replicate OD ratios by recipient (also around 1)
ggplot(cor_b1, aes(x = Drug,
                   y = slope,
                   colour = Drug)) +
  geom_boxplot() + 
  theme_classic() + 
  scale_x_discrete(position = "top") + 
  labs(x = '',
       y = "Replicate OD ratio",
       title = "Replicate concordance by recipient in batch 1") + 
  theme(legend.position = 'none',
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

# check replicate correlation within and across the two batches
# for each donor drug

allrep = lapply(batch2_annot, rep_corr)
allrep = lapply(allrep, function (x) transform(x, slope = x$OD_r1/x$OD_r2))
allrep = lapply(allrep, function (x) transform(x, ratio_AUC = x$AUCnorm_r1/x$AUCnorm_r2))
cor_b2 = plyr::ldply(allrep, data.frame, .id = "plate")


label_b2 = cor_b2 %>% group_by(donor) %>%
  dplyr::summarise(pos = quantile(slope, 0.6))

OD_rep_b2 = ggplot(cor_b2, aes(x = donor,
                   y = slope,
                   colour = donor)) +
  geom_boxplot() + 
  geom_text(data = label_b2,
            aes(label = donor,
                y = pos)) + 
  theme_classic() + 
  scale_x_discrete(position = "top") + 
  ylim(c(0,3)) + 
  labs(x = '',
       y = "Replicate OD ratio",
       title = "Replicate concordance by donor in batch 2") + 
  theme(legend.position = 'none',
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

ggsave(OD_rep_b2, 
       filename = paste0(figdir, "ODrep-batch2.pdf"),
       width = 12, height = 8)



# stratify by time point
rep_tp_b1 = ggplot(cor_b1,
       aes(x = factor(Time_points),
           y = slope,
           color = factor(Time_points))) + 
  geom_boxplot() +
  scale_color_viridis(discrete = T) +
  labs(x = "Time points",
       y = "Replicate OD ratio",
       title = "Replicate concordance split by run in batch 1") + 
  facet_wrap(~ run) + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

ggsave(rep_tp_b1, filename = paste0(figdir, "repTp-batch1.pdf"))


# facet by position
rep_pos_b1 = ggplot(cor_b1,
                   aes(x = factor(position),
                       y = slope,
                       color = factor(position))) + 
  geom_boxplot() +
  scale_color_viridis(discrete = T) +
  labs(x = "Position",
       y = "Replicate OD ratio",
       title = "Replicate concordance split by position and time point in batch 1") + 
  facet_wrap(~ Time_points) + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))


# do the same for batch 2
rep_tp_b2 = ggplot(cor_b2,
                   aes(x = factor(Time_points),
                       y = slope,
                       color = factor(Time_points))) + 
  geom_boxplot() +
  scale_color_viridis(discrete = T) +
  labs(x = "Time points",
       y = "Replicate OD ratio",
       title = "Replicate concordance split by run in batch 2") + 
  facet_wrap(~ run) + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

ggsave(rep_tp_b2, filename = paste0(figdir, "repTp-batch2.pdf"))




label_b1 = cor_b1 %>% group_by(donor) %>%
  dplyr::summarise(pos = quantile(ratio_AUC, 0.55))


auc_conc_b1 = ggplot(cor_b1,
       aes(x = donor,
           y = ratio_AUC,
           colour = donor)) +
  ylim(c(0,2)) + 
  geom_boxplot() + 
  geom_text(data = label_b1,
            aes(label = donor,
                y = pos)) + 
  labs(x = '',
       y = "Replicate AUCnorm_bug_rob ratio",
       title = "Replicate concordance (norm AUC ratios) by donor in batch 1") + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(auc_conc_b1,
       filename = paste0(figdir, "No-bsub-AUC_replicates-batch1.pdf"),
       width = 12,
       height = 8)


OD.vs.AUC = filter(cor_b1, ratio_AUC > 2) %>%
  group_by(donor, Drug, donor_conc, Concentration.x) %>%
  dplyr::summarise(slope = mean(slope),
            ratio_AUC = mean(ratio_AUC))

auc_outliers = ggplot(OD.vs.AUC, aes(x = slope,
                      y = ratio_AUC)) + 
  geom_point() + ylim(c(2,10)) + 
  labs(x = "Mean replicate OD ratio",
       y = "Replicate AUC (norm.) ratio",
       title = "AUCnorm replicate ratio outliers") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(auc_outliers, filename = paste0(figdir, "AUC-outliers-batch1.pdf"))
# what are the donors in these discrepant replicates?

ggplot(OD.vs.AUC, aes(x = donor,
                      fill = donor)) +
  geom_bar() + theme_bw() + 
  labs(x = "Donor",
       y = "Count",
       title = "Donors with discrepant AUC ratios") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5))


# what are some of the outlier drugs
ggplot(OD.vs.AUC, aes(x = Drug,
                      fill = Drug)) +
  geom_bar() + theme_bw() + 
  labs(x = "Recipient",
       y = "Count",
       title = "Recipients with discrepant AUC ratios") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5))


amx = cor_b1 %>% 
  dplyr::filter(Drug == "CTX" & donor == "AMX" & 
                  donor_conc == 2) %>%
  dplyr::distinct(donor, Drug, 
                  Drug_Flag_Conc, donor_conc, slope,
                  Time_points, OD_r1, OD_r2,ratio_AUC,
                  OD_no1st_r1, OD_no1st_r2)


ggplot(filter(amx, Drug_Flag_Conc == "CTX_03_0.625"), aes(x = Time_points,
                y = OD_r1)) +
  geom_point() + 
  geom_point(data = filter(amx, Drug_Flag_Conc == "CTX_03_0.625"), aes(x = Time_points,
                                                                y = OD_r2, color ='red')) +
  ylim(c(0,0.5))


DescTools::AUC(filter(amx, Drug_Flag_Conc == "CTX_03_0.625")$Time_points,
               filter(amx, Drug_Flag_Conc == "CTX_03_0.625")$OD_no1st_r2)

## Check growth curves for drugs with bad AUC ratio
## check for example batch 2 data
## select the drugs for which the spread of AUC ratio
## is worse than for the OD ratio (i.e. when background
## subtraction leads to lower replicate concordance)

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

batch2 <- list.files(paste0(home_dir, "batch2"))

# process batch 2 data
batch2_annot = lapply(paste0(home_dir, "batch2/", batch2), read.csv, header=T,
                      check.names = FALSE)
names(batch2_annot) <- file_path_sans_ext(batch2)

tsb2 = batch2_annot[grep("TSB", names(batch2_annot))]
tsb2 = compute_plate_stats(annotated = tsb2)
tsb2_df = plyr::ldply(tsb2, .id = "plate")

# some of the drugs we can look at are
# AUR, CCCP, CLI, CLR, CTC, DOX, IPM, LZD, NVB, RIF
# which have a particularly bad AUC ratio
dr = c("AUR", "CCCP", "CLI", "CLR", "CTC",
       "DOX", "IPM", "LZD", "NVB", "RIF")
tsb2_df = filter(tsb2_df, Drug %in% dr)


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

allrep = lapply(tsb2, rep_corr)
allrep = lapply(allrep, function (x) transform(x, slope = x$OD_r1/x$OD_r2))
allrep = lapply(allrep, function (x) transform(x, ratio_AUC = x$AUCnorm_r1/x$AUCnorm_r2))
cor_b1 = plyr::ldply(allrep, data.frame, .id = "plate")
cor_b1 = dplyr::filter(cor_b1, Drug != "BUG")
cor_b1 = filter(cor_b1, Drug %in% dr)


bad.rep = bind_rows(filter(cor_b1, ratio_AUC < 0.5),
                    filter(cor_b1, ratio_AUC > 1.5) )
bad.rep = distinct(bad.rep, plate, Drug, Concentration.x)


# take one observation
pdf(file = paste0(figdir, "check-bgsub.pdf"),
    onefile = T)
for (i in 1:nrow(bad.rep)) {
  
  df = filter(tsb2_df, Drug == bad.rep[i,]$Drug &
                plate == bad.rep[i,]$plate & Concentration == bad.rep[i,]$Concentration.x)
  
  od.pl = ggplot(df,
                 aes(x = Time_h,  y = OD, colour = factor(Replicate))) + 
    geom_point() + geom_line() + 
    theme_bw() + 
    labs(x = "Time point [hours]",
         y = "OD",
         colour = "Replicate") + 
    theme(aspect.ratio = 1,
          legend.position = "none") + 
    ylim(c(0, max(df$OD)))
  
  
  odno1.pl = ggplot(df,
                    aes(x = Time_h,  y = OD_no1st, colour = factor(Replicate))) + 
    geom_point() + geom_line() + 
    theme_bw() + 
    labs(x = "Time point [hours]",
         y = "OD - OD(t1)",
         colour = "Replicate") + 
    theme(legend.position = "none",
          aspect.ratio = 1) + 
    ylim(c(0, max(df$OD)))
  
  grpl = grid.arrange(od.pl, odno1.pl, ncol = 2, top = paste(df$Drug[1],
                                                             df$Concentration[1], "uM",
                                                             "Plate:", df$plate[1]))
  
  print(grpl)

}
dev.off()



# plot for one of the normal growing wells
df = filter(tsb2_df, Drug == "AUR" & plate == "TSB_1_K_1_10_16_02" & Concentration == 0.0625)

od.pl = ggplot(df,
               aes(x = Time_h,  y = OD, colour = factor(Replicate))) + 
  geom_point() + geom_line() + 
  theme_bw() + 
  labs(x = "Time point [hours]",
       y = "OD",
       colour = "Replicate") + 
  theme(aspect.ratio = 1,
        legend.position = "none") + 
  ylim(c(0, max(df$OD)))


odno1.pl = ggplot(df,
                  aes(x = Time_h,  y = OD_no1st, colour = factor(Replicate))) + 
  geom_point() + geom_line() + 
  theme_bw() + 
  labs(x = "Time point [hours]",
       y = "OD - OD(t1)",
       colour = "Replicate") + 
  theme(legend.position = "none",
        aspect.ratio = 1) + 
  ylim(c(0, max(df$OD)))

grpl = grid.arrange(od.pl, odno1.pl, ncol = 2, top = paste(df$Drug[1],
                                                           df$Concentration[1], "uM",
                                                           "Plate:", df$plate[1]))



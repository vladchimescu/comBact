#installing/loading needed packages
pacman::p_load(dplyr, reshape2, plyr, stringr, gtools, Rmisc, data.table, RColorBrewer, zoo, hexbin, tools, directlabels, grid, ggpubr, gridExtra,
               readr, ggplot2, stringi, ggrepel, viridis, scales, tidyr, sigr, lemon, here, smoothmest)


#here are gonna be the extracted reads (to be annotated)
toannotate_dir <- "~/Documents/embl/screenK/screenK1_annotated_reads/all_annotated_AUC_scaled_25_5_18/"

#list them
setwd(toannotate_dir)
list_reads <- list.files(toannotate_dir)
names(list_reads) <- file_path_sans_ext(list_reads)

#read them
scr <- lapply(list_reads, read.table, header=T, sep=",")
names(scr) <- file_path_sans_ext(list_reads)

#make df with only raw value, median and sd of the plate
scr <- lapply(scr, function(x) {
  x <- x[ , c("Time_h", "Time_points", "OD", "OD_no1st", "Replicate", "AUC_well", "AUC_well_ODno1st", "z_score", 
              "Drug", "Drug_Flag_Conc", "Concentration", "flag_fixed", "donor", "donor_conc", "plate_name", "position", "run", "batch",
              "wellID")]
  return(x) }
)




#calculate robust mean of BUG wells and normalise for that
robustMean <- function(x){
  if(length(x) == 1){return(x)}
  else if (unique(x) == 0) {return (mean(x))}
  else{
    return(smhuber(x)$mu)
  } 
}


plate_mean_bug <- df %>% dplyr::group_by(Time_points) %>%
  dplyr::summarise(robmean_plate_bug = robustMean(OD_no1st))

df = inner_join(df, plate_mean_bug, by = "Time_points") %>% mutate(ODno1st_normmeanbug = OD_no1st/mean_plate_bug,
                                                                   ODno1st_normRobmeanbug = OD_no1st/robmean_plate_bug,
                                                                   ODno1st_normmedbug = OD_no1st/med_plate_bug)


# df = scr$AMX_1_K_1_4_01_01
# df = filter(df, Drug == "BUG")
# #df = group_by(df, Time_points)
# df1 = subset(df, Time_points == 1)
# smhuber(df1$OD_no1st)
# #this doesnt work: doesnt know how to deal with all zeros

#smhuber cannot deal with 0s --> OD at time 1 is 0
#using AUCs
#to make all lighter get rid of OD
# scr <- lapply(scr, function(x) {
#   x <- x[ , c("Time_h", "Time_points", "Replicate", "AUC_well", "AUC_well_ODno1st",
#               "Drug", "Drug_Flag_Conc", "Concentration", "flag_fixed", "donor", "donor_conc", "plate_name", "position", "run", "batch",
#               "wellID")]
#   return(x) }
# )


POC_bugwells <- function (df) {
  plate_mean_bug <- dplyr::filter(df, Drug == "BUG") %>%
    dplyr::group_by(Time_points) %>%
    dplyr::summarise(mean_AUC_bug_shift = mean(AUC_well_ODshifted), med_AUC_bug_shift = median(AUC_well_ODshifted),
                     mean_AUC_bug_rob_shift = robustMean(AUC_well_ODshifted))
  
  df = inner_join(df, plate_mean_bug, by = "Time_points") %>% mutate(AUCnorm_bug_shift = AUC_well_ODshifted/mean_AUC_bug_shift,
                                                                     AUCnorm_bug_rob_shift = AUC_well_ODshifted/mean_AUC_bug_rob_shift,
                                                                     AUCnorm_bug_med_shift = AUC_well_ODshifted/med_AUC_bug_shift)
  return(df)
}

scr = lapply(scr, POC_bugwells)

#ran this for AUC scaled, shifted, and 1st tp sub

#save these files
# home_dir = "~/Documents/EMBL/combinatorials/screen_K/screenK1_annotated_reads/all_annotated_robustmean_6_18/"
# for (i in seq_along(scr)) {
#   filename = paste0(home_dir, names(scr)[i], ".csv")
#   write.csv(scr[[i]], filename, row.names=FALSE, append = TRUE)
# }

robmean <- lapply(scr, function(x) {
  x <- x[ , c("Time_h", "Time_points", "OD", "OD_no1st", "Replicate", "AUC_well", "AUC_well_ODno1st", "z_score", 
              "Drug", "Drug_Flag_Conc", "flag_fixed", "donor", "donor_conc", "plate_name", "position", "run", "batch",
              "wellID", "AUCnorm_bug", "AUCnorm_bug_rob", "AUCnorm_bug_med",
              "AUCnorm_bug_scal", "AUCnorm_bug_rob_scal", "AUCnorm_bug_med_scal",
              "AUCnorm_bug_shift", "AUCnorm_bug_med_shift", "AUCnorm_bug_rob_shift")]
  return(x) }
)

time6 = plyr::ldply(lapply(robmean, function(x) {
  x[x$Time_points == 6,]
}), 'data.frame')

pdf("~/Documents/EMBL/combinatorials/screen_K/screenK1_plots/robustmean_6_18/boxplot_plates_AUC_robmean_scaled_shifted.pdf", width=32, height=10)

ordered.data <- time6[order(time6$run),]
time6$plate_name <- factor(time6$plate_name, levels = unique(ordered.data$plate_name))

#OD, OD bg sub and z score
ggplot (time6, aes(x= plate_name, y = OD, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("raw OD") +
  scale_y_continuous(limits = c(min(time6$OD), max(time6$OD))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = OD_no1st, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("OD bg sub") +
  scale_y_continuous(limits = c(min(time6$OD_no1st), max(time6$OD_no1st))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = z_score, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("z score") +
  scale_y_continuous(limits = c(min(time6$z_score), max(time6$z_score))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

#AUC raw and bg sub
ggplot (time6, aes(x= plate_name, y = AUC_well, colour = factor(run))) +
  geom_boxplot() + 
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC wells raw") +
  scale_y_continuous(limits = c(min(time6$AUC_well), max(time6$AUC_well))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUC_well_ODno1st, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC calculated from OD after subtracting 1st tp") +
  scale_y_continuous(limits = c(min(time6$AUC_well_ODno1st), max(time6$AUC_well_ODno1st))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

#AUC norm
ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_med, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD no1st normalised by median of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_med), max(time6$AUCnorm_bug_med))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD no1st normalised by mean of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug), max(time6$AUCnorm_bug))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_rob, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD no1st normalised by robust mean of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_rob), max(time6$AUCnorm_bug_rob))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

#AUC norm shifted (min sub)
ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_med_shift, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD min sub normalised by median of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_med_shift), max(time6$AUCnorm_bug_med_shift))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_shift, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD min sub normalised by mean of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_shift), max(time6$AUCnorm_bug_shift))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_rob_shift, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD min sub normalised by robust mean of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_rob_shift), max(time6$AUCnorm_bug_rob_shift))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

#AUC norm scaled
ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_med_scal, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC scaled normalised by median of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_med_scal), max(time6$AUCnorm_bug_med_scal))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_scal, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC scaled normalised by mean of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_scal), max(time6$AUCnorm_bug_scal))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_rob_scal, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC scaled normalised by robust mean of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_rob_scal), max(time6$AUCnorm_bug_rob_scal))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))


#--------------------ordered by donor
time6_bydonor = time6[order(time6$donor), ]
time6$plate_name <- factor(time6$plate_name, levels = unique(time6_bydonor$plate_name))

#OD, OD bg sub and z score
ggplot (time6, aes(x= plate_name, y = OD, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("raw OD by donor") +
  scale_y_continuous(limits = c(min(time6$OD), max(time6$OD))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = OD_no1st, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("OD bg sub by donor") +
  scale_y_continuous(limits = c(min(time6$OD_no1st), max(time6$OD_no1st))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = z_score, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("z score by donor") +
  scale_y_continuous(limits = c(min(time6$z_score), max(time6$z_score))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

#AUC raw and bg sub
ggplot (time6, aes(x= plate_name, y = AUC_well, colour = factor(run))) +
  geom_boxplot() + 
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC wells raw by donor") +
  scale_y_continuous(limits = c(min(time6$AUC_well), max(time6$AUC_well))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUC_well_ODno1st, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC calculated from OD after subtracting 1st tp by donor") +
  scale_y_continuous(limits = c(min(time6$AUC_well_ODno1st), max(time6$AUC_well_ODno1st))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

#AUC norm
ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_med, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD no1st normalised by median of BUG wells by donor") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_med), max(time6$AUCnorm_bug_med))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD no1st normalised by mean of BUG wells by donor") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug), max(time6$AUCnorm_bug))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_rob, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD no1st normalised by robust mean of BUG wells by donor") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_rob), max(time6$AUCnorm_bug_rob))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

#AUC norm shifted (min sub)
ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_med_shift, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD min sub normalised by median of BUG wells by donor") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_med_shift), max(time6$AUCnorm_bug_med_shift))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_shift, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD min sub normalised by mean of BUG wells by donor") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_shift), max(time6$AUCnorm_bug_shift))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_rob_shift, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD min sub normalised by robust mean of BUG wells by donor") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_rob_shift), max(time6$AUCnorm_bug_rob_shift))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

#AUC norm scaled
ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_med_scal, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC scaled normalised by median of BUG wells by donor") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_med_scal), max(time6$AUCnorm_bug_med_scal))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_scal, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC scaled normalised by mean of BUG wells by donor") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_scal), max(time6$AUCnorm_bug_scal))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_rob_scal, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC scaled normalised by robust mean of BUG wells by donor") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_rob_scal), max(time6$AUCnorm_bug_rob_scal))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

dev.off()
  
 
  
#--------------------only tsb plates----------------------
time6 = subset(time6, donor == "TSB")

pdf("~/Documents/EMBL/combinatorials/screen_K/screenK1_plots/robustmean_6_18/boxplot_plates_AUC_robmean_scaled_shifted_TSB.pdf", width=32, height=10)

ordered.data <- time6[order(time6$run),]
time6$plate_name <- factor(time6$plate_name, levels = unique(ordered.data$plate_name))

#OD, OD bg sub and z score
ggplot (time6, aes(x= plate_name, y = OD, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("raw OD") +
  scale_y_continuous(limits = c(min(time6$OD), max(time6$OD))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=19, angle=60, hjust=1),
        axis.text.y = element_text(size=19),
        plot.title = element_text(hjust = 0.5,size=28))

ggplot (time6, aes(x= plate_name, y = OD_no1st, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("OD bg sub") +
  scale_y_continuous(limits = c(min(time6$OD_no1st), max(time6$OD_no1st))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=19, angle=60, hjust=1),
        axis.text.y = element_text(size=19),
        plot.title = element_text(hjust = 0.5,size=28))

ggplot (time6, aes(x= plate_name, y = z_score, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("z score") +
  scale_y_continuous(limits = c(min(time6$z_score), max(time6$z_score))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=19, angle=60, hjust=1),
        axis.text.y = element_text(size=19),
        plot.title = element_text(hjust = 0.5,size=28))

#AUC raw and bg sub
ggplot (time6, aes(x= plate_name, y = AUC_well, colour = factor(run))) +
  geom_boxplot() + 
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC wells raw") +
  scale_y_continuous(limits = c(min(time6$AUC_well), max(time6$AUC_well))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=19, angle=60, hjust=1),
        axis.text.y = element_text(size=19),
        plot.title = element_text(hjust = 0.5,size=28))

ggplot (time6, aes(x= plate_name, y = AUC_well_ODno1st, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC calculated from OD after subtracting 1st tp") +
  scale_y_continuous(limits = c(min(time6$AUC_well_ODno1st), max(time6$AUC_well_ODno1st))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=19, angle=60, hjust=1),
        axis.text.y = element_text(size=19),
        plot.title = element_text(hjust = 0.5,size=28))

#AUC norm
ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_med, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD no1st normalised by median of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_med), max(time6$AUCnorm_bug_med))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=19, angle=60, hjust=1),
        axis.text.y = element_text(size=19),
        plot.title = element_text(hjust = 0.5,size=28))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD no1st normalised by mean of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug), max(time6$AUCnorm_bug))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=19, angle=60, hjust=1),
        axis.text.y = element_text(size=19),
        plot.title = element_text(hjust = 0.5,size=28))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_rob, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD no1st normalised by robust mean of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_rob), max(time6$AUCnorm_bug_rob))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=19, angle=60, hjust=1),
        axis.text.y = element_text(size=19),
        plot.title = element_text(hjust = 0.5,size=28))

#AUC norm shifted (min sub)
ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_med_shift, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD min sub normalised by median of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_med_shift), max(time6$AUCnorm_bug_med_shift))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=19, angle=60, hjust=1),
        axis.text.y = element_text(size=19),
        plot.title = element_text(hjust = 0.5,size=28))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_shift, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD min sub normalised by mean of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_shift), max(time6$AUCnorm_bug_shift))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=19, angle=60, hjust=1),
        axis.text.y = element_text(size=19),
        plot.title = element_text(hjust = 0.5,size=28))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_rob_shift, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD min sub normalised by robust mean of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_rob_shift), max(time6$AUCnorm_bug_rob_shift))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=19, angle=60, hjust=1),
        axis.text.y = element_text(size=19),
        plot.title = element_text(hjust = 0.5,size=28))

#AUC norm scaled
ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_med_scal, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC scaled normalised by median of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_med_scal), max(time6$AUCnorm_bug_med_scal))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=19, angle=60, hjust=1),
        axis.text.y = element_text(size=19),
        plot.title = element_text(hjust = 0.5,size=28))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_scal, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC scaled normalised by mean of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_scal), max(time6$AUCnorm_bug_scal))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=19, angle=60, hjust=1),
        axis.text.y = element_text(size=19),
        plot.title = element_text(hjust = 0.5,size=28))

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_rob_scal, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC scaled normalised by robust mean of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_rob_scal), max(time6$AUCnorm_bug_rob_scal))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=19, angle=60, hjust=1),
        axis.text.y = element_text(size=19),
        plot.title = element_text(hjust = 0.5,size=28))

dev.off()
  
  
  
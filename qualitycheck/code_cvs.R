pacman::p_load(dplyr, reshape2, plyr, stringr, gtools, Rmisc, data.table, RColorBrewer, zoo, hexbin, tools, directlabels, grid, ggpubr, gridExtra,
               readr, ggplot2, stringi, ggrepel, viridis, scales, tidyr, sigr, lemon, here, smoothmest,EnvStats)


#here are gonna be the extracted reads (to be annotated)
toannotate_dir <- "/Users/Elisabetta/Documents/EMBL/combinatorials/screen_K/screenK1_annotated_reads/all_annotated_robustmean_6_18/"

#list them
setwd(toannotate_dir)
list_reads <- list.files(toannotate_dir)
names(list_reads) <- file_path_sans_ext(list_reads)

#read them
scr <- lapply(list_reads, read.table, header=T, sep=",")
names(scr) <- file_path_sans_ext(list_reads)

#make df with only raw value, median and sd of the plate
# scr <- lapply(scr, function(x) {
#   x <- x[ , c("Time_h", "Time_points", "OD", "OD_no1st", "Replicate", "AUC_well", "AUC_well_ODno1st", "z_score", 
#               "Drug", "Drug_Flag_Conc", "Concentration", "flag_fixed", "donor", "donor_conc", "plate_name", "position", "run", "batch",
#               "recflag", "recconc", "donorconc", "combid", "recid", "donrec", "fitrec", "med_AUC_plate_notnorm", "med_AUC_plate_norm_bug",
#               "rowID", "colID", "wellID")]
#   return(x) }
# )

scr <- lapply(scr, function(x) {
  x <- x[ , c("Time_h", "Time_points", "Replicate", "AUC_well_ODno1st", "AUCnorm_bug_rob",
              "Drug", "Drug_Flag_Conc", "Concentration", "flag_fixed",
              "donor", "donor_conc", "plate_name", "position", "run", "batch",
              "donorconc", "med_AUC_plate_notnorm", "med_AUC_plate_norm_bug")]
  return(x) }
)

time6 = plyr::ldply(lapply(scr, function(x) {
  x[x$Time_points == 6,]
}), 'data.frame')

##---------------plot cvs per plate and per well across plates---------------
# cvs = time6 %>% dplyr::group_by(plate_name) %>%
#   dplyr::summarise(cvs_AUC_plate_bgsub = cv(AUC_well_ODno1st), cvs_AUCrobmeannorm = cv(AUCnorm_bug_rob)) %>% ungroup()
# time6 = inner_join(time6, cvs, by= "plate_name")

cvs = time6 %>% dplyr::group_by(plate_name) %>%
  dplyr::summarise(cvs_AUC_plate_bgsub = cv(AUC_well_ODno1st), cv_AUC_robmean = cv(AUCnorm_bug_rob)) %>% ungroup()
time6 = inner_join(time6, cvs, by= "plate_name")

# time6$slope = time6$cvs_AUC_plate_raw / time6$cvs_AUCrobmeannorm
# time6 = time6[ , c("cvs_AUC_plate_raw", "cvs_AUCrobmeannorm", "slope", "plate_name", "position",
#                          "run", "batch", "donor", "donor_conc", "donorconc")]

ggplot(time6, aes(x=cvs_AUC_plate_bgsub, fill=as.factor(batch))) +
  geom_density(alpha = 0.3, adjust = 2)
# xlim(c(0.55, 1)) +
# ggtitle("R scores distribution per shaking regimen (runs 1-13: shaking 50 rpm until run 11)") +
# # annotate(x = .65, y = 9,
#           label = "Real range: 0.82-0.96","text", size=4) +
#geom_vline(data=cdat, aes(xintercept=shake_mean,  colour=as.factor(shake)),
# linetype="dashed", size=1)

hist(time6$cvs_AUC_plate_bgsub)


pdf("~/Documents/EMBL/combinatorials/screen_K/screenK1_plots/cvs_AUC_bgsub.pdf", width=35, height=10)

ordered.data <- time6[order(time6$donor),]
time6$plate_name <- factor(time6$plate_name, levels = unique(ordered.data$plate_name))

ggplot(time6, aes(x = plate_name, y = cvs_AUC_plate_bgsub, fill = factor(run))) +
  geom_bar(stat = "identity") +
  # geom_errorbar(aes(ymin = mean_auce - sd_auce, ymax = mean_auce + sd_auce), colour="black", size = 1, width=.3) +
  # labs(x = "blank/wt/mutant and drug concentration", y = "Mean and sd of AUC (no fit)", fill = "Concentration") +
  # labs(title = title_plot) +
  theme(axis.text.x = element_text(size = 10, angle=60, hjust=1),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 29, angle = 90), 
        axis.title.x = element_text(size=29),
        #plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
        plot.margin = unit(c(3, 3, 3, 3), "cm"), 
        legend.title = element_text(colour="black", size=26,face="bold"),
        legend.text = element_text(colour="black", size=22),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines'))

#scale_y_continuous(breaks = scales::pretty_breaks(n = 20))

dev.off()




time6 = time6[ , c("donor", "donor_conc", "plate_name", "position", "run", "batch",
                   "donorconc", "cvs_AUC_plate_bgsub", "cv_AUC_robmean")]
time6 = unique(time6)

pdf("~/Documents/EMBL/combinatorials/screen_K/screenK1_plots/screen_cvs.pdf", width=15, height=15) #onefile=TRUE

#if you want 96 panels per page (12x8: good format for titration plates) use the following
#pages = ceiling(length(unique(annotated_read$Drug_Flag_Conc)) / 54)

ggplot(time6, aes(x=cvs_AUC_plate_bgsub, y=cv_AUC_robmean, colour = as.factor(run))) +
  geom_point(aes(size = 1), alpha = 0.7) +
  geom_text_repel(data=subset(time6, cv_AUC_robmean > 2.5 & cvs_AUC_plate_bgsub > 2.5) , aes(label= donorconc), size = 8) +
  coord_fixed() +
  geom_smooth(method = "lm", se=F, color="black", formula = y ~ x, linetype="dashed", size = 0.5) +
  # annotate(x = min(df$cvs_OD_plate) + 0.05, y = max(df$cvs_AUCrobmeannorm) - 0.05,
  #          label = paste0("R = ", round(unique(df$ODr_corr), digits = 2)),"text", size=8) +
  geom_abline(slope=1, size = 0.3) +
  labs(title = "CVs per plate", x="CV AUC plate raw", y="CV AUC normalised for robust mean", colour = "Run") +
  # scale_x_continuous(breaks = seq(0, 0.45, by = 0.025), limits = c(0, 0.45)) +
  # scale_y_continuous(breaks = seq(0, 0.45, by = 0.025), limits = c(0, 0.45)) +
  theme(axis.text.x = element_text(size = 22, angle=60, hjust=1),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 29, angle = 90), 
        axis.title.x = element_text(size=29),
        #plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
        plot.margin = unit(c(3, 3, 3, 3), "cm"), 
        legend.title = element_text(colour="black", size=26,face="bold"),
        legend.text = element_text(colour="black", size=22),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines'))

dev.off()    

#write down the cvs
# write.table (time6,
#              "/Users/Elisabetta/Documents/EMBL/combinatorials/screen_K/good_dfs/cvs.tsv",
#              sep = "\t")


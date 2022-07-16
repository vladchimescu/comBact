argv = commandArgs(trailingOnly = TRUE)

# argv[1] = data/ directory
# argv[2] = batch num
# argv[3] = strain

#argv = c("/Volumes/gitlab/parallel/pneumo-main/data/screen_P", 10, "P")

#installing/loading needed packages
pacman::p_load(dplyr, reshape2, plyr, stringr, gtools, Rmisc,
               data.table, RColorBrewer, zoo, hexbin, tools,
               directlabels, grid, ggpubr, gridExtra,
               readr, ggplot2, stringi, ggrepel, viridis,
               scales, tidyr, sigr, lemon, here, smoothmest,EnvStats)

params = rjson::fromJSON(file = paste0("params/batch", argv[2], ".json" ))
if ("excludetp" %in% names(params)) excludetp = params$excludetp
# trim AUC curves at this time point
trimtp = params$trimtp

library(ggforce)

strain = paste0("screen", argv[3])

script_dir = "~/Documents/embl/gitlab/drug_comb_screen/dataprep/"
if(!dir.exists(script_dir)) script_dir = "dataprep/"

figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/", strain, "/batch", argv[2], "/QC-plots/")
  dir.create(figdir, recursive = T)
}


home_dir = file.path(argv[1], 
                     paste0(strain, "_annotated_reads"),
                     paste0("batch", argv[2], "/"))

stats_dir = file.path(argv[1], 
                      paste0(strain, "_stats"),
                      paste0("batch", argv[2], "-stats", "/"))
if(!dir.exists(stats_dir)) {
  dir.create(stats_dir)
}

# source the script with core functions
source(paste0(script_dir, "code_screen_core_fun.R"))

#========LOAD ANNOTATED DATA=============================

# load annotated data
list_all_plates <- list.files(home_dir)

#read them
annotated <- lapply(paste0(home_dir, list_all_plates), read.csv, header=T,
                    check.names = FALSE)
names(annotated) <- file_path_sans_ext(list_all_plates)


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
# unlist the 'annotated'

max_hour = ifelse(argv[3] == "P", 5, 7.8)
# this has to be added just to make sure that run 70 gets trimmed at the right time point
annotated = lapply(annotated, function(x) subset(x, x$Time_h < max_hour))

if(exists("excludetp")) {
  annotated = lapply(annotated, function(x) {
    if(!unique(x$run) %in% c(82,83)) subset(x, !x$Time_points %in% excludetp)
    else x
  })
}

annotated_unl <- data.table::rbindlist(annotated)

stats_all <- stats_whole(annotated_unl,
                         output_dir = stats_dir)


if(argv[3] == "P") {
  # estimate background from individual plates
  # and plot the curves used for plate-level background estimation
  pdf(paste0(figdir, "bgsub_estimation_curves.pdf"), width=13, height=10)
  
  for (i in 1:length(annotated)) {
    title_plot <- names(annotated[i])  
    
    df <- as.data.frame (annotated[[i]])
    
    # remove spiking wells
    # assume within the first 4 time points
    # the curve is monotonically increasing
    non.spiking = dplyr::group_by(df, Well) %>%
      dplyr::filter(OD[Time_points == 2] > OD[Time_points == 1]) %>%
      dplyr::filter(OD[Time_points == 3] > OD[Time_points == 1]) %>%
      dplyr::filter(OD[Time_points == 4] > OD[Time_points == 1]) 
    
    bg_est = ungroup(non.spiking) %>%
      group_by(Time_points) %>%
      mutate(OD = median(OD))
    if(length(unique(as.character(bg_est$Well))) > 5) {
      p = ggplot(non.spiking, aes(x=Time_h, y=OD, group=Well)) + 
        #shape = Well) +
        geom_line(color = "#b3b3b3") +
        # scale_y_continuous(limits = c(min(non.spiking$OD),
        #                               max(non.spiking$OD))) +
        #guides(colour=guide_legend(title="Recipient_well")) + 
        ggtitle(title_plot) +
        geom_point(data=bg_est, color='red') +
        geom_line(data=bg_est, color='red') + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size =22))
      #scale_color_manual(values = rep(getPalette(18), each = 2), name= "Recip_well" )
      print(p)
    }
    
  }      
  
  dev.off()
  
  pdf(paste0(figdir, "spiking_curves.pdf"), width=13, height=10)
  for (i in 1:length(annotated)) {
    title_plot <- names(annotated[i])  
    
    df <- as.data.frame (annotated[[i]])

    # remove spiking wells
    # assume within the first 4 time points
    # the curve is monotonically increasing
    non.spiking = dplyr::group_by(df, Well) %>%
      dplyr::filter(OD[Time_points == 2] > OD[Time_points == 1]) %>%
      dplyr::filter(OD[Time_points == 3] > OD[Time_points == 1]) %>%
      dplyr::filter(OD[Time_points == 4] > OD[Time_points == 1]) 
    
    # median of non-spiking wells as background estimate for spiking wells
    bg_est = ungroup(non.spiking) %>%
      group_by(Time_points) %>%
      mutate(OD = median(OD))
    
    spiking = dplyr::anti_join(df, non.spiking, by="Well")
    if(nrow(spiking)) {
      OD_range = stats::IQR(dplyr::filter(spiking, Time_points ==1)$OD, na.rm = T)
      if(OD_range > 0.05 & !is.na(OD_range)) {
        p = ggplot(spiking, aes(x=Time_h, y=OD, group=Well)) + 
          #shape = Well) +
          geom_line(color = "#b3b3b3") +
          # scale_y_continuous(limits = c(min(non.spiking$OD),
          #                               max(non.spiking$OD))) +
          #guides(colour=guide_legend(title="Recipient_well")) + 
          ggtitle(title_plot) +
          geom_point(data=bg_est, color='red') +
          geom_line(data=bg_est, color='red') + 
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size =22))
        #scale_color_manual(values = rep(getPalette(18), each = 2), name= "Recip_well" )
        print(p)
      }
    }
    
  }      
  
  dev.off()
  
  
  non.spiking = dplyr::group_by(annotated_unl, plate_name, Well)%>%
    dplyr::filter(OD[Time_points == 2] > OD[Time_points == 1]) %>%
    dplyr::filter(OD[Time_points == 3] > OD[Time_points == 1]) %>%
    dplyr::filter(OD[Time_points == 4] > OD[Time_points == 1])
  
  spiking = dplyr::anti_join(annotated_unl, non.spiking, by=c("plate_name", "Well"))
  
  spiking = dplyr::filter(spiking, Time_points == 1) %>%
    dplyr::group_by(plate_name) %>%
    dplyr::mutate(OD_range = stats::IQR(OD, na.rm = T)) %>%
    dplyr::filter(OD_range > 0.05 & !is.na(OD_range)) %>%
    ungroup()
  
  non.spiking = dplyr::anti_join(annotated_unl, spiking,
                                 by = c("plate_name", "Well"))
  
  spike.summary = bind_rows(dplyr::mutate(ungroup(non.spiking), bubble="non spiking"),
            dplyr::mutate(spiking, bubble="spiking")) %>%
    dplyr::distinct(plate_name, Well, position, bubble, run)
  
  p = ggplot(spike.summary, aes(x = factor(position), fill=bubble))+
    geom_bar()+ theme_bw()+
    facet_wrap(~ run) + 
    xlab("Position") + ylab("Count")
  
  ggsave(p, filename = paste0(figdir, "spike-vs-position.pdf"),
         width = 10)
  
  spike.runs = inner_join(group_by(spike.summary, run, bubble) %>%
               dplyr::summarise(n=n()),
             group_by(spike.summary, run) %>%
               dplyr::summarise(ntot=n()),
             by = "run") %>%
    mutate(proportion = n/ntot)
  p = ggplot(spike.runs,
         aes(x = factor(run), y = proportion, fill = bubble))+
    geom_bar(stat = 'identity') +
    labs(x = "Run", y = "Proportion of spiking wells", fill = "") +
    theme_bw()
  ggsave(p, filename = paste0(figdir, "spike-proportion-by-run.pdf"),
         width = 10)
  
  spiking = dplyr::distinct(spiking, plate_name, Well, .keep_all=T)
  spike.wells = dplyr::group_by(spiking, Well, Well_letter, Well_number, position) %>%
    dplyr::summarise(n=n())
  
  p = ggplot(spike.wells, aes (Well_number, Well_letter)) +
    geom_tile(aes(fill = n), colour = "white", size = 0.25) +
    #labs(title = ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ position) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_fill_gradient(low="white", high = "blue") +
    labs(x="", y="") +
    #geom_hex () +
    #maintains aspect ratio
    coord_fixed() +
    #set a base size for all fonts
    theme_grey(base_size=8)+
    #theme options
    theme(
      #bold font for both axis text
      axis.text=element_text(face="bold"),
      axis.text.x = element_text(size=7, angle=60, hjust=1), 
      #set thickness of axis ticks
      axis.ticks=element_line(size=0.4),
      #remove plot background
      plot.background=element_blank(),
      strip.text = element_text(size=15),
      #remove plot border
      #panel.border=element_blank())
      panel.spacing = unit(2, "lines"))
  
  ggsave(p, filename = paste0(figdir, "spike-vs-plate.pdf"),
         width = 10, height = 8)
  

  
  spike.wells = dplyr::group_by(spiking, Well, Well_letter, Well_number, run) %>%
    dplyr::summarise(n=n())
  
  p = ggplot(spike.wells, aes (Well_number, Well_letter)) +
    geom_tile(aes(fill = n), colour = "white", size = 0.25) +
    #labs(title = ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ run) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_fill_gradient(low="white", high = "blue") +
    labs(x="", y="") +
    #geom_hex () +
    #maintains aspect ratio
    coord_fixed() +
    #set a base size for all fonts
    theme_grey(base_size=8)+
    #theme options
    theme(
      #bold font for both axis text
      axis.text=element_text(face="bold"),
      axis.text.x = element_text(size=7, angle=60, hjust=1), 
      #set thickness of axis ticks
      axis.ticks=element_line(size=0.4),
      #remove plot background
      plot.background=element_blank(),
      strip.text = element_text(size=15),
      #remove plot border
      #panel.border=element_blank())
      panel.spacing = unit(2, "lines"))
  
  ggsave(p, filename = paste0(figdir, "spike-vs-run.pdf"),
         width = 10, height = 8)
  
  
  spike.wells = dplyr::group_by(spiking, Well, Well_letter, Well_number, plate_name) %>%
    dplyr::summarise(spike=n()) %>%
    dplyr::mutate(spike = ifelse(spike > 0, T, F))
  
  pdf(paste0(figdir,"heatmap_plates_spikes.pdf"), width=26, height=18)
  pages = ceiling(length(unique(spike.wells$plate_name))/12)
  n_last = length(unique(spike.wells$plate_name)) - floor(length(unique(spike.wells$plate_name)) / 12) * 12
  for (i in seq_len(pages)) {
    nc = 4; nr = 3
    page <-  ggplot(spike.wells, aes (Well_number, Well_letter)) +
      geom_tile(aes(fill = spike), colour = "white", size = 0.25) +
      #labs(title = ) +
      theme(plot.title = element_text(hjust = 0.5)) +
      facet_wrap_paginate(~ plate_name, ncol= ifelse(i == pages & n_last < nc, n_last, nc),
                          nrow = ifelse(i == pages & n_last < nc, 1, nr), page=i) +
      scale_y_discrete(expand=c(0,0)) +
      scale_x_discrete(expand=c(0,0)) +
      labs(x="", y="") +
      coord_fixed() +
      #set a base size for all fonts
      theme_grey(base_size=8)+
      #theme options
      theme(
        #bold font for both axis text
        axis.text=element_text(face="bold"),
        axis.text.x = element_text(size=7, angle=60, hjust=1), 
        #set thickness of axis ticks
        axis.ticks=element_line(size=0.4),
        #remove plot background
        plot.background=element_blank(),
        strip.text = element_text(size=15),
        panel.spacing = unit(2, "lines"))
    
    print(page)   
  }
  
  dev.off()
  
  
  pdf(paste0(figdir, "spike_individual_plates.pdf"), width=60, height=20)
  p = ggplot(group_by(spike.wells, plate_name) %>%
    dplyr::summarise(n=n()), aes(x = plate_name, y = n, fill = plate_name) )+
    geom_bar(stat = 'identity') +
    theme_bw() + 
    labs(x = "", y = "Number of spiking wells") + 
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0))+
    theme(axis.text.x = element_text(size = 20, angle=90, hjust=1),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 29, angle = 90), 
          axis.title.x = element_text(size=29),
          plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
          plot.margin = unit(c(3, 3, 3, 3), "cm"),
          legend.position = "none") 
  print(p)
  dev.off()
  
  
}


#----------------------extract TSB plates to assess plain behaviour of recipient wells----------------------------
TSB = annotated[grep("TSB", names(annotated))]

drvec = unique(unlist(lapply(TSB,function(df) unique(df$Drug_Flag_Conc))))
excl = c(grep("BUG", drvec, value = T),
         grep("NOT", drvec, value = T),
         grep("ONE", drvec, value = T) )
#order them so that it plots drugs nicely in alphabetical order
drvec.red =  sort(setdiff(drvec, excl) )

# max(unlist(lapply(TSB,function(df) max(df$OD))))
# min(unlist(lapply(TSB,function(df) min(df$OD))))

pdf(paste0(figdir, "TSB_plates_drugsasrec.pdf"),
    onefile = T, width=13, height=10)
for (drug in drvec.red) {
  recip = lapply(TSB, function(x) {
    x[x$Drug_Flag_Conc == drug,]
  })
  #drug="PEN_02_0.039"
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  p = ggplot(plyr::ldply(recip, 'data.frame'), 
             aes(x=Time_h, y=OD_no1st, colour = as.factor(position), shape = as.factor(run))) +
    geom_point() +
    geom_line(aes(group=plate_wellID, colour = as.factor(position))) +
    ggtitle(drug) + theme_bw() +
    scale_y_continuous(limits = c(0.0, min(c(0.5, max(plyr::ldply(TSB)$OD_no1st))))) +
    theme(plot.title = element_text(hjust = 0.5, size =22)) +
    scale_color_manual(values = getPalette(18), name= "Position" )
  #scale_shape_manual(values = c(19,17), name="Recip_well")
  print(p)
  
}

dev.off()


pdf(paste0(figdir, "TSB_plates_drugsasrec_rawOD.pdf"),
    onefile = T, width=13, height=10)
for (drug in drvec.red) {
  recip = lapply(TSB, function(x) {
    x[x$Drug_Flag_Conc == drug,]
  })
  #drug="PEN_02_0.039"
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  p = ggplot(plyr::ldply(recip, 'data.frame'), 
             aes(x=Time_h, y=OD, colour = as.factor(position), shape = as.factor(run))) +
    geom_point() +
    geom_line(aes(group=plate_wellID, colour = as.factor(position))) +
    ggtitle(drug) + theme_bw() +
    scale_y_continuous(limits = c(0.0, min(c(0.6, max(plyr::ldply(TSB)$OD))))) +
    theme(plot.title = element_text(hjust = 0.5, size =22)) +
    scale_color_manual(values = getPalette(18), name= "Position" )
  #scale_shape_manual(values = c(19,17), name="Recip_well")
  print(p)
  
}

dev.off()


#plot medians of all plates facetted by run
pdf(paste0(figdir, "median_plates_bg_platename_wrap.pdf"),
    onefile = T, width=18, height=18)

pages = ceiling(length(unique(stats_all$run))/2)

for (i in seq_len(pages)) {
  page <- ggplot(stats_all, aes(x=Time_points, y=med_plate, group = plate_name, colour = plate_name)) +
    geom_point() + geom_line() +
    #facet_wrap_paginate(Drug ~ Flag_Conc, ncol= 9, nrow = 6, page = i) +
    facet_wrap_paginate( ~ run , ncol = 1, nrow = 2, page = i) +
    scale_y_continuous(limits = c(min(stats_all$med_plate), max(stats_all$med_plate))) +
    scale_colour_discrete(guide = 'none') +
    geom_dl(aes(label = plate_name), method = list(dl.combine("last.points"), 'last.bumpup', cex = 0.8, rot=c(30))) +
    ggtitle("Plate medians bg sub by run") +
    theme_bw() +
    theme(plot.margin = unit(c(1, 1.5, 1, 1.5), "cm")) +
    theme(plot.title = element_text(hjust = 0.5, size =22))

  print(page)
}

dev.off()


#plot medians all at once
pdf(paste0(figdir, "median_plates_altogether.pdf"),
    onefile = T, width=20, height=12)

p <- ggplot(stats_all, aes(x=Time_points, y=med_plate, group = plate_name, colour = run)) +
  geom_point() + geom_line() +
  #facet_wrap_paginate(Drug ~ Flag_Conc, ncol= 9, nrow = 6, page = i) +
  scale_y_continuous(limits = c(min(stats_all$med_plate), max(stats_all$med_plate))) +
  #scale_colour_discrete(guide = 'none') +
  geom_dl(aes(label = plate_name), method = list(dl.combine("last.points"), 'last.bumpup', cex = 0.8, rot=c(30))) +
  ggtitle("Plate medians bg sub altogether") +
  theme_bw() +
  #theme(plot.margin = unit(c(1, 1.5, 1, 1.5), "cm")) +
  theme(plot.title = element_text(hjust = 0.5, size =22))

print(p)

dev.off()

#------------------------------plot control wells------------------------------------
# plot ONE wells (drug as single donor)

#to have them separated for each plate
pdf(paste0(figdir, "donor_controls_wells.pdf"), width=13, height=10)
for (i in 1:length(annotated)) {
  title_plot <- names(annotated[i])  
  
  df <- as.data.frame (annotated[[i]])
  one <- df[df$Drug == "ONE",]
  #df = mutate(df, idwell = paste(.id, Well))
  
  p = ggplot(one, aes(x=Time_h, y=OD_no1st, colour = Well)) + 
    #shape = Well) +
    geom_point() + geom_line() +
    scale_y_continuous(limits = c(min(one$OD_no1st),
                                  max(one$OD_no1st))) +
    #guides(colour=guide_legend(title="Recipient_well")) + 
    ggtitle(title_plot) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size =22))
  #scale_color_manual(values = rep(getPalette(18), each = 2), name= "Recip_well" )
  print(p)
}      
dev.off()


pdf(paste0(figdir, "donor_controls_wells_rawOD.pdf"), width=13, height=10)
for (i in 1:length(annotated)) {
  title_plot <- names(annotated[i])  
  
  df <- as.data.frame (annotated[[i]])
  one <- df[df$Drug == "ONE",]
  #df = mutate(df, idwell = paste(.id, Well))
  
  p = ggplot(one, aes(x=Time_h, y=OD, colour = Well)) + 
    #shape = Well) +
    geom_point() + geom_line() +
    scale_y_continuous(limits = c(min(one$OD),
                                  max(one$OD))) +
    #guides(colour=guide_legend(title="Recipient_well")) + 
    ggtitle(title_plot) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size =22))
  #scale_color_manual(values = rep(getPalette(18), each = 2), name= "Recip_well" )
  print(p)
}      
dev.off()


#to have them separated for each plate
pdf(paste0(figdir, "donor_controls_wells-by-time-points.pdf"), width=13, height=10)

for (i in 1:length(annotated)) {
  title_plot <- names(annotated[i])  
  
  df <- as.data.frame (annotated[[i]])
  one <- df[df$Drug == "ONE",]
  #df = mutate(df, idwell = paste(.id, Well))
  
  p = ggplot(one, aes(x=Time_points, y=OD_no1st, colour = Well)) + 
    #shape = Well) +
    geom_point() + geom_line() +
    scale_y_continuous(limits = c(min(one$OD_no1st),
                                  max(one$OD_no1st))) +
    #guides(colour=guide_legend(title="Recipient_well")) + 
    ggtitle(title_plot) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size =22))
  #scale_color_manual(values = rep(getPalette(18), each = 2), name= "Recip_well" )
  print(p)
}      

dev.off()


pdf(paste0(figdir, "donor_controls_wells-by-time-points_rawOD.pdf"), width=13, height=10)

for (i in 1:length(annotated)) {
  title_plot <- names(annotated[i])  
  
  df <- as.data.frame (annotated[[i]])
  one <- df[df$Drug == "ONE",]
  #df = mutate(df, idwell = paste(.id, Well))
  
  p = ggplot(one, aes(x=Time_points, y=OD, colour = Well)) + 
    #shape = Well) +
    geom_point() + geom_line() +
    scale_y_continuous(limits = c(min(one$OD),
                                  max(one$OD))) +
    #guides(colour=guide_legend(title="Recipient_well")) + 
    ggtitle(title_plot) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size =22))
  #scale_color_manual(values = rep(getPalette(18), each = 2), name= "Recip_well" )
  print(p)
}      

dev.off()


#BUG controls all in one plot (from all different plates)
#--> to do this put in annotated only one strain, otherwise too many lines in one plot
#also plots according to donor conc

pdf(paste0(figdir, "bug_wells_once.pdf"), width=13, height=10, onefile=TRUE)
recip = lapply(annotated, function(x) {
  dplyr::filter(x, Drug == "BUG" & OD_no1st < 1.5)
})

getPalette = colorRampPalette(brewer.pal(9, "Set1"))

p = ggplot(plyr::ldply(recip, 'data.frame'), 
           aes(x=Time_h, y=OD_no1st,
               group=interaction(plate_name,Well), colour = factor(run))) +
  geom_point() + geom_line() +
  ggtitle("cells control wells - by run") + theme_bw() +
  scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD_no1st)))),
                                max(unlist(lapply(recip,function(df) max(df$OD_no1st)))))) +
  guides(colour=guide_legend(title="Run")) + 
  theme(plot.title = element_text(hjust = 0.5, size =22))

print(p)

p = ggplot(plyr::ldply(recip, 'data.frame'), 
           aes(x=Time_h, y=OD_no1st, group=interaction(plate_name,Well), colour = factor(donor_conc))) +
  geom_point() + geom_line() +
  ggtitle("cells control wells - by donor concentration") + theme_bw() +
  scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD_no1st)))),
                                max(unlist(lapply(recip,function(df) max(df$OD_no1st)))))) +
  guides(colour=guide_legend(title="Donor conc")) + 
  theme(plot.title = element_text(hjust = 0.5, size =22))

print(p)

#---------- highlighting different donor concentrations-------

p = ggplot(plyr::ldply(recip, 'data.frame') %>%
             subset(donor_conc == 1),
           aes(x=Time_h, y=OD_no1st, group=interaction(plate_name, Well), colour = factor(run))) +
  
  geom_point(aes(shape = factor(donor_conc))) + 
  geom_line(aes(linetype= factor(donor_conc))) +
  ggtitle("cells control wells - donor conc = highest - by run") + theme_bw() +
  scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD_no1st)))),
                                max(unlist(lapply(recip,function(df) max(df$OD_no1st)))))) +
  theme(plot.title = element_text(hjust = 0.5, size =22))

print(p)

p = ggplot(plyr::ldply(recip, 'data.frame') %>%
             subset(donor_conc == 2),
           aes(x=Time_h, y=OD_no1st, group=interaction(plate_name, Well), colour = factor(run))) +
  geom_point(aes(shape = factor(donor_conc))) + geom_line(aes(linetype= factor(donor_conc))) +
  ggtitle("cells control wells - donor conc = intermediate - by run") + theme_bw() +
  scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD_no1st)))),
                                max(unlist(lapply(recip,function(df) max(df$OD_no1st)))))) +
  theme(plot.title = element_text(hjust = 0.5, size =22))

print(p)

p = ggplot(plyr::ldply(recip, 'data.frame') %>%
             subset(donor_conc == 3),
           aes(x=Time_h, y=OD_no1st, group=interaction(plate_name, Well), colour = factor(run))) +
  geom_point(aes(shape = factor(donor_conc))) + geom_line(aes(linetype= factor(donor_conc))) +
  ggtitle("cells control wells - donor conc = lowest - by run") + theme_bw() +
  scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD_no1st)))),
                                max(unlist(lapply(recip,function(df) max(df$OD_no1st)))))) +
  theme(plot.title = element_text(hjust = 0.5, size =22))

print(p)

#----------------highlighting different positions---------------

p = ggplot(plyr::ldply(recip, 'data.frame'),
           aes(x=Time_h, y=OD_no1st, group=interaction(plate_name, Well), colour = factor(position))) +
  geom_point(aes(shape = factor(donor_conc))) + geom_line(aes(linetype= factor(donor_conc))) +
  ggtitle("cells control wells - position") + theme_bw() +
  scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD_no1st)))),
                                max(unlist(lapply(recip,function(df) max(df$OD_no1st)))))) +
  theme(plot.title = element_text(hjust = 0.5, size =22))

print(p)

dev.off()


pdf(paste0(figdir, "bug_wells_once-by-time-points.pdf"), width=13, height=10, onefile=TRUE)

getPalette = colorRampPalette(brewer.pal(9, "Set1"))

p = ggplot(plyr::ldply(recip, 'data.frame'), 
           aes(x=Time_points, y=OD_no1st,
               group=interaction(plate_name,Well), colour = factor(run))) +
  geom_point() + geom_line() +
  ggtitle("cells control wells - by run") + theme_bw() +
  scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD_no1st)))),
                                max(unlist(lapply(recip,function(df) max(df$OD_no1st)))))) +
  guides(colour=guide_legend(title="Run")) + 
  theme(plot.title = element_text(hjust = 0.5, size =22))

print(p)

p = ggplot(plyr::ldply(recip, 'data.frame'), 
           aes(x=Time_points, y=OD_no1st, group=interaction(plate_name,Well), colour = factor(donor_conc))) +
  geom_point() + geom_line() +
  ggtitle("cells control wells - by donor concentration") + theme_bw() +
  scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD_no1st)))),
                                max(unlist(lapply(recip,function(df) max(df$OD_no1st)))))) +
  guides(colour=guide_legend(title="Donor conc")) + 
  theme(plot.title = element_text(hjust = 0.5, size =22))

print(p)

#---------- highlighting different donor concentrations-------

p = ggplot(plyr::ldply(recip, 'data.frame') %>%
             subset(donor_conc == 1),
           aes(x=Time_points, y=OD_no1st, group=interaction(plate_name, Well), colour = factor(run))) +
  
  geom_point(aes(shape = factor(donor_conc))) + 
  geom_line(aes(linetype= factor(donor_conc))) +
  ggtitle("cells control wells - donor conc = highest - by run") + theme_bw() +
  scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD_no1st)))),
                                max(unlist(lapply(recip,function(df) max(df$OD_no1st)))))) +
  theme(plot.title = element_text(hjust = 0.5, size =22))

print(p)

p = ggplot(plyr::ldply(recip, 'data.frame') %>%
             subset(donor_conc == 2),
           aes(x=Time_points, y=OD_no1st, group=interaction(plate_name, Well), colour = factor(run))) +
  geom_point(aes(shape = factor(donor_conc))) + geom_line(aes(linetype= factor(donor_conc))) +
  ggtitle("cells control wells - donor conc = intermediate - by run") + theme_bw() +
  scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD_no1st)))),
                                max(unlist(lapply(recip,function(df) max(df$OD_no1st)))))) +
  theme(plot.title = element_text(hjust = 0.5, size =22))

print(p)

p = ggplot(plyr::ldply(recip, 'data.frame') %>%
             subset(donor_conc == 3),
           aes(x=Time_points, y=OD_no1st, group=interaction(plate_name, Well), colour = factor(run))) +
  geom_point(aes(shape = factor(donor_conc))) + geom_line(aes(linetype= factor(donor_conc))) +
  ggtitle("cells control wells - donor conc = lowest - by run") + theme_bw() +
  scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD_no1st)))),
                                max(unlist(lapply(recip,function(df) max(df$OD_no1st)))))) +
  theme(plot.title = element_text(hjust = 0.5, size =22))

print(p)

#----------------highlighting different positions---------------

p = ggplot(plyr::ldply(recip, 'data.frame'),
           aes(x=Time_points, y=OD_no1st, group=interaction(plate_name, Well), colour = factor(position))) +
  geom_point(aes(shape = factor(donor_conc))) + geom_line(aes(linetype= factor(donor_conc))) +
  ggtitle("cells control wells - position") + theme_bw() +
  scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD_no1st)))),
                                max(unlist(lapply(recip,function(df) max(df$OD_no1st)))))) +
  theme(plot.title = element_text(hjust = 0.5, size =22))

print(p)

dev.off()


# plot control wells for each plate individually
pdf(paste0(figdir, "bug_wells_by_plate.pdf"), width=13, height=10, onefile=TRUE)

for (i in 1:length(annotated)) {
  title_plot <- names(annotated[i])  
  
  df <- as.data.frame (annotated[[i]])
  one <- df[df$Drug == "BUG",]
  #df = mutate(df, idwell = paste(.id, Well))
  
  p = ggplot(one, aes(x=Time_h, y=OD_no1st, colour = Well)) + 
    #shape = Well) +
    geom_point() + geom_line() +
    scale_y_continuous(limits = c(min(one$OD_no1st),
                                  max(one$OD_no1st))) +
    #guides(colour=guide_legend(title="Recipient_well")) + 
    ggtitle(title_plot) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size =22))
  #scale_color_manual(values = rep(getPalette(18), each = 2), name= "Recip_well" )
  print(p)
}      

dev.off()



pdf(paste0(figdir, "bug_wells_by_plate_rawOD.pdf"), width=13, height=10, onefile=TRUE)

for (i in 1:length(annotated)) {
  title_plot <- names(annotated[i])  
  
  df <- as.data.frame (annotated[[i]])
  one <- df[df$Drug == "BUG",]
  #df = mutate(df, idwell = paste(.id, Well))
  
  p = ggplot(one, aes(x=Time_h, y=OD, colour = Well)) + 
    #shape = Well) +
    geom_point() + geom_line() +
    scale_y_continuous(limits = c(min(one$OD),
                                  max(one$OD))) +
    #guides(colour=guide_legend(title="Recipient_well")) + 
    ggtitle(title_plot) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size =22))
  #scale_color_manual(values = rep(getPalette(18), each = 2), name= "Recip_well" )
  print(p)
}      

dev.off()

# data frame with control wells
recip_df = plyr::ldply(recip, 'data.frame')
#-------------the runs separated----------
pdf(paste0(figdir, "bug_wells_once_stackingornot.pdf"),
    width=13, height=10, onefile=TRUE)

for (r in unique(recip_df$run)) {
  p = ggplot(recip_df %>%
               subset(run == r),
             aes(x=Time_points, y=OD_no1st, 
                 group=interaction(plate_name,Well),
                 colour = factor(donor_conc))) +
    geom_point() + geom_line() +
    ggtitle(paste("cells control wells - run", r,
                  " - by donor concentration")) +
    theme_bw() +
    scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD_no1st)))),
                                  max(unlist(lapply(recip,function(df) max(df$OD_no1st)))))) +
    geom_dl(aes(label = position),
            method = list(dl.trans(x = x + 0.2), "last.points", 'last.bumpup', cex = 0.9, rot=c(30))) +
    guides(colour=guide_legend(title="Donor conc")) + 
    theme(plot.title = element_text(hjust = 0.5, size =22))
  
  print(p)
  
}

dev.off()


pdf(paste0(figdir, "bug_wells_once_stackingornot_rawOD.pdf"),
    width=13, height=10, onefile=TRUE)

for (r in unique(recip_df$run)) {
  p = ggplot(recip_df %>%
               subset(run == r),
             aes(x=Time_points, y=OD, 
                 group=interaction(plate_name,Well),
                 colour = factor(donor_conc))) +
    geom_point() + geom_line() +
    ggtitle(paste("cells control wells - run", r,
                  " - by donor concentration")) +
    theme_bw() +
    scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD)))),
                                  max(unlist(lapply(recip,function(df) max(df$OD)))))) +
    geom_dl(aes(label = position),
            method = list(dl.trans(x = x + 0.2), "last.points", 'last.bumpup', cex = 0.9, rot=c(30))) +
    guides(colour=guide_legend(title="Donor conc")) + 
    theme(plot.title = element_text(hjust = 0.5, size =22))
  
  print(p)
  
}

dev.off()


pdf(paste0(figdir, "bugOD_by_well_run.pdf"))
p = ggplot(recip_df, aes(x = Well, y = OD, fill = Well)) +
  geom_boxplot() +
  facet_wrap(~ run) +
  theme_bw() +
  labs(x = "")
print(p)
dev.off()

#median vs BUG wells for each plate
## HERE
# pdf(paste0(figdir, "median_vs_bugs.pdf"),
#     onefile = T, width=20, height=12)
# ordered.data <- annotated_unl[order(annotated_unl$run),]
# annotated_unl$plate_name <- factor(annotated_unl$plate_name, levels = unique(ordered.data$plate_name))
# 
# pages = ceiling(length(unique(annotated_unl$plate_name)) / 30)
# 
# for (i in seq_len(pages)) {
#   page <- ggplot(annotated_unl, aes(x = Time_points, y = med_plate_bug, colour=factor(run))) + 
#     geom_point() + geom_line() +
#     geom_line(aes(x=Time_points, y=med_plate_bug), colour = "red") +
#     #geom_errorbar(aes(ymin=mean_plate-sd, ymax=mean_plate+sd), colour="black", width=.3) +
#     scale_y_continuous(limits = c(min(annotated_unl$med_plate_bug), max(annotated_unl$med_plate_bug))) +
#     #guides(colour=guide_legend(title="Recipient_well")) + 
#     theme(plot.title = element_text(hjust = 0.5, size =22)) +
#     ggtitle("Plate medians all plates vs BUG wells medians") +
#     facet_wrap_paginate( ~ plate_name, ncol = 6 , nrow = 5, page = i) 
#   
#   print(page)
#   
# }
# 
# dev.off()


#--------histogram of BUG controls with their median highlighted at time6

pdf(paste0(figdir, "hist_BUG_wells.pdf"), width=60, height=20)
# select Time_point with the smallest difference to the final time point (max_hour)
time6 = lapply(annotated, function(x) {
  dplyr::group_by(x, Well) %>%
    dplyr::mutate(time_diff = abs(Time_h - max_hour) ) %>%
    dplyr::filter(time_diff <= min(time_diff)) %>%
    dplyr::select(-time_diff)
})
time6 <- data.table::rbindlist(time6)
ordered.data <- time6[order(time6$run),]
time6$plate_name <- factor(time6$plate_name, levels = unique(ordered.data$plate_name))

ggplot(filter(time6, Drug == "BUG"), 
       aes(x = plate_name, y = mean_plate_bug,
           fill = factor(donor), color = factor(plate_name))) +
  geom_bar(position=position_dodge(), stat = "identity") +
  geom_errorbar(aes(ymin = min_plate_bug, ymax = max_plate_bug),
                colour="black", size = 1, width=.3, position=position_dodge(.9)) +
  labs(x = "BUG controls per plate", y = paste("Mean, min and max of OD at tp =", params$trimtp), fill = "plate") +
  labs(title = "BUG controls per plate: min max mean") +
  theme_bw() + 
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(size = 20, angle=90, hjust=1),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 29, angle = 90), 
        axis.title.x = element_text(size=29),
        plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
        plot.margin = unit(c(3, 3, 3, 3), "cm"),
        legend.position = "none") 
#scale_y_continuous(breaks = scales::pretty_breaks(n = 20))

dev.off()



#------- plot NOT controls (just medium)
pdf(paste0(figdir, "not_wells_once.pdf"), width=13, height=10, onefile=TRUE)

recip = lapply(annotated, function(x) {
  x[x$Drug == "NOT",]
})

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
recip_df = plyr::ldply(recip, 'data.frame')

p = ggplot(recip_df, 
           aes(x=Time_points, y=OD_no1st,
               group=interaction(plate_name,Well), colour = factor(run))) +
  geom_point() + geom_line() +
  ggtitle("medium control wells - by run") + theme_bw() +
  scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD_no1st)))),
                                max(unlist(lapply(recip,function(df) max(df$OD_no1st)))))) +
  guides(colour=guide_legend(title="Run")) + 
  theme(plot.title = element_text(hjust = 0.5, size =22))

print(p)


for (r in unique(recip_df$run)) {
  p = ggplot(recip_df %>%
               subset(run == r),
             aes(x=Time_points, y=OD_no1st, group=interaction(plate_name,Well), colour = factor(donor_conc))) +
    geom_point() + geom_line() +
    ggtitle(paste("medium control wells - run", r, " - by donor concentration")) + 
    theme_bw() +
    scale_y_continuous(limits = c(min(unlist(lapply(recip,function(df) min(df$OD_no1st)))),
                                  max(unlist(lapply(recip,function(df) max(df$OD_no1st)))))) +
    geom_dl(aes(label = plate_name), method = list(dl.trans(x = x - 2), "last.points", 'last.bumpup', cex = 0.9, rot=c(30))) +
    guides(colour=guide_legend(title="Donor conc")) + 
    theme(plot.title = element_text(hjust = 0.5, size =22))
  
  print(p)
}

dev.off()

#------------------plot heatmap for 6th time point

#chunk to make the average per quadrant across all the plates (to check if quadrant effect) - OD is raw
# med_quadrants = dplyr::group_by (annotated_unl, wellID, quadr, Time_points, run, Well_letter, Well_number) %>%
#   dplyr::summarise (med_plate = median(OD), sd=sd(OD), se = sd(OD)/sqrt(length(OD)),
#                     mean_plate = mean(OD))
# 
# pdf(paste0(figdir, "heatmap_plates_quadrants_alltp.pdf"), width=26, height=18)
# 
# #8 is number of time points
# for (i in 1:8) {
#   
#   page <-  ggplot(med_quadrants, aes (Well_number, Well_letter)) +
#     geom_tile(aes(fill = med_plate), colour = "white", size = 0.25) +
#     #labs(title = ) +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     facet_wrap_paginate(Time_points ~ quadr, ncol= 2, nrow = 2, page=i) +
#     scale_fill_continuous(limits = c(min(med_quadrants$med_plate),max(med_quadrants$med_plate))) +
#     scale_y_discrete(expand=c(0,0)) +
#     scale_x_discrete(expand=c(0,0)) +
#     scale_fill_gradient(low="white", high = "blue") +
#     labs(x="", y="") +
#     #geom_hex () +
#     #maintains aspect ratio
#     coord_fixed() +
#     #set a base size for all fonts
#     theme_grey(base_size=8)+
#     #theme options
#     theme(
#       #bold font for both axis text
#       axis.text=element_text(face="bold"),
#       axis.text.x = element_text(size=7, angle=60, hjust=1), 
#       #set thickness of axis ticks
#       axis.ticks=element_line(size=0.4),
#       #remove plot background
#       plot.background=element_blank(),
#       strip.text = element_text(size=15),
#       #remove plot border
#       #panel.border=element_blank())
#       panel.spacing = unit(2, "lines"))
#   
#   print(page)   
# }
# 
# dev.off()

#----heatmap of whole plate--------

time6 = lapply(annotated, function(x) {
  dplyr::group_by(x, Well) %>%
    dplyr::mutate(time_diff = abs(Time_h - max_hour) ) %>%
    dplyr::filter(time_diff <= min(time_diff)) %>%
    dplyr::select(-time_diff)
})
time6 <- data.table::rbindlist(time6)

#reverse order of letters to get A-P in the heatmap as in the 384 well
time6$Well_letter <- factor(as.character(time6$Well_letter),levels=rev(levels(time6$Well_letter)))

pdf(paste0(figdir, "heatmap_plates_", max_hour, "hr_ODno1st.pdf"), width=26, height=18)

pages = ceiling(length(unique(time6$plate_name))/9)

for (i in seq_len(pages)) {
  page <-  ggplot(time6, aes (Well_number, Well_letter)) +
    geom_tile(aes(fill = OD_no1st), colour = "white", size = 0.25) +
    #labs(title = ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap_paginate(~ plate_name, ncol= 3, nrow = 3, page=i) +
    scale_fill_continuous(limits = c(min(time6$OD_no1st),max(time6$OD_no1st))) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_fill_gradient(low="white", high = "blue") +
    labs(x="", y="") +
    #geom_hex () +
    #maintains aspect ratio
    coord_fixed() +
    #set a base size for all fonts
    theme_grey(base_size=8)+
    #theme options
    theme(
      #bold font for both axis text
      axis.text=element_text(face="bold"),
      axis.text.x = element_text(size=7, angle=60, hjust=1), 
      #set thickness of axis ticks
      axis.ticks=element_line(size=0.4),
      #remove plot background
      plot.background=element_blank(),
      strip.text = element_text(size=15),
      #remove plot border
      #panel.border=element_blank())
      panel.spacing = unit(2, "lines"))
  
  print(page)   
}

dev.off()

#raw OD
pdf(paste0(figdir,"heatmap_plates_",max_hour, "hr_ODraw.pdf"), width=26, height=18)
pages = ceiling(length(unique(time6$plate_name))/9)

for (i in seq_len(pages)) {
  page <-  ggplot(time6, aes (Well_number, Well_letter)) +
    geom_tile(aes(fill = OD), colour = "white", size = 0.25) +
    #labs(title = ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap_paginate(~ plate_name, ncol= 3, nrow = 3, page=i) +
    scale_fill_continuous(limits = c(min(time6$OD),max(time6$OD))) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_fill_gradient(low="white", high = "blue") +
    labs(x="", y="") +
    #geom_hex () +
    #maintains aspect ratio
    coord_fixed() +
    #set a base size for all fonts
    theme_grey(base_size=8)+
    #theme options
    theme(
      #bold font for both axis text
      axis.text=element_text(face="bold"),
      axis.text.x = element_text(size=7, angle=60, hjust=1), 
      #set thickness of axis ticks
      axis.ticks=element_line(size=0.4),
      #remove plot background
      plot.background=element_blank(),
      strip.text = element_text(size=15),
      #remove plot border
      #panel.border=element_blank())
      panel.spacing = unit(2, "lines"))
  
  print(page)   
}

dev.off()


#do also heatmap of AUCs normalised and not
annotated_unl$Well_letter <- factor(as.character(annotated_unl$Well_letter),levels=rev(levels(annotated_unl$Well_letter)))
pdf(paste0(figdir, "heatmap_plates_AUC_well_bug_good.pdf"), width=26, height=18)

pages = ceiling(length(unique(annotated_unl$plate_name))/9)

for (i in seq_len(pages)) {
  page <-  ggplot(annotated_unl, aes (Well_number, Well_letter)) +
    geom_tile(aes(fill = AUC_norm_bug), colour = "white", size = 0.25) +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap_paginate(~ plate_name, ncol= 3, nrow = 3, page=i) +
    scale_fill_continuous(limits = c(min(annotated_unl$AUC_norm_bug),max(annotated_unl$AUC_norm_bug))) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_fill_gradient(low="white", high = "blue") +
    labs(x="", y="") +
    #geom_hex () +
    #maintains aspect ratio
    coord_fixed() +
    #set a base size for all fonts
    theme_grey(base_size=8)+
    #theme options
    theme(
      #bold font for both axis text
      axis.text=element_text(face="bold"),
      axis.text.x = element_text(size=7, angle=60, hjust=1), 
      #set thickness of axis ticks
      axis.ticks=element_line(size=0.4),
      #remove plot background
      plot.background=element_blank(),
      strip.text = element_text(size=15),
      #remove plot border
      #panel.border=element_blank())
      panel.spacing = unit(2, "lines"))
  
  print(page)   
}

dev.off()


annotated_unl$Well_letter <- factor(as.character(annotated_unl$Well_letter),levels=rev(levels(annotated_unl$Well_letter)))
pdf(paste0(figdir, "heatmap_plates_AUC_well_top_good.pdf"), width=26, height=18)

pages = ceiling(length(unique(annotated_unl$plate_name))/9)

for (i in seq_len(pages)) {
  page <-  ggplot(annotated_unl, aes (Well_number, Well_letter)) +
    geom_tile(aes(fill = AUC_norm_top), colour = "white", size = 0.25) +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap_paginate(~ plate_name, ncol= 3, nrow = 3, page=i) +
    scale_fill_continuous(limits = c(min(annotated_unl$AUC_norm_top),max(annotated_unl$AUC_norm_top))) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_fill_gradient(low="white", high = "blue") +
    labs(x="", y="") +
    #geom_hex () +
    #maintains aspect ratio
    coord_fixed() +
    #set a base size for all fonts
    theme_grey(base_size=8)+
    #theme options
    theme(
      #bold font for both axis text
      axis.text=element_text(face="bold"),
      axis.text.x = element_text(size=7, angle=60, hjust=1), 
      #set thickness of axis ticks
      axis.ticks=element_line(size=0.4),
      #remove plot background
      plot.background=element_blank(),
      strip.text = element_text(size=15),
      #remove plot border
      #panel.border=element_blank())
      panel.spacing = unit(2, "lines"))
  
  print(page)   
}

dev.off()


#------------boxplots of plates at 6th tp----------
pdf(paste0(figdir, "boxplot_plates.pdf"), width=30, height=10)

#time6 = subset(annotated_unl, Time_points == trimtp)

ordered.data <- time6[order(time6$run),]
time6$plate_name <- factor(time6$plate_name, levels = unique(ordered.data$plate_name))

ggplot (time6, aes(x= plate_name, y = OD, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("raw OD") +
  scale_y_continuous(limits = c(min(time6$OD), max(time6$OD))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))
# facet_wrap( ~ run, nrow = 5, ncol = 2)

ggplot (time6, aes(x= plate_name, y = OD_no1st, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("OD bg sub") +
  scale_y_continuous(limits = c(min(time6$OD_no1st), max(time6$OD_no1st))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUC_well_ODno1st, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC calculated from OD after trimming 1st tp") +
  scale_y_continuous(limits = c(min(time6$AUC_well_ODno1st), max(time6$AUC_well_ODno1st))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))


ggplot (time6, aes(x= plate_name, y = AUC_well_top, colour = factor(run))) +
  geom_boxplot() + 
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC normalised by median of top 4% wells") +
  #scale_y_continuous(limits = c(min(time6$AUC_well_top), max(time6$AUC_well_top))) +
  scale_y_continuous(limits = c(1.2, 6600)) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))


ggplot (time6, aes(x= plate_name, y = AUC_norm_bug, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC normalised by median of BUG wells") +
  #scale_y_continuous(limits = c(min(time6$AUC_norm_bug), max(time6$AUC_norm_bug))) +
  scale_y_continuous(limits = c(0.0, 2.4)) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))


#ordered by donor
time6_bydonor = time6[order(time6$donor), ]
time6$plate_name <- factor(time6$plate_name, levels = unique(time6_bydonor$plate_name))

ggplot (time6, aes(x= plate_name, y = OD, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("raw OD - comparing donors") +
  scale_y_continuous(limits = c(min(time6$OD), max(time6$OD))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = OD_no1st, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("OD bg sub - comparing donors") +
  scale_y_continuous(limits = c(min(time6$OD_no1st), max(time6$OD_no1st))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUC_well_ODno1st, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC calculated from OD after trimming 1st tp - comparing donors") +
  scale_y_continuous(limits = c(min(time6$AUC_well_ODno1st), max(time6$AUC_well_ODno1st))) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))

ggplot (time6, aes(x= plate_name, y = AUC_norm_bug, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC normalised by median of BUG wells - comparing donors") +
  #scale_y_continuous(limits = c(min(time6$AUC_norm_bug), max(time6$AUC_norm_bug))) +
  scale_y_continuous(limits = c(0.0, 2.4)) +
  guides(geom_point=guide_legend(title="Control wells", labels = c("bugs", "medium", "single donor"))) +
  theme(axis.text.x = element_text(size=6, angle=60, hjust=1))


dev.off()


#warning: this takes lots of time
# pdf(paste0(figdir, "well_QC_AUCnormbug.pdf"), width=30, height=25)
# 
# time6$donorconc <- paste0(time6$donor, "_", time6$donor_conc)
# #order by values of AUC
# ordered.data <- time6[order(time6$AUC_norm_bug,time6$Drug_Flag_Conc_Rep), ]
# time6$plate_name <- factor(time6$plate_name, levels = unique(ordered.data$plate_name))
# 
# for (i in 1:96){
#   page <- ggplot(time6, 
#                  aes(x=plate_name, y=AUC_norm_bug, colour = factor(donor))) +
#     geom_point() +
#     geom_line(aes(x= plate_name, y = AUC_norm_bug, group=1)) +
#     geom_text_repel(data=time6, aes(label=donorconc)) +
#     theme(axis.text.x = element_text(size=7, angle=45, hjust=1)) +
#     ggtitle ("Wells across all plates AUC norm bug") +
#     scale_y_continuous(limits = c(-0.1, 6)) +
#     facet_wrap_paginate ( ~ Drug_Flag_Conc_Rep, nrow = 2, ncol = 2, page =i) +
#     theme(strip.text = element_text(size=14))
#   #theme(plot.title = element_text(hjust = 0.5, size =22))
#   #scale_color_manual(values = rep(getPalette(18), each = 2), name= "Recip_well" )
#   print (page)
# }
# 
# dev.off()



#------check replicate correlation-------
#replicate correlation

rep_corr = function (x) {
  x = x[(x$Drug != "ONE" | x$Drug != "NOT"), ]
  rep1 = subset(x, Replicate == 1)
  rep2 = subset(x, Replicate == 2)
  rep1 = dplyr::rename (rep1, OD_r1 = OD, OD_no1st_r1 = OD_no1st, AUC_norm_bug_r1 = AUC_norm_bug, AUC_well_ODno1st_r1 = AUC_well_ODno1st)
  rep2 = dplyr::rename (rep2, OD_r2 = OD, OD_no1st_r2 = OD_no1st, AUC_norm_bug_r2 = AUC_norm_bug, AUC_well_ODno1st_r2 = AUC_well_ODno1st)
  rep1$Replicate = NULL
  rep2$Replicate = NULL
  rep1$wellID = NULL
  rep2$wellID = NULL
  allrep = merge(rep1, rep2, by = c("Drug_Flag_Conc", "Time_points", "position", "flag_fixed", "batch", "run", "Drug", "plate_name", "donor", "donor_conc"))
  return(allrep)
}
# fix the replicate correlation plot
allrep = lapply(annotated, rep_corr)

#this is to plot
#allrep = lapply(allrep, function (x) transform(x, slope = x$AUC_norm_bug_r1/x$AUC_norm_bug_r1))

#plot replicate correlation first for tsb plates to check recipient well intrinsic behaviour
#then for all plates
#unlist

allrep_unl = plyr::ldply(allrep, 'data.frame' )

#extract one plate with clearly bad correlation to check this
onlyrs = allrep_unl[ , c("Time_points", "OD_r1", "OD_r2", "OD_no1st_r1", "OD_no1st_r2", "AUC_well_ODno1st_r1",
                         "AUC_well_ODno1st_r2", "AUC_norm_bug_r1", "AUC_norm_bug_r2", "plate_name", "position",
                         "run", "batch", "donor", "donor_conc", "Drug", "Drug_Flag_Conc", "flag_fixed")]
#get r per plate
onlyrs$ODr_corr = NULL
## NEW
onlyrs = unique(onlyrs)
maxtp = max(onlyrs$Time_points)
onlyrs = subset(onlyrs, Time_points == min(c(maxtp, trimtp)))

rcorr = onlyrs %>% dplyr::group_by(plate_name) %>%
  dplyr::summarise(ODr_corr=cor(OD_no1st_r1,OD_no1st_r2)) %>% ungroup()
onlyrs = inner_join(onlyrs, rcorr, by= "plate_name")

plot_reps <- function (df, plot_dir) {
  
  pdf(paste0(plot_dir, "rep_corr_OD_allplates_scalesfree.pdf"), width=36, height=25) #onefile=TRUE
  pages = ceiling(length(unique(df$plate_name)) / 15)
  #if you want 96 panels per page (12x8: good format for titration plates) use the following
  #pages = ceiling(length(unique(annotated_read$Drug_Flag_Conc)) / 54)
  
  for (i in seq_len(pages)){
    page <- ggplot(df, aes(x=OD_no1st_r1, y=OD_no1st_r2)) +
      geom_point(aes(size = 0.8, colour = as.factor(flag_fixed)), alpha = 0.7) +
      geom_text_repel(data=subset(df, abs(OD_no1st_r1 - OD_no1st_r2) > 0.25),
                      aes(label=Drug, colour = as.factor(flag_fixed)), size = 8) +
      facet_wrap_paginate(run ~ plate_name, ncol= 5, nrow = 3, page = i) +
      geom_point(data=subset(df, Drug == "BUG"), colour="black", shape = 8, size = 4) +
      
      geom_smooth(method = "lm", se=F, color="black", formula = y ~ x, linetype="dashed", size = 0.5) +
      ggpubr::stat_cor(method = "pearson", label.y = max(df$OD_no1st_r2) - 0.05, 
                       label.x = min(df$OD_no1st_r1), size = 14) +
      # annotate(x = min(df$AUC_norm_bug_r1) + 0.05, y = max(df$AUC_norm_bug_r1) - 0.05,
      #          label = paste0("R = ", round(unique(df$ODr_corr), digits = 2)),"text", size=8) +
      geom_abline(slope=1, size = 0.3) +
      labs(title = "Replicate correlation", x="OD_no1st rep1", y="OD_no1st rep2", colour = "Drug") +
      #coord_capped_cart(bottom='both', left='both') +
      coord_fixed() +
      # scale_x_continuous(breaks = seq(0, 0.45, by = 0.025), limits = c(0, 0.45)) +
      # scale_y_continuous(breaks = seq(0, 0.45, by = 0.025), limits = c(0, 0.45)) +
      theme(plot.title = element_text(hjust = 0.5,size=18),
            panel.spacing = unit(1, "lines"),
            axis.text.x = element_text(size=14, angle=60, hjust=1),
            axis.text.y = element_text(size=14),
            legend.position = 'none',
            strip.text = element_text(size=26),
            aspect.ratio = 1)
    
    print(page)
    
  }
  dev.off()
}  

#in here all of them, also tsb plates so it's ok
plot_reps(onlyrs, plot_dir = figdir)

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

#plot and make lighter
robmean <- lapply(annotated, function(x) {
  x <- x[ , c("Time_h", "Time_points", "OD", "OD_no1st", "Replicate","AUC_well_ODno1st",
              "Drug", "Drug_Flag_Conc", "flag_fixed", "donor", "donor_conc",
              "plate_name", "position", "run", "batch",
              "wellID", "AUCnorm_bug_mean", "AUCnorm_bug_rob", "AUCnorm_bug_med")]
  return(x) }
)

time6 = lapply(annotated, function(x) {
  dplyr::group_by(x, Well) %>%
    dplyr::mutate(time_diff = abs(Time_h - max_hour) ) %>%
    dplyr::filter(time_diff <= min(time_diff)) %>%
    dplyr::select(-time_diff)
})
time6 <- data.table::rbindlist(time6)


#remove runs 15 and 18 from time 6 otherwise it screws the scale of plots
time6 = subset(time6, run != 15 & run != 18)


pdf(paste0(figdir, "boxplot_plates_AUC_robmean.pdf"), width=32, height=10)

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

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_mean, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD no1st normalised by mean of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_mean), max(time6$AUCnorm_bug_mean))) +
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

#--------------------ordered by donor
time6_bydonor = time6[order(time6$donor), ]
time6$plate_name <- factor(time6$plate_name, levels = unique(time6_bydonor$plate_name))

#OD, OD bg sub
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

#AUC raw

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

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_mean, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD no1st normalised by mean of BUG wells by donor") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_mean), max(time6$AUCnorm_bug_mean))) +
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

dev.off()



#--------------------only tsb plates----------------------
time6 = subset(time6, donor == "TSB")

pdf(paste0(figdir, "boxplot_plates_AUC_robmean_onlyTSB_allruns.pdf"), width=32, height=10)

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

#AUC bg sub
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

ggplot (time6, aes(x= plate_name, y = AUCnorm_bug_mean, colour = factor(run))) +
  geom_boxplot() +
  geom_point(data=subset(time6, Drug == "BUG"), colour="black") +
  geom_point(data=subset(time6, Drug == "NOT"), colour="red", shape = 8) +
  geom_point(data=subset(time6, Drug == "ONE"), colour="blue", shape = 2) +
  ggtitle("AUC from OD no1st normalised by mean of BUG wells") +
  scale_y_continuous(limits = c(min(time6$AUCnorm_bug_mean), max(time6$AUCnorm_bug_mean))) +
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


dev.off()


## Overshooting drug analysis
## April 1, 2019
## V. Kim

# Single Drug Fitness of the two strains
argv = commandArgs(trailingOnly = TRUE)
# argv[1] = strain 1 data directory
# argv[2] = strain 2 data directory
# argv[3] = strain 3 data directory
# argv = c("~/Documents/embl/screen_K", "~/Documents/embl/screen_M", "~/Documents/embl/screen_D")

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


get.OD.data <- function(home_dir) {
  batches = grep("batch", list.files(home_dir), value = T)
  
  alldata = lapply(batches, function(b) {
    # load annotated data
    batch <- list.files(paste0(home_dir, b))
    # process batch 1 data
    batch_annot = lapply(paste0(home_dir, b, "/", batch), read.csv, header=T,
                         check.names = FALSE)
    names(batch_annot) <- file_path_sans_ext(batch)
    batch_annot = plyr::ldply(batch_annot, data.frame, .id = "plate")
    batch_annot
  })
  
  names(alldata) = batches
  alldata_df = plyr::ldply(alldata, data.frame, .id=NULL)
  
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

# maximum well fitness in plates with high growth
maxfitwell = 0.9


data.1 = get.fitness.data(home_dir = paste0(argv[1], 
                                            "/screen",
                                            strain_1,
                                            "_annotated_reads/"))
# filter out dead plates
plate_med = get_plate_median(data.1)
data.1 = inner_join(data.1, plate_med, by = "plate")

# subset overgrowing wells
data.1.overgr = dplyr::filter(data.1, donor != "TSB" & Drug != "BUG") %>%
  dplyr::filter(AUCnorm_bug_rob >= maxfitwell)


overgr.heat = dplyr::group_by(data.1.overgr, donor, Drug) %>%
  dplyr::summarise(n=n())


filler = expand.grid(donor = unique(data.1$donor), 
                     Drug = unique(data.1$Drug))
overgr.heat = dplyr::right_join(overgr.heat, filler)
overgr.heat = dplyr::mutate(overgr.heat, n= ifelse(is.na(n), 0, n))


heatmat = as.data.frame(tidyr::spread(overgr.heat, Drug, n))
rownames(heatmat) = heatmat$donor


library(pheatmap)

png(paste0(figdir, "overgrowing-heatmap", strain_1,".png"),
    width = 11, height = 10, units = "in", res = 300)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(heatmat[,-1], color = inferno(50),
         main = paste("Number of overgrowing (>1) combination wells in strain", strain_1))
setHook("grid.newpage", NULL, "replace")
grid.text("Recipient", y=-0.07, gp=gpar(fontsize=16))
grid.text("Donor", x=-0.07, rot=90, gp=gpar(fontsize=16))
dev.off()

ph1 = pheatmap(heatmat[,-1], color = inferno(50),
               main = paste("Number of overgrowing (>0.8) combination wells in strain", strain_1))

plot(ph1$tree_row)
rowcl = cutree(ph1$tree_row, h=250)

don_over_1 = names(rowcl[rowcl == 1])

plot(ph1$tree_col)
colcl = cutree(ph1$tree_col, h=100)
rec_over_1 = names(colcl[colcl == 1])

#oversh1 = expand.grid(donor = don_over_1, Drug = rec_over_1)

oversh1 = dplyr::arrange(overgr.heat, desc(n)) 
oversh1_df = as.data.frame(table(c(as.character(oversh1$donor), as.character(oversh1$Drug))))


# AND DO THE same for strain M

data.2 = get.fitness.data(home_dir = paste0(argv[2], 
                                            "/screen",
                                            strain_2,
                                            "_annotated_reads/"))

# dead plates already filtered out
plate_med = get_plate_median(data.2)

data.2 = inner_join(data.2, plate_med, by = "plate")

# subset overgrowing wells
data.2.overgr = dplyr::filter(data.2, donor != "TSB" & Drug != "BUG") %>%
  dplyr::filter(AUCnorm_bug_rob >= maxfitwell)


overgr.heat = dplyr::group_by(data.2.overgr, donor, Drug) %>%
  dplyr::summarise(n=n())


filler = expand.grid(donor = unique(data.2$donor), 
                     Drug = unique(data.2$Drug))
overgr.heat = dplyr::right_join(overgr.heat, filler)
overgr.heat = dplyr::mutate(overgr.heat, n= ifelse(is.na(n), 0, n))


heatmat = as.data.frame(tidyr::spread(overgr.heat, Drug, n))
rownames(heatmat) = heatmat$donor


png(paste0(figdir, "overgrowing-heatmap", strain_2,".png"),
    width = 11, height = 10, units = 'in', res=300)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(heatmat[,-1], color = inferno(50),
         main = paste("Number of overgrowing (>0.8) combination wells in strain", strain_2))
setHook("grid.newpage", NULL, "replace")
grid.text("Recipient", y=-0.07, gp=gpar(fontsize=16))
grid.text("Donor", x=-0.07, rot=90, gp=gpar(fontsize=16))
dev.off()


ph2 = pheatmap(heatmat[,-1], color = inferno(50),
               main = paste("Number of overgrowing (>0.8) combination wells in strain", strain_2))

plot(ph2$tree_row)
rowcl = cutree(ph2$tree_row, h=150)

don_over_2 = names(rowcl[rowcl %in% c(2,3)])

plot(ph2$tree_col)
colcl = cutree(ph2$tree_col, h=100)
rec_over_2 = names(colcl[colcl == 1])

#oversh2 = expand.grid(donor = don_over_2, Drug = rec_over_2)
oversh2 = dplyr::arrange(overgr.heat, desc(n)) 
oversh2_df = as.data.frame(table(c(as.character(oversh2$donor), as.character(oversh2$Drug))))

colnames(oversh1_df) = c("drug", "n")
colnames(oversh2_df) = c("drug", "n")

oversh1_df$strain = "K"
oversh2_df$strain = "M"

oversh_KM = bind_rows(oversh1_df, oversh2_df) %>%
  arrange(strain, desc(n))

write.table(oversh_KM,
            file = "~/Desktop/overshooting-KM.tsv",
            sep = "\t", quote = F, row.names = F)

common.oversh = inner_join(oversh1, oversh2)

ov.comb1 = anti_join(oversh1, oversh2)
ov.comb2 = anti_join(oversh2, oversh1)



plate_maps1 = dplyr::group_by(data.1.overgr, donor, Drug, Well) %>%
  dplyr::summarise(n=n())


plate_maps1 = dplyr::mutate(plate_maps1, 
                            row = gsub("^([A-Z])([0-9]+)", "\\1", Well),
                            col = gsub("^([A-Z])([0-9]+)", "\\2", Well))

plate_maps1 = inner_join(plate_maps1, oversh1)

pl_ef1 = ggplot(plate_maps1,
                aes(x = factor(col, levels = 1:24),
                    y = factor(row, levels = rev(LETTERS[1:16])), fill = n)) +
  geom_tile(colour = "white", size = 0.25) +
  scale_fill_gradient(low="white", high = "blue") +
  labs(x = "", y = "",
       title = paste("Number of overgrowing cases split by donor in", strain_1),
       fill = "Num. overgrowing") + 
  facet_wrap(~ donor) + 
  coord_fixed() + 
  theme_grey(base_size=8)+
  #theme options
  theme(
    axis.text=element_text(face="bold"),
    axis.text.x = element_text(size=7, angle=60, hjust=1), 
    panel.background = element_blank(),
    axis.ticks=element_line(size=0.4),
    plot.background=element_blank(),
    strip.text = element_text(size=10),
    panel.spacing = unit(2, "lines"),
    plot.title = element_text(hjust = 0.5, size=12))

ggsave(pl_ef1,
       width = 12, height = 10,
       filename = paste0(figdir, "plate-effect-overgrowing", strain_1, ".png"))


plate_maps2 = dplyr::group_by(data.2.overgr, donor, Drug, Well) %>%
  dplyr::summarise(n=n())
plate_maps2 = dplyr::mutate(plate_maps2, 
                            row = gsub("^([A-Z])([0-9]+)", "\\1", Well),
                            col = gsub("^([A-Z])([0-9]+)", "\\2", Well))
plate_maps2 = inner_join(plate_maps2, oversh2)

pl_ef2 = ggplot(plate_maps2,
                aes(x = factor(col, levels = 1:24),
                    y = factor(row, levels = rev(LETTERS[1:16])), fill = n)) +
  geom_tile(colour = "white", size = 0.25) +
  scale_fill_gradient(low="white", high = "blue") +
  labs(x = "", y = "",
       title = paste("Number of overgrowing cases split by donor in", strain_2),
       fill = "Num. overgrowing") + 
  facet_wrap(~ donor) + 
  coord_fixed() + 
  theme_grey(base_size=8)+
  #theme options
  theme(
    axis.text=element_text(face="bold"),
    axis.text.x = element_text(size=7, angle=60, hjust=1), 
    panel.background = element_blank(),
    axis.ticks=element_line(size=0.4),
    plot.background=element_blank(),
    strip.text = element_text(size=10),
    panel.spacing = unit(2, "lines"),
    plot.title = element_text(hjust = 0.5, size=12))

ggsave(pl_ef2,
       width = 12, height = 10,
       filename = paste0(figdir, "plate-effect-overgrowing", strain_2, ".png"))




# plot growth curves for strain K
OD.1 = get.OD.data(home_dir = paste0(argv[1], 
                                            "/screen",
                                            strain_1,
                                            "_annotated_reads/"))


pdf(file = paste0(figdir, "growth-curves-overshooting-K.pdf"),
    width = 12, height = 12)
for (i in 1:50) {
  don = oversh1$donor[i]
  rec = oversh1$Drug[i]
  
  gr.curves = ggplot(filter(OD.1, Drug == rec & donor == don ), 
              aes(x = Time_h, y = OD_no1st, 
                  colour = factor(Replicate))) + 
    geom_point() + 
    geom_line() + facet_wrap(~ plate + Drug_Flag_Conc) +
    labs(colour = "Replicate") +
    theme_bw() +
    theme(aspect.ratio = 1)
  
  print(gr.curves)
  
}
dev.off()



# plot growth curves for strain M
OD.2 = get.OD.data(home_dir = paste0(argv[2], 
                                     "/screen",
                                     strain_2,
                                     "_annotated_reads/"))


pdf(file = paste0(figdir, "growth-curves-overshooting-M.pdf"),
    width = 12, height = 12)
for (i in 1:50) {
  don = oversh2$donor[i]
  rec = oversh2$Drug[i]
  
  gr.curves = ggplot(filter(OD.2, Drug == rec & donor == don ), 
                     aes(x = Time_h, y = OD_no1st, 
                         colour = factor(Replicate))) + 
    geom_point() + 
    geom_line() + facet_wrap(~ plate + Drug_Flag_Conc) +
    labs(colour = "Replicate") +
    theme_bw() +
    theme(aspect.ratio = 1)
  
  print(gr.curves)
  
}
dev.off()






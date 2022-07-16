## Thresholds for filtering
## V. Kim
## April 12, 2019
argv = commandArgs(trailingOnly = TRUE)
# argv = c("~/Documents/embl/screen_K", "~/Documents/embl/screen_M",
# "~/Documents/embl/screen_D")
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

data.2 = get.fitness.data(home_dir = paste0(argv[2], 
                                            "/screen",
                                            strain_2,
                                            "_annotated_reads/"))
trimtp = 9
data.3 = get.fitness.data(home_dir = paste0(argv[3], 
                                            "/screen",
                                            strain_3,
                                            "_annotated_reads/"))


rep_cor1 = get_cor(data = data.1, strain = strain_1)
rep_cor2 = get_cor(data = data.2, strain = strain_2)
rep_cor3 = get_cor(data = data.3, strain = strain_3)

rep_cor = bind_rows(rep_cor1, rep_cor2, 
                    filter(rep_cor3, correl > 0.3))


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


plate_med = bind_rows(get_plate_median(data.1),
          get_plate_median(data.2),
          get_plate_median(data.3))

low_cor_plates = filter(rep_cor, correl < 0.7)
low_cor_plates = inner_join(low_cor_plates, plate_med)


low_fitness_plates = filter(plate_med, med_plate < 0.1 & iqr_plate < 0.1)
low_fitness_plates = inner_join(low_fitness_plates, rep_cor) %>%
  filter(correl >= 0.7)

## Generate OD data for the Shiny app
## V. Kim
## Dec. 13, 2019


argv = commandArgs(trailingOnly = TRUE)
# argv[1] = data/ directory
# argv[2] = strain
# argv = c("/g/huber/users/vkim/gitlab/parallel/staph-devel/data/screen_M", "M")
pacman::p_load(dplyr, reshape2, plyr, stringr, gtools, Rmisc,
               data.table, RColorBrewer, zoo, hexbin, tools,
               directlabels, grid, ggpubr, gridExtra,
               readr, ggplot2, stringi, ggrepel, viridis,
               scales, tidyr, sigr, lemon, here, smoothmest,EnvStats)

strain = paste0("screen", argv[2])

script_dir = "~/Documents/embl/gitlab/drug_comb_screen/dataprep/"
if(!dir.exists(script_dir)) script_dir = "dataprep/"

figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/wilcoxon/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/", strain, "/signif/")
}

home_dir = "~/Documents/embl/screen_M/screenM_annotated_reads/"
if(!dir.exists(home_dir)) {
  home_dir = file.path(argv[1], 
                       paste0(strain, "_annotated_reads/"))
}

wilc_dir = "~/Documents/embl/screen_M/wilcoxon/"
if(!dir.exists(wilc_dir)) {
  wilc_dir = file.path(argv[1], "wilcoxon/")
}

# source the script with core functions
source(paste0(script_dir, "code_screen_core_fun.R"))
# downstream processing of screen results
source(paste0(script_dir, "downstream_fun.R"))

#========LOAD ANNOTATED DATA=============================
druglegend = read.table(file.path(argv[1], "legend_gramnegpos.txt"),
                        sep = "\t", header = T,
                        stringsAsFactors = F)

batches = grep("batch", list.files(home_dir), value = T)

alldata = lapply(batches, function(b) {
  # load annotated data
  batch <- list.files(paste0(home_dir, b))
  batch_annot = lapply(paste0(home_dir, b, "/", batch), read.csv, header=T,
                       check.names = FALSE)
  names(batch_annot) <- file_path_sans_ext(batch)
  
  batch_annot = plyr::ldply(batch_annot, data.frame, .id = "plate")
  batch_annot = dplyr::select(batch_annot, plate, donor, donor_conc, Drug, Concentration, Replicate,
                              batch, biol_rep, Time_points, Time_h, OD, OD_no1st)
  batch_annot = dplyr::mutate(batch_annot,  donor = plyr::mapvalues(donor, from = druglegend$code,
                                                         to = druglegend$Drug),
                          Drug = plyr::mapvalues(Drug, from = druglegend$code,
                                               to = druglegend$Drug))
  batch_annot
})

names(alldata) = batches
alldata_df = plyr::ldply(alldata, data.frame, .id=NULL)

saveRDS(alldata_df, file = paste0(argv[1], "/forshiny/combOD.rds"))



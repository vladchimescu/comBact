argv = commandArgs(trailingOnly = TRUE)

# argv[1] = data/ directory
# argv[2] = batch num
# argv[3] = strain


#installing/loading needed packages
pacman::p_load(dplyr, reshape2, plyr, stringr, gtools, Rmisc,
               data.table, RColorBrewer, zoo, hexbin, tools,
               directlabels, grid, ggpubr, gridExtra,
               readr, ggplot2, stringi, ggrepel, viridis,
               scales, tidyr, sigr, lemon, here, smoothmest,EnvStats)


#===============SET ALL INPUT AND OUTPUT DIRECTORIES============================================
#
strain = paste0("screen", argv[3])

script_dir = "~/Documents/embl/gitlab/drug_comb_screen/dataprep/"
if(!dir.exists(script_dir)) script_dir = "dataprep/"

toextract_dir <- "~/Documents/embl/screen_K/screenK_biomek_files/batch3/"
if(!dir.exists(toextract_dir)) toextract_dir = file.path(argv[1],
                                                         paste0(strain, "_biomek_files"),
                                                         paste0("batch", argv[2], "/"))

toannotate_dir <- "~/Documents/embl/screen_K/screenK_extract_reads/to_annotate/batch3/"
if(!dir.exists(toannotate_dir)) {
  toannotate_dir = file.path(argv[1],
                             paste0(strain, "_extract_reads/to_annotate"),
                             paste0("batch", argv[2], "/"))
  dir.create(toannotate_dir)
}

home_dir = "~/Documents/embl/screen_K/screenK_annotated_reads/batch3/"
if(!dir.exists(home_dir)) {
  home_dir = file.path(argv[1], 
                       paste0(strain, "_annotated_reads"),
                       paste0("batch", argv[2], "/"))
  dir.create(home_dir)
}

stats_dir = "~/Documents/embl/screen_K/screenK_stats/batch3-stats/"
if(!dir.exists(stats_dir)) {
  stats_dir = file.path(argv[1], 
                        paste0(strain, "_stats"),
                        paste0("batch", argv[2], "-stats", "/"))
  dir.create(stats_dir)
}

keyplate_dir = "~/Documents/embl/screen_K/"
if(!dir.exists(keyplate_dir)) {
  keyplate_dir = paste0(argv[1], "/")
}

# source the script with core functions
source(paste0(script_dir, "code_screen_core_fun.R"))
#============================IF READS AREN'T ANNOTATED RUN THIS====================================
#extract them from robot files
#list them

list_biomek <- list.files(toextract_dir)
names(list_biomek) <- file_path_sans_ext(list_biomek)

#extract reads
screen_reads <- lapply(paste0(toextract_dir, list_biomek), 
                       Ext_384plates_Eli, 
                       output_dir = toannotate_dir,
                       barcode = FALSE)


#list them
list_reads <- list.files(toannotate_dir)

library(readxl)
name_map = read_xlsx(paste0(keyplate_dir, "key_plates_biomek.xlsx"), sheet=1)
plate_names <- file_path_sans_ext(list_reads)

modify_names = c(grep("^[A-Z]+_[0-9]+$", plate_names),
                 grep("^[A-Z]+_[0-9]+_[A-Z]_[0-9]+_[0-9]+$", plate_names))

plate_names[modify_names] = plyr::mapvalues(plate_names[modify_names],
                                            from = name_map$biomek_name, to = name_map$correct,
                                            warn_missing = F)

add_batch = grep("^[A-Z]+_[0-9]+_[A-Z]_[0-9]+_[0-9]+_[0-9]+$", plate_names)
if(length(add_batch) > 0) {
  plate_names[add_batch] = paste0(plate_names[add_batch], "_",
                                  str_pad(argv[2], pad = 0, width = 2))
}

#read them
screen <- lapply(paste0(toannotate_dir, list_reads),
                 read.table, header=T, sep="\t")
names(screen) <- plate_names


library(growthcurver)


df = screen$U1_2_K_1_24_17_03
gc_out = SummarizeGrowthByPlate(df, plot_fit = T, plot_file ="~/Desktop/growthcurves_U1.pdf")


mapdir = "~/Documents/embl/screen_K"
if(!dir.exists(mapdir)) mapdir = argv[1]
mapname = paste0("screen_map_batch", argv[2], ".txt")

#===========load the plate map===========
map <- read.table(file.path(mapdir, mapname),
                  header =TRUE, sep="\t")

gc_annot = inner_join(dplyr::rename(gc_out, Well = sample),
           map)


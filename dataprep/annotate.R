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

toextract_dir = file.path(argv[1],
                          paste0(strain, "_biomek_files"),
                          paste0("batch", argv[2], "/"))

toannotate_dir = file.path(argv[1],
                           paste0(strain, "_extract_reads/to_annotate"),
                           paste0("batch", argv[2], "/"))
if(!dir.exists(toannotate_dir)) dir.create(toannotate_dir)

home_dir = file.path(argv[1], 
                     paste0(strain, "_annotated_reads"),
                     paste0("batch", argv[2], "/"))
if(!dir.exists(home_dir)) dir.create(home_dir)

stats_dir = file.path(argv[1], 
                      paste0(strain, "_stats"),
                      paste0("batch", argv[2], "-stats", "/"))
if(!dir.exists(stats_dir)) dir.create(stats_dir)

keyplate_dir = paste0(argv[1], "/")


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

#melt them
myfun <- function(x){
  melt(x,id.vars = "Time_h", variable.name = "wells", value.name= "OD")
}

out <- lapply(screen,myfun)

#fix names
colnames = c("Time_h","Well","OD")
out <- lapply(out, setNames, colnames)

#-----check some stuff-----
#get the max OD for each drug plate in the list
lapply(out,function(df) max(df$OD))

# hist(unlist(lapply(out,function(df) max(df$OD))), 30,
#      xlab = "Max well OD")
# #get the range OD for each plate
# lapply(out,function(df) range(df$OD))
# #get the mean OD
# lapply(out,function(df) mean(df$OD))
# #get the max OD for all drug plates at once
# max(unlist(lapply(out,function(df) max(df$OD))))

# #get the min OD at time point = 0 (convenient, u just look at the beginning of the reads file and check if it's the right one)
# OD_at_0 <- lapply(out, function(x) { x[x$Time_h == 0.000000 ,] })
# min_init_OD <- lapply(OD_at_0, function (df) min (df$OD))

# #CHANGE THIS ACCORDINGLY
# output = file("~/Documents/EMBL/combinatorials/min_init_ODs/min_init_OD_screen_K_1.txt", "a")
# min_init_OD <- cbind.data.frame(min_init_OD) %>%
#         reshape2::melt (value.name = "min_start_OD", variable.name= 'plate_name')
# write.table(min_init_OD, file = output, sep = "\t", row.names = FALSE, append=T)
# close(output)

mapdir = argv[1]
mapname = paste0("screen_map_batch", argv[2], ".txt")

#===========load the plate map===========
map <- read.table(file.path(mapdir, mapname),
                  header =TRUE, sep="\t")
#===fix the map for subsequent plotting====
#add leading 0 before digits 1-9 in Conc_Flag (otherwise it thinks 2 > 10)
flag_fixed <- sprintf("%02d", as.numeric(map$Conc_flag))
#put the fixed column back in the map
map <- cbind(map, flag_fixed)
#remove old flag
map$Conc_flag = NULL
#round decimals in concentrations to avoid long labels
map$Concentration <- round(as.numeric(map$Concentration), digits = 4)
#create columns "Flag + real conc" and "Drug_flag_conc" with this digit fixed
map$Flag_Conc <- paste0(map$flag_fixed,"_",map$Concentration)
map$Drug_Flag_Conc <- paste0(map$Drug,"_",map$Flag_Conc)
map$quadr <- sapply(strsplit(as.character(map$quadrant.id), "_"), '[', 2)
map$quadrant.id <- NULL
map$colID <- sprintf("%02d", as.numeric(map$Well_number))
map$rowID = paste0(map$Well_letter,0)
map$wellID = paste0(map$rowID, map$colID)

#==============match reads with the map=============
annotated <- list()

for (i in 1:length(out)) {
  #df <- as.data.frame (out[[i]])
  df = dplyr::distinct(out[[i]])
  annotated [[i]] <- inner_join(df, map, by="Well")
}

names(annotated) <- names(out)

#check if all wells are covered
length(unique(annotated[[1]]$Well))
length(unique(annotated[[1]]$Time_h))
#this should make 384
length(annotated[[1]]$OD)/length(unique(out[[1]]$Time_h))


#==============================adding useful variables=============================
#get a column for each df with the plate id (i.e. the name of the df in the list)
#these both work
# annotated2 <- setNames ( lapply (names(annotated), function(x) {
#   new_fr <- annotated[[x]]
#   new_fr[["plate_name"]] <- x
#   return(new_fr)
# }), names(annotated))

annotated <- sapply (names(annotated), function(x) {
  new_fr <- annotated[[x]]
  new_fr[["plate_name"]] <- x
  return(new_fr)
}, simplify = FALSE, USE.NAMES = TRUE)

#make column with run number
annotated <- lapply(annotated,
                    function (x) transform(x,
                                           run = sapply(strsplit(as.character(x$plate_name), "_"), '[', 5)))

#this was for the plain recipient run
#tst <- lapply(tst, function (x) transform (x, run = rep("1", nrow(x))))

#make column with donor name
annotated <- lapply(annotated, function (x) transform(x, donor = sapply(strsplit(as.character(x$plate_name), "_"), '[', 1)))

#make column with donor concentration
annotated <- lapply(annotated, function (x) transform(x, donor_conc = sapply(strsplit(as.character(x$plate_name), "_"), '[', 2)))
#make it with appending 0 (useful later when comparing with rec concentrations)
annotated <- lapply(annotated, function (x) transform(x, donor_conc_fixed = sprintf("%02d", as.numeric(x$donor_conc))))

#make column with strain
annotated <- lapply(annotated, function (x) transform(x, strain = sapply(strsplit(as.character(x$plate_name), "_"),'[', 3)))

#make column with plate position
annotated <- lapply(annotated, function (x) transform(x, position = sapply(strsplit(as.character(x$plate_name), "_"),'[', 6)))

#make column with biological replicate
annotated <- lapply(annotated, function (x) transform(x, biol_rep = sapply(strsplit(as.character(x$plate_name), "_"),'[', 4)))

#make column with recipient plate batch number
annotated <- lapply(annotated, function (x) transform(x, batch = sapply(strsplit(as.character(x$plate_name), "_"),'[', 7)))

#make column with donor concentration fixed
annotated <- lapply(annotated, function (x) transform(x, donor_conc_fixed = sprintf("%02d", as.numeric(x$donor_conc))))

#make column with recipient concentration fixed -- THIS CHECK IF NEEDED
annotated <- lapply(annotated, function (x) transform(x, flag_fixed = sprintf("%02d", as.numeric(x$flag_fixed))))

#column that corresponds to unique well id (drug_flag_conc + rep)
annotated <- lapply(annotated, function (x) transform(x, Drug_Flag_Conc_Rep = paste0(x$Drug_Flag_Conc, "_", x$Replicate)))

#drop column Well_id
annotated <- lapply(annotated, function(x) { x["Well_id"] <- NULL; x })
annotated <- lapply(annotated, function(x) { x["X3847WP_id"] <- NULL; x })

#add time points
add_tps <- function (df) {
  v <- c(1:length(unique(df$Time_h)))
  df$Time_points <- rep (v)
  return(df)}
annotated <- lapply (annotated, add_tps)

#column plate_tp (useful for plotting bug vs median plate)
annotated <- lapply(annotated, function (x) transform(x, plate_tp = paste0(x$plate_name, "_", x$Time_points)))

annotated <- lapply(annotated, function (x) transform(x, plate_wellID = paste0(x$plate_name, "_", x$wellID)))

#=========================background subtraction===============
# for Staph strains smooth ODs between time points 4 and 5
if(strain %in% c("screenK", "screenM")) {
  annotated <- lapply(annotated, smooth_OD)
}

# this function subtracts the first time point 
# if there are any negative values -> set to 0
#do the bg sub
annotated <- lapply (annotated, bg_sub_mod)
# for Pneumo estimate background from non-spiking curves of the same plate
if (argv[3] == "P") {
  annotated = lapply(annotated, bgsub_estimate)
}

#write the annotated reads as csv files --> this overwrites!!! that's why I have a dedicated folder
for (i in seq_along(annotated)) {
  filename = paste0(home_dir, names(annotated)[i], ".csv")
  write.csv(annotated[[i]], filename, row.names=FALSE, append = TRUE)
}






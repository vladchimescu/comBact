## script for checkerboard assay plotting
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
batches = grep("batch", list.files(home_dir), value = T)

alldata = lapply(batches, function(b) {
  # load annotated data
  batch <- list.files(paste0(home_dir, b))
  batch_annot = lapply(paste0(home_dir, b, "/", batch), read.csv, header=T,
                       check.names = FALSE)
  names(batch_annot) <- file_path_sans_ext(batch)
  
  params = rjson::fromJSON(file = paste0("params/",b, ".json" ))
  batch_annot = compute_plate_stats_OD(annotated = batch_annot,
                                    params = params)
  batch_annot = plyr::ldply(batch_annot, data.frame, .id = "plate")
  batch_annot
})

names(alldata) = batches

alldata_df = plyr::ldply(alldata, data.frame, .id=NULL)

donfit = lapply(unique(alldata_df$batch), function(b) {
  batchdata = filter(alldata_df, batch == b)
  epsbatch = estimate_single_fit(batchdata, norm = "L2", mode = "donor")
  epsbatch
})
names(donfit) = unique(alldata_df$batch)
donfit = plyr::ldply(donfit, .id = "batch")
donfit = tidyr::separate(donfit, donconc,
                         c("donor", "donor_conc", "biol_rep"), sep="_")

recfit = lapply(unique(alldata_df$batch), function(b) {
  batchdata = filter(alldata_df, batch == b)
  epsbatch = estimate_single_fit(batchdata, norm = "L2", mode = "recipient")
  epsbatch
})
names(recfit) = unique(alldata_df$batch)
recfit = plyr::ldply(recfit, .id = "batch")
recfit = tidyr::separate(recfit, recconcrep,
                         c("Drug", "Concentration", "Replicate"), sep="_")

alldata_masked = lapply(unique(alldata_df$batch), function(b) {
  batch_data = filter(alldata_df, batch == b)
  params = rjson::fromJSON(file = paste0("params/batch",b, ".json" ))
  mask_excluded_wells(data = batch_data, params = params)
  
})

alldata_clean = plyr::ldply(alldata_masked)

# load the strongest associations
strong.inter = read.table(file = paste0(wilc_dir, "binomial-strongassoc.tsv"),
                          stringsAsFactors = F, header = T)

strong.inter = arrange(strong.inter, padj, desc(abs(medpbliss)))

params = rjson::fromJSON(file = "params.json")
maxfit = params$maxfit

# load real donor concentrations
donconc = read.table(file.path(argv[1], paste0("donconc", argv[2], ".txt")),
                     sep = "\t", header = T, stringsAsFactors = F)

donfit = inner_join(mutate(donfit, donor_conc = as.integer(donor_conc)), donconc)

## function with no replicate averaging
get.chboard.noav <- function(d1, d2) {
  d1.df = filter(alldata_clean, donor == d1 & Drug == d2)
  
  if(nrow(d1.df)) {
    chboard = distinct(d1.df, plate, Drug, Concentration, Replicate, .keep_all = T)
    chboard = inner_join(chboard, donconc)
    # join with tsb plates
    tsb.ch = filter(recfit, Drug %in% c(d2) & batch %in% unique(chboard$batch))
    tsb.ch = distinct(tsb.ch, batch, Drug, Concentration, Replicate, .keep_all = T)
    tsb.ch = mutate(tsb.ch, donconc = paste0(d1, "_0"),
                    recconc = paste(Drug, Concentration, sep="_"))
    tsb.ch = mutate(tsb.ch, donor = gsub("(.+)(_[0-9])", '\\1', donconc))
    
    # biol_rep is not meaningful for TSB plates, hence add missing biol_rep
    
    tsb.ch = inner_join(mutate(tsb.ch, batch = as.integer(as.character(batch))),
                        distinct(d1.df, batch, biol_rep))
    # use real concentrations
    chboard = mutate(chboard,
                     recconc = paste(Drug, Concentration, sep="_"),
                     donconc = paste(donor, donor_real_conc, sep="_"))
    don.ch = mutate(inner_join(mutate(donfit, batch = as.integer(as.character(batch)),
                                      biol_rep = as.integer(biol_rep)),
                               distinct(d1.df, batch, biol_rep, donor, donor_conc, run)),
                    recconc = paste0(d2, "_0"),
                    donconc = paste(donor, donor_real_conc, sep="_"))
    don.ch = bind_rows(mutate(don.ch, Replicate = 1),
                       mutate(don.ch, Replicate = 2))
    tsb.ch = mutate(tsb.ch, Replicate = as.integer(Replicate))
    
    don.ch = mutate(don.ch, batch = as.integer(as.character(batch)),
                    biol_rep = as.integer(biol_rep))
    
    
    chboard = chboard[,colnames(chboard) %in% intersect(colnames(tsb.ch), colnames(don.ch))]
    chboard = bind_rows(chboard, dplyr::select(tsb.ch, -c(Drug, Concentration)), 
                        dplyr::select(don.ch, -c(donor_conc, donor_real_conc, run)) )
    
    chboard = mutate(chboard, 
                     mode = ifelse(donor == d1, paste(d1, "as donor"),
                                   paste(d1, "as recipient")),
                     batch = paste("Batch", batch))
    
    ## add a BUG well combination (OD = 1)
    bug = distinct(ungroup(chboard), batch, biol_rep, .keep_all = T)
    bug = ungroup(bug) %>% mutate(donconc = paste0(d1, "_0"),
                                  recconc = paste0(d2, "_0"),
                                  ODnorm_bug_rob = 1)
    bug =bind_rows(mutate(bug, Replicate = 1),
                   mutate(bug, Replicate = 2))
    chboard = bind_rows(chboard, bug)
    
    
    chboard
  } else {
    data.frame()
  }
}

# checkerboards without averaging across replicates
pdf(paste0(figdir, "binomial-checkerboard.pdf"),
    onefile = T, width = 12, height = 10)
for (i in 1:nrow(strong.inter)) {
  d1 = strong.inter$don[i]
  d2 = strong.inter$rec[i]
  
  chboard.av = bind_rows(get.chboard.noav(d1 = d1, d2 = d2),
                         get.chboard.noav(d1 = d2, d2 = d1))
  
  # only complete checkerboards
  only.compl = distinct(chboard.av, batch, mode, biol_rep, donconc, recconc) %>%
    group_by(batch, mode, biol_rep) %>%
    dplyr::summarise(n=n()) %>%
    filter(n > 8)
  
  lims = c(min(chboard.av$ODnorm_bug_rob),
           ifelse(max(chboard.av$ODnorm_bug_rob) < 2, max(chboard.av$ODnorm_bug_rob),
                  2))
  
  cbd.df = inner_join(chboard.av, only.compl)
  cbd.df = mutate(cbd.df, biol_rep = paste("biol_rep", biol_rep))
  
  ch.plot = ggplot(cbd.df, aes(x = factor(recconc,
                                          levels = gtools::mixedsort(levels(factor(recconc)))),
                               y = factor(donconc,
                                          levels = gtools::mixedsort(levels(factor(donconc)))),
                               fill = ODnorm_bug_rob)) + 
    geom_tile() + 
    labs(x = "",
         y = "",
         fill = "Normalized OD",
         title = paste("Combination checkerboard", d1, "+", d2,
                       strong.inter$type[i],
                       "padj =", format(strong.inter$padj[i], scientific = T, digits = 2),
                       ", medpbliss =", format(strong.inter$medpbliss[i], digits=2))) + 
    
    scale_fill_gradient(low = "white", high = ifelse(strong.inter$type[i] == "synergy",
                                                     "#4c928e", "#a44157"),
                        limits=lims) + 
    theme_classic() + 
    facet_wrap(~ batch + biol_rep + mode + Replicate, scales = "free") + 
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1,
          axis.text.x = element_text(angle = 90),
          plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  
  print(ch.plot)
  
}

dev.off()


## ROC curves for the final benchmarking set
## September 2, 2020
## V. Kim

library(dplyr)
library(ggplot2)

call_interactions <- function(eps_stats, multpTest_cutoff, strong_cutoff) {
  eps_stats = dplyr::mutate(eps_stats, type = ifelse(padj < multpTest_cutoff,
                                                     ifelse(efsize < -strong_cutoff, 'S', 
                                                            ifelse(efsize > strong_cutoff, 'A', 'N')), 'N'))
  
  eps_stats = dplyr::mutate(eps_stats,
                            error = ifelse(type %in% c("A", "S"),
                                                    ifelse(benchmarking == type, "TP", "FP"),
                                                    ifelse(benchmarking == "N", "TN", "FN")))
  tpvec = table(eps_stats$error)
  
  tpr = ifelse(is.na(tpvec["TP"]), 0,
               ifelse(multpTest_cutoff == 1, 1, tpvec["TP"] / (tpvec["TP"] + tpvec["FN"])))
  fpr = ifelse(is.na(tpvec["FP"]), 0, 
               ifelse(multpTest_cutoff == 1, 1, tpvec["FP"] / (tpvec["FP"] + tpvec["TN"])))
  
  prec = ifelse(is.na(tpvec["TP"]), 0,
                ifelse(multpTest_cutoff == 1, (tpvec["TP"] + tpvec["FN"]) / sum(tpvec), tpvec["TP"] / (tpvec["TP"] + tpvec["FP"])))
  
  rec = ifelse(is.na(tpvec["TP"]), 0,
               ifelse(multpTest_cutoff == 1, 1, tpvec["TP"] / (tpvec["TP"] + tpvec["FN"])))
  
  data.frame(tpr = tpr,
             fpr = fpr,
             prec = prec,
             rec =rec,
             fdr = multpTest_cutoff,
             row.names = NULL)
  
}

simple_auc <- function(TPR, FPR){
  # inputs already sorted, best scores first 
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}

get_roc <- function(strong_cutoff) {
  eli_roc = lapply(seq(0,1,length.out = 200), function(x) call_interactions(bench, multpTest_cutoff = x,
                                                                            strong_cutoff = strong_cutoff))
  roc_df_eli = plyr::ldply(eli_roc, data.frame)
  auc = with(roc_df_eli, simple_auc(tpr, fpr))
  roc_df_eli$thresh = paste0("efsize = ", 
                             format(strong_cutoff, digits = 2),
                             ", AUCROC = ",
                             format(auc, digits = 3))
  roc_df_eli
}

get_tpfp_weakcons <- function(bench, multpTest_cutoff,
                              strong_cutoff, weak_cutoff) {
  newman = filter(screen_df, strain == 'Newman')
  newman =  mutate(newman, type = ifelse(padj < multpTest_cutoff,
                                         ifelse(efsize < -strong_cutoff, 'S', 
                                                ifelse(efsize > strong_cutoff, 'A', 'N')), 'N'))
  
  dsm = filter(screen_df, strain == 'DSM20231')
  dsm =  mutate(dsm, type = ifelse(padj < multpTest_cutoff,
                                   ifelse(efsize < -strong_cutoff, 'S', 
                                          ifelse(efsize > strong_cutoff, 'A', 'N')), 'N'))
  KM = inner_join(newman, dsm, by="comb") %>%
    filter(comb %in% unique(bench$comb)) %>%
    mutate(type.x = ifelse(type.x == "N" & type.y != "N",
                           ifelse(abs(efsize.x) > weak_cutoff, type.y, type.x), type.x)) %>%
    mutate(type.y = ifelse(type.y == "N" & type.x != "N",
                           ifelse(abs(efsize.y) > weak_cutoff, type.x, type.y), type.y))
  
  K = select(KM, comb, type.x, strain.x) %>%
    rename(type = type.x, strain=strain.x)
  
  K = inner_join(K, select(bench, comb, strain, benchmarking))
  
  M = select(KM, comb, type.y, strain.y) %>%
    rename(type = type.y, strain=strain.y)
  M = inner_join(M, select(bench, comb, strain, benchmarking))
  
  bsub = filter(screen_df, strain == "B. subtilis")
  bsub =  mutate(bsub, type = ifelse(padj < multpTest_cutoff,
                                     ifelse(efsize < -strong_cutoff, 'S', 
                                            ifelse(efsize > strong_cutoff, 'A', 'N')), 'N'))
  
  D = inner_join(select(bsub, comb, type, strain),
                 select(bench, comb, strain, benchmarking))
  
  bench = bind_rows(K, M, D)
  
  bench = dplyr::mutate(bench,
                        error = ifelse(type %in% c("A", "S"),
                                       ifelse(benchmarking == type, "TP", "FP"),
                                       ifelse(benchmarking == "N", "TN", "FN")))
  tpvec = table(bench$error)
  
  tpr = ifelse(is.na(tpvec["TP"]), 0,
               ifelse(multpTest_cutoff == 1, 1, tpvec["TP"] / (tpvec["TP"] + tpvec["FN"])))
  fpr = ifelse(is.na(tpvec["FP"]), 0, 
               ifelse(multpTest_cutoff == 1, 1, tpvec["FP"] / (tpvec["FP"] + tpvec["TN"])))
  
  prec = ifelse(is.na(tpvec["TP"]), 0,
                ifelse(multpTest_cutoff == 1, (tpvec["TP"] + tpvec["FN"]) / sum(tpvec), tpvec["TP"] / (tpvec["TP"] + tpvec["FP"])))
  
  rec = ifelse(is.na(tpvec["TP"]), 0,
               ifelse(multpTest_cutoff == 1, 1, tpvec["TP"] / (tpvec["TP"] + tpvec["FN"])))
  
  data.frame(tpr = tpr,
             fpr = fpr,
             prec = prec,
             rec =rec,
             fdr = multpTest_cutoff,
             row.names = NULL)
}

datadir = "~/Documents/embl/gitlab/drug_comb_screen/data/"
figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/"

# mount /g/huber drive or change directories to local
dirs = c('/Volumes/gitlab/parallel/staph-main/data/screen_K/wilcoxon/',
         '/Volumes/gitlab/parallel/staph-main/data/screen_M/wilcoxon/',
         '/Volumes/gitlab/parallel/bacillus-main/data/screen_D/wilcoxon/')
data_screen = lapply(dirs, function(x) read.csv(paste0(x, 'pvals-all.txt'), sep = '\t'))
names(data_screen) = c("Newman", "DSM20231", "B. subtilis")
screen_df = plyr::ldply(data_screen, .id = 'strain')


# first load the benchmarking set without human-targeting drugs
bench = read.csv(paste0(datadir, "benchmarked_310820_noHTDs.txt"),
                 sep = "\t")
# make strain names compatible with the ones in the files on the server
bench = mutate(bench, strain = plyr::mapvalues(strain,
                                               from = c("Bsubtilis",
                                                        "SADSM20231",
                                                        "SANewman"),
                                               to = c("B. subtilis",
                                                      "DSM20231",
                                                      "Newman")))
# make 'comb' id's compatible
bench = mutate(bench, Drug.x = gsub("(.+)_(.+)", "\\1", comb),
               Drug.y = gsub("(.+)_(.+)", "\\2", comb))
bench$comb = purrr::map2_chr(as.character(bench$`Drug.x`),
                                    as.character(bench$`Drug.y`),
                                    ~ paste(sort(c(.x, .y)), collapse = "_"))
bench = select(bench, comb, strain, benchmarking)

bench = inner_join(bench, screen_df,
           by=c("comb", "strain"))
bench = mutate(bench, benchmarking = plyr::mapvalues(benchmarking,
                                        from = c("synergy", "neutral", "antagonism"),
                                        to = c("S", "N", "A")))

combs_nohuman = unique(bench$comb)

# compute the ROC curve
roc_list = lapply(seq(0.07, 0.14, by=0.01), get_roc)
roc_df = plyr::ldply(roc_list)
roc_df = mutate(roc_df, youden = tpr - fpr)

roc_df$thresh[which.max(roc_df$youden)]
roc_df$fdr[which.max(roc_df$youden)]

point_sel = call_interactions(bench, multpTest_cutoff = 0.05,
                              strong_cutoff = 0.1)

thresh_label = unique(filter(roc_df, grepl('efsize = 0.1,', thresh))$thresh)
point_sel$thresh = thresh_label

weak_cons = get_tpfp_weakcons(bench, multpTest_cutoff = 0.05,
                              strong_cutoff = 0.1, weak_cutoff = 0.08)
weak_cons$thresh = thresh_label

p = ggplot(roc_df, aes(x = fpr, y = tpr, color = thresh)) + 
  geom_line(alpha = 0.6, size=1.3) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "False positive rate",
       y = "True positive rate",
       color = "",
       title = "ROC curve, Gram-positive validation set")+
  geom_abline(slope = 1, color = 'blue', linetype = "dotted")+
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8,0.2)) +
  geom_point(data = point_sel, show.legend = F, size=2) +
  geom_point(data = weak_cons, shape=8, show.legend = F, size=2) +
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("FDR", format(fdr, digits = 2))),
                           nudge_y=0.01,
                           nudge_x = -0.04,
                           segment.color = NA,
                           show.legend = F) +
  ggrepel::geom_text_repel(data = weak_cons,
                           aes(label = "Weak: 0.08"),
                           nudge_y = 0.02,
                           nudge_x = -0.02,
                           segment.color = NA,
                           show.legend = F) 
ggsave(p, filename = paste0(figdir, "ROC-Elisabetta-benchmarking-noHTDs.pdf"),
       width = 8)


# make precision-recall plot
roc_df = group_by(roc_df, thresh) %>%
  mutate(prec = ifelse(prec == 0, 
                       prec[fpr == min( fpr[fpr!=min(fpr)] )],
                       prec),
         f1score = 2*prec*rec / (prec + rec))

g = ggplot(roc_df, aes(x = rec, y = prec, color = thresh))+
  geom_line(alpha = 0.5)+
  geom_point(alpha=0.5)+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  labs(x = "Recall",
       y = "Precision",
       color ="",
       title = "Precision-recall, Gram-positive validation set")  +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_point(data = weak_cons, shape = 8, size = 3, show.legend = F) +
  geom_point(data = point_sel, shape = 8, size = 3, show.legend = F) + 
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("q <", format(fdr, digits = 2))),
                           segment.color = NA,
                           nudge_y = 0.03,
                           nudge_x = 0.04,
                           show.legend = F) +
  ggrepel::geom_text_repel(data = weak_cons,
                           aes(label = "Weak: 0.08"),
                           segment.colour = NA,
                           nudge_y = 0,
                           nudge_x = 0.08,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "precrecall-Elisabetta-benchmarking-noHTDs.pdf"),
       width = 9)

# add weak thresholds for all effect size thresholds
weak_labels = distinct(roc_df, thresh)
weak_labels$efsize = seq(0.07, 0.14, by=0.01)

weak_cons = lapply(weak_labels$efsize, function(x) {
  get_tpfp_weakcons(bench, multpTest_cutoff = 0.05,
                    strong_cutoff = x,
                    weak_cutoff = x-0.02)
})
names(weak_cons) = weak_labels$efsize
weak_cons = plyr::ldply(weak_cons, .id='efsize')

weak_cons = inner_join(mutate(weak_cons, efsize = as.character(efsize)),
                       mutate(weak_labels, efsize = as.character(efsize)))
p = ggplot(roc_df, aes(x = fpr, y = tpr, color = thresh)) + 
  geom_line(alpha = 0.6, size=1.3) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "False positive rate",
       y = "True positive rate",
       color = "",
       title = "ROC curve, Gram-positive validation set")+
  geom_abline(slope = 1, color = 'blue', linetype = "dotted")+
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8,0.2)) +
  geom_point(data = point_sel, show.legend = F, size=2) +
  geom_point(data = weak_cons, shape=8, show.legend = F, size=2) +
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("FDR", format(fdr, digits = 2))),
                           nudge_y=0.01,
                           nudge_x = -0.04,
                           segment.color = NA,
                           show.legend = F)
ggsave(p, filename = paste0(figdir, "ROC-benchmarking-noHTDs-allweak-thresholds.pdf"),
       width = 8)

g = ggplot(roc_df, aes(x = rec, y = prec, color = thresh))+
  geom_line(alpha = 0.5)+
  geom_point(alpha=0.5)+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  labs(x = "Recall",
       y = "Precision",
       color ="",
       title = "Precision-recall, Gram-positive validation set")  +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.2,0.2) )+
  geom_point(data = weak_cons, shape =8, size = 3, show.legend = F) +
  geom_point(data = point_sel, size = 3, show.legend = F) + 
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("q <", format(fdr, digits = 2))),
                           segment.color = NA,
                           nudge_y = 0.03,
                           nudge_x = 0.02,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "precrecall-benchmarking-noHTDs-allweak-thresholds.pdf"),
       width = 9)

# save the data frames to plot ROC and PR curves
# save(roc_df, point_sel, weak_cons, file = '~/Desktop/ROC-PR-noHDT-Bliss.Rdata')


# now benchmarking with human drugs
bench = read.csv(paste0(datadir, "benchmarked_310820.txt"),
                 sep = "\t")
# make strain names compatible with the ones in the files on the server
bench = mutate(bench, strain = plyr::mapvalues(strain,
                                               from = c("Bsubtilis",
                                                        "SADSM20231",
                                                        "SANewman"),
                                               to = c("B. subtilis",
                                                      "DSM20231",
                                                      "Newman")))
# make 'comb' id's compatible
bench = mutate(bench, Drug.x = gsub("(.+)_(.+)", "\\1", comb),
               Drug.y = gsub("(.+)_(.+)", "\\2", comb))
bench$comb = purrr::map2_chr(as.character(bench$`Drug.x`),
                             as.character(bench$`Drug.y`),
                             ~ paste(sort(c(.x, .y)), collapse = "_"))
bench = select(bench, comb, strain, benchmarking)

bench = inner_join(bench, screen_df,
                   by=c("comb", "strain"))
bench = mutate(bench, benchmarking = plyr::mapvalues(benchmarking,
                                                     from = c("synergy", "neutral", "antagonism"),
                                                     to = c("S", "N", "A")))

# compute the ROC curve
roc_list = lapply(seq(0.07, 0.14, by=0.01), get_roc)
roc_df = plyr::ldply(roc_list)
roc_df = mutate(roc_df, youden = tpr - fpr)

roc_df$thresh[which.max(roc_df$youden)]
roc_df$fdr[which.max(roc_df$youden)]

point_sel = call_interactions(bench, multpTest_cutoff = 0.05,
                              strong_cutoff = 0.1) 
thresh_label = unique(filter(roc_df, grepl('efsize = 0.1,', thresh))$thresh)
point_sel$thresh = thresh_label

# this function has to be tweaked for the set with human-targeting drugs
get_tpfp_weakcons <- function(bench, multpTest_cutoff,
                              strong_cutoff, weak_cutoff) {
  newman = filter(screen_df, strain == 'Newman')
  newman =  mutate(newman, type = ifelse(padj < multpTest_cutoff,
                                         ifelse(efsize < -strong_cutoff, 'S', 
                                                ifelse(efsize > strong_cutoff, 'A', 'N')), 'N'))
  
  dsm = filter(screen_df, strain == 'DSM20231')
  dsm =  mutate(dsm, type = ifelse(padj < multpTest_cutoff,
                                   ifelse(efsize < -strong_cutoff, 'S', 
                                          ifelse(efsize > strong_cutoff, 'A', 'N')), 'N'))
  KM = inner_join(newman, dsm, by="comb") %>%
    filter(comb %in% unique(bench$comb)) %>%
    mutate(type.x = ifelse(type.x == "N" & type.y != "N",
                           ifelse(abs(efsize.x) > weak_cutoff, type.y, type.x), type.x)) %>%
    mutate(type.y = ifelse(type.y == "N" & type.x != "N",
                           ifelse(abs(efsize.y) > weak_cutoff, type.x, type.y), type.y))
  
  K = select(KM, comb, type.x, strain.x) %>%
    rename(type = type.x, strain=strain.x)
  
  K = inner_join(K, select(bench, comb, strain, benchmarking))
  
  M = select(KM, comb, type.y, strain.y) %>%
    rename(type = type.y, strain=strain.y)
  M = inner_join(M, select(bench, comb, strain, benchmarking))
  
  bsub = filter(screen_df, strain == "B. subtilis")
  bsub =  mutate(bsub, type = ifelse(padj < multpTest_cutoff,
                                     ifelse(efsize < -strong_cutoff, 'S', 
                                            ifelse(efsize > strong_cutoff, 'A', 'N')), 'N'))
  
  D = inner_join(select(bsub, comb, type, strain),
                 select(bench, comb, strain, benchmarking))
  
  comb_missing = anti_join(bench, bind_rows(K, M, D), by=c("comb", "strain")) %>%
    mutate(type = ifelse(padj < multpTest_cutoff,
                         ifelse(efsize < -strong_cutoff, 'S', 
                                ifelse(efsize > strong_cutoff, 'A', 'N')), 'N')) %>%
    select(comb, strain, type) %>%
    inner_join(select(bench, comb, strain, benchmarking))
  
  bench = bind_rows(K, M, D, comb_missing)
  
  bench = dplyr::mutate(bench,
                        error = ifelse(type %in% c("A", "S"),
                                       ifelse(benchmarking == type, "TP", "FP"),
                                       ifelse(benchmarking == "N", "TN", "FN")))
  tpvec = table(bench$error)
  
  tpr = ifelse(is.na(tpvec["TP"]), 0,
               ifelse(multpTest_cutoff == 1, 1, tpvec["TP"] / (tpvec["TP"] + tpvec["FN"])))
  fpr = ifelse(is.na(tpvec["FP"]), 0, 
               ifelse(multpTest_cutoff == 1, 1, tpvec["FP"] / (tpvec["FP"] + tpvec["TN"])))
  
  prec = ifelse(is.na(tpvec["TP"]), 0,
                ifelse(multpTest_cutoff == 1, (tpvec["TP"] + tpvec["FN"]) / sum(tpvec), tpvec["TP"] / (tpvec["TP"] + tpvec["FP"])))
  
  rec = ifelse(is.na(tpvec["TP"]), 0,
               ifelse(multpTest_cutoff == 1, 1, tpvec["TP"] / (tpvec["TP"] + tpvec["FN"])))
  
  data.frame(tpr = tpr,
             fpr = fpr,
             prec = prec,
             rec =rec,
             fdr = multpTest_cutoff,
             row.names = NULL)
}

weak_cons = get_tpfp_weakcons(bench, multpTest_cutoff = 0.05,
                              strong_cutoff = 0.1, weak_cutoff = 0.08)
weak_cons$thresh = thresh_label

p = ggplot(roc_df, aes(x = fpr, y = tpr, color = thresh)) + 
  geom_line(alpha = 0.6, size=1.3) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "False positive rate",
       y = "True positive rate",
       color = "",
       title = "ROC curve, Gram-positive validation set")+
  geom_abline(slope = 1, color = 'blue', linetype = "dotted")+
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8,0.2)) +
  geom_point(data = point_sel, show.legend = F, size=2) +
  geom_point(data = weak_cons, shape=8, show.legend = F, size=2) +
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("FDR", format(fdr, digits = 2))),
                           nudge_y=0.02,
                           nudge_x = -0.04,
                           segment.color = NA,
                           show.legend = F) +
  ggrepel::geom_text_repel(data = weak_cons,
                           aes(label = "Weak: 0.08"),
                           nudge_y = 0.02,
                           segment.color = NA,
                           show.legend = F) 
ggsave(p, filename = paste0(figdir, "ROC-Elisabetta-benchmarking-withHTDs.pdf"),
       width = 8)


# make precision-recall plot
roc_df = group_by(roc_df, thresh) %>%
  mutate(prec = ifelse(prec == 0, 
                       prec[fpr == min( fpr[fpr!=min(fpr)] )],
                       prec),
         f1score = 2*prec*rec / (prec + rec))

g = ggplot(roc_df, aes(x = rec, y = prec, color = thresh))+
  geom_line(alpha = 0.5)+
  geom_point(alpha=0.5)+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  labs(x = "Recall",
       y = "Precision",
       color ="",
       title = "Precision-recall, Gram-positive validation set")  +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_point(data = weak_cons, shape = 8, size = 3, show.legend = F) +
  geom_point(data = point_sel, shape = 8, size = 3, show.legend = F) + 
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("q <", format(fdr, digits = 2))),
                           segment.color = NA,
                           nudge_y = 0.03,
                           nudge_x = 0.02,
                           show.legend = F) +
  ggrepel::geom_text_repel(data = weak_cons,
                           aes(label = "Weak: 0.08"),
                           segment.colour = NA,
                           nudge_y = 0,
                           nudge_x = 0.08,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "precrecall-Elisabetta-benchmarking-withHTDs.pdf"),
       width = 9)

# add weak thresholds for all effect size thresholds
weak_labels = distinct(roc_df, thresh)
weak_labels$efsize = seq(0.07, 0.14, by=0.01)

weak_cons = lapply(weak_labels$efsize, function(x) {
  get_tpfp_weakcons(bench, multpTest_cutoff = 0.05,
                    strong_cutoff = x,
                    weak_cutoff = x-0.02)
})
names(weak_cons) = weak_labels$efsize
weak_cons = plyr::ldply(weak_cons, .id='efsize')

weak_cons = inner_join(mutate(weak_cons, efsize = as.character(efsize)),
                       mutate(weak_labels, efsize = as.character(efsize)))
p = ggplot(roc_df, aes(x = fpr, y = tpr, color = thresh)) + 
  geom_line(alpha = 0.6, size=1.3) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "False positive rate",
       y = "True positive rate",
       color = "",
       title = "ROC curve, Gram-positive validation set")+
  geom_abline(slope = 1, color = 'blue', linetype = "dotted")+
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8,0.2)) +
  geom_point(data = point_sel, show.legend = F, size=2) +
  geom_point(data = weak_cons, shape=8, show.legend = F, size=2) +
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("FDR", format(fdr, digits = 2))),
                           nudge_y=0.01,
                           nudge_x = -0.04,
                           segment.color = NA,
                           show.legend = F)
ggsave(p, filename = paste0(figdir, "ROC-benchmarking-withHTDs-allweak-thresholds.pdf"),
       width = 8)

# make precision-recall plot
roc_df = group_by(roc_df, thresh) %>%
  mutate(prec = ifelse(prec == 0, 
                       prec[fpr == min( fpr[fpr!=min(fpr)] )],
                       prec),
         f1score = 2*prec*rec / (prec + rec))

g = ggplot(roc_df, aes(x = rec, y = prec, color = thresh))+
  geom_line(alpha = 0.5)+
  geom_point(alpha=0.5)+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  labs(x = "Recall",
       y = "Precision",
       color ="",
       title = "Precision-recall, Gram-positive validation set")  +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.2,0.2) )+
  geom_point(data = weak_cons, shape =8, size = 3, show.legend = F) +
  geom_point(data = point_sel, size = 3, show.legend = F) + 
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("q <", format(fdr, digits = 2))),
                           segment.color = NA,
                           nudge_y = 0.03,
                           nudge_x = 0.02,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "precrecall-benchmarking-withHTDs-allweak-thresholds.pdf"),
       width = 9)

# for human-targeting drugs separate benchmarking (ROC and PR curves)
bench = filter(bench, !comb %in% combs_nohuman)

# compute the ROC curve
roc_list = lapply(seq(0.07, 0.14, by=0.01), get_roc)
roc_df = plyr::ldply(roc_list)
roc_df = mutate(roc_df, youden = tpr - fpr)

roc_df$thresh[which.max(roc_df$youden)]
roc_df$fdr[which.max(roc_df$youden)]

point_sel = call_interactions(bench, multpTest_cutoff = 0.05,
                              strong_cutoff = 0.1)

thresh_label = unique(filter(roc_df, grepl('efsize = 0.1,', thresh))$thresh)
point_sel$thresh = thresh_label

p = ggplot(roc_df, aes(x = fpr, y = tpr, color = thresh)) + 
  geom_line(alpha = 0.6, size=1.3) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "False positive rate",
       y = "True positive rate",
       color = "",
       title = "ROC curve, Gram-positive validation set")+
  geom_abline(slope = 1, color = 'blue', linetype = "dotted")+
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8,0.2)) +
  geom_point(data = point_sel, show.legend = F, size=2) +
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("FDR", format(fdr, digits = 2))),
                           nudge_y=0.01,
                           nudge_x = -0.04,
                           segment.color = NA,
                           show.legend = F) 
ggsave(p, filename = paste0(figdir, "ROC-benchmarking-onlyHTDs.pdf"),
       width = 8)


# make precision-recall plot
roc_df = group_by(roc_df, thresh) %>%
  mutate(prec = ifelse(prec == 0, 
                       prec[fpr == min( fpr[fpr!=min(fpr)] )],
                       prec),
         f1score = 2*prec*rec / (prec + rec))

g = ggplot(roc_df, aes(x = rec, y = prec, color = thresh))+
  geom_line(alpha = 0.5)+
  geom_point(alpha=0.5)+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  labs(x = "Recall",
       y = "Precision",
       color ="",
       title = "Precision-recall, Gram-positive validation set")  +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_point(data = point_sel, shape = 8, size = 3, show.legend = F) + 
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("q <", format(fdr, digits = 2))),
                           segment.color = NA,
                           nudge_y = 0.03,
                           nudge_x = 0.04,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "precrecall-benchmarking-onlyHTDs.pdf"),
       width = 9)

# now benchmarking with Loewe additivity model
bench = read.csv(paste0(datadir, "loeweblisscomparison.txt"),
                 sep = "\t")
# make strain names compatible with the ones in the files on the server
bench = mutate(bench, strain = plyr::mapvalues(strain,
                                               from = c("Bsubtilis",
                                                        "SADSM20231",
                                                        "SANewman"),
                                               to = c("B. subtilis",
                                                      "DSM20231",
                                                      "Newman")))
# make 'comb' id's compatible
bench = mutate(bench, Drug.x = gsub("(.+)_(.+)", "\\1", comb),
               Drug.y = gsub("(.+)_(.+)", "\\2", comb))
bench$comb = purrr::map2_chr(as.character(bench$`Drug.x`),
                             as.character(bench$`Drug.y`),
                             ~ paste(sort(c(.x, .y)), collapse = "_"))
bench = select(bench, comb, strain, type_loewe)
bench = rename(bench, benchmarking = type_loewe)
bench = inner_join(bench, screen_df,
                   by=c("comb", "strain"))
bench = mutate(bench, benchmarking = plyr::mapvalues(benchmarking,
                                                     from = c("synergy", "neutral", "antagonism"),
                                                     to = c("S", "N", "A")))

# compute the ROC curve
roc_list = lapply(seq(0.07, 0.14, by=0.01), get_roc)
roc_df = plyr::ldply(roc_list)
roc_df = mutate(roc_df, youden = tpr - fpr)

roc_df$thresh[which.max(roc_df$youden)]
roc_df$fdr[which.max(roc_df$youden)]

point_sel = call_interactions(bench, multpTest_cutoff = 0.05,
                              strong_cutoff = 0.1)

thresh_label = unique(filter(roc_df, grepl('efsize = 0.1,', thresh))$thresh)
point_sel$thresh = thresh_label

weak_cons = get_tpfp_weakcons(bench, multpTest_cutoff = 0.05,
                              strong_cutoff = 0.1, weak_cutoff = 0.08)
weak_cons$thresh = thresh_label

p = ggplot(roc_df, aes(x = fpr, y = tpr, color = thresh)) + 
  geom_line(alpha = 0.6, size=1.3) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "False positive rate",
       y = "True positive rate",
       color = "",
       title = "ROC curve, Gram-positive validation set")+
  geom_abline(slope = 1, color = 'blue', linetype = "dotted")+
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8,0.2)) +
  geom_point(data = point_sel, show.legend = F, size=2) +
  geom_point(data = weak_cons, shape=8, show.legend = F, size=2) +
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("FDR", format(fdr, digits = 2))),
                           nudge_y=0.01,
                           nudge_x = -0.04,
                           segment.color = NA,
                           show.legend = F) +
  ggrepel::geom_text_repel(data = weak_cons,
                           aes(label = "Weak: 0.08"),
                           nudge_y = 0.02,
                           nudge_x = -0.02,
                           segment.color = NA,
                           show.legend = F) 
ggsave(p, filename = paste0(figdir, "ROC-Elisabetta-benchmarking-noHTDs.pdf"),
       width = 8)


# make precision-recall plot
roc_df = group_by(roc_df, thresh) %>%
  mutate(prec = ifelse(prec == 0, 
                       prec[fpr == min( fpr[fpr!=min(fpr)] )],
                       prec),
         f1score = 2*prec*rec / (prec + rec))

g = ggplot(roc_df, aes(x = rec, y = prec, color = thresh))+
  geom_line(alpha = 0.5)+
  geom_point(alpha=0.5)+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  labs(x = "Recall",
       y = "Precision",
       color ="",
       title = "Precision-recall, Gram-positive validation set")  +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_point(data = weak_cons, shape = 8, size = 3, show.legend = F) +
  geom_point(data = point_sel, shape = 8, size = 3, show.legend = F) + 
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("q <", format(fdr, digits = 2))),
                           segment.color = NA,
                           nudge_y = 0.03,
                           nudge_x = 0.04,
                           show.legend = F) +
  ggrepel::geom_text_repel(data = weak_cons,
                           aes(label = "Weak: 0.08"),
                           segment.colour = NA,
                           nudge_y = 0,
                           nudge_x = 0.08,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "precrecall-Elisabetta-benchmarking-noHTDs.pdf"),
       width = 9)

# add weak thresholds for all effect size thresholds
weak_labels = distinct(roc_df, thresh)
weak_labels$efsize = seq(0.07, 0.14, by=0.01)

weak_cons = lapply(weak_labels$efsize, function(x) {
  get_tpfp_weakcons(bench, multpTest_cutoff = 0.05,
                    strong_cutoff = x,
                    weak_cutoff = x-0.02)
})
names(weak_cons) = weak_labels$efsize
weak_cons = plyr::ldply(weak_cons, .id='efsize')

weak_cons = inner_join(mutate(weak_cons, efsize = as.character(efsize)),
                       mutate(weak_labels, efsize = as.character(efsize)))
p = ggplot(roc_df, aes(x = fpr, y = tpr, color = thresh)) + 
  geom_line(alpha = 0.6, size=1.3) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "False positive rate",
       y = "True positive rate",
       color = "",
       title = "ROC curve, Gram-positive validation set")+
  geom_abline(slope = 1, color = 'blue', linetype = "dotted")+
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8,0.2)) +
  geom_point(data = point_sel, show.legend = F, size=2) +
  geom_point(data = weak_cons, shape=8, show.legend = F, size=2) +
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("FDR", format(fdr, digits = 2))),
                           nudge_y=0.01,
                           nudge_x = -0.04,
                           segment.color = NA,
                           show.legend = F)
ggsave(p, filename = paste0(figdir, "ROC-benchmarking-noHTDs-allweak-thresholds.pdf"),
       width = 8)

g = ggplot(roc_df, aes(x = rec, y = prec, color = thresh))+
  geom_line(alpha = 0.5)+
  geom_point(alpha=0.5)+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  labs(x = "Recall",
       y = "Precision",
       color ="",
       title = "Precision-recall, Gram-positive validation set")  +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.2,0.2) )+
  geom_point(data = weak_cons, shape =8, size = 3, show.legend = F) +
  geom_point(data = point_sel, size = 3, show.legend = F) + 
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("q <", format(fdr, digits = 2))),
                           segment.color = NA,
                           nudge_y = 0.03,
                           nudge_x = 0.02,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "precrecall-benchmarking-noHTDs-allweak-thresholds.pdf"),
       width = 9)

# save the data frames to plot ROC and PR curves
#save(roc_df, point_sel, weak_cons, file = '~/Desktop/ROC-PR-noHDT-Loewe.Rdata')

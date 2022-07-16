## try different effect size measures (such as median Bliss, quantiles)
## and compare these with medpbliss ROC and precision-recall

## V. Kim
## Dec. 12 , 2019

library(dplyr)
library(ggplot2)

call_interactions <- function(eps_stats, multpTest_cutoff, strong_cutoff) {
  eps_stats = dplyr::mutate(eps_stats, type = ifelse(padj < multpTest_cutoff,
                                                     ifelse(medpbliss < -strong_cutoff, 'S', 
                                                            ifelse(medpbliss >= strong_cutoff, 'A', 'N')), 'N'))
  
  eps_stats = dplyr::mutate(eps_stats,
                            binomial_error = ifelse(type %in% c("A", "S"),
                                                    ifelse(Benchmarking == "N", "FP", "TP"),
                                                    ifelse(Benchmarking == "N", "TN", "FN")))
  tpvec = table(eps_stats$binomial_error)
  
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

datadir = "~/Documents/embl/gitlab/drug_comb_screen/data/"
figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/"

bench = read.csv(paste0(datadir, "bench_total_final.txt"),
                 sep = "\t")
bench = rename(bench, Benchmarking =type.y) %>%
  mutate(Benchmarking = plyr::mapvalues(Benchmarking,
                                        from = c("synergy", "neutral", "antagonism"),
                                        to = c("S", "N", "A")))

get_roc <- function(strong_cutoff) {
  eli_roc = lapply(seq(0,1,length.out = 200), function(x) call_interactions(bench, multpTest_cutoff = x,
                                                                            strong_cutoff = strong_cutoff))
  roc_df_eli = plyr::ldply(eli_roc, data.frame)
  auc = with(roc_df_eli, simple_auc(tpr, fpr))
  roc_df_eli$thresh = paste0("|efsize| > ", 
                             format(strong_cutoff, digits = 2))
  roc_df_eli$AUCROC = with(roc_df_eli, simple_auc(tpr, fpr))
  roc_df_eli
}

roc_list = lapply(seq(0.01, 0.08, by=0.01), get_roc)
roc_df = plyr::ldply(roc_list)
roc_df = mutate(roc_df, youden = tpr - fpr)

roc_df$thresh[which.max(roc_df$youden)]
roc_df$fdr[which.max(roc_df$youden)]

point_sel = call_interactions(bench, multpTest_cutoff = 0.05,
                              strong_cutoff = 0.04) 
point_sel$thresh = "|efsize| > 0.04"
point_sel$efsize = "Median polish Bliss"


get_tpfp_weakcons <- function(bench, multpTest_cutoff,
                              strong_cutoff, weak_cutoff) {
  newman = read.table("~/Documents/embl/screen_K/wilcoxon/pvals-all.txt",
                      header = T, stringsAsFactors = F)
  newman$strain = "K"
  newman = mutate(newman, padj = p.adjust(pval, method = "BH"))
  newman =  mutate(newman, type = ifelse(padj < multpTest_cutoff,
                                         ifelse(medpbliss <= -strong_cutoff, 'S', 
                                                ifelse(medpbliss >= strong_cutoff, 'A', 'N')), 'N'))
  
  dsm = read.table("~/Documents/embl/screen_M/wilcoxon/pvals-all.txt",
                   header = T, stringsAsFactors = F)
  dsm$strain = "M"
  dsm = mutate(dsm, padj = p.adjust(pval, method = "BH"))
  dsm =  mutate(dsm, type = ifelse(padj < multpTest_cutoff,
                                   ifelse(medpbliss <= -strong_cutoff, 'S', 
                                          ifelse(medpbliss >= strong_cutoff, 'A', 'N')), 'N'))
  KM = inner_join(newman, dsm, by="comb") %>%
    filter(comb %in% unique(bench$uniquecomb)) %>%
    mutate(type.x = ifelse(type.x == "N" & type.y != "N",
                           ifelse(abs(medpbliss.x) > weak_cutoff, type.y, type.x), type.x)) %>%
    mutate(type.y = ifelse(type.y == "N" & type.x != "N",
                           ifelse(abs(medpbliss.y) > weak_cutoff, type.x, type.y), type.y))
  
  K = select(KM, comb, type.x, strain.x) %>%
    rename(uniquecomb = comb, type = type.x, strain=strain.x)
  
  K = inner_join(K, select(bench, uniquecomb, strain, Benchmarking))
  
  M = select(KM, comb, type.y, strain.y) %>%
    rename(uniquecomb = comb, type = type.y, strain=strain.y)
  M = inner_join(M, select(bench, uniquecomb, strain, Benchmarking))
  
  bench = dplyr::mutate(bench, type = ifelse(padj < multpTest_cutoff,
                                             ifelse(medpbliss <= -strong_cutoff, 'S', 
                                                    ifelse(medpbliss >= strong_cutoff, 'A', 'N')), 'N'))
  D = select(bench, type, uniquecomb, strain, Benchmarking) %>%
    filter(strain == "D")
  
  bench = bind_rows(K, M, D)
  
  bench = dplyr::mutate(bench,
                            binomial_error = ifelse(type %in% c("A", "S"),
                                                    ifelse(Benchmarking == "N", "FP", "TP"),
                                                    ifelse(Benchmarking == "N", "TN", "FN")))
  tpvec = table(bench$binomial_error)
  
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
                              strong_cutoff = 0.04, weak_cutoff = 0.02)
weak_cons$thresh = "|efsize| > 0.04"
weak_cons$efsize = "Median polish Bliss"

roc_df_medp = roc_df
roc_df_medp$efsize = "Median polish Bliss"


bench = rename(select(bench, -medpbliss), medpbliss = eps.med.x)
roc_list = lapply(seq(0.01, 0.08, by=0.01), get_roc)
roc_df_median = plyr::ldply(roc_list)
roc_df_median = mutate(roc_df_median, youden = tpr - fpr)
roc_df_median$efsize = "Median Bliss"

roc_df = bind_rows(roc_df_medp,
          roc_df_median)

p = ggplot(roc_df, aes(x = fpr, y = tpr, group= interaction(thresh, efsize),
                       linetype = efsize, color = thresh)) + 
  geom_line(alpha = 0.6, size=1) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "False positive rate",
       y = "True positive rate",
       color = "Threshold",
       linetype = "Effect size measure",
       title = "ROC curve, Gram-positive validation set (OD-based interactions)")+
  geom_abline(slope = 1, color = 'blue', linetype = "dotted")+
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  scale_linetype_manual(values=c("dashed", "solid"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8,0.3)) +
  geom_point(data = point_sel, show.legend = F, size=2) +
  geom_point(data = weak_cons, shape=8, show.legend = F, size=2) +
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("FDR", format(fdr, digits = 2))),
                           nudge_y=0.01,
                           nudge_x = -0.02,
                           segment.color = NA,
                           show.legend = F) +
  ggrepel::geom_text_repel(data = weak_cons,
                           aes(label = "Weak: 0.02"),
                           nudge_y = 0.02,
                           segment.color = NA,
                           show.legend = F) 
ggsave(p, filename = paste0(figdir, "ROC-benchmarking-medp-vs-median.pdf"),
       width = 8)

# make precision-recall plot
roc_df = group_by(roc_df, thresh) %>%
  mutate(prec = ifelse(prec == 0, 
                             prec[fpr == min( fpr[fpr!=min(fpr)] )],
                             prec),
         f1score = 2*prec*rec / (prec + rec))

g = ggplot(roc_df, aes(x = rec, y = prec, group= interaction(thresh, efsize),
                       linetype = efsize, color = thresh))+
  geom_line(alpha = 0.5)+
  geom_point(alpha=0.5)+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#7f78d2", "#ce2e6c",
                                "#dfcdc3", "#719192",
                                "#3c4245", "#2f416d"))+
  scale_linetype_manual(values=c("dashed", "solid"))+
  labs(x = "Recall",
       y = "Precision",
       color = "Threshold",
       linetype = "Effect size measure",
       title = "Precision-recall, Gram-positive validation set (OD-based interactions)")  +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_point(data = weak_cons, shape = 8, show.legend = F) +
  geom_point(data = point_sel, shape = 8, size = 2, show.legend = F) + 
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("q <", format(fdr, digits = 2))),
                           segment.color = NA,
                           nudge_y = 0.02,
                           nudge_x = 0.04,
                           show.legend = F) +
  ggrepel::geom_text_repel(data = weak_cons,
                           aes(label = "Weak: 0.02"),
                           segment.colour = NA,
                           nudge_y = -0.02,
                           nudge_x = 0.06,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "precrecall-benchmarking-medp-vs-median.pdf"),
       width = 9)


# now compare Median polish vs quantiles
bench = read.csv(paste0(datadir, "bench_total_final.txt"),
                 sep = "\t")
bench = rename(bench, Benchmarking =type.y) %>%
  mutate(Benchmarking = plyr::mapvalues(Benchmarking,
                                        from = c("synergy", "neutral", "antagonism"),
                                        to = c("S", "N", "A")))

bench = mutate(bench, efsize = ifelse(medpbliss < 0, eps.q1.x, eps.q3.x))
bench = rename(select(bench, -medpbliss), medpbliss = efsize)
roc_list = lapply(seq(0.08, 0.12, by=0.01), get_roc)
roc_df_quant = plyr::ldply(roc_list)
roc_df_quant = mutate(roc_df_quant, youden = tpr - fpr)
roc_df_quant$efsize = "Bliss Quartiles"

roc_df = bind_rows(filter(roc_df_medp, thresh == "|efsize| > 0.04"),
                   roc_df_quant)

quant_sel = call_interactions(bench, multpTest_cutoff = 0.05,
                              strong_cutoff = 0.1) 
quant_sel$thresh = "|efsize| > 0.1"
quant_sel$efsize = "Bliss Quartiles"


get_tpfp_weakcons <- function(bench, multpTest_cutoff,
                              strong_cutoff, weak_cutoff) {
  newman = read.table("~/Documents/embl/screen_K/wilcoxon/pvals-all.txt",
                      header = T, stringsAsFactors = F)
  newman$strain = "K"
  newman = mutate(newman, padj = p.adjust(pval, method = "BH"))
  newman = mutate(newman, efsize = ifelse(medpbliss < 0, eps.q1, eps.q3))
  newman = rename(select(newman, -medpbliss), medpbliss = efsize)
  newman =  mutate(newman, type = ifelse(padj < multpTest_cutoff,
                                         ifelse(medpbliss <= -strong_cutoff, 'S', 
                                                ifelse(medpbliss >= strong_cutoff, 'A', 'N')), 'N'))
  
  dsm = read.table("~/Documents/embl/screen_M/wilcoxon/pvals-all.txt",
                   header = T, stringsAsFactors = F)
  dsm$strain = "M"
  dsm = mutate(dsm, padj = p.adjust(pval, method = "BH"))
  dsm = mutate(dsm, efsize = ifelse(medpbliss < 0, eps.q1, eps.q3))
  dsm = rename(select(dsm, -medpbliss), medpbliss = efsize)
  dsm =  mutate(dsm, type = ifelse(padj < multpTest_cutoff,
                                   ifelse(medpbliss <= -strong_cutoff, 'S', 
                                          ifelse(medpbliss >= strong_cutoff, 'A', 'N')), 'N'))
  KM = inner_join(newman, dsm, by="comb") %>%
    filter(comb %in% unique(bench$uniquecomb)) %>%
    mutate(type.x = ifelse(type.x == "N" & type.y != "N",
                           ifelse(abs(medpbliss.x) > weak_cutoff, type.y, type.x), type.x)) %>%
    mutate(type.y = ifelse(type.y == "N" & type.x != "N",
                           ifelse(abs(medpbliss.y) > weak_cutoff, type.x, type.y), type.y))
  
  K = select(KM, comb, type.x, strain.x) %>%
    rename(uniquecomb = comb, type = type.x, strain=strain.x)
  
  K = inner_join(K, select(bench, uniquecomb, strain, Benchmarking))
  
  M = select(KM, comb, type.y, strain.y) %>%
    rename(uniquecomb = comb, type = type.y, strain=strain.y)
  M = inner_join(M, select(bench, uniquecomb, strain, Benchmarking))
  
  bench = dplyr::mutate(bench, type = ifelse(padj < multpTest_cutoff,
                                             ifelse(medpbliss <= -strong_cutoff, 'S', 
                                                    ifelse(medpbliss >= strong_cutoff, 'A', 'N')), 'N'))
  D = select(bench, type, uniquecomb, strain, Benchmarking) %>%
    filter(strain == "D")
  
  bench = bind_rows(K, M, D)
  
  bench = dplyr::mutate(bench,
                        binomial_error = ifelse(type %in% c("A", "S"),
                                                ifelse(Benchmarking == "N", "FP", "TP"),
                                                ifelse(Benchmarking == "N", "TN", "FN")))
  tpvec = table(bench$binomial_error)
  
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
weak_quant = get_tpfp_weakcons(bench, multpTest_cutoff = 0.05,
                              strong_cutoff = 0.1, weak_cutoff = 0.08)
weak_quant$thresh = "|efsize| > 0.1"
weak_quant$efsize = "Bliss Quartiles"


p = ggplot(roc_df, aes(x = fpr, y = tpr, group= interaction(thresh, efsize),
                       linetype = efsize, color = thresh)) + 
  geom_line(alpha = 0.6, size=1) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "False positive rate",
       y = "True positive rate",
       color = "Threshold",
       linetype = "Effect size measure",
       title = "ROC curve, Gram-positive validation set (OD-based interactions)")+
  geom_abline(slope = 1, color = 'blue', linetype = "dotted")+
  scale_color_manual(values = c("#ce2e6c",
                                "#5f6769", "#f8a978",
                                "#51eaea",
                                "#dfcdc3", "#71a95a"))+
  scale_linetype_manual(values=c("longdash", "solid"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8,0.3)) +
  geom_point(data = point_sel, show.legend = F, size=2) +
  geom_point(data = weak_cons, shape=8, show.legend = F, size=2) +
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("FDR", format(fdr, digits = 2))),
                           nudge_y=0.01,
                           nudge_x = -0.02,
                           segment.color = NA,
                           show.legend = F) +
  ggrepel::geom_text_repel(data = weak_cons,
                           aes(label = "Weak: 0.02"),
                           nudge_x = -0.04,
                           nudge_y = 0.02,
                           segment.color = NA,
                           show.legend = F) +
  geom_point(data = quant_sel, show.legend = F, size=2) +
  geom_point(data = weak_quant, shape=8, show.legend = F, size=2) +
  ggrepel::geom_text_repel(data = quant_sel,
                           aes(label = paste("FDR", format(fdr, digits = 2))),
                           nudge_y=0.01,
                           nudge_x = -0.02,
                           segment.color = NA,
                           show.legend = F)+
  ggrepel::geom_text_repel(data = weak_quant,
                           aes(label = "Weak: 0.08"),
                           nudge_x = 0.04,
                           segment.color = NA,
                           show.legend = F) 
ggsave(p, filename = paste0(figdir, "ROC-benchmarking-medp-vs-quartiles.pdf"),
       width = 10, height = 9)



# make precision-recall plot
roc_df = group_by(roc_df, thresh) %>%
  mutate(prec = ifelse(prec == 0, 
                       prec[fpr == min( fpr[fpr!=min(fpr)] )],
                       prec),
         f1score = 2*prec*rec / (prec + rec))

g = ggplot(roc_df, aes(x = rec, y = prec, group= interaction(thresh, efsize),
                       linetype = efsize, color = thresh))+
  geom_line(alpha = 0.5)+
  geom_point(alpha=0.5)+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_color_manual(values = c("#ce2e6c",
                                "#5f6769", "#f8a978",
                                "#51eaea",
                                "#dfcdc3", "#71a95a"))+
  scale_linetype_manual(values=c("longdash", "solid"))+
  labs(x = "Recall",
       y = "Precision",
       color = "Threshold",
       linetype = "Effect size measure",
       title = "Precision-recall, Gram-positive validation set (OD-based interactions)")  +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_point(data = point_sel, shape = 8, size = 2, show.legend = F) + 
  ggrepel::geom_text_repel(data = point_sel,
                           aes(label = paste("q <", format(fdr, digits = 2))),
                           segment.color = NA,
                           nudge_y = 0.02,
                           nudge_x = 0.04,
                           show.legend = F) +
  geom_point(data = quant_sel, shape=8, show.legend = F, size=2) +
  ggrepel::geom_text_repel(data = quant_sel,
                           aes(label = paste("q <", format(fdr, digits = 2))),
                           nudge_y=0.01,
                           nudge_x = 0.04,
                           segment.color = NA,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "precrecall-benchmarking-medp-vs-quartiles.pdf"),
       width = 9)

# write.table(format(as.data.frame(distinct(filter(roc_df, abs(fdr - 0.05) < 1e-3),
#                                           thresh, AUCROC, f1score,
#                                           tpr, fpr, prec, rec, efsize)), digits=2 ),
#             file  = "~/Desktop/quartiles_ROC_F1score.txt", col.names = T, quote = F, row.names = F, sep="\t")

## new code here
newman = read.table("~/Documents/embl/screen_K/wilcoxon/pvals-all.txt",
                    header = T, stringsAsFactors = F)
newman$strain = "K"

dsm = read.table("~/Documents/embl/screen_M/wilcoxon/pvals-all.txt",
                 header = T, stringsAsFactors = F)
dsm$strain = "M"

bsub = read.table("~/Documents/embl/screen_D/wilcoxon/pvals-all.txt",
                 header = T, stringsAsFactors = F)
bsub$strain = "D"

# reload the benchmarking file (just to make sure)
bench = read.csv(paste0(datadir, "bench_total_final.txt"),
                 sep = "\t")
bench = rename(bench, Benchmarking =type.y) %>%
  mutate(Benchmarking = plyr::mapvalues(Benchmarking,
                                        from = c("synergy", "neutral", "antagonism"),
                                        to = c("S", "N", "A")))

bench = inner_join(bind_rows(newman, dsm, bsub) %>%
  select(comb, strain, medpconc) %>%
    rename(uniquecomb = comb), bench)

g = ggplot(bench, aes(x = medpbliss, y = medpconc, color = Benchmarking))+
  geom_point(size = 3, alpha= 0.8)+
  labs(x = "Median polish Bliss",
       y = "Median polish + Concentration Effect",
       color = "Type") + 
  theme_classic()+
  scale_color_manual(values = c("#FFCC33", "#8c8c8c", "#009999")) +
  xlim(c(-0.4, 0.4))+
  ylim(c(-0.4,0.4))+
  geom_abline(slope = 1, linetype ="dotted")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(g, filename = paste0(figdir, "Eli-validation-medp-vs-medpconc.pdf"))

bench = rename(select(bench, -medpbliss), medpbliss = medpconc)

roc_list = lapply(seq(0.08, 0.12, by=0.01), get_roc)
roc_df_medpconc = plyr::ldply(roc_list)
roc_df_medpconc = mutate(roc_df_medpconc, youden = tpr - fpr)
roc_df_medpconc$efsize = "Medpolish + Conc Bliss"

roc_df = bind_rows(roc_df_medpconc, roc_df_quant)


p = ggplot(roc_df, aes(x = fpr, y = tpr, group= interaction(thresh, efsize),
                                linetype = efsize, color = thresh)) + 
  geom_line(alpha = 0.6, size=1.3) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "False positive rate",
       y = "True positive rate",
       color = "",
       title = "ROC curve, Gram-positive validation set (OD-based interactions)")+
  geom_abline(slope = 1, color = 'blue', linetype = "dotted")+
  scale_color_manual(values = c("#5f6769", "#6b591d",
                                "#ce2e6c", "#7f78d2", 
                                "#dfcdc3"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8,0.2)) 

ggsave(p, filename = paste0(figdir, "ROC-benchmarking-medpconc-vs-quartilees.pdf"),
       width = 8)


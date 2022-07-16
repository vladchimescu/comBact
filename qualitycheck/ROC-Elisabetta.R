## plot the roc curve for Eli's validation set
## V. Kim
## Nov. 7 , 2019

library(dplyr)
library(ggplot2)

call_interactions <- function(eps_stats, multpTest_cutoff, strong_cutoff) {
  eps_stats = dplyr::mutate(eps_stats, type = ifelse(padj <= multpTest_cutoff,
                                                     ifelse(medpbliss <= -strong_cutoff, 'S', 
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

g = ggplot(bench, aes(x = eps.med.x, y = eps.med.y, color = type.y))+
  geom_point(size = 3, alpha= 0.8) +
  labs(x = "Median Bliss score in the screen",
       y = "Median Bliss score in validation",
       title = "Effect size comparison (OD-based interactions)",
       color = "Benchmarking") +
  theme_classic() +
  scale_color_manual(values = c("#FFCC33", "#8c8c8c", "#009999")) +
  xlim(c(-0.4, 0.4))+
  ylim(c(-0.4,0.4))+
  geom_abline(slope = 1, linetype ="dotted")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(g, filename = paste0(figdir, "Eli-validation-efsize-OD-based.pdf"))


bench = rename(bench, Benchmarking =type.y) %>%
  mutate(Benchmarking = plyr::mapvalues(Benchmarking,
                                        from = c("synergy", "neutral", "antagonism"),
                                        to = c("S", "N", "A")))


get_roc <- function(strong_cutoff) {
  eli_roc = lapply(seq(0,1,length.out = 200), function(x) call_interactions(bench, multpTest_cutoff = x,
                                                                            strong_cutoff = strong_cutoff))
  roc_df_eli = plyr::ldply(eli_roc, data.frame)
  auc = with(roc_df_eli, simple_auc(tpr, fpr))
  roc_df_eli$thresh = paste0("medp = ", 
                             format(strong_cutoff, digits = 2),
                             ", AUCROC = ",
                             format(auc, digits = 3))
  roc_df_eli
}

roc_list = lapply(seq(0.01, 0.08, by=0.01), get_roc)
roc_df = plyr::ldply(roc_list)
roc_df = mutate(roc_df, youden = tpr - fpr)

roc_df$thresh[which.max(roc_df$youden)]
roc_df$fdr[which.max(roc_df$youden)]

point_sel = call_interactions(bench, multpTest_cutoff = 0.05,
                              strong_cutoff = 0.04) 
point_sel$thresh = "medp = 0.04, AUCROC = 0.791"


get_tpfp_weakcons <- function(bench, multpTest_cutoff,
                              strong_cutoff, weak_cutoff) {
  newman = read.table("~/Documents/embl/screen_K/wilcoxon/pvals-all.txt",
                      header = T, stringsAsFactors = F)
  newman$strain = "K"
  newman = mutate(newman, padj = p.adjust(pval, method = "BH"))
  newman =  mutate(newman, type = ifelse(padj <= multpTest_cutoff,
                                         ifelse(medpbliss <= -strong_cutoff, 'S', 
                                                ifelse(medpbliss >= strong_cutoff, 'A', 'N')), 'N'))
  
  dsm = read.table("~/Documents/embl/screen_M/wilcoxon/pvals-all.txt",
                   header = T, stringsAsFactors = F)
  dsm$strain = "M"
  dsm = mutate(dsm, padj = p.adjust(pval, method = "BH"))
  dsm =  mutate(dsm, type = ifelse(padj <= multpTest_cutoff,
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
  
  bench = dplyr::mutate(bench, type = ifelse(padj <= multpTest_cutoff,
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
weak_cons$thresh = "medp = 0.04, AUCROC = 0.791"

p = ggplot(roc_df, aes(x = fpr, y = tpr, color = thresh)) + 
  geom_line(alpha = 0.6, size=1.3) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "False positive rate",
       y = "True positive rate",
       color = "",
       title = "ROC curve, Gram-positive validation set (OD-based interactions)")+
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
                           nudge_x = -0.02,
                           segment.color = NA,
                           show.legend = F) +
  ggrepel::geom_text_repel(data = weak_cons,
                           aes(label = "Weak: 0.02"),
                           nudge_y = 0.02,
                           segment.color = NA,
                           show.legend = F) 
ggsave(p, filename = paste0(figdir, "ROC-Elisabetta-benchmarking-OD-based.pdf"),
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
ggsave(g, filename = paste0(figdir, "precrecall-Elisabetta-benchmarking-OD-based.pdf"),
       width = 9)

# roc_df  = filter(roc_df, thresh == "medp = 0.04, AUC = 0.791")
# roc_df$thresh[which.max(roc_df$f1score)]
# roc_df$fdr[which.max(roc_df$f1score)]

roc_OD_df = roc_df
# load AUC-based interactions

newman = read.table("~/Documents/embl/screen_K/wilcoxon/pvals-all-AUC.txt",
                    header = T, stringsAsFactors = F)
newman$strain = "K"
newman = mutate(newman, padj = p.adjust(pval, method = "BH"))

dsm = read.table("~/Documents/embl/screen_M/wilcoxon/pvals-all-AUC.txt",
                    header = T, stringsAsFactors = F)
dsm$strain = "M"
dsm = mutate(dsm, padj = p.adjust(pval, method = "BH"))

bsub = read.table("~/Documents/embl/screen_D/wilcoxon/pvals-all-AUC.txt",
                    header = T, stringsAsFactors = F)
bsub$strain = "D"
bsub = mutate(bsub, padj = p.adjust(pval, method = "BH"))


bench = rename(bench, comb = uniquecomb,
               eps.med = eps.med.y) %>%
  select(-c(eps.med.x, comb.x, type.x, pval, padj, noise, medpbliss, eps.mean))

bench = bind_rows(inner_join(newman, bench, by = c("comb", "strain")),
          inner_join(dsm, bench, by = c("comb", "strain")),
          inner_join(bsub, bench, by = c("comb", "strain"))
)

g = ggplot(bench, aes(x = eps.med.x, y = eps.med.y, color = Benchmarking))+
  geom_point(size = 3, alpha= 0.8) +
  labs(x = "Median Bliss score in the screen",
       y = "Median Bliss score in validation",
       title = "Effect size comparison (AUC-based interactions)",
       color = "Benchmarking") +
  theme_classic() +
  scale_color_manual(values = c("#FFCC33", "#8c8c8c", "#009999")) +
  xlim(c(-0.4, 0.4))+
  ylim(c(-0.4,0.4))+
  geom_abline(slope = 1, linetype ="dotted")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(g, filename = paste0(figdir, "Eli-validation-efsize-AUC-based.pdf"))


roc_list = lapply(seq(0.01, 0.08, by=0.01), get_roc)
roc_df = plyr::ldply(roc_list)
roc_df = mutate(roc_df, youden = tpr - fpr)

roc_df$thresh[which.max(roc_df$youden)]
roc_df$fdr[which.max(roc_df$youden)]


g = ggplot(roc_df, aes(x = fpr, y = tpr, color = thresh)) + 
  geom_line(alpha = 0.5) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "FPR", y = "TPR", title = "ROC curves in Gram-positive validation set")+
  labs(x = "False positive rate",
       y = "True positive rate",
       color = "",
       title = "ROC curve, Gram-positive validation set (AUC-based interactions)")+
  geom_abline(slope = 1, color = 'blue', linetype = "dotted")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_point(data = roc_df[which.max(roc_df$youden),], shape = 8, show.legend = F) +
  ggrepel::geom_text_repel(data = roc_df[which.max(roc_df$youden),],
                           aes(label = paste("FDR", format(fdr, digits = 2))),
                           nudge_y=0.03,
                           segment.color = NA,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "ROC-Elisabetta-benchmarking-AUC-based.pdf"),
       width = 9)


# make precision-recall plot
roc_df = group_by(roc_df, thresh) %>%
  mutate(prec = ifelse(prec == 0, 
                       prec[fpr == min( fpr[fpr!=min(fpr)] )],
                       prec),
         f1score = 2*prec*rec / (prec + rec))

roc_df$thresh[which.max(roc_df$f1score)]
roc_df$fdr[which.max(roc_df$f1score)]


g = ggplot(roc_df, aes(x = rec, y = prec, color = thresh))+
  geom_line(alpha = 0.5)+
  geom_point(alpha=0.5)+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "Recall",
       y = "Precision",
       color ="",
       title = "Precision-recall, Gram-positive validation set (AUC-based interactions)")  +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_point(data = roc_df[which.max(roc_df$f1score),], shape = 8, show.legend = F) +
  ggrepel::geom_text_repel(data = roc_df[which.max(roc_df$f1score),],
                           aes(label = paste("FDR", format(fdr, digits = 2))),
                           nudge_y = 0.02,
                           nudge_x = 0.03,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "precrecall-Elisabetta-benchmarking-AUC-based.pdf"),
       width = 9)

roc_AUC_df = roc_df

roc_OD_df$type = "OD-based"
roc_AUC_df$type = "AUC-based"

roc_df = bind_rows(roc_OD_df,
                   roc_AUC_df)

g = ggplot(roc_df, aes(x = fpr, y = tpr, color = type)) + 
  geom_point(alpha=0.2)+
  geom_smooth(alpha = 0.5, se = F) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "False positive rate",
       y = "True positive rate",
       color = "",
       title = "ROC curve, Gram-positive validation set")+
  geom_abline(slope = 1, color = 'blue', linetype = "dotted")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_point(data = mutate(weak_cons, type = "OD-based"), shape = 8, show.legend = F)+
  geom_point(data = mutate(point_sel, type = "OD-based"), shape = 8, show.legend = F)+
  ggrepel::geom_text_repel(data = mutate(point_sel, type = "OD-based"),
                           aes(label = paste("q <", format(fdr, digits = 2))),
                           segment.color = NA,
                           nudge_y = 0.02,
                           nudge_x = -0.01,
                           show.legend = F) +
  ggrepel::geom_text_repel(data = mutate(weak_cons, type = "OD-based"),
                           aes(label = "Weak: 0.02"),
                           segment.colour = NA,
                           nudge_y = 0.02,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "ROC-Elisabetta-benchmarking-AUC-vs-OD.pdf"),
       width = 9)



g = ggplot(roc_df, aes(x = rec, y = prec, color = type))+
  geom_smooth(se = F)+
  geom_point(alpha=0.5)+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "Recall",
       y = "Precision",
       color ="",
       title = "Precision-recall, Gram-positive validation set")  +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_point(data = mutate(weak_cons, type = "OD-based"), shape = 8, show.legend = F)+
  geom_point(data = mutate(point_sel, type = "OD-based"), shape = 8, show.legend = F)+
  ggrepel::geom_text_repel(data = mutate(point_sel, type = "OD-based"),
                           aes(label = paste("q <", format(fdr, digits = 2))),
                           segment.color = NA,
                           nudge_y = 0.02,
                           nudge_x = 0.04,
                           show.legend = F) +
  ggrepel::geom_text_repel(data = mutate(weak_cons, type = "OD-based"),
                           aes(label = "Weak: 0.02"),
                           segment.colour = NA,
                           nudge_y = -0.02,
                           nudge_x = 0.06,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "precrecall-Elisabetta-benchmarking-AUC-vs-OD.pdf"),
       width = 9)


# take this (optimised single-drug fitness)
roc_L2 = roc_OD_df

bench = read.csv(paste0(datadir, "bench_total_final.txt"),
                 sep = "\t")
bench = rename(bench, Benchmarking =type.y) %>%
  mutate(Benchmarking = plyr::mapvalues(Benchmarking,
                                        from = c("synergy", "neutral", "antagonism"),
                                        to = c("S", "N", "A")))

newman = read.table("~/Documents/embl/screen_K/wilcoxon/pvals-all-TSB.txt",
                    header = T, stringsAsFactors = F)
newman$strain = "K"
newman = mutate(newman, padj = p.adjust(pval, method = "BH"))

dsm = read.table("~/Documents/embl/screen_M/wilcoxon/pvals-all-TSB.txt",
                 header = T, stringsAsFactors = F)
dsm$strain = "M"
dsm = mutate(dsm, padj = p.adjust(pval, method = "BH"))

bsub = read.table("~/Documents/embl/screen_D/wilcoxon/pvals-all-TSB.txt",
                  header = T, stringsAsFactors = F)
bsub$strain = "D"
bsub = mutate(bsub, padj = p.adjust(pval, method = "BH"))


bench = rename(bench, comb = uniquecomb,
               eps.med = eps.med.y) %>%
  select(-c(eps.med.x, comb.x, type.x, pval, padj, noise, medpbliss, eps.mean))

bench = bind_rows(inner_join(newman, bench, by = c("comb", "strain")),
                  inner_join(dsm, bench, by = c("comb", "strain")),
                  inner_join(bsub, bench, by = c("comb", "strain"))
)

g = ggplot(bench, aes(x = eps.med.x, y = eps.med.y, color = Benchmarking))+
  geom_point(size = 3, alpha= 0.8) +
  labs(x = "Median Bliss score in the screen",
       y = "Median Bliss score in validation",
       title = "Effect size comparison (OD-based interactions with TSB/top single fitness)",
       color = "Benchmarking") +
  theme_classic() +
  scale_color_manual(values = c("#FFCC33", "#8c8c8c", "#009999")) +
  xlim(c(-0.4, 0.4))+
  ylim(c(-0.4,0.4))+
  geom_abline(slope = 1, linetype ="dotted")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(g, filename = paste0(figdir, "Eli-validation-efsize-OD-based-TSB-top.pdf"))


roc_list = lapply(seq(0.01, 0.08, by=0.01), get_roc)
roc_df = plyr::ldply(roc_list)
roc_df = mutate(roc_df, youden = tpr - fpr)

roc_df$thresh[which.max(roc_df$youden)]
roc_df$fdr[which.max(roc_df$youden)]

roc_df = group_by(roc_df, thresh) %>%
  mutate(prec = ifelse(prec == 0, 
                       prec[fpr == min( fpr[fpr!=min(fpr)] )],
                       prec),
         f1score = 2*prec*rec / (prec + rec))

roc_df$thresh[which.max(roc_df$f1score)]
roc_df$fdr[which.max(roc_df$f1score)]

g = ggplot(roc_df, aes(x = fpr, y = tpr, color = thresh)) + 
  geom_line(alpha = 0.5) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "FPR", y = "TPR", title = "ROC curves in Gram-positive validation set")+
  labs(x = "False positive rate",
       y = "True positive rate",
       color = "",
       title = "ROC curve, Gram-positive validation set (OD-based, TSB/top)")+
  geom_abline(slope = 1, color = 'blue', linetype = "dotted")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_point(data = roc_df[which.max(roc_df$youden),], shape = 8, show.legend = F) +
  ggrepel::geom_text_repel(data = roc_df[which.max(roc_df$youden),],
                           aes(label = paste("FDR", format(fdr, digits = 2))),
                           nudge_y=0.03,
                           segment.color = NA,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "ROC-Elisabetta-benchmarking-OD-based-TSB-top.pdf"),
       width = 9)

# make precision-recall plot
roc_df = group_by(roc_df, thresh) %>%
  mutate(prec = ifelse(prec == 0, 
                       prec[fpr == min( fpr[fpr!=min(fpr)] )],
                       prec),
         f1score = 2*prec*rec / (prec + rec))

roc_df$thresh[which.max(roc_df$f1score)]
roc_df$fdr[which.max(roc_df$f1score)]

g = ggplot(roc_df, aes(x = rec, y = prec, color = thresh))+
  geom_line(alpha = 0.5)+
  geom_point(alpha=0.5)+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "Recall",
       y = "Precision",
       color ="",
       title = "Precision-recall, Gram-positive validation set (OD-based, TSB/top)")  +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(g, filename = paste0(figdir, "precrecall-Elisabetta-benchmarking-OD-based-TSB-top.pdf"),
       width = 9)

roc_TSB = roc_df
roc_L2$type = "Optimized single drug fitness"
roc_TSB$type = "TSB / Top fitness"

roc_df = bind_rows(roc_L2,
                   roc_TSB)

g = ggplot(roc_df, aes(x = fpr, y = tpr, color = type)) + 
  geom_point(alpha=0.2)+
  geom_smooth(alpha = 0.5, se = F) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "False positive rate",
       y = "True positive rate",
       color = "",
       title = "Aggreagated ROC curves (OD-based interactions)")+
  geom_abline(slope = 1, color = 'blue', linetype = "dotted")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8, 0.2)) +
  geom_point(data = mutate(weak_cons, type = "Optimized single drug fitness"), shape = 8, show.legend = F)+
  geom_point(data = mutate(point_sel, type = "Optimized single drug fitness"), shape = 8, show.legend = F)+
  ggrepel::geom_text_repel(data = mutate(point_sel, type = "Optimized single drug fitness"),
                           aes(label = paste("q <", format(fdr, digits = 2))),
                           segment.color = NA,
                           nudge_y = 0.02,
                           nudge_x = -0.01,
                           show.legend = F) +
  ggrepel::geom_text_repel(data = mutate(weak_cons, type = "Optimized single drug fitness"),
                           aes(label = "Weak: 0.02"),
                           segment.colour = NA,
                           nudge_y = 0.02,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "ROC-Elisabetta-benchmarking-L2-vs-TSB-top.pdf"))

# precision-recall
g = ggplot(roc_df, aes(x = rec, y = prec, color = type))+
  geom_smooth(se = F)+
  geom_point(alpha=0.5)+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "Recall",
       y = "Precision",
       color ="",
       title = "Precision-recall, Gram-positive validation set")  +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_point(data = mutate(weak_cons, type = "Optimized single drug fitness"), shape = 8, show.legend = F)+
  geom_point(data = mutate(point_sel, type = "Optimized single drug fitness"), shape = 8, show.legend = F)+
  ggrepel::geom_text_repel(data = mutate(point_sel, type = "Optimized single drug fitness"),
                           aes(label = paste("q <", format(fdr, digits = 2))),
                           segment.color = NA,
                           nudge_y = 0.02,
                           nudge_x = 0.04,
                           show.legend = F) +
  ggrepel::geom_text_repel(data = mutate(weak_cons, type = "Optimized single drug fitness"),
                           aes(label = "Weak: 0.02"),
                           segment.colour = NA,
                           nudge_y = -0.02,
                           nudge_x = 0.06,
                           show.legend = F)
ggsave(g, filename = paste0(figdir, "precrecall-Elisabetta-benchmarking-L2-vs-TSB-top.pdf"),
       width = 9)






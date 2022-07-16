## Local FDR test for significance of interactions
## Sep 18, 2019
## V. Kim
require(graphics)
require(stats)
require(Hmisc)
require(zoo)
require(gplots)
require(LSD)
require(corrplot)
library(igraph)
require(plotrix)
require (purrr)
library(dplyr)
library(ggExtra)
library(locfdr)
argv = commandArgs(trailingOnly = TRUE)
# argv[1] = data/ directory
# argv[2] = strain

params = rjson::fromJSON(file = "params.json")
strain = paste0("screen", argv[2])

figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/wilcoxon/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/", strain, "/signif/")
}

wilc_dir = "~/Documents/embl/screen_M/wilcoxon/"
if(!dir.exists(wilc_dir)) {
  wilc_dir = file.path(argv[1], "wilcoxon/")
}

pval_strain = read.table(paste0(wilc_dir, "pvals-all.txt"),
                         header = T, stringsAsFactors = F)
pval_strain = dplyr::select(pval_strain, -c(pval))

# test per strain as each combination is repeated 3 times if 
# the summarised Bliss scores are considered all together
pdf(paste0(figdir, "lfdr-null-distrib.pdf"))
w_strain = locfdr(zz = pval_strain$medpbliss, bre = 200)
dev.off()

pval_strain$lfdr = w_strain$fdr

xleft = w_strain$mat[which.min(abs(w_strain$mat[,'Fdrleft'] - 0.1)), 'x']
xright = w_strain$mat[which.min(abs(w_strain$mat[,'Fdrright'] - 0.1)), 'x']
fdrleft = w_strain$mat[which.min(abs(w_strain$mat[,'Fdrleft'] - 0.1)), 'fdr']
fdrright = w_strain$mat[which.min(abs(w_strain$mat[,'Fdrright'] - 0.1)), 'fdr']

pval_strain$fdrleft = fdrleft
pval_strain$fdrright = fdrright

write.table(pval_strain,
            file = paste0(wilc_dir, "localfdr.txt"),
            sep = "\t", quote = F, row.names = F)

signif = filter(pval_strain, (lfdr < fdrleft & medpbliss < 0) | 
                  (lfdr < fdrright & medpbliss > 0))

p = ggplot(as.data.frame(w_strain$mat),
           aes(x = x,
               y = Fdrleft)) +
  geom_point(color = "#ffbb33") +
  geom_line(data = as.data.frame(w_strain$mat),
            aes(x = x, y = fdr)) +
  theme_bw() +
  geom_hline(yintercept = 0.1, color = "#0066cc", linetype = "dashed") +
  annotate("text", x = -0.2, y = 0.12, label = "FDR = 0.1") + 
  geom_vline(xintercept = xleft,
             linetype = "dotted") +
  annotate("text", x = -0.12, y = 0.24, label = paste("fdr =", 
                                                      format(fdrleft, digits = 2)))+
  geom_point(data = as.data.frame(w_strain$mat),
             aes(x = x,
                 y = Fdrright),
             color = "#595959") +
  geom_vline(xintercept = xright,
             linetype = "dotted") +
  labs(x = "Effect size",
       y = "Local (solid) and tail-area (points) FDR",
       title = paste(nrow(signif), "interactions in", argv[2])) + 
  annotate("text", x = 0.12, y = 0.33, label = paste("fdr =", 
                                                     format(fdrright, digits = 2)))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(p, filename = paste0(figdir, "locfdr-FDR.pdf"),
       width = 9, height = 7)


signif = dplyr::mutate(signif,
                       type = ifelse(medpbliss < 0, "synergy",
                                     ifelse(medpbliss > 0,"antagonism", "none")))

signif1 = mutate(signif, rec = gsub("(.+)(_.*)", "\\1", comb),
                 don = gsub("(.*_)(.+)", "\\2" , comb))

signif2 = mutate(signif, don = gsub("(.+)(_.*)", "\\1", comb),
                 rec = gsub("(.*_)(.+)", "\\2" , comb))

signif = bind_rows(signif1, signif2)


filler = expand.grid(don = signif$don, rec = signif$rec)
signif_hmap = right_join(signif, filler)

signif_hmap = distinct(signif_hmap, don, rec, type, .keep_all = T)

x.labels = sort(unique(signif_hmap$rec))
y.labels = sort(unique(signif_hmap$don))

inter_heatmap = ggplot(signif_hmap, 
                       aes(x = factor(rec, levels = x.labels),
                           y = factor(don, levels = y.labels),
                           fill = factor(type))) + 
  geom_tile(colour = "black") + theme_classic() + 
  scale_fill_manual(values = c("#C9788A", "#78A1A7"),
                    labels = c("Antagonism", "Synergy", "None"),
                    na.value = "#737373") + 
  labs(x = "Drug 1",
       y = "Drug 2",
       fill = "Interaction type",
       title = "Significant drug interactions") + 
  coord_equal() + 
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(vjust = 0.5, angle = 90),
        plot.title = element_text(hjust = 0.5))

ggsave(inter_heatmap,
       filename = paste0(figdir, "lfdr-drugheatmap-all.pdf"),
       width = 10, height = 9)



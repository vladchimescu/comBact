# Script for multiple testing correction
# significant hits are annotated using the drug classification
library(dplyr)
library(ggplot2)

argv = commandArgs(trailingOnly = TRUE)
# argv[1] = data/ directory
# argv[2] = strain
params = rjson::fromJSON(file = "params.json")
# effect size thresholds
strongsyn = params$strongsyn
strongant = params$strongant

strain = paste0("screen", argv[2])

wilc_dir = "~/Documents/embl/screen_K/wilcoxon/"
if(!dir.exists(wilc_dir)) {
  wilc_dir = file.path(argv[1], "wilcoxon/")
}

figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/wilcoxon/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/", strain, "/signif/")
  dir.create(figdir)
}

pvals =read.table(paste0(wilc_dir, "pvals-all.txt"),
                  sep = "\t",
                  header = T, stringsAsFactors = F)

# unfiltered p-value histogram
p_all = ggplot(pvals, aes(x = pval)) + 
  geom_histogram(bins = 50, fill = "#7EBEC5", colour = "black") + 
  theme_classic() + 
  labs(x = "P value (unadjusted)",
       y = "Count",
       title = "P value histogram of all interactions") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(p_all, filename = paste0(figdir, "phist_all.pdf"))


# the reason is quite trivial. Most of the interactions
# in this tail are "neutral" (neither synergistic nor antagonistic)

#--------------------STRONG ASSOCIATIONS BASED ON THE COMPLETE DISTRIBUTION------------
# get rid of them first
p_nonneut = ggplot(filter(pvals, abs(medpbliss) > 0.01), aes(x = pval)) + 
  geom_histogram(bins = 50, fill = "#F2CCC3", colour = "black") + 
  theme_classic() + 
  labs(x = "P value (unadjusted)",
       y = "Count",
       title = "P value histogram of non-neutral interactions") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(p_nonneut, filename = paste0(figdir, "phist_nonneutral.pdf"))

# pval.without.neutral = filter(pvals, synergistic | antagonistic) %>%
#   mutate(padj = p.adjust(pval, method = "BH"))


#=================SIGNIFICANT INTERACTIONS==============================
alpha = 0.05
# plot p-values against their rank
# ggplot(mutate(pvals, rank = rank(pval)), aes(x = rank, y = pval)) +
#   geom_point() +
#   geom_abline(slope = alpha / nrow(pvals)) +
#   xlim(c(0,300)) + ylim(c(0,0.01))
# pvals = mutate(pvals, padj = p.adjust(pval, method = "BH"))
signif = filter(pvals, padj <= alpha)

# here 'don' and 'conc' are simply 'Drug 1' and 'Drug 2'
# i.e. there is no inherent order implied by (donor_conc)
signif1 = mutate(signif, rec = gsub("(.+)(_.*)", "\\1", comb),
               don = gsub("(.*_)(.+)", "\\2" , comb))

signif2 = mutate(signif, don = gsub("(.+)(_.*)", "\\1", comb),
                 rec = gsub("(.*_)(.+)", "\\2" , comb))

signif = bind_rows(signif1, signif2)

# remove BUG donors and recipients
signif = filter(signif, rec != "BUG" & don != "BUG" & rec != "NOT" & don != "NOT")

signif = mutate(signif, 
                type = ifelse(efsize < -strongsyn,
                              "synergy",
                              ifelse(efsize > strongant,
                                     "antagonism", "none")))

## remove neutral combinations
signif = filter(signif, type != "none")

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
        axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5))

ggsave(inter_heatmap,
       filename = paste0(figdir, "drugheatmap-all.pdf"),
       width = 10, height = 9)


# write strong significant associations
write.table(distinct(signif, comb, .keep_all = T),
            file = paste0(wilc_dir, "strongassoc.tsv"), quote = F,
            sep = "\t", row.names = F)

# #----------------RESTRICTED ASSOCIATIONS------------------------
# p_nonneut = ggplot(filter(pvals, synrest | antrest), aes(x = pval)) + 
#   geom_histogram(bins = 50, fill = "#F2CCC3", colour = "black") + 
#   theme_classic() + 
#   labs(x = "P value (unadjusted)",
#        y = "Count",
#        title = "P value histogram of non-neutral interactions") + 
#   theme(plot.title = element_text(hjust = 0.5))
# 
# ggsave(p_nonneut, filename = paste0(figdir, "phist_rest_nonneutral.pdf"))
# 
# pval.without.neutral = filter(pvals, synrest | antrest) %>%
#   mutate(padj = p.adjust(pval, method = "BH"))
# 
# 
# alpha = 0.05
# signif.rest = filter(pval.without.neutral, padj <= alpha)
# 
# # here 'don' and 'conc' are simply 'Drug 1' and 'Drug 2'
# # i.e. there is no inherent order implied by (donor_conc)
# signif1 = mutate(signif.rest, rec = gsub("([A-Z]+)(_.*)", "\\1", comb),
#                  don = gsub("(.*_)([A-Z]+)", "\\2" , comb))
# 
# signif2 = mutate(signif.rest, don = gsub("([A-Z]+)(_.*)", "\\1", comb),
#                  rec = gsub("(.*_)([A-Z]+)", "\\2" , comb))
# 
# signif.rest = bind_rows(signif1, signif2)
# 
# # remove BUG donors and recipients
# signif.rest = filter(signif.rest, rec != "BUG" & don != "BUG" & rec != "NOT" & don != "NOT")
# 
# 
# signif.rest = mutate(signif.rest, 
#                 type = ifelse(synrest & !antrest,
#                               "synergy",
#                               ifelse(antrest & !synrest,
#                                      "antagonism", "none")))
# 
# 
# 
# rest_signif_inter = filter(signif.rest, comb %in% setdiff(signif.rest$comb, signif$comb))
# rest_signif_inter = distinct(rest_signif_inter, comb, .keep_all = T)
# stopifnot(nrow(rest_signif_inter) > 0)
# write.table(rest_signif_inter, 
#             file = paste0(wilc_dir, "restrictedsignif.tsv"), quote = F,
#             sep = "\t", row.names = F)

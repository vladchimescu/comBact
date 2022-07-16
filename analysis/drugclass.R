# Script for multiple testing correction
# significant hits are annotated using the drug classification

library(dplyr)
library(ggplot2)


## script for checkerboard assay plotting
argv = commandArgs(trailingOnly = TRUE)
# argv[1] = data/ directory
# argv[2] = strain


strain = paste0("screen", argv[2])

data_dir = "~/Documents/embl/screen_M"
if(!dir.exists(data_dir)) data_dir = argv[1]

script_dir = "~/Documents/embl/gitlab/drug_comb_screen/dataprep/"
if(!dir.exists(script_dir)) script_dir = "dataprep/"

figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/wilcoxon/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/", strain, "/signif/")
}

wilc_dir = "~/Documents/embl/screen_M/wilcoxon/"
if(!dir.exists(wilc_dir)) {
  wilc_dir = file.path(argv[1], "wilcoxon/")
}


# load the strongest associations
strong.inter = read.table(file = paste0(wilc_dir, "binomial-strongassoc.tsv"),
                          stringsAsFactors = F, header = T)

signif1 = mutate(strong.inter, rec = gsub("(.+)(_.*)", "\\1", comb),
                 don = gsub("(.*_)(.+)", "\\2" , comb))

signif2 = mutate(strong.inter, don = gsub("(.+)(_.*)", "\\1", comb),
                 rec = gsub("(.*_)(.+)", "\\2" , comb))

signif = bind_rows(signif1, signif2)

# load the drug annotation
drugannot = read.table(paste0(data_dir, "/classes_drugs-NT.txt"),
                       sep = "\t",
                       header = T, stringsAsFactors = F)
drugannot = drugannot[,-6]


signif = mutate(signif, recclass = ifelse(rec %in% drugannot$code,
                                         plyr::mapvalues(rec, from=drugannot$code,
                                                     to = drugannot$Master_Class),
                                         NA))

signif = mutate(signif, donclass = ifelse(don %in% drugannot$code,
                                          plyr::mapvalues(don, from=drugannot$code,
                                                          to = drugannot$Master_Class),
                                          NA))

signif = arrange(signif, recclass, donclass)

# set up the theme
t_heatmap =  theme(axis.ticks = element_blank(),
                   axis.text.y = element_text(size = 10, hjust=0),
                   axis.text.x = element_text(size = 10,
                                              margin = margin(t = 0),
                                              hjust = 1, angle = 45),
                   axis.title = element_text(size=10),
                   panel.background = element_blank(),
                   plot.title = element_text(hjust=0.5),
                   plot.background = element_rect(fill = "transparent",colour = NA),
                   legend.background = element_blank(),
                   legend.position = "top",
                   legend.title = element_blank(),
                   plot.margin = unit(rep(0,4), "cm"))

tmp = na.omit(signif) %>% group_by(donclass, recclass) %>%
  mutate(n = n(),
            eps.med = median(medpbliss))

lev = union(sort(unique(tmp$recclass)),
      sort(unique(tmp$donclass)) )


drugclass_int = ggplot(tmp, aes(x = factor(recclass, levels = lev),
                y = factor(donclass, levels = lev),
                colour = eps.med,
                size = n)) + 
  geom_point() + 
  guides(size = guide_legend(title = "Significant interactions"),
         colour = guide_colorbar(title = "Median effect size")) + 
  ylab("") + xlab("") + 
  ggtitle(paste("Interaction landscape of strain", argv[2])) + 
  coord_equal() + 
  scale_color_distiller(type = "div") + 
  scale_size(range = c(3,10)) + 
  t_heatmap + 
  theme(legend.title = element_text(),
        legend.background = element_blank(),
        legend.position = "left",
        legend.key = element_blank())

ggsave(drugclass_int,
       filename = paste0(figdir, "drugclass-interactions.pdf"),
       width = 7.5,
       height = 5.5)


# plot bar charts for synergies
drug_barpl = ggplot(signif,
       aes(x = factor(don, levels = rev(sort(unique(signif$don)))),
           fill = type)
) + 
  geom_bar() + 
  scale_fill_manual(values = c("#FFD23F", "#3CBBB1")) + 
  labs(x = "",
       y = "Count",
       title = paste("Significant interaction types by drug in strain", argv[2]),
       fill = "" ) + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  coord_flip() 
  

ggsave(drug_barpl,
       filename = paste0(figdir, "signif-inter-barplot.pdf"),
       width = 8,
       height = 12)


# split the dotplot for synergies and antagonisms
tmp = na.omit(signif) %>% 
  filter(type == "synergy") %>%
  group_by(donclass, recclass) %>%
  mutate(n = n(),
         eps.med = median(medpbliss))

lev = union(sort(unique(tmp$recclass)),
            sort(unique(tmp$donclass)) )


drugclass_syn = ggplot(tmp,
                       aes(x = factor(recclass, levels = lev),
                                y = factor(donclass, levels = rev(lev)),
                                colour = eps.med,
                                size = n)) + 
  geom_point() + 
  guides(size = guide_legend(title = "Significant synergistic interactions"),
         colour = guide_colorbar(title = "Median effect size")) + 
  ylab("") + xlab("") + 
  ggtitle(paste("Synergistic interactions of strain", argv[2])) + 
 # scale_y_discrete(labels=rev(labels)) + 
  #scale_x_discrete(labels=labels) + 
  coord_equal() + 
  scale_color_distiller(type = "seq", palette = "Greens") + 
  scale_size(range = c(3,10)) + 
  t_heatmap + 
  theme(legend.title = element_text(),
        legend.background = element_blank(),
        legend.position = "left",
        legend.key = element_blank())

ggsave(drugclass_syn,
       filename = paste0(figdir, "drugclass-synergies.pdf"),
       width = 7.5,
       height = 5.5)

levs = rev(sort(unique(na.omit(signif)$recclass)))

class_barpl = ggplot(na.omit(signif), 
                     aes(x = factor(recclass, levels = levs),
                   fill = factor(type))) + 
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("#FFD23F", "#3CBBB1")) + 
  labs(x = "",
       y = "Count",
       title = paste("Significant interaction types by drug class in strain", argv[2]),
       fill = "" ) + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  coord_flip()

ggsave(class_barpl,
       filename = paste0(figdir, "drugclass-barplot.pdf"),
       width = 7)

tmp = na.omit(signif) %>% 
  filter(type == "antagonism") %>%
  group_by(donclass, recclass) %>%
  mutate(n = n(),
         eps.med = median(medpbliss))

lev = union(sort(unique(tmp$recclass)),
            sort(unique(tmp$donclass)) )


drugclass_ant = ggplot(tmp,
                       aes(x = factor(recclass, levels = lev),
                           y = factor(donclass, levels = rev(lev)),
                           colour = eps.med,
                           size = n)) + 
  geom_point() + 
  guides(size = guide_legend(title = "Significant antagonistic interactions"),
         colour = guide_colorbar(title = "Median effect size")) + 
  ylab("") + xlab("") + 
  ggtitle(paste("Antagonistic interactions of strain", argv[2])) + 
  coord_equal() + 
  scale_color_distiller(type = "seq", palette = "Oranges") + 
  scale_size(range = c(3,10)) + 
  t_heatmap + 
  theme(legend.title = element_text(),
        legend.background = element_blank(),
        legend.position = "left",
        legend.key = element_blank())

ggsave(drugclass_ant,
       filename = paste0(figdir, "drugclass-antagonisms.pdf"),
       width = 7.5,
       height = 5.5)



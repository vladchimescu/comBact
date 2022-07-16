pacman::p_load(dplyr, reshape2, plyr, stringr, gtools, Rmisc, data.table, RColorBrewer, zoo, hexbin, tools, directlabels, grid, ggpubr, gridExtra,
               readr, ggplot2, stringi, ggrepel, viridis, scales, tidyr, sigr, lemon, here, smoothmest,EnvStats, scales)

#set working dir
home_dir = "/Users/Elisabetta/Documents/EMBL/combinatorials/screen_K/wilcoxon/"
setwd(home_dir)

#read p value file as computed from cluster
bug_interactions = read.table(paste0(home_dir,"P-vals_0.1.txt"), header = T, sep = "\t")

#exp_fit_threshold = 0.2
multpTest_cutoff = 0.05
interac_cutoff = 0.1
weak_interac_cutoff = 0.06
range = "both"

#correct for multiple testing
bug_interactions[2] = sapply(bug_interactions[2], FUN=p.adjust,method = "BH")
bug_interactions[4] = sapply(bug_interactions[4], FUN=p.adjust,method = "BH")
bug_interactions[6] = sapply(bug_interactions[6], FUN=p.adjust,method = "BH")

#subset significant p values
bug_inter_sign = subset (bug_interactions, Complete_p_val <= 0.05 | Antag_p_val <= 0.05 | Syn_p_val <= 0.05)

#is there any case where the complete p-value (not filtering in advance for syn or ant, i.e. exp fit < or > than real fit) is significant
#but neither of the syn and ant p values is significant?
check = subset (bug_interactions, Complete_p_val <= 0.05 & ( Antag_p_val >= 0.05 & Syn_p_val >= 0.05))
#only one case
rm(check)


#play with significant values before continuing with Ana's stuff
#merge back drug classes
bug_inter_sign$drug1 <- sapply(strsplit(as.character(bug_inter_sign$Combid), "_"), '[', 1)
bug_inter_sign$drug2 <- sapply(strsplit(as.character(bug_inter_sign$Combid), "_"), '[', 2)
#create drug legend
drug = unique(c(unique(bug_inter_sign$drug1), unique(bug_inter_sign$drug2)))
#load legend
legend = read.table("/Users/Elisabetta/Documents/EMBL/combinatorials/screen_K/interactions/drug_legend.txt", sep="\t", header = T)
legend$drug = legend$code
legend$Drug_full = legend$Drug
legend$Drug = NULL
legend$code = NULL

legend$drug1 = NULL
legend$drug2 = NULL

#how many drugs did I lose after p-value calculation? i.e. didnt lead to any significant interaction?
drug_beginning = unique(c(unique(for_wilc$donor), unique(for_wilc$Drug)))
setdiff(drug_beginning, drug)
#lost 11 drugs
rm(drug)
rm(drug_beginning)


#how many BUG controls are actually spotted as interactions??
check = subset (bug_inter_sign, drug1 == "BUG" | drug2 == "BUG")
nrow(check[check$Antag_p_val <= 0.05, ])
nrow(check[check$Syn_p_val <= 0.05, ])
nrow(check[check$Syn_p_val <= 0.05 & check$Antag_p_val <= 0.05, ])
weirdos = unique(c(unique(check$drug1), unique(check$drug2)))
rm(check)
rm(weirdos)

#check how many interactions are only significant for syn, antag or both (weird cases)
bug_inter_sign[is.na(bug_inter_sign)] <- 100

bug_inter_sign$asb <- ifelse(bug_inter_sign$Antag_p_val <= 0.05 & bug_inter_sign$Syn_p_val > 0.05, "A",
                      ifelse(bug_inter_sign$Syn_p_val <= 0.05 & bug_inter_sign$Antag_p_val > 0.05 , "S",
                        ifelse(bug_inter_sign$Antag_p_val <= 0.05 & bug_inter_sign$Syn_p_val <= 0.05, "B",
                               NA  )))[,1]

bug_inter_sign$Antag_p_val <- as.numeric(bug_inter_sign$Antag_p_val)
bug_inter_sign$Syn_p_val <- as.numeric(bug_inter_sign$Syn_p_val)
bug_inter_sign$Complete_p_val <- as.numeric(bug_inter_sign$Complete_p_val)

#get column with categories asb and their count
bug_inter_sign$value <- as.numeric(ave(bug_inter_sign$asb, bug_inter_sign$asb, FUN = length))
bug_inter_sign$labels <- paste0(round((bug_inter_sign$value / 108)  * 100, 1), "%")
bug_inter_sign$midpoint <- 108 - bug_inter_sign$value / 2


#plot this
pdf("~/Documents/EMBL/combinatorials/screen_K/interactions/pie_interaction_sign.pdf", width=10, height=10)

midpoint <- 108 - bug_inter_sign$value / 2

ggplot(bug_inter_sign, aes(x=factor(1), fill= factor(asb),label=asb), position = position_stack(vjust = 0.5)) +
  geom_bar(width=1) +
  coord_polar(theta = "y") +
  labs(x = "", y ="", title = "Interaction sign", fill = "Significant as antagonism, synergy or both") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour=guide_legend(title.position="top", 
                                        title.hjust =0.5, nrow=NULL)) +
  scale_y_continuous(breaks=midpoint, labels = bug_inter_sign$labels) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5), 
        panel.border = element_blank())
  
dev.off()



#check who are the "both"
both = subset (bug_inter_sign, asb == "B")

to_merge = cbind.data.frame (unique(c(unique(both$drug1), unique(both$drug2))))
colnames(to_merge) = "drug"
to_merge = inner_join(legend, to_merge, by = "drug")

pdf("~/Documents/EMBL/combinatorials/screen_K/interactions/twosided_interaction_classes.pdf", width=10, height=10)
ggplot(to_merge , aes(x = Master_Class, fill = factor(Master_Class))) +
  geom_bar() +
  labs(x = "Drug Class", y = "count", fill = "Drug Class") +
  labs(title = "Two-sided interactions - drug classes") +
  theme(axis.text.x = element_text(size = 12, angle=90, hjust=1),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 20, angle = 90),
      axis.title.x = element_text(size=20),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      plot.margin = unit(c(3, 3, 3, 3), "cm"),
      legend.title = element_text(colour="black", size=24,face="bold"),
      legend.text = element_text(colour="black", size=20),
      legend.key = element_rect(size = 5),
      legend.key.size = unit(2, 'lines')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

dev.off()


#look at their epsilon distribution
eps_df = cbind.data.frame(for_wilc$unique, for_wilc$eps)
colnames(eps_df) = c("Combid", "eps")

#to_merge = inner_join(eps_df, both, by = "Combid")
to_merge = inner_join(eps_df, bug_inter_sign, by = "Combid")
to_merge = subset(to_merge, Combid != "ASA_FEP")

hist(to_merge$eps)


to_merge$synorant <- ifelse(to_merge$Antag_p_val > to_merge$Syn_p_val , "antag",
                             ifelse(to_merge$Syn_p_val > to_merge$Antag_p_val , "syn",
                                    NA))


pdf("~/Documents/EMBL/combinatorials/screen_K/interactions/twosided_vs_onesided_epsilon_distribution.pdf", width=8, height=8)
ggplot(to_merge, aes(x= eps, fill = asb)) +
         geom_density(alpha=.3) +
  labs(x = "epsilon", y = "density", fill = "Sign of interaction", title = "Two-sided vs one-sided epsilon distribution") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
dev.off()

pdf("~/Documents/EMBL/combinatorials/screen_K/interactions/before_wilc_epsilon_distribution.pdf", width=7, height=7)
ggplot(for_wilc, aes(x= eps)) +
  geom_density(alpha=.3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
dev.off()

pdf("~/Documents/EMBL/combinatorials/screen_K/interactions/after_wilc_allepsilon_distribution.pdf", width=7, height=7)
ggplot(to_merge, aes(x= eps)) +
  geom_density(alpha=.3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
dev.off()

#write to merge -->Cytoscape
write.table(bug_inter_sign, "~/Documents/EMBL/combinatorials/screen_K/interactions/table.tsv",row.names=F,quote=FALSE,col.names = T, sep = "\t")


test = bug_inter_sign
test$drug = test$drug2


test = inner_join(legend, test, by = "drug")

write.table(test, "~/Documents/EMBL/combinatorials/screen_K/interactions/table_withclasses.tsv",row.names=F,quote=FALSE,col.names = T, sep = "\t")

bla = subset (to_merge, drug1 == "BUG")

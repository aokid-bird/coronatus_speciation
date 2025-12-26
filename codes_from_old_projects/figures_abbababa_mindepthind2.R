#===========================#
#### Libraries & Sources ####
#===========================#
library(tidyverse)
library(magrittr)
library(patchwork)
library(scales)
library(treeio)
library(ggtree)
library(khroma)
library(colorblindr)
source("code/dataprep.R")
source("../functions/ngs_functions.R")
source("../functions/cluster_functions.R")
source("../functions/general_functions.R")

#===================#
#### Directories ####
#===================#

resdir_angsd_abbababa.local <- "res/angsd_abbababa_mindepthind2"

lv <- c("Tomakomai", "Kushiro", "Russia", "Aomori")
pal <- color("bright")(2) %>% as.character()

df.dstat <- 
  read_csv(str_interp("${resdir_angsd_abbababa.local}/df.dstat.csv")) %>% 
  mutate(H1 = factor(H1, levels = lv),
         H2 = factor(H2, levels = lv),
         H3 = factor(H3, levels = lv)) %>% 
  arrange(H1, H2, H3) %>% 
  mutate(max = D+3*sd,
         signif = if_else(max < 0, "***", ""))

df.dstat %>% 
  filter(file == "Uncorrected") %>% 
  ggplot(data=.) +
  geom_hline(yintercept = 0.0, lty = 2)+
  geom_linerange(aes(x = pair, ymin = D-3*sd, ymax = D+3*sd), 
                 color = pal[1], lwd = 1)+
  geom_linerange(aes(x = pair, ymin = D-sd, ymax = D+sd), 
                 color = pal[2], lwd = 2)+
  geom_point(aes(x = pair, y = D), size = 3)+
  geom_text(aes(x = pair, y = D, label = signif), 
            nudge_x = 0.15, nudge_y = 0.007)+
  coord_flip()+
  theme_bw()+
  ylab("D-Statistics")+xlab("Population pair")+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))

figdir <- "figures/revise/abbababa"
mkdir(figdir)
ggsave(str_interp("${figdir}/Fig_abbababa.png"), width = 20, height = 20, units = "cm", dpi = 500)
ggsave(str_interp("${figdir}/Fig_abbababa.eps"), width = 20, height = 20, units = "cm")

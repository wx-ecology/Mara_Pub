## this script create a stacked bar figure of HMSC variance partitioning 
pal8 = c("#827580", "#e6dae4" , "#74e3c4", "#74c7e3","#429EBD", "#053F5C", "#F7AD19", "#ff7700")
#  ----------  read library ----------# 
library(ggplot2)
library(tidyverse)

df <- readRDS("./data/Hmsc_VP_R1.RDS")
df <- df$vals
df <- data.frame(variable = row.names(df), df) %>%
  pivot_longer(cols = 2:13, names_to = "species", values_to = "value") %>%
  mutate(
    variable =case_when((variable == "Random: month" | variable == "Random: plot") ~ "Random: spatiotemporal", TRUE ~ variable),
    variable = factor(variable, levels = c("Random: sample", "Random: spatiotemporal",
                                                "Forage quality","Site utilization", "Forage quantity",  "Precipitation", "Season", "Distance to border"
                                  #   "Forage-quality", "Forage-quantity", "Precipitation", "Season","Distance-to-border"
                                  )),
    species = factor(species, levels = colnames(df)),
    ) 


p <- ggplot(df, aes(fill = variable, y = value, x = species )) +
  geom_bar(position="fill", stat="identity") +
  bbplot::bbc_style()  + 
  scale_fill_manual(values = pal8) +
  theme (axis.text.x = element_text(size = 20),
         axis.text = element_text(size = 20),
         axis.title.y = element_blank(),
         legend.text = element_text(size = 20, color = "grey40"),
         legend.box.margin = margin(t = 20),
         legend.background = element_rect(
           color = "grey40", 
           fill = "grey95",
           size = .3
         ),
         strip.text = element_blank()) +
  guides(fill = guide_legend(nrow = 2, reverse = T)) +
  coord_flip()
ggsave("./figures/hmsc_vp2_R1.png", p,
        width = 20, height = 8, device = ragg::agg_png)

pal12 <- c( "#fdcce5", "#bd7ebe",  "#BF7E7E", "#ffb55a",     "#ffee65",     "#beb9db",   "#b2e061",   "#FEFEE3",    "#E9C09B",   "#7eb0d5",   "#8bd3c7" , "#fd7f6f" )
#         "Buffalo" ,  "Dikdik"  ,"Eland"  ,"Elephant", "Grants_Gazelle",  "Impala",  "Thompsons_Gazelle", "Topi",   "Waterbuck", "Wildebeest", "Zebra"  ,  "Cattle" 

library(tidyverse)
library(ggplot2)
library(lubridate)
library(bbplot)
### figure - waffle plot showing animal composition over time overlapped with precipitation ###

data.waffle <- read_rds("./data/mara_animal_compiled.rds") %>% 
  pivot_longer(cols = Cattle:Ostrich, names_to = "Species", values_to = "Count") %>%
  filter(!Species %in% c("Giraffe", "Sheep_Goats", "Ostrich")) %>%
  group_by(Yr_Mo, Species) %>%
  summarise(Count = sum(Count, na.rm = T)) %>%
  mutate(
    Species = factor(Species, levels = c( "Buffalo" ,  "Dikdik"  ,"Eland"  ,
                                          "Elephant", "Grants_Gazelle",  "Impala",  "Thompsons_Gazelle", 
                                          "Topi",   "Waterbuck", "Wildebeest", "Zebra",  "Cattle" )),
    Yr_Mo = ym (Yr_Mo))

##########################################################################################
############# stack bar plot + precipitation ######################################################
##########################################################################################
data.pr <- data %>% select(Yr_Mo, Precip) %>% distinct() %>% mutate(Yr_Mo = ym (Yr_Mo))
coeff.pr = 1/24
p_count_pr <-  ggplot() +
  geom_bar(data = data.waffle, aes(fill = Species, y = Count, x = Yr_Mo),
           color = "white", size = 0.1, position = "stack", stat="identity") +
    geom_line(data = data.pr, aes(x = Yr_Mo, y = pr / coeff.pr), 
              color = "grey40", linewidth = 1.5, alpha = 0.5, linetype="dashed") +
    scale_y_continuous(
      name = "Habitat usage (Dung count)",
      sec.axis = sec_axis(~.*coeff.pr, name="Precipitation (mm)")
    ) +
  bbplot::bbc_style()  +
  scale_fill_manual(values = pal12) +
  theme (axis.text.x = element_blank(),
         axis.text = element_text(size = 14),
         axis.title.y = element_text(size=16,
                                   color="#222222"),
         legend.text = element_text(size = 14, color = "grey40"),
         legend.box.margin = margin(t = 10),
         legend.background = element_rect(
           color = "grey40", 
           fill = "grey95",
           size = .3
         ),
         strip.text = element_blank()) +
  guides(fill = guide_legend(ncol = 6))
p_count_pr
ggsave("./figures/materials/count_precipitation.png", p_count_pr,
        width = 20, height = 8, device = ragg::agg_png)


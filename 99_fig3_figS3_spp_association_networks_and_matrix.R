## this script create a figure comparing raw correlation, HMSC derived spp covariance, abd ecoCopula derived spp network 

#  ----------  read library ----------# 
library(ggraph)
library(tidyverse)
library(igraph)
library(gridExtra)
library(circlize)
library(ggcorrplot)

# ---------- the function below is adapted from ecoCoupla package utils.R --------- #
## turn a correlation matrix into a graph object 
graph_from_cor <- function(cor){
  diag(cor)<-0
  vertex <- data.frame(name = colnames(cor))
  Edges=igraph::as_data_frame(igraph::graph_from_adjacency_matrix(cor,mode="undirected",weighted = TRUE))
  if(ncol(Edges)==3){
    colnames(Edges)[3]<- "cor"
  }
  igraph::graph_from_data_frame(Edges, directed=FALSE, vertices=vertex)
}

#  ----------  data prep ----------# 
raw_count <- read_rds("./data/mara_animal_compiled.rds") %>% 
  dplyr::select(Cattle, Wildebeest, Zebra, Thompsons_Gazelle, Impala, Topi,
                Eland, Buffalo, Grants_Gazelle, Waterbuck, Dikdik, Elephant) %>% replace(is.na(.), 0)



##########################################
## compile and plot correlation matrix ###
##########################################
raw_cor <- raw_count %>%
  as.matrix(.) %>% # 1,455 x 35. Some months/sites missing values. Dropped. 
  cor(.) 
hmsc_hi <- readRDS("./data/Hmsc_network_95_R1.RDS") 
copula_partcor <- readRDS("./data/ecoCopula_partcor_network.RDS") 

#########################################################
############# supp figure correlation matrix ############
#########################################################
a <- ggcorrplot(raw_cor, type = "lower", ggtheme = ggplot2::theme_minimal,  
                colors = c("#b2182b", "white","#2166ac"),
                lab = TRUE, digits = 2, lab_size = 3, show.legend = F, tl.cex = 10)
b <- ggcorrplot(hmsc_hi, type = "lower", ggtheme = ggplot2::theme_minimal,  
                colors = c("#b2182b", "white","#2166ac"),
                lab = TRUE, digits = 2, lab_size = 3, show.legend = F, tl.cex = 10)
c <- ggcorrplot(copula_partcor, type = "lower", ggtheme = ggplot2::theme_minimal, 
                colors = c("#b2182b", "white","#2166ac"), 
                lab = TRUE, digits = 2, lab_size = 3, show.legend = F, tl.cex = 10)

# ggsave("./figures/materials/mara_corrplot_2_R1.png", grid.arrange(a, b, c, nrow = 2),
#        width = 12, height = 10, device = ragg::agg_png)

##################################
########## plot networks #########
##################################

# turn to graph for plotting
raw_cor <- raw_cor %>% graph_from_cor(.)

copula_igraph_partcor <- copula_partcor  %>% graph_from_cor(.)

hmsc_igraph_hi <- hmsc_hi %>% graph_from_cor(.)  # posterior probability of the correlation estimate to be true is greater than 95%

# get node cooordinates so the two graphs will have the same node locations
seed = 1
Coords <- layout_with_kk(raw_cor) %>% 
  as_tibble %>%
  bind_cols(tibble(names = names(V(copula_igraph_partcor))))
# name the lay out better
Coords <- Coords %>% 
  mutate(V1 = case_when(names == "Cattle" ~ V1 - 1.7, 
                        names == "Buffalo" ~ V1 - 0.3, 
                        names == "Zebra" ~ V1 + 0.2, 
                        names == "Elephant" ~ V1 + 0.5, 
                        TRUE ~ V1),
         V2 = case_when(names == "Cattle" ~ V2 - 1, 
                        names == "Buffalo" ~ V2 + 0.2,
                        names == "Elephant" ~ V2 + 0.3, 
                        TRUE ~ V2))

a <- raw_cor %>%  
  ggraph(layout = as.matrix(Coords[,1:2])) + # see ?layout_tbl_graph_igraph
  geom_edge_fan0(aes( colour = cor, width= cor, alpha = 0.9)) +
  scale_edge_width(range = c(1, 4))+
  scale_edge_color_gradient2(low="#b2182b",mid="white",high="#2166ac")+
  geom_node_label(aes(label=name), repel = TRUE, nudge_y = 0)+
  geom_node_point(aes( size= 2))+
  theme_void() +
  labs(subtitle = "Raw correlation") +
  theme(
    legend.position = 'none',
    plot.title = element_text(color="black", size=14, face="bold.italic", vjust = 5, hjust = 0.5 ),
    plot.margin = margin(1, 1.5, 1, 1, "cm")
  )

b <- hmsc_igraph_hi %>%  
  ggraph(layout = as.matrix(Coords[,1:2])) + # see ?layout_tbl_graph_igraph
  geom_edge_fan0(aes( colour = cor, width= cor, alpha = 0.8)) +
  scale_edge_width(range = c(1, 4))+
  scale_edge_color_gradient2(low="#b2182b",mid="white",high="#2166ac")+
  geom_node_label(aes(label=name), repel = TRUE, nudge_y = 0)+
  geom_node_point(aes( size= 2))+
  theme_void() +
  labs(subtitle = "Hierarchical Modeling of Species Communities") +
  theme(
    legend.position = 'none',
    plot.title = element_text(color="black", size=14, face="bold.italic", vjust = 5, hjust = 0.5 ),
    plot.margin = margin(1, 1, 1, 1.5, "cm")
  )

c <- copula_igraph_partcor %>%  
  ggraph(layout = as.matrix(Coords[,1:2])) + # see ?layout_tbl_graph_igraph
  geom_edge_fan0(aes( colour = cor, width= cor, alpha = 0.8)) +
  scale_edge_width(range = c(1, 4))+
  scale_edge_color_gradient2(low="#b2182b",mid="white",high="#2166ac")+
  geom_node_label(aes(label=name), repel = TRUE, nudge_y = 0)+
  geom_node_point(aes( size= 2))+
  theme_void() +
  labs(subtitle = "Gaussian copula graphical models") +
  theme(
    legend.position = 'none',
    plot.title = element_text(color="black", size=14, face="bold.italic", vjust = 5, hjust = 0.5 ),
    plot.margin = margin(1, 1.5, 1, 1, "cm")
  )

# ggsave("./figures/materials/mara_graph_R1.png", grid.arrange(a, b, c, nrow = 1),
#        width = 13, height = 6, device = ragg::agg_png)

##########################################
########### descriptive stats ############
##################################ß########
# any species with negative association with cattle? 
raw_cor_tb <- raw_cor %>% as_tibble(.) %>% 
  mutate(spp2 = rownames(raw_cor)) %>% select(Cattle, spp2)
raw_cor_tb %>% filter(Cattle < 0)

hmsc_hi_tb <- hmsc_hi %>% as_tibble(.) %>% 
  mutate(spp2 = rownames(hmsc_hi)) %>% select(Cattle, spp2)
hmsc_hi_tb %>% filter(Cattle < 0)

copula_partcor_tb <- copula_partcor %>% ßas_tibble(.) %>% 
  mutate(spp2 = rownames(copula_partcor)) %>% select(Cattle, spp2)
copula_partcor_tb %>% filter(Cattle < 0)

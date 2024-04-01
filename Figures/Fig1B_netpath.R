# load packages
library(igraph)
library(intergraph)
library(dichromat)
library(scales)
library(GGally)
library(ggplot2)
library(ggrepel)
library(ggnetwork)

# read in a matrix of rg values that were significant (padj < 0.05)
# this is a modified version of the correlation table in Supplementary Table 2, except all correlations where
# padj > 0.05 are set to NA
SM<-as.matrix(read.table("matrix_p05_rg30.csv", header=TRUE, row.names=1, sep=","))

# duplicate the matrix
AM = SM

# make the diagonal equal to 0
diag(AM) = 0

# define a network graph based on the duplicated matrix
g = graph_from_adjacency_matrix(AM, 
                                mode = "lower",
                                weighted = TRUE)

# scale edges by rg
edge_attr(g)$label<-edge_attr(g)$weight
edge_attr(g)$scaled<-rescale(edge_attr(g)$label)
edge_attr(g)$weight<-abs(edge_attr(g)$weight)

# define colour palette
colours<-colorRamp(c("grey80", "darkslateblue"))

# plot the network graph
# the graph is generated anew each time the code is run. It is possible to run this
# chunk of code multiple times until a suitable graph is generated
net=g
plot<-ggnet2(net,
       size=10,
       color="grey90",
       label=TRUE,
       label.size=3,
       edge.size=edge_attr(g)$weight*9,
       edge.alpha=0.9,
       edge.color=rgb(colours(edge_attr(g)$scaled),maxColorValue=256),
       edge.label=as.character(edge_attr(g)$label),
       edge.label.fill=NA,
       edge.label.size=4
      )+
  geom_point(shape=21, fill="grey90", color="black", size=12, stroke=1)+
  geom_text(aes(label=label), size=3.5, fontface="bold")

plot

# save as .pdf
ggsave("netpath.pdf", plot)

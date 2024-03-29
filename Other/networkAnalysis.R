---
  title: "Organoid Network Analysis"
output: html_notebook
---
  
  The code here is an exploratory network analysis on the functional connectivity data (based on correlation of spike trains) obtained from organoids for The Organoid Project 2018. 



# Setting up some pre-requisites

Loading the library 

```{r}
library(igraph)
library(R.matlab)
library(ggplot2)
```

Loading the data 

```{r}
network <- readMat('D:/MECP2_2019_AD/Scripts_and_Output/S1.2.File_Conversion_Output/190830_slice1stim8_current_cSpikes_L-0.1254')

# convert from a list to a matrix 

network <- matrix(unlist(network), ncol = 60, byrow = TRUE)

```

Some extra functions for convinience

Function for calculating the node strength distribution
```{r}
# From here: https://stackoverflow.com/questions/8344565/how-do-i-calculate-weighted-degree-distributions-with-igraph-in-r
graph.strength.distribution <- function (graph, cumulative = FALSE, ...)
{
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  # graph.strength() instead of degree()
  cs <- graph.strength(graph, ...)
  hi <- hist(cs, -1:max(cs), plot = FALSE)$density
  if (!cumulative) {
    res <- hi
  }
  else {
    res <- rev(cumsum(rev(hi)))
  }
  res
}
```


# Initial visualisation 

Converting our weighted adjacency matrix to an `igraph` object 

```{r}

network <- graph_from_adjacency_matrix(network, mode = "undirected", weighted = TRUE)

```

Check some of the attributes of the `network` object

- note that the weights are stored in `E(network)$weight`

```{r}
graph_attr(network)

edge_attr(network)

vertex_attr(network)
```


Visualise the network (very rough overview)

```{r}
plot(network)
```

Let us also look at the adjacency matrix 

```{r}
netm <- get.adjacency(network, attr = "weight", sparse = F)

palf <- colorRampPalette(c("gold", "dark orange"))

# Create figure of a specific size so the plot does not look compressed 
dev.new(width = 20, height = 20)

heatmap(netm[, 60:1], Rowv = NA, Colv = NA, col = palf(100), 
        scale = "none", margins = c(2, 2))
# need to add tick marks 
# need to reduce the axis labels to 1, 10, 20... 60 
# need to create colourmap

```


Make some sensible adjustments so that the network plot provides more information in a more clear way, in particular: 
  
  - What are the main nodes; which nodes are more connected to others? 
  - We (at least I) am not interested in self-connections, so lets remove that 
- Lets also plot it a way that is easiest to look at

```{r}

# Set node size based on the number of degree (note that this discounts the weigtedness)
# V(network)$size = degree(network) # nope, doens't work yet, the network is fully connected according to them
# need to convert NA values to 0


# use a sensible layout (that is not grid), since this is for the purpose of visualising the network relation 
# rather than the spatial relationship 

l <- layout_with_fr(net.bg)

plot(network, layout = l)
```


Let's set some cut off point and plot again

```{r}

# convert NA values ot zero 
E(network)$weight[is.na(E(network)$weight)] <- 0 

# remove self loops 
network.pruned <- simplify(network, remove.loops = TRUE)

# create cut off point, and remove connections lower than the cut off
# cut.off <- mean(E(network.pruned)$weight)
cut.off <- 0.80
network.pruned <- delete_edges(network.pruned, E(network.pruned)[weight < cut.off])

# remove vertices without edges 

# network.pruned <- delete_vertices(network.pruned, degree(network.pruned) == 0)

# set node size according to degree
# V(network.pruned)$size = degree(network.pruned)

plot(network.pruned, layout = l)

```
## Visualise network in a grid (to match MEA plot)

- also need the edge weight to be proportional to the edge thickness
- layout needs to be a 8 x 8 grid with the corners removed


```{r}

yTemp = 1:8; 
# yCoord = repmat(fliplr(yTemp), 1, 8); # matlab code to be translated
yCoord = rep(rev(yTemp), 8);

xTemp = 1:8;
xCoord = rep(xTemp, each = 8);

coord = matrix(c(xCoord, yCoord), nrow = 64, ncol = 2)

# remove the 4 corners 
coord <- coord[!((coord[, 1] == 1) & (coord[, 2] == 1)), ]
coord <- coord[!((coord[, 1] == 1) & (coord[, 2] == 8)), ]
coord <- coord[!((coord[, 1] == 8) & (coord[, 2] == 1)), ]
coord <- coord[!((coord[, 1] == 8) & (coord[, 2] == 8)), ]


# make edge thickness proportional to edge weight 
# E(network.pruned)$ width <- 1
E(network.pruned)$width <- 1 + E(network)$weight * 3 

# Set node size according to node strength 
# V(network.pruned)$size = 1 + strength(network.pruned) * 1

# Set node size according to node degree
V(network.pruned)$size = 1 + igraph::degree(network.pruned) # it may be masked by the sna package

# remove vertex labels
V(network.pruned)$label <- NA

# Chane border colour 
V(network.pruned)$frame.color = "light blue" # white, black, or grey (or NA)

# save as eps 
postscript("/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/figures/paper_figures_first_draft/figure_4/igraph_network_plots/networkplot_threshold0.8_degree1_noBorder_edgeWeight_veretexDegree_greyborder.eps", width = 500, height = 500)

# png("/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/figures/paper_figures_first_draft/figure_4/igraph_network_plots/networkplot_threshold0.8_degree1_noBorder_edgeWeight_veretexDegree.png", width = 30, height = 20, units = 'cm', res = 300)

# plot(network.pruned, layout = layout_on_grid)
plot(network.pruned, layout = coord)

# legend for edges 

legend(x = 1.2, y = 0.2, c("0.8", "0.9", "1"), col=c("grey","grey", "grey"), lwd = c(1, 2.25, 4), box.col = "white", title = "Correlation")

# Legend for vertices : Node strength
# legend(x = 1.2, y = -0.25, c("0.3", "0.6", "1"), col = c("orange", "orange", "orange"), pch = 19, box.col = "white", title = "Normalized strength", pt.cex = c(1.5, 3, 5), y.intersp = 1.85, x.intersp = 1.5)
# pch = 19 means filled circle

# Legend for vertices: Node degree (threshold = 0.75)
# legend(x = 1.2, y = -0.25, c("5", "10", "20"), col = c("black", "black", "black"), pch = 21, pt.bg = "orange", box.col = "white", title = "Node degree", pt.cex = c(1.25, 2.5, 5), y.intersp = 1.85, x.intersp = 1.5)

# Legend for vertices: Node degree (threshold = 0.80)
legend(x = 1.2, y = -0.25, c("1", "5", "10"), col = c("black", "black", "black"), pch = 21, pt.bg = "orange", box.col = "white", title = "Node degree", pt.cex = c(0.75, 2, 4), y.intersp = 1.85, x.intersp = 1.5)

dev.off()



# library(GGally)
# library(ggnetwork)

# ggnet2(network, edge.size = "weight")


```

Attempt to make a legend

```{r}
# plot(network.pruned, layout = coord)
library(ggnetwork)
# save as eps 
postscript("/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/figures/paper_figures_first_draft/figure_4/igraph_network_plots/networkplot_legend_attempt1.eps", width = 500, height = 500)

ggplot(data = ggnetwork(network.pruned, layout = coord),
aes(x, y, xend = xend, yend = yend)) +
geom_edges(aes(size = width), color = "grey") +
geom_nodes(aes(size = size), color = "orange") +
scale_size_continuous(name = " ", 
breaks = c(2, 5, 7, 2, 5, 7),
guide = guide_legend(override.aes = list(
shape = c(rep(16, 3), c(rep(NA, 3))),
linetype = c(rep(0, 3), c(rep(1, 3)))
)
)
) + theme_blank()

dev.off()

```



# Some basic summary statistics 

## Node strength distribution 

Node strength is the sum of the weighted edges

```{r}
hist(strength(network, loops = FALSE))
```

## Edge statistics

Edge density

```{r}
edge_density(network.pruned, loops = F)
```

Edge weight distribution 

```{r}
hist(E(network.pruned)$weight) # this gives the ferquency 

# lets also plot the proportion 

h = hist(E(network.pruned)$weight)
h$density = h$counts / sum(h$counts) * 100 
plot(h, freq = FALSE, ylab = "Proportion (%)")

# A cumulative distribution plot 
plot(h$breaks[2:length(h$breaks)], cumsum(h$density))
```



# Hubs and authorities 

```{r}
hs <- hub_score(network.pruned)$vector # note that weights are used by default

l <- layout_in_circle(network.pruned)
plot(network.pruned, vertex.size = hs * 10, main = "Hubs", layout=layout_with_lgl)

```

# Centrality analysis

Closeness
```{r}
closeness(network.pruned)
```

Eigencentrality 

```{r}
eigen_centrality(network.pruned, directed = F)
```



# Subgroups and communities 


# Rich club topology 

`brainGraph` provides a function to claculate the normalised rich-club coefficient

```{r}
rich_club_norm(network.pruned, N = 100)
```


# Interactive plot 

```{r}
tkid <- tkplot(network.pruned, vertex.size = hs * 10, vertex.color = "orange", layout = layout_with_lgl) 
```



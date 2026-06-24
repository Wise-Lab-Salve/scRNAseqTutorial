## Validation

### scType
In order to validate or check your annotations, there are publically available tools. We will show you how to use scType to confirm cell types.

### Setup scType Tool

scType needs quite a few additional packages, you can install and open as below.

```
packages <- c(
  "dplyr", "Seurat", "HGNChelper", "openxlsx",
  "ggraph", "igraph", "tidyverse", "data.tree"
)

install.packages(setdiff(packages, rownames(installed.packages())))

# Load packages
lapply(packages, library, character.only = TRUE)

install.packages("remotes")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R")
```

## scType on their github page has R files that can help with the analysis. We use the source command to load these premade scripts into our R

```
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
source("https://raw.githubusercontent.com/kris-nader/sc-type/master/R/sctype_wrapper.R"); 
```
There is no documentation for each command used below, instead the developer has made a tutorial. 

[Click here to visit the documentation for scType](https://github.com/IanevskiAleksandr/sc-type/blob/master/README.md) 


## Provide a tissue type your data belongs to
e.g. Immune system, Pancreas, Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
```
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Immune system" 
gs_list <- gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") 

```

# Compare new identities with scType

We will rename the object so that we do not edit our current object with proper labels. 
```
seurat_object <- so_v1
```

Prepare Gene Sets from scType with the `gene_sets_prepare` command.
```
gs_list <- gene_sets_prepare(db_, tissue)
```

Check that our R data is v5
```
seurat_package <- isFALSE('counts' %in% names(attributes(seurat_object[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package, "v5", "v4")))
```

Extract scaled scRNA-seq counts data as it will need to switch format into a matrix
```
scRNAseqData_scaled <- if (seurat_package) as.matrix(seurat_object[["RNA"]]$scale.data) else as.matrix(seurat_object[["RNA"]]@scale.data)
```
Run scType with the `sctype_score` command.
```
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
```

Now that scType was run, we need to compare the cell labels it created to our clusters and labels.<br>
scType provides a rather long piece of code to extract the clusters and merge with its output.<br> 
We grabbed this code from the scType website https://github.com/IanevskiAleksandr/sc-type and copied it here. 

```
CL_Results <- do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl){
       es.max.cl = sort(rowSums(es.max[,rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
       head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_object@meta.data$seurat_clusters==cl)), 10)
   }))
sctype_scores <- CL_Results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])
```

Overlay the identified cell types on UMAP plot
```
so_sctype <- run_sctype(seurat_object,known_tissue_type="Immune system",custom_marker_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",name="sctype_classification",plot=TRUE)
```
Continue overlay and comparison of clusters
```
CL_Results <- CL_Results[order(CL_Results$cluster),]; edges = CL_Results; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL
nodes_lvl1 <- sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
```
Update colors for each cluster
```
ccolss <-  c("aquamarine3", "darkseagreen2", "cornflowerblue", "pink2", "goldenrod1", "aquamarine", "skyblue1","lightgoldenrod1", "thistle2","palevioletred2","azure3","orchid1","navajowhite2","salmon1","darkslateblue","lightcyan3")
```
## Graphing
Now that the clusters are compared and we have colors that we would like to use, we will make a graph from scType commands.
```
for (i in 1:length(unique(CL_Results$cluster))){
  dt_tmp = CL_Results[CL_Results$cluster == unique(CL_Results$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes <- rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db <- openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)

nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

nodes_v2 <- nodes %>% distinct() ###Can not have cluster column of nodes duplicated
nodes_v2 <- nodes %>% distinct(cluster, .keep_all = T)

mygraph <- graph_from_data_frame(edges, vertices=nodes_v2)

gggr <- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("black"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
```

Use this code for the same graph but without cluster labels on each bubble
``` 
gggr <- ggraph(mygraph, layout = 'circlepack', weight = I(ncells)) + 
  geom_node_circle(aes(filter = ord == 1, fill = I("#F5F5F5"), colour = I("#D3D3D3")), alpha = 0.9) + 
  geom_node_circle(aes(filter = ord == 2, fill = I(Colour), colour = I("#D3D3D3")), alpha = 0.9) +
  geom_node_text(aes(
    filter = ord == 2,
    label = shortName,
    colour = I("black"),
    repel = !1,
    parse = TRUE,
    size = I(log(ncells, 25) * 1.5)
  )) +
  theme_void()
```

Saving gggr `+` umap to your files
```
pdf(file = paste(plot.dir,"ScTypeGraphwithclusters_v1.pdf", sep = ""), height=20, width=10)
DimPlot(seurat_object, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss)+ gggr
dev.off()
```

Finally, we can open this graph and compare the output.

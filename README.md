Introduction:
-------------

Computes Shannon entropy for gene expression across cell types or
conditions to find cell type or condition specific gene expression
patterns genomic regions.

Installation:
-------------

Install the package from Github.

     devtools::install_github("nashchem/Shantropy")

Quick start:
------------

Load the package into R session.

     library(Shantropy)

The enrichmotifpairR package includes an example dataset called
example\_data.

    data(example_data)

Compute cell type/condition specific expression.

    cell_type_data <- compute_ShannonEntropy(data = data, exp_thr = 1, 
                                             cell_specificity_thr = 1)

Plot the cell type/condition specific expression data. Install/load the
required packages for visualization.

    library(pheatmap)
    library(RColorBrewer)

Plotting version 1

    data2 <- cell_type_data[, 2:4]
    pheatmap::pheatmap(log2(data2), cluster_rows = FALSE, cluster_cols = FALSE, 
                       scale = "none", show_rownames = FALSE, show_colnames = TRUE, 
                       color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11,
                       name = "RdYlBu")))(50))
    dev.copy2pdf(file="Shannon_entropy_based_cell_type_expression_heatmap_v1.pdf", 
                 width = 3, height = 4)

    path_to_fig <- '/Users/njayavel/Downloads/Shantropy/data/Shannon_entropy_based_cell_type_expression_heatmap_v1.png'
    knitr::include_graphics(path_to_fig)

<img src="/Users/njayavel/Downloads/Shantropy/data/Shannon_entropy_based_cell_type_expression_heatmap_v1.png" width="900" />

Plotting version 2

    mycolor <- colorRampPalette(c("#4462DB", "Yellow", "Red"))(100)
    data2 <- cell_type_data[, 2:4]
    pheatmap::pheatmap(log2(data2 + 1), cluster_rows = FALSE, cluster_cols = FALSE, 
                       scale = "none", show_rownames = FALSE, show_colnames = TRUE, 
                       color = mycolor)
    dev.copy2pdf(file="Shannon_entropy_based_cell_type_expression_heatmap_v2.pdf", 
                 width = 3, height = 3)

    knitr::include_graphics('/Users/njayavel/Downloads/Shantropy/data/Shannon_entropy_based_cell_type_expression_heatmap_v2.png')

<img src="/Users/njayavel/Downloads/Shantropy/data/Shannon_entropy_based_cell_type_expression_heatmap_v2.png" width="900" />

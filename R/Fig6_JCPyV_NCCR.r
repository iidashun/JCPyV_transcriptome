# 2024.09.12 created by Shun Iida
# This program visualizes the noncoding control region (NCCR) of the JCPyV genome


########## load library ##########
library(ggplot2)
library(gggenes)


########## function ##########
# Function to create NCCR plot
create_NCCR_plot <- function(data, output_file) {
  plot <- ggplot(data,
                 aes(xmin = start, xmax = end, y = molecule, fill = col)) +
    
    geom_segment(data = subset(data, type == 'spacer'),
                 aes(x = start, xend = end, y = molecule, yend = molecule),
                 size = 0, color = "#FFFFFF") +
    
    geom_segment(data = subset(data, type == 'lineg'),
                 aes(x = start, xend = end, y = molecule, yend = molecule),
                 size = 10, color = "#B3DE69") +
    
    geom_segment(data = subset(data, type == 'liner'),
                 aes(x = start, xend = end, y = molecule, yend = molecule),
                 size = 10, color = "#FB8072") +
    
    geom_segment(data = subset(data, type == 'lined'),
                 aes(x = start, xend = end, y = molecule, yend = molecule),
                 size = 1.5, linetype = "dotted", color = "#000000") +
    
    geom_segment(data = subset(data, type == 'linev'),
                 aes(x = start, xend = end, y = moleculeup, yend = moleculedown),
                 size = 1.5, linetype = "dotted", color = "#000000") +
    
    geom_segment(data = subset(data, type == 'mutation'),
                 aes(x = start, xend = end, y = molecule, yend = molecule),
                 size = 10, color = "#0000FF") +

    geom_gene_arrow(data = subset(data, type == 'box'),
                    aes(xmin = start, xmax = end, y = molecule, fill = col),
                    arrowhead_width = grid::unit(0, "mm"),
                    arrowhead_height = grid::unit(0, "mm"),
                    arrow_body_height = grid::unit(9, "mm")) +

    geom_gene_label(data = subset(data, type == 'boxlabel'),
                    aes(xmin = start, xmax = end, y = molecule, label = gene),
                    color = "#000000",
                    align = "centre",
                    height = grid::unit(9, "mm"),
                    grow = TRUE,
                    reflow = TRUE,
                    fontface = "bold") +
    
    geom_gene_label(data = subset(data, type == 'positionlabel'),
                    aes(xmin = start, xmax = end, y = molecule, label = gene),
                    color = "#000000",
                    align = "centre",
                    height = grid::unit(9, "mm"),
                    grow = TRUE,
                    reflow = TRUE) +

    geom_gene_label(data = subset(data, type == 'name'),
                    aes(xmin = start, xmax = end, y = molecule, label = gene),
                    color = "#000000",
                    align = "left",
                    height = grid::unit(9, "mm"),
                    grow = TRUE,
                    reflow = TRUE,
                    fontface = "bold") +

    theme_minimal() +
    scale_fill_identity() +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )
  
  ggsave(output_file, plot = plot, height = 25, width = 15)
}


########## plot ##########
JCPyV_NCCR <- read.csv("inputs/JCPyV/csv/JCPyV_NCCR.csv", header = TRUE, sep = ",")
create_NCCR_plot(JCPyV_NCCR, 'outputs/figs/fig6/JCPyV_NCCR.tiff')

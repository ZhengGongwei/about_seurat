library(ggplot2)
library(ggsci)
library(RColorBrewer)

plot_vln_Corgenes <- function(seurat_obj,
                         assay,
                         slot,
                         var_group,
                         vec_features,
                         vec_group_colors=NULL,
                         f_color = colorRampPalette(brewer.pal(n=11, name="RdYlBu")),
                         flip = T,
                         do_plot = F,
                         pt.size = 0,
                         feature_fontface = "bold.italic",
                         fontsize_axistext_x=12,
                         fontsize_axistext_y=12,
                         aspect.ratio =NULL,
                         figurenames = NULL,
                         width = 7,
                         height = 7
) {

  #=============prepare group and colors==================
  seurat_obj_tmp = seurat_obj
  Idents(seurat_obj_tmp) <- var_group
  levels(x = seurat_obj_tmp) = sort(unique(seurat_obj_tmp@meta.data[[var_group]]), decreasing = if (flip) T else F)

  if (is.null(vec_group_colors)) {
    n_group <- length(levels(x = seurat_obj_tmp))
    vec_group_colors <- f_color(n_group)
    names(vec_group_colors) <- levels(x = seurat_obj_tmp)
  }

  #=============generate plot list==================
  # produces a list of rows of violin plots, one per feature
  list_plot <- VlnPlot(object=seurat_obj_tmp,
                       assay=assay,
                       features = vec_features,
                       pt.size = pt.size,
                       cols=vec_group_colors,
                       sort = F,
                       #group.by = var_group,
                       same.y.lims = F,
                       slot=slot,
                       log = F,
                       combine = F,
                       flip=F)

  names(list_plot) <- vec_features

  if (is.null(aspect.ratio)) {
    aspect.ratio = 1.5*length(vec_group_colors)/length(vec_features)
    message(paste0("Using aspect ratio ", aspect.ratio))
  }

  list_plot_flip <- lapply(1:length(list_plot), function(i) {
    plot_tmp=list_plot[[i]]
    if (flip) {
      plot_tmp = plot_tmp +
        coord_flip()
    }
    plot_tmp = plot_tmp  +
      theme(
        plot.title = element_text(face = feature_fontface, size=fontsize_axistext_y),
        axis.text.y= element_blank(),
        plot.margin = margin(b=1, unit="cm"))
    #margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)

    if (i==1) {
      plot_tmp <- plot_tmp +
        theme(
          axis.text.y=element_text(
            hjust=1,
            vjust=0,
            size=fontsize_axistext_x,
            angle=30
          )
        )
    }


    plot_tmp <- plot_tmp +
      theme(
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(
          size= fontsize_axistext_x,
          angle=30),
        aspect.ratio=aspect.ratio,
        legend.position="none")
  })

  p <- patchwork::wrap_plots(... = list_plot_flip,
                             ncol = length(list_plot_flip)
  )

  if (do_plot) print(p)
  ## save the pdf figure
  if(!is.null(figurenames)){
    pdf(file = figurenames, width = width, height = height)
    print(p)
    dev.off()
  }
  #return(p)
}


top5genes <- c("IL7R","CD14","LYZ","MS4A1","CD8A","FCGR3A","MS4A7","GNLY","NKG7","FCER1A","PPBP")

plot_vln_Corgenes(seurat_obj=pbmc3k,
             assay="RNA", slot="data",
             var_group="seurat_annotations",# 细胞cluster注释列
             vec_features=top5genes,
             vec_group_colors= pal_d3(alpha =0.5)(10),
             do_plot = T
             )


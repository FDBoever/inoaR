#vis-func

# Function to tidy up ggplots
fdb_style <- function(aspect.ratio=1) {
  ggplot2::theme_classic() + ggplot2::theme(
    aspect.ratio = aspect.ratio,
    plot.title = ggplot2::element_text(hjust = 0.5),
    axis.title = ggplot2::element_text(size = 11, colour = '#000000'),
    axis.text = ggplot2::element_text(size = 10, colour = '#000000'),
    legend.justification = c(1, 1),
    legend.key.width = unit(0.25, 'cm'),
    legend.key.height = unit(0.55, 'cm'),
    legend.text = ggplot2::element_text(size = 10),
    legend.title = ggplot2::element_text(size = 11),
    panel.background = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size  = 11,  hjust = 0)
  )
}


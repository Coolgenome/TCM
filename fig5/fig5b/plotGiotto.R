#'--------------------------------------------------------------
#' filename : plotGiotto.R
#' Date : 2022-08-15
#' contributor : Yanshuo Chu
#' function: plotGiotto
#'--------------------------------------------------------------

print('<==== plotGiotto.R ====>')

suppressMessages({
    library(optparse)
    library(tidyverse)
    library(Giotto)
})

option_list = list(
    make_option(c("-f","--fovi"),
                type = 'integer',
                default = 23,
                help = 'integer',
                metavar = 'integer'),
    make_option(c("--xmin"),
                type = 'integer',
                default = 50,
                help = 'width',
                metavar = 'integer'),
    make_option(c("--xmax"),
                type = 'integer',
                default = 50,
                help = 'height',
                metavar = 'integer'),
    make_option(c("--ymin"),
                type = 'integer',
                default = 50,
                help = 'width',
                metavar = 'integer'),
    make_option(c("--ymax"),
                type = 'integer',
                default = 50,
                help = 'height',
                metavar = 'integer'),
    make_option(c("-o","--out"),
                type = 'character',
                help = 'out',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

instrs = createGiottoInstructions(python_path = "/rsrch3/home/genomic_med/ychu2/miniconda3/bin/python3",
                                  show_plot = FALSE,  
                                  return_plot = TRUE,
                                  save_plot = TRUE,
                                  save_dir = dirname(opt$out),
                                  plot_format = 'png',
                                  dpi = 200,
                                  height = 9,
                                  width = 9)

load("/rsrch3/scratch/genomic_med/ychu2/projects/p1review/CosMx/data/Giotto/Processed Data Giotto Object/SMI_Giotto_Object.RData")
gem@instructions <- instrs

metadata <- as_tibble(gem@cell_metadata$rna)
subset_cell_IDs = metadata$cell_ID[metadata$Run_Tissue_name == "Lung5_Rep2"]
gem52=subsetGiotto(gem, cell_ids = subset_cell_IDs)

fovLocGroupNum <- c(28, 20, 32, 32, 32, 30, 20, 45)
fLGN_cum <- cumsum(fovLocGroupNum)

i = 3
print(gem52@parameters$`0_fov`$fov_position.fov[(fLGN_cum[i]+1):fLGN_cum[i+1]])
print(gem52@parameters$`0_fov`$fov_position.x[(fLGN_cum[i]+1):fLGN_cum[i+1]])
print(gem52@parameters$`0_fov`$fov_position.y[(fLGN_cum[i]+1):fLGN_cum[i+1]])

ids <- gem52@parameters$`0_fov`$fov_position.fov[(fLGN_cum[i]+1):fLGN_cum[i+1]]
xs <- gem52@parameters$`0_fov`$fov_position.x[(fLGN_cum[i]+1):fLGN_cum[i+1]]
ys <- gem52@parameters$`0_fov`$fov_position.y[(fLGN_cum[i]+1):fLGN_cum[i+1]]
ids <- ids[1:30]
xs <- xs[1:30]
ys <- ys[1:30]

g <- spatPlot(gem52, point_size = 0.3) +
    geom_rect(aes(xmin = xs[opt$fovi] + opt$xmin / 5472 * (3.940 - 2.955),
                  xmax = xs[opt$fovi] + opt$xmax / 5472 * (3.940 - 2.955),
                  ymin = ys[opt$fovi] - opt$ymin / 3972 * (19.661 - 18.967),
                  ymax = ys[opt$fovi] - opt$ymax / 3972 * (19.661 - 18.967)),
              fill = NA,
              color = "red",
              size = 1) +
    geom_segment(aes(x = 0.5, xend = 1.5, y = -19, yend = -19),
                 color = "red",
                 size = 1) +
    geom_segment(aes(x = 0.5, xend = 0.5, y = -19, yend = -18.95),
                 color = "red",
                 size = 1) +
    geom_segment(aes(x = 1.5, xend = 1.5, y = -19, yend = -18.95),
                 color = "red",
                 size = 1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))
ggsave(opt$out, g)

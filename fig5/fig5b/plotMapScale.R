#'--------------------------------------------------------------
#' filename : plotMapScale.R
#' Date : 2022-08-15
#' contributor : Yanshuo Chu
#' function: plotMapScale
#'--------------------------------------------------------------

print('<==== plotMapScale.R ====>')

suppressMessages({
    library(optparse)
    library(imager)
    library(tidyverse)
})

option_list = list(
    make_option(c("-d","--data"),
                type = 'character',
                help = 'data.rds',
                metavar = 'character'),
    make_option(c("-w","--width"),
                type = 'integer',
                default = 50,
                help = 'width',
                metavar = 'integer'),
    make_option(c("-h","--height"),
                type = 'integer',
                default = 50,
                help = 'height',
                metavar = 'integer'),
    make_option(c("-o","--out"),
                type = 'character',
                help = 'out',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list, add_help_option = F);
opt = parse_args(opt_parser);

im <- load.image(opt$data)
jpeg(opt$out, width = opt$width, height = opt$height, units = "px")
par(mar = rep(0, 4))
plot(im, axes = F, xaxs="i", yaxs="i", xlim = c(1, opt$width), ylim =c(opt$height, 1))
segments(10, opt$height - 10, 288, opt$height - 10, col = "white", lwd = 2)
segments(10, opt$height - 10, 10, opt$height - 15, col = "white", lwd = 2)
segments(288, opt$height - 10, 288, opt$height - 15, col = "white", lwd = 2)
dev.off()

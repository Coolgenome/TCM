#--------------------------------------------------------------
# filename : merge.R
# Date : 2022-02-16
# contributor : Yanshuo Chu
# function: merge
#--------------------------------------------------------------

print('<==== merge.R ====>')

suppressMessages({library(optparse)
    library(ggplot2)
    library(readr)
    library(rjson)
    library(Seurat)})

print('----loading cellranger outputs----')

option_list = list(
    make_option(c("-d","--data"),
                type = 'character',
                help = 'dataFolder',
                metavar = 'character'),

    make_option(c("-o",'--out'),
                type = 'character',
                default = 'cellranger.rds',
                help = 'result file name [default = %default]',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

print(paste0("Loading data from ",opt$data))
samples = list.dirs(path = opt$data,recursive = F)
samples = basename(samples)
if(length(samples) < 1) stop(paste0('Failed to find data in ',opt$data))
obj.list = list()
for (ids in seq_along(samples)) {
    ds = samples[ids]
    print(ds)
    tenx.data = Read10X(file.path(opt$data, ds))
    tenx0 = CreateSeuratObject(counts = tenx.data, min.cells = 3,
                               min.features = 200,
                               project = ds)
    obj.list[[length(obj.list)+1]] = tenx0
}

if(length(samples) > 1) {
    combined.data = merge(x = obj.list[[1]],y = obj.list[-1],add.cell.ids = samples)
} else {
    combined.data = obj.list[[1]]
}

##Human: MT; Mouse: mt
mito.features = grep(pattern = '^MT-|^mt-', x = rownames(x = combined.data), value = T)
percent.mito <- Matrix::colSums(x = GetAssayData(object = combined.data, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = combined.data, slot = 'counts'))
combined.data[['percent.mito']] = percent.mito

print("saving output")
saveRDS(combined.data, file = file.path(opt$out, "merged.rds"))
print('----end----')


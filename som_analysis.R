### start oposSOM pipeline ###
library(oposSOM)

metacell_list <- readRDS("metacell_calc/metacell_seurat.rds")

seurat <- metacell_list[[1]]
metacell.data <- metacell_list[[2]]

env <- opossom.new(list(dataset.name = "skin_metacells",
                        database.host = "https://apr2022.archive.ensembl.org/",
                        dim.1stLvlSom = 60,
                        standard.spot.modules = "overexpression",
                        spot.coresize.modules = 2,
                        spot.threshold.modules = 0.92 ) )
# Load input data
env$indata <- metacell.data
# Define sample groups
env$group.labels <- toupper(substr(sapply( strsplit(colnames(env$indata), "_" ), "[", 3 ), 1, 1))
o <- order( env$group.labels)
env$group.labels <- env$group.labels[o]
env$indata <- env$indata[,o]
# execute
env <- opossom.run(env)

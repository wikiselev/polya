# arguments from the command line
args <- commandArgs(TRUE)

source(paste0(getwd(), "/", "../resources/ann.R"))

data.folder <- paste0(getwd(), "/", "../../data/map/", args[1], "/")
files <- list.files(path = data.folder, recursive = TRUE, include.dirs = FALSE, 
    full.names = TRUE)

# this is much faster than rbind!!! (using data.table library)
data <- data.table()
for(i in files) {
    data <- rbindlist(list(data, readRDS(i)))
}

saveRDS(data, paste0(data.folder, args[1], ".RDS"))

# for(i in files) {
# 	file.remove(i)
# }

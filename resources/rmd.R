library(knitr)
library(slidify)
library(markdown)

args <- commandArgs(TRUE)

setwd(paste0("../report/", args[2]))
unlink(args[1], recursive = TRUE)
dir.create(args[1])

rep.name <- paste0("report-", args[1])
file.copy(from = "report.Rmd", to = paste0(rep.name, ".Rmd"), overwrite = TRUE)

# tangle the R code first
purl(paste0(rep.name, ".Rmd"), documentation = 2)
# knit the document
arg <- args[1]
dir.create(paste0("figure-", arg))
knit(paste0(rep.name, ".Rmd"))
file.remove(paste0(rep.name, ".Rmd"))

######
# slidify functions

# I've changed the author function, so that now it does not open index.Rmd and
# don't stop the script therefore

copy_dir <- function(from, to){
  # if (!(file.exists(to))){
    # dir.create(to, recursive = TRUE)
    message('Copying files to ', to, '...')
    file.copy(list.files(from, full.names = T), to, recursive = TRUE)
  # }
}
######

setwd(args[1])

scaffold = system.file('skeleton', package = 'slidify')
copy_dir(scaffold, ".")

file.copy(from = paste0("../", rep.name, ".md"), to = "script.md", overwrite = TRUE)
file.remove(paste0("../", rep.name, ".md"))
markdownToHTML("script.md", "report.html",  stylesheet = "~/slidify_css/markdown.css")
file.copy(from = paste0("../", rep.name, ".R"), to = "script.R", overwrite = TRUE)
file.remove(paste0("../", rep.name, ".R"))
file.copy(from = paste0("../", rep.name, ".Rmd"), to = "script.Rmd", overwrite = TRUE)
file.remove("index.Rmd")
dir.create("figure")
files <- list.files(paste0("../figure-", args[1], "/figure/"), full.names = TRUE)
file.copy(from=files, to="figure/", recursive = TRUE, overwrite = TRUE)
unlink(paste0("../figure-", args[1], "/"), recursive = TRUE)

# slidify the knitted script
slidify("script.md")
file.rename(from = "script.txt", to = "slides.html")
file.copy("~/slidify_css/default.css",
	"libraries/frameworks/io2012/css/default.css", overwrite = TRUE)
file.copy("~/slidify_css/slidify.css",
	"libraries/frameworks/io2012/css/slidify.css", overwrite = TRUE)

unlink(paste0("~/public_html/misc/polya/", args[2], "/", args[1]), recursive = TRUE)
dir.create(paste0("~/public_html/misc/polya/", args[2], "/", args[1]))
files <- list.files(full.names = TRUE)
file.copy(from=files, to=paste0("~/public_html/misc/polya/", args[2], "/", args[1], "/"), recursive = TRUE, overwrite = TRUE)

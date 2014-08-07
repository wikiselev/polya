# load my general functions
source(paste0(getwd(), "/", "../resources/ann.R"))
out.folder <- "report/"
# pick up protein files according to iscen
prot.folder <- "../../data/prot/processed/"
prot.files <- list.files(path = prot.folder)
prot.names <- sapply(strsplit(prot.files, "\\."), "[[", 1)

magnify_scores <- function(data) {
    setnames(data, c("coord", "sum.score"), c("coord", "score"))
    data_ctrl <- data.table(coord = seq(-150, 150), score = 0)
    setkey(data, coord)
    data <- data[data_ctrl]
    data[is.na(score), score := as.integer(0)]
    
    scores <- data$score
    scores1 <- scores
    
    # for 5 nucleotide cluster peaks
    #######
    for(i in 3:(length(scores) - 2)) {
            scores1[i] <- scores[i] + scores[i - 1] + scores[i + 1] +
                    scores[i - 2] + scores[i + 2]
    }
    
    scores1[1] <- 0
    scores1[2] <- 0
    scores1[length(scores1)] <- 0
    scores1[length(scores1) - 1] <- 0
    #######
    
    data[,score:=scores1]
    data <- data[score != 0, list(coord,score)]
    
    return(data)
}

prot_data_make_up <- function(data) {
  setkeyv(data, c("coord", "cond", "gene.id", "css.id", "usage.let", "prot"))
  p.data <- data[,list(sum.score = sum(score)), by = key(data)]

  p.data.magn <- data.table()
  for(con in unique(p.data[,cond])) {
    for(p.name in prot.names) {
      dat <- 
      p.data[cond == con & gene.id == 1 & prot == p.name,
      list(coord, sum.score)]
      if(dim(dat)[1] != 0) {
        dat <- magnify_scores(dat)
        dat[,cond := con]
        dat[,gene.id := 1]
        dat[,css.id := 1]
        dat[,usage.let := "us"]
        dat[,prot := p.name]
        p.data.magn <- rbindlist(list(p.data.magn, dat))
      }
      for(c.ind in 1:2) {
        for(us in c("us_un", "un_us", "in_in")) {
          dat <- p.data[cond == con & gene.id == 2 & css.id == c.ind &
          usage.let == us & prot == p.name, list(coord, sum.score)]
          if(dim(dat)[1] != 0) {
            dat <- magnify_scores(dat)
            dat[,cond := con]
            dat[,gene.id := 2]
            dat[,css.id := c.ind]
            dat[,usage.let := us]
            dat[,prot := p.name]
            p.data.magn <- rbindlist(list(p.data.magn, dat))
          }
        }
      }
    }
  }

  setkeyv(data, c("cond", "gene.id", "css.id", "usage.let", "prot"))
  ann.data <- data[,list(gene.num = length(unique(gene_id))), by = key(data)]
  return(list(p.data = p.data.magn, ann.data = ann.data))
}

prot_plot1 <- function(d) {
  p <- ggplot(as.data.frame(d), aes(x = coord, y = score, color = cond)) +
          geom_point(size = 0.5) +
          stat_smooth(span = 0.2, se = FALSE, size = 1.3) +
          facet_wrap( ~ prot, ncol = 4, scales = "free_y") +
          # geom_text(aes(label = as.character(gene.num)), data = a,
          #     hjust=0, show_guide = FALSE, size = 3) +
          labs(x = "Coordinate, [bp]", 
               y = "Profile ~ read number")
  return(p)
}

prot_plot2 <- function(d) {
  p <- ggplot(as.data.frame(d), aes(x = coord, y = score, color = as.factor(css.id))) +
          geom_point(size = 0.5) +
          stat_smooth(span = 0.2, se = FALSE, size = 1.3) +
          facet_wrap( ~ prot, ncol = 4, scales = "free_y") +
          # geom_text(aes(label = as.character(gene.num)), data = a,
          #     hjust=0, show_guide = FALSE, size = 3) +
          labs(x = "Coordinate, [bp]", 
               y = "Profile ~ read number")
  return(p)
}

plot_separated_prot_a_seq <- function(prot.name)
{
  d <- p.data$p.data[prot == prot.name]
  p <- ggplot(as.data.frame(d[gene.id == 2]), aes(x = coord, y = score, color = as.factor(css.id))) +
          geom_point(size = 0.9) +
          stat_smooth(span = 0.2, se = FALSE, size = 1.0) +
          facet_grid(cond ~ usage.let) +
          # geom_text(aes(label = as.character(gene.num)), data = a,
          #     hjust=0, show_guide = FALSE, size = 3) +
          labs(x = "Coordinate, [bp]", 
               y = "Profile ~ read number")
  if(arg == "a-seq") {
    pdf(paste0(out.folder, prot.name, "a-seq-2.pdf"), w = 10, h = 4)
  }
  if(arg == "atlas") {
    pdf(paste0(out.folder, prot.name, "atlas-2.pdf"), w = 10, h = 12)
  }
  print(p)
  dev.off()

  p <- ggplot(as.data.frame(d[gene.id == 1]), aes(x = coord, y = score)) +
          geom_point(size = 0.9) +
          stat_smooth(span = 0.2, se = FALSE, size = 1.0) +
          facet_grid(cond ~ .) +
          # geom_text(aes(label = as.character(gene.num)), data = a,
          #     hjust=0, show_guide = FALSE, size = 3) +
          labs(x = "Coordinate, [bp]", 
               y = "Profile ~ read number")
  if(arg == "a-seq") {
    pdf(paste0(out.folder, prot.name, "a-seq-1.pdf"), w = 3, h = 4)
  }
  if(arg == "atlas") {
    pdf(paste0(out.folder, prot.name, "atlas-1.pdf"), w = 3, h = 12)
  }
  print(p)
  dev.off()
}

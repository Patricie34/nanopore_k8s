library(BioCircos)
library(dplyr)
library(ggplot2)
library(htmlwidgets)


bgcol = "#FFF9F2"
poscol = "#B0E66A"
negcol = "#E6675C"
neutcol = "#75C1E6"

indels = function(INS_fil,DEL_fil ) {
  INSdf = data.frame(Chr = INS_fil$V1,
                     Pos = INS_fil$V2,
                     Type = rep('INS', length(INS_fil$V1)))
 
  DELdf = data.frame(Chr = DEL_fil$V1,
                     Pos = DEL_fil$V2,
                     Type = rep('DEL', length(DEL_fil$V1)))
  
  data = rbind(INSdf, DELdf)
  
  tracks = BioCircosTracklist() + BioCircosBackgroundTrack(
    "bars_background",
    fillColors = bgcol,
    minRadius = 0.55,
    maxRadius = 0.85
  )
  
  for (chr in unique(data$Chr)) {
    Chrdata = data %>% filter(Chr == chr)
    Chrdata$Type = factor(Chrdata$Type, levels = c("INS", "DEL"))
    
    if (nrow(Chrdata)<2) { #One datapoint causes problems
      next
    }
    p = ggplot(Chrdata, aes(x = Pos, fill = Type)) +
      geom_density(
        adjust = 1,
        colour = NA,
        position = "fill",
        alpha = 1
      )
    pg <- ggplot_build(p)
    
      if (sum(Chrdata$Type=='INS')){
        tracks = tracks +  BioCircosBarTrack(
      paste0("Chr", chr),
      chromosome = chr,
      starts = pg$data[[1]]$x[pg$data[[1]]$group ==1],
      ends = pg$data[[1]]$x[pg$data[[1]]$group == 1] + rep(diff(pg$data[[1]]$x[pg$data[[1]]$group == 1])[1], length(pg$data[[1]]$x)),
      values = pg$data[[1]]$y[pg$data[[1]]$group ==1],
      # labels = sprintf("INS: %.0f%%", 100 * (1 - pg$data[[1]]$y[pg$data[[1]]$group == 2])),
      labels = sprintf("INS: %.0f%%", 100 *  pg$data[[1]]$y[pg$data[[1]]$group == 1]),
      color = poscol,
      #green
      range = c(0, 1),
      minRadius = 0.55,
      maxRadius = 0.85
        ) }
    
    if (sum(Chrdata$Type=='DEL')){
    tracks = tracks +  BioCircosBarTrack(
        paste0("Chr", chr),
        chromosome = chr,
        starts = pg$data[[1]]$x[pg$data[[1]]$group == 2],
        ends = pg$data[[1]]$x[pg$data[[1]]$group == 2] + rep(diff(pg$data[[1]]$x[pg$data[[1]]$group == 2])[1], length(pg$data[[1]]$x)),
        values = pg$data[[1]]$y[pg$data[[1]]$group == 2],
        labels = sprintf("DEL: %.0f%%", 100 * pg$data[[1]]$y[pg$data[[1]]$group == 2]),
        color = negcol,
        #red
        range = c(0, 1),
        minRadius = 0.55,
        maxRadius = 0.85
      )}
  }
  return (tracks)
}

# START
#CNV_ref <- read.table("hg19_lens.txt", header = FALSE, sep = "")

args = commandArgs(trailingOnly=TRUE)
sorted_cns <- args[1]
vars_tsv <- args[2]
name <- args[3]
grch_lens <- args[4]

CNV_ref <- read.table(grch_lens, header = FALSE, sep = "")
# CNV_ref <- read.table("grch38_lens.txt", header = FALSE, sep = "")
# 
# CNVs <- read.table("signatures_BRNO1727/BRNO1727.sorted.cns", header = TRUE, sep = "")
CNVs <- read.table(sorted_cns, header = TRUE, sep = "")

CNVs_fil = CNVs %>%
  filter(!grepl('GL|M|K', chromosome))
CNVs_fil$log2[CNVs_fil$log2 < -2] = -2
CNVs_fil$log2[CNVs_fil$log2 > 2] = 2


INS <- read.table("signatures/ins.bed", header = FALSE, sep = "\t")
INS_fil = INS %>%
  filter(!grepl('GL|M|K', V1)) %>%
  filter(V5 > 3)

DEL <- read.table("signatures/del.bed", header = FALSE, sep = "\t")
DEL_fil = DEL %>%
  filter(!grepl('GL|M|K', V1)) %>%
  filter(V5 > 3)

INV <- read.table("signatures/inv.bed", header = FALSE, sep = "\t")
INV_fil = INV %>%
  filter(!grepl('GL|M|K', V1)) %>%
  filter(V5 > 3)

# BNDs <- read.table("BRNO1727.slimmed.tsv", header = TRUE, sep = "")
BNDs <- read.table(vars_tsv, header = TRUE, sep = "")


BNDs_fil = BNDs %>% mutate(distance = as.integer(distance)) %>%
                            filter(Support_nonUnique > 3,
                           is.na(distance) | distance > 100000,
                           !grepl('GL|M|K', ChrTo),
                           !grepl('GL|M|K', ChrFrom))

tracklist = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0, maxRadius = 0.55, borderSize = 0,fillColors = bgcol)
BND_labels = sprintf( "Chr%s:%s / Chr%s:%s \n  Support = %s", BNDs_fil$ChrFrom, BNDs_fil$PosFrom, BNDs_fil$ChrTo,
                      BNDs_fil$PosTo, BNDs_fil$Support_nonUnique)
tracklist_BNDs = tracklist
maxsup = max(BNDs_fil$Support_nonUnique)
minsup = min(BNDs_fil$Support_nonUnique)


for (i in unique(BNDs_fil$Support_nonUnique)) {
  mask = BNDs_fil$Support_nonUnique == i
  tracklist_BNDs = tracklist_BNDs + BioCircosLinkTrack(
    'myLinkTrack',
    BNDs_fil$ChrFrom[mask],
    BNDs_fil$PosFrom[mask],
    BNDs_fil$PosFrom[mask],
    BNDs_fil$ChrTo[mask],
    BNDs_fil$PosTo[mask],
    BNDs_fil$PosTo[mask],
    maxRadius = 0.55,
    labels = BND_labels[mask],
    displayLabel = FALSE,
    color = neutcol,
    width = paste0((i-minsup)/maxsup/5+0.05, "em")

  )
}


# Create CNV track
CNVtracks =  BioCircosCNVTrack('cnv_normal_line', as.character(CNV_ref$V1), CNV_ref$V2, CNV_ref$V3, 0, color = "#D3D3D3",range = c(-2, 2),
                               minRadius = 0.85,maxRadius = 1) +
  BioCircosCNVTrack('cnv_track', as.character(CNVs_fil$chromosome), CNVs_fil$start, CNVs_fil$end, CNVs_fil$log2, 
                    color = neutcol, labels = CNVs_fil$gene, range = c(-2, 2), 
                    minRadius = 0.85, maxRadius = 1) +
  BioCircosBackgroundTrack("arcs_background", fillColors =  bgcol,
                           minRadius = 0.85, maxRadius = 1)
# Add background
indelTrack = indels(INS_fil,DEL_fil)
tracks = CNVtracks  + tracklist_BNDs + indelTrack

circ = BioCircos(tracks, genomeTicksDisplay = F, genomeFillColor = rep("white", 24), genomeLabelDy = -35, chrPad = 0.01, genomeBorderColor = F, width = "1080px", height = "720px")
saveWidget(circ, paste0(name,"_circos.html"), selfcontained = T)

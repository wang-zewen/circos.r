library(circlize)
library(grid)
library(ComplexHeatmap)
setwd("d://de/other/sun/circos/") ## working directory
seq_stat <- read.delim("4.length",stringsAsFactors = FALSE) ## chromosome length file
sample1 <- read.delim("control4.snp",stringsAsFactors = FALSE) ## snp position file of control group
sample2 <- read.delim('gdm4.snp',stringsAsFactors = FALSE) ## snp position file of experimental group
sample_list <- function(dat){
  dat <- list(dat[which(dat$alt_type == "SNP:A>T|T>A"),],
              dat[which(dat$alt_type == "SNP:A>G|T>C"),],
              dat[which(dat$alt_type == "SNP:A>C|T>G"),],
              dat[which(dat$alt_type == "SNP:G>A|C>T"),],
              dat[which(dat$alt_type == "SNP:G>T|C>A"),],
              dat[which(dat$alt_type == "SNP:G>C|C>G"),])
  return(dat)
}
sample1_list <- sample_list(sample1)
sample2_list <- sample_list(sample2)
pdf('4.pdf',width = 8,height = 8) ## out put pdf filename
circle_size = unit(1,'snpc')
circos.par(gap.degree = 3,start.degree = 90)
circos.genomicInitialize(seq_stat,plotType = c('axis','labels'),major.by = 250000,track.height = 0.05)
circos.genomicTrackPlotRegion(
  seq_stat,track.height=0.05,stack = TRUE,bg.border=NA,
  panel.fun = function(region,value,...){
    circos.genomicRect(region,value,col='#8292b4',border = NA,...)
  }
)

color_assign <- c('#BC80BD', '#FDB462', '#80B1D3', '#FB8072', '#8DD3C7', '#FFFFB3')  ##
circos.genomicTrackPlotRegion(
  sample1_list,ylim=c(1,36),track.height=0.12,bg.border='black',bg.lwd=0.4,bg.col="#f3d8ae",
  panel.fun = function(region,value,...){
    if (nrow(sample1_list[[getI(...)]]) > 0){
    circos.genomicPoints(region,value,pch=16,cex=0.5,col = color_assign[getI(...)],...)}
    circos.yaxis(labels.cex = 0.2,lwd=0.1,tick.length = convert_x(0.15,'mm'))
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.text(mean(xlim),mean(ylim),'Control',cex=1.7,col='black',facing = 'inside',niceFacing = TRUE )
  }
)

circos.genomicTrackPlotRegion(
  sample2_list,ylim=c(1,18),track.height=0.12,bg.border='black',bg.lwd=0.4,bg.col="#65c1b3",
  panel.fun = function(region,value,...){
    if (nrow(sample2_list[[getI(...)]]) > 0){
    circos.genomicPoints(region,value,pch=16,cex=0.5,col = color_assign[getI(...)],...)}
    circos.yaxis(labels.cex = 0.2,lwd=0.1,tick.length = convert_x(0.15,'mm'))
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    ylim
    circos.text(mean(xlim),mean(ylim),'GDM',cex=1.7,col='black',facing = 'inside',niceFacing = TRUE )
  }
)

snv_legend <- Legend(
  at = c(1,2,3,4,5,6),
  labels = c(' SNP: A>T|T>A', ' SNP: A>G|T>C', ' SNP: A>C|T>G', ' SNP: G>A|C>T', ' SNP: G>T|C>A', ' SNP: G>C|C>G'),
  labels_gp = gpar(fontsize=6),title='Variance type',title_gp=gpar(fontsize=7),
  grid_height = unit(0.4,'cm'),grid_width=unit(0.4,'cm'),type='points',background=NA,
  legend_gp = gpar(col= color_assign)
)

pushViewport(viewport(x = 0.5, y = 0.5))
grid.draw(snv_legend)
upViewport()

circos.clear()
dev.off()

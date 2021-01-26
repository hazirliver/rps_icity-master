if (!require('VennDiagram')) install.packages('VennDiagram'); library('VennDiagram')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('nVennR')) install.packages('nVennR'); library('nVennR')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('stringr')) install.packages('stringr'); library('stringr')



RT <- fread("./RT/icity_pfam_RT.csv")
TR <- fread("./TR/icity_pfam_TR.csv")
VR <- fread("./VR/icity_pfam_VR.csv")


############################
# Построение диаграммы Венна
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
pdf(NULL)
par(bg = NA)
venn.diagram(
  x = list(
    unique(RT$pfam) , 
    unique(TR$pfam)  , 
    unique(VR$pfam) 
  ),
  category.names = c("RT" , "TR" , "VR"),
  filename = './Venn/venn_nice.png',
  output = TRUE ,
  imagetype="png" ,
  height = 900 , 
  width = 900 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 0.25,
  col=c("#4d5156", '#4d5156', '#4d5156'),
  fill = c(alpha("#1B9E77",0.75), alpha('#D95F02',0.75), alpha('#7570B3',0.75)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#1B9E77", '#D95F02', '#7570B3'),
  rotation = 1,
  lty =2, 
  label.col = "#4d5156"
)
dev.off()
##############################
Venn <- plotVenn(list("RT"=unique(RT$pfam), "TR"=unique(TR$pfam), "VR"=unique(VR$pfam)), 
                 outFile = "./Venn/venn.svg")
intersects <- listVennRegions(Venn)
intersect_names <- c()
for(i in names(intersects)){
  write.table(intersects[[i]], paste0("./Venn/",i,".csv"), quote = F, row.names = F, col.names = F)
  i_name <- str_extract_all(i, "\\([^()]+\\)")[[1]]
  i_name <- substring(i_name, 2, nchar(i_name)-1)
  intersect_names <- c(intersect_names,i_name)
}

pfam_intersect <- stack(setNames(intersects, intersect_names))

all_gis <- rbind(cbind(RT, "RT"),
                 cbind(TR, "TR"),
                 cbind(VR, "VR"))
tmp <- merge(all_gis, pfam_intersect, by.x = "pfam", by.y = "values")
to_out <- tmp[,c("gi", "ind")]
write.table(to_out, "./Venn/gi_intersection.tsv", quote = F, row.names = F, sep = "\t")




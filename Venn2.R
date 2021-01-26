if (!require('VennDiagram')) install.packages('VennDiagram'); library('VennDiagram')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('nVennR')) install.packages('nVennR'); library('nVennR')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('stringr')) install.packages('stringr'); library('stringr')



RT_icity <- fread("./RT/icity_pfam_RT.csv")
TR_icity <- fread("./TR/icity_pfam_TR.csv")
VR_icity <- fread("./VR/icity_pfam_VR.csv")


RT_relevance <- fread("./RT/relevance_pfam_named_RT.csv")
TR_relevance <- fread("./TR/relevance_pfam_named_TR.csv")
VR_relevance <- fread("./VR/relevance_pfam_named_VR.csv")

RT_good_clusters <- RT_relevance[RT_relevance$median_dist_to_baits < 5 & RT_relevance$icity > 0.6 & RT_relevance$e_size_baits != -1]
TR_good_clusters <- TR_relevance[TR_relevance$median_dist_to_baits < 5 & TR_relevance$icity > 0.6 & TR_relevance$e_size_baits != -1]
VR_good_clusters <- VR_relevance[VR_relevance$median_dist_to_baits < 5 & VR_relevance$icity > 0.6 & VR_relevance$e_size_baits != -1]


#
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
pdf(NULL)
par(bg = NA)
venn.diagram(
  x = list(
    RT_good_clusters$Cluster_ID , 
    TR_good_clusters$Cluster_ID  , 
    VR_good_clusters$Cluster_ID 
  ),
  category.names = c("RT" , "TR" , "VR"),
  filename = './Venn/venn_nice_good_dist_icity.png',
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
#

Venn2 <- plotVenn(list("RT"=RT_good_clusters$Cluster_ID, "TR"=TR_good_clusters$Cluster_ID, 
                      "VR"=VR_good_clusters$Cluster_ID), 
                 outFile = "./Venn/venn_good_dist_icity.svg")
intersects2 <- listVennRegions(Venn2)
intersect_names2 <- c()
for(i in names(intersects2)){
  write.table(intersects2[[i]], paste0("./Venn/",i,"good_dist_icity.csv"), quote = F, row.names = F, col.names = F)
  i_name <- str_extract_all(i, "\\([^()]+\\)")[[1]]
  i_name <- substring(i_name, 2, nchar(i_name)-1)
  intersect_names2 <- c(intersect_names2,i_name)
}

pfam_intersect2 <- stack(setNames(intersects2, intersect_names2))

all_gis2 <- rbind(cbind(RT_icity, "RT"),
                 cbind(TR_icity, "TR"),
                 cbind(VR_icity, "VR"))
tmp2 <- merge(all_gis2, pfam_intersect2, by.x = "pfam", by.y = "values")
to_out2 <- tmp2[,c("gi", "ind")]
to_out2$gi <- as.character(to_out2$gi)

write.table(to_out2, "./Venn/gi_intersection_good_dist_icity.tsv", quote = F, row.names = F, sep = "\t")
write.table(to_out2$gi, "./Venn/gi_intersection_good_dist_icity_only_gis.tsv", quote = F, row.names = F, sep = "\t", col.names = F)


######
rps_vs_immun_ig_lectin <- fread("rps_blast_results_good_subset_vs_immun_ig_lectin.tsv", sep = "\t", )
library(tidyr)
library(dplyr)

y <- separate(rps_vs_immun_ig_lectin, V1, into = c("id1", "gi", "id2", "gb"), sep = "\\|", remove = T)
y$id1 <- NULL; y$id2 <- NULL

rps_with_intersections <- merge(y, to_out2, by.x = "gi", by.y = "gi")
rps_with_intersections_unique_gis <- rps_with_intersections[ , .SD[which.min(V3)], by = gi]

barplot(table(rps_with_intersections_unique_gis$ind))


pfams <- colsplit(rps_with_intersections_unique_gis$V2,",",c("pfam","description"))
pfams_table <- table(pfams$pfam)
pfams_more_5 <- pfams[pfams$pfam %in% names(pfams_table[pfams_table > 5]), ]

barplot(table(pfams_more_5$pfam))


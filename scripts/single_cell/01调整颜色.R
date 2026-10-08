###调整配色
library(RColorBrewer) 
library(viridis)
library(wesanderson)
library(ggsci)
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


#################
###########  col_vector ;  cols = col_vector
# ######  p3 = DimPlot(scRNA, group.by="ferr_group", label=F, label.size=4.5, 
# reduction='umap',cols = col_vector)+theme(title = element_blank())
###########################



pie(rep(1,n), col=sample(col_vector, n))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pie(rep(6,n), col=sample(color, n))
col_vector
col_vector =c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest1"), wes_palette("Cavalcanti1"), wes_palette("GrandBudapest2"), wes_palette("FantasticFox1"))
pal <- wes_palette("Zissou1", 10, type = "continuous")
pal2 <- wes_palette("Zissou1", 5, type = "continuous")
pal[1:10]


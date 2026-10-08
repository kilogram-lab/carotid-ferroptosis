# 从CellMarker下载细胞注释表格
# CellMarker：httLOCAL_PROJECT_PATH
# cell_marker <- read.csv("CellMarker.csv" )
cell_marker <- fread(file = "CellMarker.csv")
colnames(cell_marker)[c(5,6)] <- c("cell_name" , "cell_marker")


cell_marker <- cell_marker[,c(5, 6)] %>%
  arrange(cell_marker)

















?separate_rows
cell_marker <- cell_marker %>%
  separate_rows(Cell.Marker, sep = ", ") %>%
  na.omit() %>%
  distinct() %>% 
  arrange(Cell.Type) 
# View(cell_marker)
cell_marker$Cell.Marker <- str_trim(cell_marker$Cell.Marker, "left") 
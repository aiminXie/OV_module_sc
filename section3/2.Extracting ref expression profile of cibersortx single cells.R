#Some cells were randomly selected to construct the expression profile
cellType <- as.character(mergeData_immune$cellType)
names(cellType) <- colnames(mergeData_immune)
cell_list <- split(colnames(mergeData_immune),f = cellType)
sampleCell <- lapply(cell_list,function(x){
    if(length(x)<400){
        subsample <- x
    }else{
        subsample <- sample(x = x,size = 400,replace = F)
    }
 return(subsample)
})
sampleCell <- do.call(c,sampleCell)
cellType_sub <- cellType[sampleCell]
data <- GetAssayData(subset(mergeData_immune,cells = sampleCell),assay='RNA',slot='counts')
df <- t(t(data) / colSums(data) )* 1e6
df <- as.data.frame(x = df)

cellMarkers <- readRDS(file = './output/1.sangleCell_pre/OV_mergeData_immune_cellMarkers.RDS')
cellMarkers <- subset(cellMarkers,!grepl(pattern = '^RP[LS]',gene))
cellMarkers <- subset(cellMarkers,!grepl(pattern = '^MT-',gene))

#output the cibersortx single cell ref expression profile
df <- df[,sampleCell]
df <- df[unique(cellMarkers$gene),]
colnames(df) <- cellType_sub
x<- rownames(df)
rownames(df)<-NULL
final <- cbind(x,df)
colnames(final)[1]<- 'GeneSymbol'
y <- colnames(final)
y <- sub(pattern = '\\.[0-9]+',replacement = '',y)
colnames(final)<-NULL
final <- rbind(y,final)
write.table(x = final,file = './output/1.sangleCell_pre/mergeData_cibersort_singleCellRef.txt',sep = '\t',append = F,row.names = F,col.names = F)
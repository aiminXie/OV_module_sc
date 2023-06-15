#####
## gene ID convert for Exp matirx
###
exp_idchange <- function(data,idRaw='ENSEMBL',idTo='SYMBOL'){
  #id <- c("SYMBOL", "ENTREZID","ENSEMBL")
  #id_change <- setdiff(id,idRaw)
    
  #library(clusterProfiler)
  library(AnnoProbe)
  geneID <- rownames(data)
  #geneID.table=select(org.Hs.eg.db,keys=geneID,columns = id_change, keytype=idRaw)
  geneID.table <- annoGene(IDs = geneID,ID_type =  idRaw)[,c(idRaw,idTo)]


  exp <- data.frame(probe_id=rownames(data),data)
  exprSet_symbol <- merge.data.frame(geneID.table, exp, by.x = idRaw, by.y = "probe_id",all=F)
  exprSet_symbol <- aggregate(x = exprSet_symbol[,-c(1,2)],by = list(exprSet_symbol[,idTo]), FUN = mean)
  res_mRNA <- as.matrix(exprSet_symbol[,-1])
  rownames(res_mRNA) <- exprSet_symbol[,1]
  res_mRNA
}


#############################
#Pairwise Fisher’s exact test 
#############################

make.Pairwise.Tables <- function(geneSets1,geneSets2,background)
{
 mem1 <- do.call(cbind,lapply(geneSets1,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))
 mem2 <- do.call(cbind,lapply(geneSets2,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))

 d <- t(mem1) %*% mem2;
 b <- abs(t(mem1) %*% (mem2-1))
 c <- abs(t(mem1-1) %*% (mem2))
 a <- t(mem1-1) %*% (mem2-1);

 ij <- do.call(rbind,lapply(1:length(geneSets1),function(i,j) cbind(rep(i,length(j)),j),j = 1:length(geneSets2)))

 pairwise.tables <- lapply(1:nrow(ij),function(i,ij,a,b,c,d) as.table(matrix(c(a[ij[i,1],ij[i,2]],b[ij[i,1],ij[i,2]],c[ij[i,1],ij[i,2]],d[ij[i,1],ij[i,2]]),nrow = 2)),ij = ij,a = a,b = b,c = c,d = d)
 names(pairwise.tables) <- paste(names(geneSets1)[ij[,1]],"~vs~",names(geneSets2)[ij[,2]],sep = "")
 return(pairwise.tables)
}
                           
do.FisherExactTest <- function(table.count,N_bg = NULL)
{
 if (is.null(N_bg)) N_bg = sum(rowSums(table.count))
 
 out <- fisher.test(x = table.count,or = 1,alternative = "greater")
 odds.ratio <- out$estimate
 p.value <- out$p.value;
 geneSet1.count <- rowSums(table.count)[2]
 geneSet2.count <- colSums(table.count)[2]
 expected.count <- geneSet1.count/N_bg * geneSet2.count
 overlap.count <- table.count[2,2];
 fold.change <- overlap.count/expected.count
 
 out <- c(N_bg,geneSet1.count,geneSet2.count,expected.count,overlap.count,fold.change,odds.ratio,p.value)
 names(out) <- c("Background","set1_size","set2_size","expected.overlap","actual.overlap","enrichment.foldchange","odds.ratio","FET_pvalue")
 return(out)
}

require(parallel)
require(foreach)
require(iterators)
perform.AllPairs.FET <- function(geneSets1,geneSets2,background,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
{
 pairwise.tables <- make.Pairwise.Tables(geneSets1,geneSets2,background)
 
 if (do.multicore)
 {
  #registerDoMC(n.cores)
  set.parallel.backend(n.cores)
  fact <- do.call(c,lapply(1:n.cores,function(i,dn) rep(i,dn),dn = ceiling(length(pairwise.tables)/n.cores)))
  fact <- factor(fact[1:length(pairwise.tables)])
  split.tables <- lapply(split(1:length(pairwise.tables),fact),function(i,obj) obj[i],obj = pairwise.tables);rm(fact)
  output <- foreach (tbl = split.tables,.combine = 'c') %dopar% {
            out <- lapply(tbl,do.FisherExactTest) 
			return(out)
  }
 
 }else{
  output <- lapply(pairwise.tables,do.FisherExactTest)
 }
 # collect outputs
 output <- data.frame(set1_Name = gsub("~vs~(.*)$","",names(output)),set2_Name = gsub("^(.*)~vs~","",names(output)),as.data.frame(do.call(rbind,output)))
 int.element <- mapply(FUN = function(x,y) paste(intersect(x,y),collapse = ","),x = geneSets1[as.character(output[[1]])],y = geneSets2[as.character(output[[2]])],SIMPLIFY = F)
 int.element <- do.call(c,int.element)
 if (adjust.FET.pvalue) output <- cbind.data.frame(output,data.frame(corrected.FET.pvalue = p.adjust(output$FET_pvalue,"bonferroni"),intersects = int.element))
 
 return(output)
}   

############
##
############
##
data_explor.discrete2discrete <- function(data,core_dimension,inter_dimension,plot=F,cor_class=unique(data[,core_dimension]),remove=c("","[Not Available]","[Unknown]","[Not Evaluated]"),position = c('full','stack')){
  data <- na.omit(data[,c(core_dimension,inter_dimension)])
  data <- data[!(data[,inter_dimension] %in% remove),]
  test.table <- table(data)
  if(ncol(test.table)<=1 | nrow(test.table)<=1){
    return(NA)
  }
  #print(test.table)
  if(nrow(test.table)>2){
    test.res <- fisher.test(test.table,simulate.p.value = TRUE, B = 1e5)$p.value
  }else{
    test.res <- fisher.test(test.table)$p.value
  }
  if(plot){
    library(ggplot2)
    #pdf(path)
    data <- as.data.frame(test.table,stringsAsFactors = F)
    data$sum <- rowSums(test.table)
    data$p <- data$Freq/data$sum
    
    plot2 <- ggplot(data=data,aes(data[,core_dimension],weight = Freq))+geom_bar(aes(fill=data[,inter_dimension]),position=position) + 
      ggtitle(paste(inter_dimension,"(P.value: ",round(test.res,digits = 4),")",sep = "")) + geom_text(aes(label = ..count..),stat="count",position = position) +
      xlab(core_dimension) + ylab("percentage") +
      guides(fill=guide_legend(title=inter_dimension))
    #print(plot2)
    #dev.off()
	res <- list(plot=plot2 ,p.value=test.res)
  }else{
    res <- test.res
  }
  return(res)
}

#' @@ 探索离散维度与一个连续维度的关联
#' @param  data data.frame 提供数据用于绘图 core_dimension,inter_dimension两个参数是data中对应的列
#' @param  core_dimension 字符 核心关注的维度（data中的列）
#' @param  inter_dimension 字符 与核心维度相交互的维度（data中的列）
#' @param  plot 逻辑 是否绘制图像
#' @param  path 字符 绘制图像时对应的输出路径
#' @returnType 数值
#' @return 统计检验的P值
#' 
#' @author aimin
##
##
data_explor.discrete2continuous <- function(data,core_dimension,inter_dimension,plot=F){
  if(nrow(data)<2){
    return(list(plot=NULL ,p.value=NULL))
  }
  data <- na.omit(data[,c(inter_dimension,core_dimension)])
  
  names(data) <- c(inter_dimension,core_dimension)
  data.list <- split(data[,1],f = data[,2])
  if(length(data.list)==2){
    test.res <- wilcox.test(x = data.list[[1]],y = data.list[[2]])$p.value
  }else{
    index.combn <- combn(names(data.list),m = 2)
    #test.res <- kruskal.test(x=data[,2],g=as.factor(data[,1]))
    test.p <- apply(index.combn,2,function(x){
      wilcox.test(x = data.list[[x[1]]],y = data.list[[x[2]]])$p.value
    })
    test.res <- min(test.p)
  }
  
mytheme <- theme_bw()+ theme(axis.text.x=element_blank(), 
                 axis.text.y=element_text(size=12), 
                 axis.title=element_text(size = 13), 
                 legend.text=element_text(size=12),
                 legend.title=element_text(size=12),
                 axis.line = element_line(size=0.7), 
                 panel.border = element_blank(),
                 panel.grid = element_blank())
  if(plot){
    library(ggplot2)
    library(ggsignif)
    library(ggbeeswarm)
    #pdf(path)
    if(length(data.list)>2 && length(which(test.p<0.05))>=1){
      plot2<-ggplot(data=data,aes(data[,core_dimension],data[,inter_dimension],fill=data[,core_dimension])) +geom_violin(trim=FALSE,color="white") + geom_boxplot(width=0.2,position=position_dodge(0.9)) + 
        labs(inter_dimension) + 
        geom_signif(comparisons =  combn(c(names(data.list)),2,simplify = F)[which(test.p<0.05)], test = wilcox.test, step_increase = 0.2) +
        ggtitle(inter_dimension) + xlab(core_dimension) + ylab(inter_dimension) +scale_fill_jco()+
        guides(col=guide_legend(title=core_dimension))+mytheme
    }else{
      plot2<-ggplot(data=data,aes(data[,core_dimension],data[,inter_dimension],fill=data[,core_dimension])) + geom_violin(trim=FALSE,color="white") + geom_boxplot(width=0.2,position=position_dodge(0.9)) + 
        labs(inter_dimension) + 
        geom_signif(comparisons =  combn(c(names(data.list)),2,simplify = F), test = wilcox.test, step_increase = 0.2) +
        ggtitle(inter_dimension) + xlab(core_dimension) + ylab(inter_dimension) +scale_fill_jco()+
        guides(col=guide_legend(title=core_dimension))+mytheme
    }
    
    #print(plot2)
    return(list(plot=plot2 ,p.value=test.res))
    #dev.off()
  }
  else{
    return(list(plot=NULL ,p.value=test.res))
  }
}

#'@@散点图展示两个连续变量之间的关系
#' @param  data data.frame 提供数据用于绘图 core_dimension,inter_dimension两个参数是data中对应的列
#' @param  x_dimension 字符 核心关注的维度1（data中的列名）
#' @param  y_dimension 字符 核心关注的维度2（data中的列名）
#' @param  plot 逻辑 是否绘制图像
#' @param  color_dimension 字符/因子/连续的数值 额外维度，在图中用不同颜色区分（data中的列名）
#' @param  alpha_dimension 数值 额外维度，在图中用不同的透明度区分（data中的列名）
#' @param  cor_test 逻辑 是否对两个连续变量进行相关性检验，如果TRUE，则进行相关性检验，返回一个list，形如list(plot=plot2,cor=res)，plot为ggplot绘图对象，cor为相关性检验的结果，如果FALSE，则只返回ggplot绘图对象
#' @param  cor_method 可选的字符组 相关性检验使用的方法，从"pearson", "kendall", "spearman"中选取
#' @param  linear_fitting 逻辑 是否进行曲线拟合
#' @param  Smoothing_method 可选的字符组 曲线拟合所使用的方法，从"auto", "lm", "glm", "gam", "loess"中选取；For method = "auto" the smoothing method is chosen based on the size of the largest group (across all panels). loess() is used for less than 1,000 observations; 
#                                         otherwise mgcv::gam() is used with formula = y ~ s(x, bs = "cs"). Somewhat anecdotally, loess gives a better appearance, but is O(N^2) in memory, so does not work for larger datasets.
#                                         If you have fewer than 1,000 observations but want to use the same gam() model that method = "auto" would use, then set method = "gam", formula = y ~ s(x, bs = "cs").
#' @param  fittingByColor 逻辑 绘制拟合曲线的时候是否对散点按照不同的color_dimension进行分组之后再拟合，如果TRUE，则将得到多条拟合曲线；如果FALSE，则只拟合总体的曲线

data_explor.continuous2continuous <- function(data,x_dimension,y_dimension,
                                              color_dimension = NULL,
											  shape_dimension=NULL,
                                              cor_test = T,cor_method=c("pearson", "kendall", "spearman"),plot=F,
                                              linear_fitting = F,Smoothing_method =c("auto", "lm", "glm", "gam", "loess"),fittingByColor=F, addHist = F){
  
colNames <- c(x_dimension,y_dimension, color_dimension, shape_dimension)
 colIndex <- which(!sapply(list(x_dimension,y_dimension, color_dimension, shape_dimension),is.null))
    
 data <- na.omit(data[,colNames])
 colnames(data) <- c("x","y","color_var",'shape_var')[colIndex]
    if(!is.null(color_dimension) && is.character(data[,"color_var"])){
      data[,"color_var"] <- as.factor(data[,"color_var"])
    }
   if(!is.null(shape_dimension) && is.character(data[,"shape_var"])){
      data[,"shape_var"] <- as.factor(data[,"shape_var"])
    }
  if(cor_test){
    res <- cor.test( ~ y + x, data=data, method=cor_method[1])
  }
  if(plot){
    library(ggplot2)
    library(ggsignif)
    library(ggbeeswarm)
    #pdf(path)
	if(is.null(color_dimension)){
	plot2<-ggplot(data=data,aes(x=x,y=y)) + geom_point() +
        ggtitle(paste(x_dimension," VS ",y_dimension,sep="")) + xlab(x_dimension) + ylab(y_dimension)
	}else{
	if(is.factor(data[,"color_var"]) & !is.null(shape_dimension)){
      plot2 <-ggplot(data=data,aes(x=x,y=y)) + geom_point(aes(color=color_var,shape = shape_var)) + 
        ggtitle(paste(x_dimension," VS ",y_dimension,sep="")) + xlab(x_dimension) + ylab(y_dimension) +
        guides(col=guide_legend(title=color_dimension), shape = guide_legend(title=shape_dimension))
      
    }else if(is.factor(data[,"color_var"])){
        plot2<-ggplot(data=data,aes(x=x,y=y)) + geom_point(aes(color=color_var)) + 
          ggtitle(paste(x_dimension," VS ",y_dimension,sep="")) + xlab(x_dimension) + ylab(y_dimension) + 
          guides(col=guide_legend(title=color_dimension))
      }else if(is.numeric(data[,"color_var"])){
        plot2<-ggplot(data=data,aes(x=x,y=y)) + geom_point(aes(color=color_var)) + 
          ggtitle(paste(x_dimension," VS ",y_dimension,sep="")) + xlab(x_dimension) + ylab(y_dimension) + scale_color_gradient(low="red") +
          guides(col=guide_legend(title=color_dimension))
      }
      

	}
    
     #guides(col=guide_legend(title=core_dimension))
    if(linear_fitting){
      if(fittingByColor){
        plot2 <- plot2 + geom_smooth(aes(color=color_var),method = Smoothing_method[1])
      }else{
        plot2 <- plot2 + geom_smooth(method = Smoothing_method[1])
      }
    }
    if(cor_test){
      plot2 <- plot2 + geom_text(aes(x= max(x),y=max(y),label=paste0("cor =" ,round(res$estimate,3),"; Pvalue = ", format(res$p.value,digits = 3,scientific = T))),vjust = "inward", hjust = "inward")
    
      }
	 if(addHist){
	  if(res$estimate>=0){
	   position <- c(0,1)
	  }else{
	    position <- c(0,0)
	  }
	 plot2<- plot2+theme(legend.justification = position , legend.position = position, legend.key=element_blank(), legend.background = element_blank())
	 plot2 <- ggExtra::ggMarginal(plot2, type = "histogram")
	 } 
    
	 
    return(list(plot=plot2,cor=res))
    
    #dev.off()
  }else{
    return(res)
  }
  
}


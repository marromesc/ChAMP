# this function needs a lot of cleaning but should work

plotDMRs <- function(Anno, max=35, only_promoter=F){
  
  # define genes colors
  genes <- as.character(Anno[,'gene'])
  genes[genes == ''] <- 'no_gene'
  genes <- unique(genes)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_gene = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_gene <- sample(col_gene, length(genes))
  names(col_gene) <- genes
  
  # define features colors
  col_feature <- c('#b3cde3', '#a6dba0', '#c2a5cf', '#008837','#4dac26', '#8c96c6', "#88419d")
  names(col_feature) <-  c("1stExon", "3'UTR", "5'UTR", "Body", "ExonBnd", "TSS1500", "TSS200")
  
  if(isTRUE(only_promoter)){
    col_feature <- col_feature[names(col_feature)%in%c('1stExon', "5'UTR", 'TSS1500', 'TSS200')]
  }
  
  # define cgi colors
  col_cgi <- brewer.pal(n = 4, name = "Set2")[1:4]
  names(col_cgi) <- c("shore", "opensea", "island", "shelf")
  
  
  for(i in 1:max){
    CpGs <- rownames(Anno[Anno$DMRindex == rownames(myDMR[[1]])[i],])
    CpGs <- CpGs[order(Anno[CpGs,'MAPINFO'])]
    
    genes <- as.character(Anno[CpGs,'gene'])
    genes[genes == ''] <- 'no_gene'
    
    select <- Anno[Anno$DMRindex==paste("DMR",i,sep="_"),]
    select <- select[order(select$MAPINFO),]
    
    Group <- split(as.data.frame(t(myLoad[rownames(select),])),as.character(p))
    
    G <- lapply(Group,function(x) t(x))
    G <- lapply(G,function(h) data.frame(Sample=rep(colnames(h),each=nrow(h)),ID=rep(rownames(h),ncol(h)),pos=rep(as.numeric(as.factor(select$MAPINFO)),ncol(h)),Value=as.vector(h)))
    for(j in 1:length(G)) G[[j]] <- data.frame(G[[j]],pheno=names(G)[j])
    X <- do.call(rbind,G)
    X <- cbind(X,showtext=paste("cpgID:",X$ID," Sample:",X$Sample,sep=""))
    
    group_col <- group_color
    names(group_col) <- c(ref.group, group2)
    
    Fit <- lapply(G,function(h) data.frame(pos=unique(h$pos),ID=unique(h$ID),Mean=aggregate(h$Value,by=list(h$pos),mean)[,2]))
    for(j in 1:length(Fit)) Fit[[j]] <- data.frame(Fit[[j]],pheno=paste(names(Fit)[j]))
    df <- data.frame(x = unlist(lapply(Fit, "[[", "pos")),
                     y = unlist(lapply(Fit, "[[", "Mean")),
                     ID = unlist(lapply(Fit, "[[", "ID")),
                     cut = unlist(lapply(Fit, "[[", "pheno")))
    
    ha <- rowAnnotation(Groups = p, col = list(Groups = group_col))
    hr <- HeatmapAnnotation(Genes = genes, CGI = Anno[CpGs, 'cgi'], Feature = Anno[CpGs, 'feature'], col = list(Genes = col_gene, CGI = col_cgi, Feature = col_feature))
    ta <- HeatmapAnnotation('average beta' = anno_lines(cbind(df$y[1:length(CpGs)],df$y[(length(CpGs)+1):nrow(df)]), ylim = c(0, 1), gp = gpar(col = rev(group_color)), height = unit(2, "cm"), pch = c(1, 16)))
    
    h <- Heatmap(t(myLoad[CpGs,]),
                 name='Beta',
                 show_row_names=F,
                 row_title='Samples',
                 cluster_rows = F, 
                 cluster_columns = F, 
                 col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
                 left_annotation = ha,
                 bottom_annotation = hr,
                 top_annotation = ta, 
                 column_title = paste0('DMR ', i),
                 # row_names_gp = grid::gpar(fontsize = 2)
    )
    print(h)
  }
}
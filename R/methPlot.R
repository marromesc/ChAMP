# this function needs a lot of cleaning but should work

methplot <- function(df, myDMR_annotated = NULL, annot = T, cap = c(-3.5,3.5) ){
  # Create input for plot
  
  data <- data.frame(location = paste0(df[,'seqnames'], ":", as.integer(df[,'start']), "-", as.integer(df[,'end'])), chr = df[,'seqnames'], start = df[,'start'], end = df[,'end'], value = df[,'value'], pos = df[,'direction'] )
  
  bin.size <- 1000
  
  bins <- QDNAseq::getBinAnnotations(binSize=bin.size)
  bins <- Biobase::pData(bins)
  
  max.chr1 <- which.max(bins[bins$chromosome == 1,"end"])
  max.chr2 <- which.max(bins[bins$chromosome == 2,"end"]) + max.chr1
  max.chr3 <- which.max(bins[bins$chromosome == 3,"end"]) + max.chr2
  max.chr4 <- which.max(bins[bins$chromosome == 4,"end"]) + max.chr3
  max.chr5 <- which.max(bins[bins$chromosome == 5,"end"]) + max.chr4
  max.chr6 <- which.max(bins[bins$chromosome == 6,"end"]) + max.chr5
  max.chr7 <- which.max(bins[bins$chromosome == 7,"end"]) + max.chr6
  max.chr8 <- which.max(bins[bins$chromosome == 8,"end"]) + max.chr7
  max.chr9 <- which.max(bins[bins$chromosome == 9,"end"]) + max.chr8
  max.chr10 <- which.max(bins[bins$chromosome == 10,"end"]) + max.chr9
  max.chr11 <- which.max(bins[bins$chromosome == 11,"end"]) + max.chr10
  max.chr12 <- which.max(bins[bins$chromosome == 12,"end"]) + max.chr11
  max.chr13 <- which.max(bins[bins$chromosome == 13,"end"]) + max.chr12
  max.chr14 <- which.max(bins[bins$chromosome == 14,"end"]) + max.chr13
  max.chr15 <- which.max(bins[bins$chromosome == 15,"end"]) + max.chr14
  max.chr16 <- which.max(bins[bins$chromosome == 16,"end"]) + max.chr15
  max.chr17 <- which.max(bins[bins$chromosome == 17,"end"]) + max.chr16
  max.chr18 <- which.max(bins[bins$chromosome == 18,"end"]) + max.chr17
  max.chr19 <- which.max(bins[bins$chromosome == 19,"end"]) + max.chr18
  max.chr20 <- which.max(bins[bins$chromosome == 20,"end"]) + max.chr19
  max.chr21 <- which.max(bins[bins$chromosome == 21,"end"]) + max.chr20
  max.chr22 <- which.max(bins[bins$chromosome == 22,"end"]) + max.chr21
  max.chr23 <- which.max(bins[bins$chromosome == "X","end"]) + max.chr22
  
  binchrend <- c(max.chr1, max.chr2, max.chr3, max.chr4, max.chr5, max.chr6, max.chr7, max.chr8, max.chr9, max.chr10, max.chr11, max.chr12, max.chr13, max.chr14, max.chr15, max.chr16, max.chr17, max.chr18, max.chr19, max.chr20, max.chr21, max.chr22, max.chr23)
  
  mid.chr1 <- which.max(bins[bins$chromosome == 1,"end"])/2
  mid.chr2 <- which.max(bins[bins$chromosome == 2,"end"])/2 + max.chr1
  mid.chr3 <- which.max(bins[bins$chromosome == 3,"end"])/2 + max.chr2
  mid.chr4 <- which.max(bins[bins$chromosome == 4,"end"])/2 + max.chr3
  mid.chr5 <- which.max(bins[bins$chromosome == 5,"end"])/2 + max.chr4
  mid.chr6 <- which.max(bins[bins$chromosome == 6,"end"])/2 + max.chr5
  mid.chr7 <- which.max(bins[bins$chromosome == 7,"end"])/2 + max.chr6
  mid.chr8 <- which.max(bins[bins$chromosome == 8,"end"])/2 + max.chr7
  mid.chr9 <- which.max(bins[bins$chromosome == 9,"end"])/2 + max.chr8
  mid.chr10 <- which.max(bins[bins$chromosome == 10,"end"])/2 + max.chr9
  mid.chr11 <- which.max(bins[bins$chromosome == 11,"end"])/2 + max.chr10
  mid.chr12 <- which.max(bins[bins$chromosome == 12,"end"])/2 + max.chr11
  mid.chr13 <- which.max(bins[bins$chromosome == 13,"end"])/2 + max.chr12
  mid.chr14 <- which.max(bins[bins$chromosome == 14,"end"])/2 + max.chr13
  mid.chr15 <- which.max(bins[bins$chromosome == 15,"end"])/2 + max.chr14
  mid.chr16 <- which.max(bins[bins$chromosome == 16,"end"])/2 + max.chr15
  mid.chr17 <- which.max(bins[bins$chromosome == 17,"end"])/2 + max.chr16
  mid.chr18 <- which.max(bins[bins$chromosome == 18,"end"])/2 + max.chr17
  mid.chr19 <- which.max(bins[bins$chromosome == 19,"end"])/2 + max.chr18
  mid.chr20 <- which.max(bins[bins$chromosome == 20,"end"])/2 + max.chr19
  mid.chr21 <- which.max(bins[bins$chromosome == 21,"end"])/2 + max.chr20
  mid.chr22 <- which.max(bins[bins$chromosome == 22,"end"])/2 + max.chr21
  mid.chr23 <- which.max(bins[bins$chromosome == "X","end"])/2 + max.chr22
  
  binchrmiddle <- c(mid.chr1, mid.chr2, mid.chr3, mid.chr4, mid.chr5, mid.chr6, mid.chr7, mid.chr8, mid.chr9, mid.chr10, mid.chr11, mid.chr12, mid.chr13, mid.chr14, mid.chr15, mid.chr16, mid.chr17, mid.chr18, mid.chr19, mid.chr20, mid.chr21, mid.chr22, mid.chr23)
  
  data$bin <- NA
  for (i in 1:nrow(data)){
    data$bin[i] <- which(bins$chromosome == gsub('chr','',data$chr[i]) & bins$start <= data$start[i] & bins$end >= data$end[i])
  }
  
  # data[data$chr == 1, "bin"] <- data$start[data$chr == 1]
  # data[data$chr == 2, "bin"] <- data$start[data$chr == 2] + max.chr1
  # data[data$chr == 3, "bin"] <- data$start[data$chr == 3] + max.chr2
  # data[data$chr == 4, "bin"] <- data$start[data$chr == 4] + max.chr3
  # data[data$chr == 5, "bin"] <- data$start[data$chr == 5] + max.chr4
  # data[data$chr == 6, "bin"] <- data$start[data$chr == 6] + max.chr5
  # data[data$chr == 7, "bin"] <- data$start[data$chr == 7] + max.chr6
  # data[data$chr == 8, "bin"] <- data$start[data$chr == 8] + max.chr7
  # data[data$chr == 9, "bin"] <- data$start[data$chr == 9] + max.chr8
  # data[data$chr == 10, "bin"] <- data$start[data$chr == 10] + max.chr9
  # data[data$chr == 11, "bin"] <- data$start[data$chr == 11] + max.chr10
  # data[data$chr == 12, "bin"] <- data$start[data$chr == 12] + max.chr11
  # data[data$chr == 13, "bin"] <- data$start[data$chr == 13] + max.chr12
  # data[data$chr == 14, "bin"] <- data$start[data$chr == 14] + max.chr13
  # data[data$chr == 15, "bin"] <- data$start[data$chr == 15] + max.chr14
  # data[data$chr == 16, "bin"] <- data$start[data$chr == 16] + max.chr15
  # data[data$chr == 17, "bin"] <- data$start[data$chr == 17] + max.chr16
  # data[data$chr == 18, "bin"] <- data$start[data$chr == 18] + max.chr17
  # data[data$chr == 19, "bin"] <- data$start[data$chr == 19] + max.chr18
  # data[data$chr == 20, "bin"] <- data$start[data$chr == 20] + max.chr19
  # data[data$chr == 21, "bin"] <- data$start[data$chr == 21] + max.chr20
  # data[data$chr == 22, "bin"] <- data$start[data$chr == 22] + max.chr21
  # data[data$chr == "X", "bin"] <- data$start[data$chr == "X"] + max.chr22
  # 
  template <- data
  template <- template[
    with(template, order(chr, start)),
  ]
  
  p <- ggplot2::ggplot(template, aes(x=bin, y=value, color = pos)) +
    scale_y_continuous(name = "Meth Coef", limits = c(cap[1],cap[2]), breaks = seq(cap[1], cap[2], 0.5), expand=c(0.0000001,0.0000001) ) +
    scale_x_continuous(name = "Chromosome", 
                       #limits = c(0, tail(binchrend,1)), 
                       breaks = binchrmiddle, labels = c(1:22, "X"), expand = c(0,0)) +
    
    geom_vline(xintercept = binchrend, color = "#666666", linetype = "solid", size=0.25) +
    geom_hline(yintercept = 0, color = '#666666', size = 1.3) +
    geom_bar(stat='identity', width = 3) +
    scale_color_manual(values=c("red", "blue")) +
    
    theme_classic() +
    theme(
      axis.line = element_line(color='black'),
      axis.ticks = element_line(color='black'),
      axis.text.y = element_text(color='black', size = 12),
      axis.text.x = element_text(color="black", size = 12),
      axis.title.y = element_text(color = "black", face = "bold", size = 15),
      axis.title.x=element_text(color = "black", face = "bold", size = 15),
      panel.border = element_rect(colour = "black", fill=NA), 
      plot.title = element_text(color="black", face="bold", size = 15)
    ) 
  
  if (isTRUE(annot) & !is.null(myDMR_annotated)){
    template$gene <- NA
    for(i in 1:nrow(template)){
      template$gene[i] <- paste(unique(myDMR_annotated[myDMR_annotated$CHR == as.character(template$chr[i]) & myDMR_annotated$MAPINFO >= as.numeric(template$start[i]) & myDMR_annotated$MAPINFO <= as.numeric(template$end[i]),'gene']), collapse = ',')
    }
    
    p <- p + annotate("text", x = template$bin[template$pos == 'hypermethylated'], y = template$value[template$pos == 'hypermethylated']+(nchar(template$gene[template$pos == 'hypermethylated'])*0.02), label = template$gene[template$pos == 'hypermethylated'], size = 3, angle = 90)+
      annotate("text", x = template$bin[template$pos == 'hypomethylated'], y = template$value[template$pos == 'hypomethylated']-(nchar(template$gene[template$pos == 'hypomethylated'])*0.02), label = template$gene[template$pos == 'hypomethylated'], size = 3, angle = 90)
  }
  
  return(p)
}
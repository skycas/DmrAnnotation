### readFeatures ###

# Function: Calculate the promoter, 5'utr, exon/cds, intron, 3'utr of genes
# Author: Fuhui Xiao
# Notes: A little slow, however, you need only one time to calculate the coordinate of genes 

setGeneric("readFeatures", function(files, upstream=2000, downstream=0) standardGeneric("readFeatures"))

setMethod("readFeatures", signature(files = "character"), function(files, upstream, downstream){
  message("Reading gtf file ... \r")
  gtf <- rtracklayer::import(files,format = "gtf")
  gtf <- as.data.frame(gtf)
  gtf <- gtf[,-c(8,9,10,12,13)]
  
  all_transcripts <- unique(gtf$transcript_id)
  coding_transcript <- unique(gtf[gtf$type == "CDS",]$transcript_id)
  noncoding_transcript <- all_transcripts[!(all_transcripts %in% coding_transcript)]
  
  r_coding <- lapply(coding_transcript, function(x) {
    tr <- gtf[gtf$transcript_id == x,]
    if (tr$strand == "+") {
      cds <- tr[grepl("CDS",tr$type),]
      exon <- tr[grepl("exon",tr$type),]
      tss <- min(exon$start)
      tes <- max(exon$end)
      
      ####
      promot_start <- tss - 1 - upstream
      promot_end <- tss - 1
      promot_width <- promot_end - promot_start
      utr5_start <- min(exon$start)
      utr5_end <- min(cds$start) - 1
      utr5_width <- utr5_end - utr5_start + 1
      utr3_start <- max(cds$end) + 1
      utr3_end <- max(exon$end)
      utr3_width <- utr3_end - utr3_start + 1
      if (dim(cds)[1] == 1) {
        promoter <- data.frame(seqnames = unique(exon$seqnames), start = promot_start, end = promot_end, width = promot_width, strand = "+", source = "refGene", type = "promoter", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        utr5 <- data.frame(seqnames = unique(exon$seqnames), start = utr5_start, end = utr5_end, width = utr5_width, strand = "+", source = "refGene", type = "utr5", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        utr3 <- data.frame(seqnames = unique(exon$seqnames), start = utr3_start, end = utr3_end, width = utr3_width, strand = "+", source = "refGene", type = "utr3", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        gene_anno = rbind(promoter,utr5,cds,utr3)
        gene_anno = gene_anno[order(gene_anno$start),]
      }else{
        intr_start <- cds$end + 1
        intr_start <- intr_start[-length(intr_start)]
        intr_end <- cds$start - 1
        intr_end <- intr_end[-1]
        intr_width <- intr_end - intr_start + 1

        promoter <- data.frame(seqnames = unique(exon$seqnames), start = promot_start, end = promot_end, width = promot_width, strand = "+", source = "refGene", type = "promoter", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        intron <- data.frame(seqnames = unique(exon$seqnames), start = intr_start, end = intr_end, width = intr_width, strand = "+", source = "refGene", type = "intron", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        utr5 <- data.frame(seqnames = unique(exon$seqnames), start = utr5_start, end = utr5_end, width = utr5_width, strand = "+", source = "refGene", type = "utr5", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        utr3 <- data.frame(seqnames = unique(exon$seqnames), start = utr3_start, end = utr3_end, width = utr3_width, strand = "+", source = "refGene", type = "utr3", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
      
        gene_anno = rbind(promoter,utr5,cds,intron,utr3)
        gene_anno = gene_anno[order(gene_anno$start),]
      }
        
    }else{
      cds <- tr[grepl("CDS",tr$type),]
      exon <- tr[grepl("exon",tr$type),]
      tss <- max(exon$end)
      tes <- min(exon$start)
      
      ####
      promot_start <- tss + 1
      promot_end <- tss + 1 + upstream
      promot_width <- promot_end - promot_start
      utr5_start <- max(cds$end) + 1
      utr5_end <- max(exon$end) 
      utr5_width <- utr5_end - utr5_start + 1
      utr3_start <- min(exon$start)
      utr3_end <- min(cds$start) - 1
      utr3_width <- utr3_end - utr3_start + 1
      if (dim(cds)[1] == 1) {
        promoter <- data.frame(seqnames = unique(exon$seqnames), start = promot_start, end = promot_end, width = promot_width, strand = "-", source = "refGene", type = "promoter", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        utr5 <- data.frame(seqnames = unique(exon$seqnames), start = utr5_start, end = utr5_end, width = utr5_width, strand = "-", source = "refGene", type = "utr5", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        utr3 <- data.frame(seqnames = unique(exon$seqnames), start = utr3_start, end = utr3_end, width = utr3_width, strand = "-", source = "refGene", type = "utr3", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        gene_anno = rbind(promoter,utr5,cds,utr3)
        gene_anno = gene_anno[order(gene_anno$start),]
      }else{
        intr_start <- cds$end + 1
        intr_start <- intr_start[-length(intr_start)]
        intr_end <- cds$start - 1
        intr_end <- intr_end[-1]
        intr_width <- intr_end - intr_start + 1

        promoter <- data.frame(seqnames = unique(exon$seqnames), start = promot_start, end = promot_end, width = promot_width, strand = "-", source = "refGene", type = "promoter", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        intron <- data.frame(seqnames = unique(exon$seqnames), start = intr_start, end = intr_end, width = intr_width, strand = "-", source = "refGene", type = "intron", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        utr5 <- data.frame(seqnames = unique(exon$seqnames), start = utr5_start, end = utr5_end, width = utr5_width, strand = "-", source = "refGene", type = "utr5", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        utr3 <- data.frame(seqnames = unique(exon$seqnames), start = utr3_start, end = utr3_end, width = utr3_width, strand = "-", source = "refGene", type = "utr3", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
      
        gene_anno = rbind(promoter,utr5,cds,intron,utr3)
        gene_anno = gene_anno[order(gene_anno$start),]
      }
    }
    gene_anno = cbind(gene_anno,"coding")
    colnames(gene_anno)[10] = "coding/noncoding"
    return(gene_anno)
    })
    
  r_noncoding <- lapply(noncoding_transcript,function(x){
    tr <- gtf[gtf$transcript_id == x,]
    if (unique(tr$strand) == "+") {
      exon <- tr[grepl("exon",tr$type),]
      tss <- min(exon$start)
      tes <- max(exon$end)
      promot_start <- tss - 1 - upstream
      promot_end <- tss - 1
      promot_width <- promot_end - promot_start
      if (dim(exon)[1] == 1) {
        promoter <- data.frame(seqnames = unique(exon$seqnames), start = promot_start, end = promot_end, width = promot_width, strand = "+", source = "refGene", type = "promoter", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        gene_anno = rbind(promoter,exon)
        gene_anno = gene_anno[order(gene_anno$start),]
      }else{
        intr_start <- exon$end + 1
        intr_start <- intr_start[-length(intr_start)]
        intr_end <- exon$start - 1
        intr_end <- intr_end[-1]
        intr_width <- intr_end - intr_start + 1
        promoter <- data.frame(seqnames = unique(exon$seqnames), start = promot_start, end = promot_end, width = promot_width, strand = "+", source = "refGene", type = "promoter", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        intron <- data.frame(seqnames = unique(exon$seqnames), start = intr_start, end = intr_end, width = intr_width, strand = "+", source = "refGene", type = "intron", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
      
        gene_anno = rbind(promoter,exon,intron)
        gene_anno = gene_anno[order(gene_anno$start),]
      }
      
    }else{
      exon <- tr[grepl("exon",tr$type),]
      tss <- max(exon$end)
      tes <- min(exon$start)
      promot_start <- tss + 1
      promot_end <- tss + 1 + upstream
      promot_width <- promot_end - promot_start
      if (dim(exon)[1] == 1) {
        promoter <- data.frame(seqnames = unique(exon$seqnames), start = promot_start, end = promot_end, width = promot_width, strand = "-", source = "refGene", type = "promoter", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        gene_anno = rbind(promoter,exon)
        gene_anno = gene_anno[order(gene_anno$start),] 
      }else{
        intr_start <- exon$end + 1
        intr_start <- intr_start[-length(intr_start)]
        intr_end <- exon$start - 1
        intr_end <- intr_end[-1]
        intr_width <- intr_end - intr_start + 1
        promoter <- data.frame(seqnames = unique(exon$seqnames), start = promot_start, end = promot_end, width = promot_width, strand = "-", source = "refGene", type = "promoter", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
        intron <- data.frame(seqnames = unique(exon$seqnames), start = intr_start, end = intr_end, width = intr_width, strand = "-", source = "refGene", type = "intron", transcript_id = unique(exon$transcript_id), gene_name = unique(exon$gene_name))
      
        gene_anno = rbind(promoter,exon,intron)
        gene_anno = gene_anno[order(gene_anno$start),] 
      }
    }
    gene_anno = cbind(gene_anno,"noncoding")
    colnames(gene_anno)[10] = "coding/noncoding"
    return(gene_anno)
  })
  
  r = c(r_coding,r_noncoding)
  return(r)
  
})

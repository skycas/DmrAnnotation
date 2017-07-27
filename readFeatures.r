### readFeatures ###

setGeneric("readFeatures", function(files, upstream=2000, downstream=0) standardGeneric("readFeatures"))

setMethod("readFeatures", signature(files = "character"), function(files, upstream, downstream){
  message("Reading gtf file ... \r")
  gtf <- rtracklayer::import(files,format = "gtf")
  gtfx <- as.data.frame(gtf)
  
  message("Transcripts Features, Please wait patiently ... \r")
  transcript_id <- unique(gtf$transcript_id)
  bed <- suppressWarnings(lapply(transcript_id, function(x) {gtf[gtf$transcript_id == x,]}))
  names(bed) <- transcript_id
  
  message("Exon coordinates ... \r")
  bed <- lapply(bed, function(x){x[x$type == "exon",]})
  exons <- gtfx[gtfx$type == "exon", ]
  
  message("Intron coordinates ... \r")
  intron <- lapply(bed, function(x) {
    x = as.data.frame(x)
    n = dim(x)[1]
    if (n == 1) {
      tmp <- c()
    }else{
      t1 = x$end
      t1 = t1[-n]
      t2 = x$start
      t2 = t2[-1]
      start = t1 + 1
      end = t2 - 1
      width = end - start + 1
      tmp = x[-n,]
      tmp$start = start
      tmp$end = end
      tmp$width = width
      tmp$type = "intron"
      tmp$exon_id = "NA"
    }
    return(tmp)
  })
  intron = intron[-which(intron == "NULL")]
  
  message("Promoter coordinates ... \r")
  promoter <- lapply(bed,function(x){
    x = as.data.frame(x)
    if (x$strand == "+") {
      start = min(x$start) - upstream + downstream
      end = min(x$start) - 1
    }else{
      start = max(x$end) + 1
      end = max(x$end) + upstream - downstream
    }
    tmp = x[1,]
    tmp$start = start
    tmp$end = end
    tmp$width = upstream
    tmp$type = "promoter"
    tmp$exon_id = "NA"
    tmp$exon_number = "NA"
    return(tmp)
  })
  
  introns <- do.call(rbind,intron)
  promoters <- do.call(rbind,promoter)
  features <- suppressWarnings(rbind(promoters,exons,introns))
  
})

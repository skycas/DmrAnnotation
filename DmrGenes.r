
dmrfile <- read.table(files)
dmr <- with(dmrfile,Granges(chr,IRanges(start,end),strand=strand,pvalue=pvalue,qvalue=qvalue,meth.diff=meth.diff)

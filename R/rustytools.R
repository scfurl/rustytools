#' Get consensus
#' @description It will take a string and print the coordinates of every repeated character substring made up of either "ACGT".
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRaanges
#' @importFrom S4Vectors DataFrame
#' @export
#'


get_consensus<-function(str, add_int){
  getconsensus(str, add_int)
}


find_cuts<-function(len, n){
  cut<-ceiling(len/n)
  c(seq(1, len, by = cut), len)
}


split_string<-function(string, n){
  cuts<-find_cuts(nchar(string), n)
  string<-sapply(1:(length(cuts)-1), function(n){
    substring(string, cuts[n], cuts[n+1])
  })
  return(list(cuts, string))
}

#' Get consensus from fasta
#' @description This function will process a fasta file that has been processed using a tool such as samtools consensus.  The function will find all regions (subsequences) of the fasta that have sequence (A, C, T, G) and output a granges object containing the features and their sequences.
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRaanges
#' @importFrom S4Vectors DataFrame
#' @export
#'

get_consensus<-function(fasta, splits = 10000, cores=1, genome="hg38"){
  if(!file.exists(fasta)) {stop("Cannot find fasta")}
  message(paste0("Reading fasta: ", fasta, "..."))
  fa<-read.fasta(file = fasta,
                 seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE,
                 set.attributes = TRUE)
  message("Fasta ingestion complete")
  message("Building sequence info")
  seqinfo<-GenomeInfoDb::Seqinfo(names(fa), seqlengths=sapply(1:length(fa), function(n) nchar(fa[[n]][1])), isCircular=NA, genome=genome)
  # x<-fa[[1]][1]
  allres<-lapply(1:length(fa), function(i){
    message(paste0("Processing ", nchar(fa[[i]][1]), " bases from ", names(fa)[i]))
    xs<-split_string(fa[[i]][1], splits)
    res<-pbmcapply::pbmclapply(1:length(xs[[1]]), function(j){
      getconsensus( xs[[2]][j], xs[[1]][j])
    }, mc.cores = cores)
    if(length(res)>0){
      df<-t(data.frame(strsplit(unlist(res), "_")))
      GenomicRanges::GRanges(seqnames = rep(names(fa)[i], nrow(df)),
                             ranges = IRanges::IRanges(start = as.numeric(df[,1]), end = as.numeric(df[,2])),
                             mcols = S4Vectors::DataFrame(sequence=df[,3]), seqinfo = seqinfo)
    }
  })
  # suppressWarnings(unlist(as(allres, "GRangesList")))
  unlist(as(allres, "GRangesList"))
}

#
#
#
#
#
#
#
#
#
#
#
#
#
# substring(x, 14447574, 14447592)
# substring(x, 2244091, 2244116)
# #
# # Start: 14447574; End: 14447584; Seq: TGCTTCTATCA
# # Start: 14447586; End: 14447592; Seq: GGTCTGT
#
# substring(x, 832501, 832522)
# # "832502_832521_AGACCATGAAGGTCCATGCT"

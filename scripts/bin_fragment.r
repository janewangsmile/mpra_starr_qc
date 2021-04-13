suppressPackageStartupMessages({
    library(Rsamtools)
    library(GenomicAlignments)
    library(GenomicFiles)
    library(optparse)
    library(data.table)
    library(optparse)
    library(rtracklayer)
    library(GenomicRanges)
    library(RColorBrewer)
    library(doParallel)
    library(argparse)
});

# create parser object and add parser arguments
parser <- ArgumentParser()

parser$add_argument("-i", "--input", nargs="+", help="Input count tables at fragment-level")
#parser$add_argument("-p", "--cores", default=4, help="Number of cores to use in multicore processing")
parser$add_argument("--csaw", required=F, action="store_true", default=FALSE, help="Specify if CSAW package will be used for enhancer call")
parser$add_argument("--starrpeaker", required=F, action='store_true', default=FALSE, help="Specify if STARRPeaker package will be used for enhancer call")
parser$add_argument("--bin", required=F, action='store_true', default=FALSE, help="Specify if fragments need to be binned")
parser$add_argument('-s', "--size", required=F, default=50, type="integer", help="Specify bin size")
parser$add_argument('--o', required=F, action="store_true", default=FALSE, help="Specify if orientation needs to be separated")
parser$add_argument("-o", "--outfile", required=F, help="Output filename")
parser$add_argument("--minFragmentLength", required=F, default = -1, type="integer", help="Specify the minimum size of fragments")
parser$add_argument("--maxFragmentLength", required=F, default = -1, type="integer", help="Specify the maximum size of fragments")


# this function will be used to merge counts from replicates
merge_counts = function( x, y, name ) {
    if(length(x)) {
        hits = findOverlaps(x, y, type="equal");
        mcols(x)[hits@from, name] = mcols(x)[hits@from, name] + mcols(y)[hits@to, name];
        y = y[-hits@to,];
        for( cn in colnames(mcols(x))) {
            if( !cn %in% colnames(mcols(y)) ) {
                mcols(y)[,cn] = 0;
            }
        }
        mcols(y) = mcols(y)[,colnames(mcols(x))];
        colnames(mcols(y)) = colnames(mcols(x));
    }
    return(append(x, y));
}

# this function will merge (add) metadata columns for
# identical reads. used to bin data (if desired)
collapse_reads = function( x ) {
  u = x[!duplicated(x)];
  x = x[ duplicated(x)];
  gc();
  for( cn in 1:ncol(mcols(u)) ) {
    z = x[ mcols(x)[,cn]>0 ];
    hits = findOverlaps( u, z, type="equal" );
    if( any(duplicated(hits@from)) ) {
        hits = aggregate( mcols(z)[hits@to,cn] ~ hits@from, FUN='sum' );
        mcols(u)[hits[,1],cn] = mcols(u)[hits[,1],cn] + hits[,2];
    } else {
        mcols(u)[hits@from,cn] = mcols(u)[hits@from,cn] + mcols(z)[hits@to,cn];
    }
  }
  return(u);
}

bin_fragments = function(count, binSize=50){
    
    start(count) = binSize * as.integer(start(count)/ binSize)+1
    end(count) = binSize * as.integer(end(count)/ binSize)+1
    count = collapse_reads(count)
    return(count)
    
}

args <- parser$parse_args()

# Merging count from replicates 
frag_count = GRanges()

for( i in 1:length(args$input) ){
    x = read.csv(args$input[i], sep='\t')
    x = x[, c('seqnames', 'start', 'end', 'count')]
    
    # set sample name
    sample = strsplit(args$input[i], "/")[[1]]
    sample = strsplit(sample[length(sample)], ".csv")[[1]]
    
    # Convert dataframe to GRanges
    x = makeGRangesFromDataFrame(x, 
                               keep.extra.columns=TRUE,
                               ignore.strand=FALSE,
                               seqinfo=NULL,seqnames.field=c("seqnames", "seqname","chromosome", "chrom","Chr", "chromosome_name","seqid"),
                               start.field="Start",
                               end.field=c("End", "stop"),
                               strand.field="Strand",
                               starts.in.df.are.0based=FALSE)
    colnames(mcols(x)) = c(sample)
    
    # Filter fragments based on size
    if(args$minFragmentLength != -1 & args$maxFragmentLength != -1){
        x = x[width(x) >= args$minFragmentLength & width(x) <= args$maxFragmentLength]
    }
    
    mcols(frag_count)[,sample] = 0;
    frag_count = merge_counts(frag_count, x, sample);

}

# Bin fragments

if(args$bin){
    binned_frag_count = bin_fragments(frag_count, binSize=args$size)
#     # Save
#     write.table(binned_frag_count, file=paste0(args$outfile), quote=FALSE, sep='\t');
    
    # Convert binned fragments count GRanges to dataframe for each replicate 
    for( i in 1:length(colnames(mcols(binned_frag_count)))){
        sample = colnames(mcols(binned_frag_count))[i]
        seqlib = binned_frag_count$sample
        
        seqlib = as.data.frame(seqlib)
        cols_to_add = data.frame(name=paste0(seqlib$seqnames,"_", seqlib$start,"_",seqlib$end), 
                         score = pmin(seqlib$count,1000), 
                         barcode = '.')
        
        seqlib_df = cbind(seqlib, cols_to_add)
        col_order = c('seqnames','start','end','name','score','strand','count','barcode')
        seqlib_df = selib_df[, col_order]
        
        # Save
        write.table(seqlib_df, file=paste0(args$outfile[i]), quote=FALSE, sep='\t');
    }
}





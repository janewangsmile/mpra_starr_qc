#!/usr/bin/env Rscript

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

parser$add_argument("-i", "--input", nargs="+", help="BAM files with reads")
#parser$add_argument("-p", "--cores", default=4, help="Number of cores to use in multicore processing")
parser$add_argument("--UMI", required=F, action="store_true", default=FALSE, help="Specify if UMI/barcode present for fragments")
parser$add_argument("-o", "--outfile", required=F, help="Output filename to create bed file.")

args <- parser$parse_args()


registerDoParallel(cores=args$cores);

# Shorthand for "pretty number" formatting
pn = function(value) {
    prettyNum(value, big.mark=",")
}

# Shorthand to print the full argument list
msgout = function(...) {
    write(paste(...), stdout());
}

if (!(args$UMI)){
    bfilters = ScanBamParam(mapqFilter=10, flag=scanBamFlag(isSecondaryAlignment=F));
} else {
    bfilters = ScanBamParam(mapqFilter=10, flag=scanBamFlag(isSecondaryAlignment=F, isDuplicate=F));
}

count_reads = function(reads) {
    uniq = unique(reads);
    # sum over duplicates to get a count for each unique 5'/3' end
    uniq$count = countOverlaps( uniq, reads, type="equal" );
    return( uniq );
}

yield.bam = function(X) {
    y = GRanges( readGAlignmentPairs(X, use.names=F, param=bfilters ));
    return(y);
}

map.bam = function(X) {
    return(X);
}

reduce.bam = function(x, y) {
    x = append(x, y);
    # print the number of readpairs processed
    msgout(pn(length(x)), 'mapped human reads');
    return(x);
}

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


ctFile = args$input;
msgout( "Processing ", ctFile );
infile = BamFile(ctFile, yieldSize=1 * 10^6, asMates=T );
aligned = reduceByYield( infile, yield.bam, map.bam, reduce.bam, parallel=F );

#msgout(pn(length(aligned)), 'mapped reads');

# compute coverage from identical reads => 'count' column
seqlib = count_reads(aligned);
print(seqlib)

# convert the GRanges object to dataframe and add columns corresponding to name, score and barcode information as per the common file format needed

seqlib = as.data.frame(seqlib)
cols_to_add = data.frame(name=paste0(seqlib$seqnames,"_", seqlib$start,"_",seqlib$end), 
                         score = pmin(seqlib$count,1000), 
                         barcode = '.')

seqlib_df = cbind(seqlib, cols_to_add)
col_order = c('seqnames','start','end','name','score','strand','count','barcode')
seqlib_df = seqlib_df[, col_order]

# write the dataframe into output format required

write.table(seqlib_df, file=paste0(args$outfile), quote=FALSE, sep='\t');

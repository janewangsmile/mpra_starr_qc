#!/usr/bin/env Rscript
library(tidyverse,quietly = T)
library(GGally,quietly = T)
theme_set(theme_bw())


main_wrap = function(raw_DNA_folder, raw_RNA_folder, threshold.reads=-Inf, 
                     method="MPRA", output.name) {
  
  # consider using 'argparse' for better input handing
  
  if(method!="MPRA") stop("No current implementation for non-MPRA experiments!")
  
  ### glossary:
  # read - numbers in column 7, usually number of raw sequencing reads
  # unit - column 8 (if provided, currently ignored), UMIs/regions/barcodes
  # element - column 4 (aka sequence, oligo, insert, tile...)
  
  # find files
  if(missing("output.name")) stop("Argument 'output.name' is missing!")
  raw_DNA_files = list.files(raw_DNA_folder, pattern = '.bed(.gz)?',full.names = T)
  raw_RNA_files = list.files(raw_RNA_folder, pattern = '.bed(.gz)?',full.names = T)
  
  ### check number of samples
  DNA_len = length(raw_DNA_files)
  RNA_len = length(raw_RNA_files)
  if(RNA_len==0 | DNA_len ==0) stop("Couldn't find DNA or RNA files (BED or BED.GZ)")
  if(RNA_len==DNA_len) {
    message("Found ", DNA_len, " DNA and RNA replicates.") } else {
      if(RNA_len>DNA_len & DNA_len==1) {
        message("Found only one DNA reference for ", RNA_len, " RNA replicates.") } else {
          stop("Found ", list.files(raw_DNA_folder, pattern = '.bed(.gz)?'),
               " DNA samples and ", list.files(raw_RNA_folder, pattern = '.bed(.gz)?'), 
               " RNA samples. That's unacceptable, aborting.")
        }
    }
  
  ### make output folder
  output.dir = file.path("report",output.name)
  dir.create(output.dir,showWarnings = F,recursive = T)
  # round(file.size(raw_DNA_files[replicate])/1024^2,1) # file size in MB, if needed
  
  ### extract stats
  stats = list(); ratios=list()
  recalculate_threshold = is.infinite(threshold.reads)
  for(replicate in seq_len(RNA_len)) {
    
    message("Replicate ",replicate)
    
    # for special case of DNA_len==1 and RNA_len>1, do not read DNA multiple times
    if(!(DNA_len==1 & exists("already_read_dna"))) {
      
      if(method == "STARR") {
        # unwritten function to bin STARR regions
        # ideally, should only change the 'name' (4) column
        # consider binning without reading whole table into memory
        STARR_handler(raw_DNA_files[replicate])
        DNA_file = name_of_STARR_handler_file
      } else {
        # for MPRA just read in the raw_mpra format
        DNA_file = raw_DNA_files[replicate]
      }
      
      # consider data.table::fread for large files
      dna = read_tsv(DNA_file,
                     col_names = c('chr','start','end','name','score','strand',
                                   'DNA','barcode'), 
                     # skip header, if it's there
                     skip = as.integer(grepl("start|end|strand",
                                             readLines(raw_DNA_files[replicate],n = 1))),
                     col_types = "ciicicic") %>% 
        group_by(name) %>% 
        summarise(DNA.reads = sum(DNA),
                  DNA.units = n())
      
      ### thresholding: currently only total DNA reads, not units (ie UMI, regions etc)
      # if no threshold provided, use Q1-1.5*IQR
      if(recalculate_threshold) {
        threshold.reads = as.integer(10^(quantile(log10(dna$DNA.reads), probs = c(0.25)) - 
                                     1.5*IQR(log10(dna$DNA.reads))))
      }
      dna_filter = dna %>% filter(DNA.reads>=threshold.reads)
      already_read_dna = T
      
      ### histogram of DNA reads per element
      hist = dna %>%
        select(DNA.reads,name) %>% # units currently ignored
        gather('key','value',-name) %>% # rewrite as pivot_longer
        ggplot(aes(value))+
        geom_histogram(bins=50) +
        geom_vline(data = data.frame(xintercept = c(threshold.reads),
                                     key = c("DNA.reads")), 
                    aes(xintercept=xintercept), color="red") +
        facet_grid(~key) +
        scale_x_log10() +
        labs(title = paste0("DNA replicate ", replicate))+
        xlab('DNA reads or units per element')
      ggsave(file.path(output.dir, paste0("DNA_rep",replicate,"_hist_plot.png")),
             hist, w=4,h=2)
      
    }
    
    rna = read_tsv(raw_RNA_files[replicate],
                   col_names = c('chr','start','end','name','score','strand',
                                 'RNA','barcode'), 
                   # skip header, if it's there
                   skip = as.integer(grepl("start|end|strand",
                                           readLines(raw_RNA_files[replicate],n = 1))),
                   col_types = "ciicicic") %>% 
      group_by(name) %>% 
      summarise(RNA.reads = sum(RNA),
                RNA.units =n())
    
    ### calculate ratios (assuming "DNA is right", ie left_join)
    ratios[[replicate]] = dna_filter %>%
      left_join(rna,by="name") %>%
      mutate(DNA.norm = log2( (DNA.reads+1)/sum(DNA.reads)*1e6 ),
             RNA.norm = log2( (RNA.reads+1)/sum(RNA.reads)*1e6 ),
             log.ratio = RNA.norm - DNA.norm,
             replicate)
    
    ### histogram of RNA reads per element
    hist2 = ratios[[replicate]] %>%
      select(RNA.reads, name) %>% 
      replace_na(list(RNA.reads=1)) %>% 
      gather('key','value',-name) %>% 
      ggplot(aes(value)) +
      geom_histogram(bins=50) +
      facet_grid(~key) +
      scale_x_log10() +
      labs(title = paste0("RNA replicate ", replicate))+
      xlab('RNA reads or units per element')
    ggsave(file.path(output.dir, paste0("RNA_rep",replicate,"_hist_plot.png")),
           hist2, w=4,h=2)
    
    stats[[replicate]] = tibble(output.name, 
                                replicate, 
                                dna_read_threshold = threshold.reads,
                                saturation_rna = round(1-sum(rna$RNA.units)/sum(rna$RNA.reads),2), 
                                saturation_dna = round(1-sum(dna$DNA.units)/sum(dna$DNA.reads),2), 
                                dna_elements_total = nrow(dna), # number of elements
                                dna_elements_filtered = nrow(dna_filter), # number of elements with reads above threshold
                                dna_reads_total = sum(dna$DNA.reads), # number of reads
                                dna_reads_filtered = sum(dna_filter$DNA.reads), # numbers of reads in filtered elements
                                dna_units_filtered = sum(dna_filter$DNA.units), # number of units, e.g. UMIs
                                med_units_per_elem_filtered = median(dna_filter$DNA.units))
  }
  
  # write result to files
  write_tsv(bind_rows(stats),file.path(output.dir,"complete_stats.tsv"))
  complete_ratios = bind_rows(ratios) %>% 
    mutate(label = paste0("rep",replicate))
  write_tsv(complete_ratios,gzfile(file.path(output.dir,"complete_ratios.tsv.gz")))
  

  ### make correlation plots (RNA,DNA,ratios)
  # suppressing warnings because of pairwise incomplete observations
  pl1 = complete_ratios %>% 
    select(label,RNA.norm,name) %>% 
    spread(label,RNA.norm) %>% 
    select(-name) %>% 
    ggpairs(title = "DNA RPM correlation")
  suppressWarnings(ggsave(file.path(output.dir, "cor_RNA.png"),pl1,w=7,h=7))
  
  pl2 = complete_ratios %>% 
    select(label,DNA.norm,name) %>% 
    spread(label,DNA.norm) %>% 
    select(-name) %>% 
    ggpairs(title = "RNA RPM correlation")
  suppressWarnings(ggsave(file.path(output.dir, "cor_DNA.png"),pl2,w=7,h=7))
  
  pl3 = complete_ratios %>% 
    select(label,log.ratio,name) %>% 
    spread(label,log.ratio) %>% 
    select(-name) %>% 
    ggpairs(title = "ratio correlation")
  suppressWarnings(ggsave(file.path(output.dir, "cor_ratio.png"),pl3,w=7,h=7))
  
}  



### run from script:
# main_wrap(raw_DNA_folder = c("data/processed/hepg2/DNA/"),
#           raw_RNA_folder = c("data/processed/hepg2/RNA/"),
#           output.name = "hepg2_normal")
# main_wrap(raw_DNA_folder = c("data/processed/hepg2/DNA/"),
#           raw_RNA_folder = c("data/processed/hepg2/RNA/"),
#           output.name = "hepg2_normal", method = "STARR") # not implemented
# main_wrap(raw_DNA_folder = c("data/processed/k562/DNA/"),
#           raw_RNA_folder = c("data/processed/k562/RNA/"),
#           output.name = "k562_normal")

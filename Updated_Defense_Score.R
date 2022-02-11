library(tidyverse)
library(readxl)
library(genbankr)
library(Biostrings)
library(rhmmer)


# Neighrborhood Defense Score ---------------------------------------------
neighborhood_defense_function <-
  function(neighborhood.file,
           upstream.window.size,
           downstream.window.size,
           output.loc) {
    defense_pfams <-
      read_xlsx(
        ##path to Supplemental Table S4
      ) %>%
      select(Family, 'Annotation' = 'Annotation / Comments', 'System' = 'Defense system type') %>%
      filter(!str_detect(Family, 'COG')) %>%
      dplyr::mutate(Family = str_replace(Family, 'pfam', 'PF'))
    
    neighborhood.hits <- read_csv(neighborhood.file) %>%
      mutate(family = str_replace_all(family, '-', ', ')) %>%
      mutate(family = str_replace_all(family, 'PF', 'pfam')) %>%
      mutate(family = strsplit(as.character(family), ", ")) %>%
      unnest(family) %>%
      dplyr::select(
        id,
        gene_key,
        family,
        ipro_family,
        desc,
        accession,
        rel_start,
        rel_stop,
        num,
        sort_key,
        taxon_id,
        direction
      ) %>%
      # unique() %>%
      group_by(gene_key) %>%
      arrange(gene_key) %>%
      mutate(rel.diff = rel_stop + lag(rel_stop))
    
    target.location <- neighborhood.hits %>%
      dplyr::slice(which.max(rel.diff > 0))
    
    ind <-
      which(neighborhood.hits$sort_key %in% target.location$sort_key)
    
    upstream.window <- 1:upstream.window.size
    downstream.window <- 1:downstream.window.size
    sliced.hits <-
      sort(Reduce(c, c(
        lapply(ind, function(x)
          x - upstream.window),
        ind,
        lapply(ind, function(x)
          x + downstream.window)
      )))
    replaced.sliced.hits <-
      replace(sliced.hits, which(sliced.hits <= 0), NA)
    cleared.sliced.hits <-
      replaced.sliced.hits[!is.na(replaced.sliced.hits)]
    
    ##removing pfam that is not presence in pfam dataset
    pfam_counts <- neighborhood.hits %>%
      dplyr::slice(cleared.sliced.hits) %>%
      dplyr::select(
        'Acc_Num' = id,
        'Family' = family,
        desc,
        accession,
        sort_key,
        rel.diff,
        direction
      ) %>%
      mutate(Family = str_replace(Family, 'pfam', 'PF'))
    
    ##sample of 50 genomes from list of genes with DnmA homologues
    sampled_genome_acc <- sample(unique(pfam_counts$Acc_Num), 50)
    write_tsv(
      tibble(acc_num = sampled_genome_acc),
      paste0(output.loc, 'Bs_Sample_Acc.txt'),
      col_names = FALSE
    )
    
    #if pfam is in defense system pfam list, make column with 1
    hits <- pfam_counts %>%
      inner_join(defense_pfams, 'Family') %>%
      mutate(Defense_System = 1)
    
    #if pfam is not in defense system pfam list, make column with 0
    misses <- pfam_counts %>%
      anti_join(hits, 'Family') %>%
      mutate(Defense_System = 0)
    
    #merge pfam hits and misses tables and create column 'Observed'
    dnmA_hits <<- bind_rows(hits, misses) %>%
      arrange(Acc_Num, sort_key) %>%
      mutate(Type = 'Observed') %>%
      write_csv(paste0(output.loc, 'DnmA_Neighborhood_DS.csv'))
  }
#neighborhood file <- path to excel sheet from SSN output
neighborhood.file #<-
 ##path to Supplemental Table S5 (convert to CSV format first)

#window sizes = number (interger) of genes upstream and downstream of target
upstream.window.size #<- ##number of genes upstream
downstream.window.size #<- ##number of genes downstream

##output.loc = output filepath
output.loc #<-
  ##filepath to output directory
  
##run function with given inputs 
neighborhood_defense_function(neighborhood.file,
                              upstream.window.size,
                              downstream.window.size,
                              output.loc)



##Extract FASTA nucleotide sequences using sample of accession numbers from the neighborhood files ("Bs_Acc_Num.txt")
##Annotate FASTA nucleotide sequences using prokka (https://github.com/tseemann/prokka)

# Samlpe Neighborhood Generation ------------------------------------------

sample_neighborhoods_function <-
  function(genbank.filepath, output.loc) {
    ##read in genbank filepaths from Prokka annotated genomes
    genbank.files <- list.files(
      genbank.filepath,
      pattern = '*.gbk',
      all.files = TRUE,
      full.names = TRUE,
      recursive = TRUE
    )
    
    ##function to read in Genbank file using genbankR::parseGenBank and converting the features to a tibble
    genbank.function <- function(var) {
      x.dummy <- genbankr::parseGenBank(var)
      as_tibble(as.data.frame(x.dummy$FEATURES[!grepl(pattern = 'GRanges', x.dummy$FEATURES)]))
    }
    
    ##read in and appened genbank files for all genomes as data frame
    random.genomes.genbank <-
      suppressMessages(tibble(filename = genbank.files)) %>%
      mutate(Acc_Num = map(filename, function(x)
        sub(".*prokka_files/ *(.*)*/.*", "\\1", x))) %>%
      unnest(cols = c(Acc_Num)) %>%
      # mutate(Acc_Num = str_remove(Acc_Num, '_example.gbk')) %>%
      mutate(file_contents = map(filename, genbank.function)) %>%
      select(-filename) %>%
      unnest(file_contents) %>%
      unnest(2) %>%
      dplyr::rename('Gene_Name' = locus_tag) %>%
      dplyr::filter(type == 'CDS')
    write_csv(random.genomes.genbank,
              paste0(output.loc, 'DnmA_Sampled_Neighborhoods.csv'))
    
    ##Create 20 gene windows from random genes
    sampled_neighborhoods <-
      random.genomes.genbank %>%
      dplyr::select(Acc_Num, start, end, strand, Gene_Name, product, translation) %>%
      group_by(Acc_Num) %>%
      mutate(index = row_number(), start2 = sample(index, 1)) %>%
      filter(index >= start2 & index <= (start2 + 20)) %>%
      filter(n() >= 20) %>%
      ungroup()
    
    ##write out AA FASTA sequences for hmmscan analysis
    sampled_neighborhoods_seq <-
      AAStringSet(sampled_neighborhoods$translation)
    names(sampled_neighborhoods_seq) <- sampled_neighborhoods$Gene_Name
    writeXStringSet(sampled_neighborhoods_seq,
                    paste0(output.loc, 'Sample_Neighborhoods_Seq.fa'))
  }

genbank.filepath #<-
  ##path to prokka files for sample genomes
output.loc <-
  ##filepath to output directory
sample_neighborhoods_function(genbank.filepath, output.loc)




##Run amino acid sequences of all proteins in sample genome neighborhoods through hmmscan (hmmer package) using pfam-A database (http://pfam.xfam.org)


# Compiling DnmA and Random Neighborhood Scores ---------------------------

parse_hmmscan_function <-
  function(hmmer.output,
           aa.seqs,
           random.neighborhood.genbank,
           dnmA.neighborhood.genbank,
           output.loc) {
    ##read in hmmscan tblout using rhmmer package
    hmmer_output <- suppressMessages(read_tblout(hmmer.output)) %>%
      group_by(query_name) %>%
      slice_min(order_by = sequence_evalue, n = 1) %>%
      arrange(query_name) %>%
      dplyr::select(1, 2, 3, 'Annotation' = 19) %>%
      dplyr::rename(Family = 'domain_accession')
    
    sampled_neighborhoods_seq <- readAAStringSet(aa.seqs)
    random_neighborhood_genbank <- read_csv(random.neighborhood.genbank)
    dnmA_neighborhood_genbank <- read_csv(dnmA.neighborhood.genbank)
    ##which sequences were not found during hhmscan using the pfam database
    missing <-
      setdiff(names(sampled_neighborhoods_seq), hmmer_output$query_name)
    
    ##create table for missing genes with "None" in domain name and domain accession and append to hmmeroutput
    hmmer_unknown <-
      suppressMessages(
        tibble(
          domain_name = 'None',
          Family = 'None',
          query_name = missing,
          Annotation = 'None'
        )
      ) %>%
      bind_rows(hmmer_output) %>%
      arrange(query_name) %>%
      dplyr::rename('Gene_Name' = query_name)
    
    ##join PFAM ids to sample genomes and calculate number of Defense Associated PFAMs
    sample.data <-
      inner_join(random_neighborhood_genbank, hmmer_unknown, c('Gene_Name')) %>%
      dplyr::select(Acc_Num,
                    start,
                    end,
                    strand,
                    Gene_Name,
                    product,
                    domain_name,
                    Family,
                    Annotation) %>%
      mutate(Family = str_sub(Family, start = 1, end = 7)) %>%
      distinct() %>%
      mutate(Defense_Systems = case_when(
        Family %in% defense_pfams$Family ~ 1,!(Family %in% defense_pfams$Family) ~ 0
      )) %>%
      group_by(Acc_Num) %>%
      dplyr::summarize(
        size = n(),
        score = sum(Defense_Systems),
        perc = score / size
      ) %>% mutate(Type = 'Expected')
    
    ##combine defense associated PFAMs in Expected and Observed datasets
    pooled_data <- dnmA_neighborhood_genbank %>%
      group_by(Acc_Num, Type) %>%
      dplyr::summarize(
        size = n(),
        score = sum(Defense_System),
        perc = score / size
      ) %>% filter(size == 21) %>% bind_rows(sample.data)
    
    write_csv(pooled_data, paste0(output.loc,'DnmA_Random_Defense_Score.csv'))
    
    top_5_pfam <- dnmA_neighborhood_genbank %>%
      group_by(desc) %>%
      filter(Defense_System == 1) %>%
      dplyr::summarize(freq = n()) %>%
      arrange(desc(freq)) %>%
      mutate(perc = freq / sum(freq)) %>%
      dplyr::slice(1:5)
    
    write_csv(top_5_pfam, paste0(output.loc,'DnmA_Top5_Adjacent_PFAM.csv'))
  }

hmmer.output #<-
  ##path to --tblout hmmscan output file
aa.seqs <-
  ##path to 'Sample_Neighborhoods_Seq.fa' (will be located in output.loc)
random.neighborhood.genbank <-
  ##path to 'DnmA_Sampled_Neighborhoods.csv' (will be located in output.loc)
dnmA.neighborhood.genbank <-
  ##path to ''DnmA_Neighborhood_DS.csv' (will be located in output.loc)
output.loc #<- 
  ##filepath to output directory



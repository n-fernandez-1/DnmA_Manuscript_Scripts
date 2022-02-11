library(tidyverse)
library(readxl)
library(genbankr)
library(Biostrings)
library(rhmmer)
source('/Volumes/Macintosh HD/Users/nicofernandez/Rscripts/Rhmmer_source.R')
source('/Volumes/Macintosh HD/Users/nicofernandez/Rscripts/std.theme.R')


neighborhood_defense_function <-
  function(neighborhood.file,
           upstream.window.size,
           downstream.window.size,
           output.loc) {
    defense_pfams <-
      read_xlsx(
        '/Volumes/GoogleDrive/My Drive/Simmons - Bioinformatics/Coincidence/Defense_Pfams.xlsx'
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
neighborhood.file <-
  '/Volumes/GoogleDrive/My Drive/Simmons - Bioinformatics/Coincidence/Coincidence_input/Finished/Misc/DnmA_neighbors.csv'

#window sizes = number (interger) of genes upstream and downstream of target
upstream.window.size <- 10
downstream.window.size <- 10

##output.loc = output filepath
output.loc <-
  '/Volumes/GoogleDrive/My Drive/Simmons - Bioinformatics/Coincidence/Coincidence_output/02102022/'

neighborhood_defense_function(neighborhood.file,
                              upstream.window.size,
                              downstream.window.size,
                              output.loc)

##Random Genomic Island Generations
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

genbank.filepath <-
  '/Volumes/Macintosh HD/Users/nicofernandez/Documents/Genbanks_DnmA_Sample/prokka_files'
output.loc <-
  '/Volumes/GoogleDrive/My Drive/Simmons - Bioinformatics/Coincidence/Coincidence_output/02102022/'
sample_neighborhoods_function(genbank.filepath, output.loc)

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
    
    ggplot(pooled_data, aes(score, group = Type)) +
      geom_histogram(
        aes(
          y = (..count..) / sum(..count..),
          fill = Type,
          group = Type
        ),
        binwidth = 1,
        color = 'black',
        position = position_dodge(1)
      ) +
      scale_x_continuous(breaks = seq(0, 9, 1)) +
      scale_y_continuous(expand = c(0, 0)) +
      std.theme()
    # geom_bar(aes(y = (..count..)/sum(..count..), fill = Type), position = position_dodge(0.8))
    
    write_csv(pooled_data, paste0(output.loc, 'Defense_Scores_DnmA.csv'))
    
    summary_by_type <- pooled_data %>%
      group_by(Type, score) %>%
      dplyr::summarize(count = n()) %>%
      group_by(Type) %>%
      mutate(Perc = count / sum(count))
    
    dnmA_neighborhood_plot <- ggplot(summary_by_type, aes(score, Perc)) +
      geom_col(
        aes(fill = Type),
        alpha = 0.5,
        position = position_dodge(.9, preserve = 'single'),
        color = 'black'
      ) +
      scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, .61),
        breaks = seq(0, .6, .1),
        name = 'Portion of Neighborhoods'
      ) +
      scale_x_continuous(name = '# of Defense Associated Protein \n Families in Neighborhood') +
      # scale_fill_manual(values = c('grey77','grey22'))+
      scale_fill_discrete(
        name = "Data Type",
        labels = c("Random", 'DnmA'),
        type = c('grey77', 'grey22')
      ) +
      std.theme()
    
    ggsave(filename = 'DnmA_Neighborhood_Plot.jpeg',
           plot = dnmA_neighborhood_plot,device = 'jpeg',
           path =  output.loc,
           height = 10 * 1.618, 
           width = 14,
           units = 'cm',dpi = 'retina')
  }

hmmer.output <-
  '/Volumes/GoogleDrive/My Drive/Simmons - Bioinformatics/Coincidence/Coincidence_output/02102022/DnmA_Sample_HMMout.txt'
aa.seqs <-
  '/Volumes/GoogleDrive/My Drive/Simmons - Bioinformatics/Coincidence/Coincidence_output/02102022/Sample_Neighborhoods_Seq.fa'
random.neighborhood.genbank <-
  '/Volumes/GoogleDrive/My Drive/Simmons - Bioinformatics/Coincidence/Coincidence_output/02102022/DnmA_Sampled_Neighborhoods.csv'
dnmA.neighborhood.genbank <-
  '/Volumes/GoogleDrive/My Drive/Simmons - Bioinformatics/Coincidence/Coincidence_output/02102022/DnmA_Neighborhood_DS.csv'

##top_5
top_5_pfam <- dnmA_neighborhood_genbank %>%
  group_by(desc) %>%
  filter(Defense_System == 1) %>%
  dplyr::summarize(freq = n()) %>%
  arrange(desc(freq)) %>%
  mutate(perc = freq / sum(freq)) %>%
  dplyr::slice(1:5)

top_5_pfam_plot <- dnmA_neighborhood_genbank %>%
  dplyr::filter(n() >= 10) %>%
  mutate(rel.diff = as.integer(rel.diff)) %>%
  mutate(position = case_when(rel.diff < 0 ~ 'Upstream',
                              rel.diff > 0 ~ 'Downstream')) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(rel.diff)) %>%
  group_by(Acc_Num, position) %>%
  mutate(rel.pos = case_when(
    position == 'Upstream' ~ -1 * rev(1:n()),
    position == 'Downstream' ~ 1 * 1:n()
  )) %>%
  filter(Defense_System == 1) %>%
  mutate(direction = case_when(
    direction == 'complement' ~ 'Reverse',
    direction == 'normal' ~ 'Forward'
  )) %>%
  mutate(rel.pos = case_when(
    direction == 'Reverse' ~ rel.pos / -1,
    direction == 'Forward' ~ rel.pos / 1
  )) %>%
  group_by(desc) %>%
  mutate(freq = n()) %>%
  # dplyr::filter(dense_rank(plyr::desc(freq)) %in% 1:5) %>%
  mutate(desc_reordered = factor(desc, levels = unique(top_5_pfam$desc))) %>%
  dplyr::select(-gene_key,-accession,-sort_key,-System) %>%
  drop_na() %>%
  pivot_longer(cols = rel.pos,
               names_to = 'Data_Type',
               values_to = 'Position') %>%
  group_by(desc_reordered, Data_Type) %>%
  dplyr::count(Position) %>%
  mutate(perc = n / sum(n)) %>%
  dplyr::ungroup() %>% distinct()

top_5_pfam_plot %>% view()

ggplot(top_5_pfam_plot, aes(Position, perc)) +
  # geom_area(aes(fill = Data_Type), position = 'identity')+
  geom_bar(
    fill = 'black',
    color = 'black',
    stat = 'identity',
    position = 'identity'
  ) +
  # geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = FALSE, aes(fill = desc))+
  facet_grid(desc_reordered ~ ., drop = TRUE) + std.theme() + scale_x_continuous(
    expand = c(0, 0),
    limit = c(-5.5, 9.5),
    breaks = seq(-5, 9, 1)
  ) + scale_y_continuous(
    name = 'Proportion',
    expand = c(0, 0),
    breaks = seq(0, .8, .2),
    limits = c(0, 0.85)
  ) +
  theme(strip.text = element_blank())
  

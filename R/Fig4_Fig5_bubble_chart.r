# 2024.08.27 created by Shun Iida
# This program visualizes read counts of JCPyV wraparound amplicon sequecing by bubble chart


########## load library ##########

library(plyr)
library(tidyverse)
library(scales)
library(pheatmap)
library(patchwork)
library(RColorBrewer)
library(gggenes)
library(reshape2)
library(ggrepel)
library(hash)
library(ggridges)
library(UpSetR)
library(svglite)
library(ggplot2)
library(reshape2)

mutate <- dplyr::mutate


########## constants ##########

# JCPyV
genome_lengths <- c("JCPyV" = 5130)

jc_genes <- tibble(
  molecule = "JC",
  gene = c("Agnoprotein", "VP1","VP2", "VP3", "ST","LT", "LT"),
  start = c(277, 1469, 526, 883, 4495, 2603, 4771),
  end = c(492, 2533, 1560, 1560, 5013, 4426, 5013),
  strand = c("forward", "forward", "forward", "forward", "reverse", "reverse", "reverse"),
  direction = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
  exon = c(1, 1, 1, 1, 1, 1, 2)
)

virus_gene_fill = c("Agnoprotein" = "#3cb371",
                    "VP1" = "#00429d",
                    "VP2" = "#73a2c6",
                    "VP3" = "#a5d5d8",
                    "ST" = "#f4777f",
                    "LT" = "#93003a",
                    "MT" = "#9c6644")

span_colors <- c("exon" = '#e63946',
                 "intron" = "#FFDDD6")

# wraparound transcripts class
class_WA_VP1_IMR <- c("IMR-L3", "IMR-L14", "IMR-L2", "IMR-L10", "IMR-L92", "IMR-L32", "IMR-L67")
class_WA_VP1_293 <- c("293-L2", "293-L7", "293-L60", "293-L3", "293-L16", "293-L27", "293-L52")
class_WA_VP23_IMR <- c("IMR-L7", "IMR-L37", "IMR-L102", "IMR-L5", "IMR-L22", "IMR-L46", "IMR-L48")
class_WA_VP23_293 <- c("293-L10", "293-L9", "293-L14", "293-L38")
class_SuperT_IMR <- c("IMR-E14", "IMR-E15", "IMR-E16")
class_SuperT_293 <- c("293-E17", "293-E16")

# introns of SuperT transcripts NOT detected in IMR-32/293
appendix_SuperT <- tibble(
  introns = c("9403_4426", 
    "9403_4426\n14533_9556\n19663_14686\n24793_19816", 
    "9403_4426\n14533_9556\n19663_14686\n24793_19816\n29923_24946", 
    "9403_4426\n14533_9556\n19663_14686\n24793_19816\n29923_24946\n35053_30076\n40183_35206"), 
  PCR_predicted_length = c(288, 747, 900, 1206)
)

# mapped read counts determined by samtools
read_count_PML_VP1 <- tibble(
  case_ID = c(1, 2, 3, 4, 5, 6, 7, 8),
  read_count = c(32377, 66188, 42751, 43564, 26326, 39034, 33194, 39913)
)
read_count_PML_VP23 <- tibble(
  case_ID = c(1, 2, 3, 4, 5, 6, 7, 8),
  read_count = c(115812, 159327, 103464, 79037, 100833, 80118, 29178, 163899)
)
read_count_PML_SuperT <- tibble(
  case_ID = c(1, 2, 3, 4, 5, 6, 7, 8),
  read_count = c(67296, 155899, 110033, 289425, 42767, 128892, 21884, 76030)
)


########## functions ##########

format_dRNAseq_inputs <- function(in_tbl) {
  # inverse the start and the end of transcripts from early region
  in_tbl <- in_tbl %>%
    mutate(newtx_start = ifelse(strand == "-", tx_end, tx_start),
           newtx_end = ifelse(strand == "-", tx_start, tx_end),
           tx_start = newtx_start,
           tx_end = newtx_end) %>%
    select(-newtx_start, -newtx_end) 
  # inverse the start and the end of introns/exons of transcripts from early region
  in_tbl <- in_tbl %>%
    mutate(new_start = ifelse(strand == "-", end, start),
           new_end = ifelse(strand == "-", start, end),
           start = new_start,
           end = new_end) %>%
    select(-new_start, -new_end)
  # overwrite tx class of unspliced transcripts
  in_tbl <- in_tbl %>%
    mutate(tx_class = ifelse(tx_class == "u-", -1, tx_class)) %>%
    mutate(tx_class = ifelse(tx_class == "u+", -2, tx_class)) %>%
    mutate(tx_class = as.numeric(tx_class))
  in_tbl
}


aggregate_PCR <- function(PCR_tbl, genome_length, PCR_strand, list_IMR, list_293, list_WA_late, list_WA_early, pos_delimiter = "_", intron_delimiter = "\n") {

  # Definition of late and early regions
  half_genome <- round(genome_length/2)
  late_region <- c(
   seq(1, half_genome),
   seq(genome_length, 3*half_genome),
   seq(2*genome_length, 5*half_genome),
   seq(3*genome_length, 7*half_genome),
   seq(4*genome_length, 9*half_genome),
   seq(5*genome_length, 11*half_genome)
  )
  early_region <- c(
   seq(half_genome, genome_length),
   seq(3*half_genome, 2*genome_length),
   seq(5*half_genome, 3*genome_length),
   seq(7*half_genome, 4*genome_length),
   seq(9*half_genome, 5*genome_length),
   seq(11*half_genome, 6*genome_length)
  )

  # Making wraparound transcript lists
  PCR_tbl_class <- PCR_tbl %>%
   filter(span_type == "intron") %>%
   select(start, end, strand, tx_class, tx_class_count) %>%
   distinct() %>%
   arrange(tx_class)
  
  PCR_tbl_class_pass <- PCR_tbl_class %>%
   filter(
     if (PCR_strand == "late") {
       (strand == "+" & ((start %in% late_region) | (end %in% late_region)))
     } else if (PCR_strand == "early") {
       (strand == "-" & ((start %in% early_region) | (end %in% early_region)))
     } else {
       FALSE
     }
   )

  passing_tx_classes <- PCR_tbl_class_pass$tx_class

  PCR_tbl <- PCR_tbl %>%
    filter(tx_class %in% passing_tx_classes) %>%
    filter(span_type == "intron") %>%
    select(start, end, tx_class, tx_class_count) %>%
    distinct() %>%
    arrange(-tx_class_count) %>%
    mutate(junc = paste(start, end, sep = pos_delimiter)) %>%
    group_by(tx_class, tx_class_count) %>%
    summarise(introns = paste(junc, collapse = intron_delimiter)) %>%
    ungroup()
  
  PCR_tbl$ID_IMR <- NA
  PCR_tbl$ID_293 <- NA
  PCR_tbl$description <- NA
  PCR_tbl$WA <- NA

  for (i in 1:nrow(PCR_tbl)) {
    ID_IMR <- list_IMR$ID[list_IMR$introns == PCR_tbl$introns[i]]
    description <- list_IMR$description[list_IMR$introns == PCR_tbl$introns[i]]
    ID_293 <- list_293$ID[list_293$introns == PCR_tbl$introns[i]]
    description <- list_293$description[list_293$introns == PCR_tbl$introns[i]]
    len_ID_IMR <- length(ID_IMR)
      if (len_ID_IMR > 0) {
        PCR_tbl$ID_IMR[i] <- ID_IMR[len_ID_IMR]
      }
    len_ID_293 <- length(ID_293)
      if (len_ID_293 > 0) {
        PCR_tbl$ID_293[i] <- ID_293[len_ID_293]
      }
    len_description <- length(description)
      if (len_description > 0) {
        PCR_tbl$description[i] <- description[len_description]
      }
    if (PCR_tbl$introns[i] %in% list_WA_late$introns) {
        PCR_tbl$WA[i] <- "WA"
      } else if (PCR_tbl$introns[i] %in% list_WA_early$introns) {
        PCR_tbl$WA[i] <- "WA"
      } else {
        FALSE
      }
  }
  return(PCR_tbl)
}


calc_exon_length <- function(PCR_tbl, tx_class_list) {
  PCR_tbl <- PCR_tbl %>%
   filter(tx_class %in% tx_class_list$tx_class) %>%
   filter(span_type == "exon") %>%
   mutate(exon_distance = abs(start - end)) %>%
   group_by(name) %>%
   mutate(mapped_exon_length = sum(exon_distance)) %>%
   select(name, strand, tx_class, tx_class_count, mapped_exon_length) %>%
   distinct()
}


calc_exon_length_mean <- function(PCR_tbl) {
  PCR_tbl <- PCR_tbl %>%
   group_by(tx_class) %>%
   mutate(mean_exon_length = round(mean(mapped_exon_length), 0)) %>%
   mutate(sd_exon_length = round(sd(mapped_exon_length), 1)) %>%
   ungroup() %>%
   select(tx_class, mean_exon_length, sd_exon_length) %>%
   distinct() %>%
   arrange(tx_class)
}


combine_list_length <- function(list_WA, predicted_length_WA, list_length_mean) {
  predicted_length_IMR <- predicted_length_WA %>%
    select(ID_IMR, predicted_length) %>%
    filter(ID_IMR != "NA") %>%
    rename("predicted_length_IMR" = predicted_length)
  predicted_length_293 <- predicted_length_WA %>%
    select(ID_293, predicted_length) %>%
    filter(ID_293 != "NA") %>%
    rename("predicted_length_293" = predicted_length)

  list_WA <- dplyr::left_join(list_WA, predicted_length_IMR, by = "ID_IMR")
  list_WA <- dplyr::left_join(list_WA, predicted_length_293, by = "ID_293")
  
  list_WA <- list_WA %>%
    mutate(predicted_length = case_when(
      !is.na(predicted_length_IMR) ~ predicted_length_IMR,
      !is.na(predicted_length_293) ~ predicted_length_293,
      TRUE ~ NA_real_
    )) %>%
    select(-predicted_length_IMR, -predicted_length_293)
  
  list_WA <- dplyr::left_join(list_WA, list_length_mean, by = "tx_class")
}


integrate_PML_samples <- function(result_tbl_left, caseNo_left, result_tbl_right, caseNo_right, read_count_tbl) {
  
  read_count_left <- read_count_tbl %>%
    filter(case_ID == caseNo_left) %>%
    pull(read_count)
  read_count_right <- read_count_tbl %>%
    filter(case_ID == caseNo_right) %>%
    pull(read_count)
  col_name_left = paste("read_count_Case", caseNo_left, sep="")
  col_name_right = paste("read_count_Case", caseNo_right, sep="")

  if(caseNo_left != 0) {
    result_tbl_left <- result_tbl_left %>%
    select(predicted_length, introns, ID_IMR, ID_293, tx_class_count) %>%
    rename(product_length = predicted_length) %>%
    mutate(tx_class_proportion = round(100*tx_class_count/read_count_left, 2)) %>%
    mutate(!!col_name_left := paste(tx_class_count, "\n(", tx_class_proportion,")", sep="")) %>%
    select(product_length, introns, ID_IMR, ID_293, !!sym(col_name_left))
  }

  result_tbl_right <- result_tbl_right %>%
    select(predicted_length, introns, ID_IMR, ID_293, tx_class_count) %>%
    rename(product_length = predicted_length) %>%
    mutate(tx_class_proportion = round(100*tx_class_count/read_count_right, 2)) %>%
    mutate(!!col_name_right := paste(tx_class_count, "\n(", tx_class_proportion,")", sep="")) %>%
    select(product_length, introns, ID_IMR, ID_293, !!sym(col_name_right))
  
  result_tbl <- dplyr::full_join(result_tbl_left, result_tbl_right, by = c("product_length", "introns", "ID_IMR", "ID_293")) %>%
  arrange(-product_length)
  return(result_tbl)
}


integrate_PML_bubble_samples <- function(result_tbl_left, caseNo_left, result_tbl_right, caseNo_right) {
  
  col_name_left = paste("read_count_Case", caseNo_left, sep="")
  col_name_right = paste("read_count_Case", caseNo_right, sep="")

  if(caseNo_left != 0) {
    result_tbl_left <- result_tbl_left %>%
    select(predicted_length, introns, ID_IMR, ID_293, tx_class_count) %>%
    rename(product_length = predicted_length) %>%
    mutate(!!col_name_left := paste(tx_class_count)) %>%
    select(product_length, introns, ID_IMR, ID_293, !!sym(col_name_left))
  }

  result_tbl_right <- result_tbl_right %>%
    select(predicted_length, introns, ID_IMR, ID_293, tx_class_count) %>%
    rename(product_length = predicted_length) %>%
    mutate(!!col_name_right := paste(tx_class_count)) %>%
    select(product_length, introns, ID_IMR, ID_293, !!sym(col_name_right))
  
  result_tbl <- dplyr::full_join(result_tbl_left, result_tbl_right, by = c("product_length", "introns", "ID_IMR", "ID_293")) %>%
  arrange(-product_length)
  return(result_tbl)
}


########## importing data ##########

# VP1 PCR
JCPyV_PCR_VP1_PML_1_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP1_PML_1_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_VP1_PML_2_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP1_PML_2_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_VP1_PML_3_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP1_PML_3_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_VP1_PML_4_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP1_PML_4_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_VP1_PML_5_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP1_PML_5_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_VP1_PML_6_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP1_PML_6_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_VP1_PML_7_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP1_PML_7_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_VP1_PML_8_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP1_PML_8_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()

# VP2/3 PCR
JCPyV_PCR_VP23_PML_1_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP23_PML_1_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_VP23_PML_2_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP23_PML_2_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_VP23_PML_3_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP23_PML_3_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_VP23_PML_4_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP23_PML_4_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_VP23_PML_5_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP23_PML_5_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_VP23_PML_6_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP23_PML_6_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_VP23_PML_7_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP23_PML_7_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_VP23_PML_8_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_VP23_PML_8_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()

# SuperT PCR
JCPyV_PCR_SuperT_PML_1_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_SuperT_PML_1_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_SuperT_PML_2_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_SuperT_PML_2_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_SuperT_PML_3_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_SuperT_PML_3_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_SuperT_PML_4_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_SuperT_PML_4_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_SuperT_PML_5_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_SuperT_PML_5_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_SuperT_PML_6_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_SuperT_PML_6_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_SuperT_PML_7_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_SuperT_PML_7_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()
JCPyV_PCR_SuperT_PML_8_raw <- read_delim("inputs/JCPyV/PML/JCV_WA-PCR_SuperT_PML_8_spans.txt", delim = "\t") %>%
  format_dRNAseq_inputs()

# transcript list
all_list_IMR32 <- read_delim("inputs/JCPyV/csv/JCPyV_IMR32_transcript_all.csv", delim = ",")
all_list_293 <- read_delim("inputs/JCPyV/csv/JCPyV_293_transcript_all.csv", delim = ",")
all_WA_list_late <- read_delim("inputs/JCPyV/csv/JCPyV_WA_late_list.csv", delim = ",")
all_WA_list_early <- read_delim("inputs/JCPyV/csv/JCPyV_WA_SuperT_list.csv", delim = ",")
WA_PCR_predicted_length <- read_delim("inputs/JCPyV/csv/WA_PCR_predicted_length.csv", delim = ",")


########## making lists ##########

JCPyV_PCR_VP1_PML_1_list <- aggregate_PCR(JCPyV_PCR_VP1_PML_1_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_VP1_PML_2_list <- aggregate_PCR(JCPyV_PCR_VP1_PML_2_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_VP1_PML_3_list <- aggregate_PCR(JCPyV_PCR_VP1_PML_3_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_VP1_PML_4_list <- aggregate_PCR(JCPyV_PCR_VP1_PML_4_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_VP1_PML_5_list <- aggregate_PCR(JCPyV_PCR_VP1_PML_5_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_VP1_PML_6_list <- aggregate_PCR(JCPyV_PCR_VP1_PML_6_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_VP1_PML_7_list <- aggregate_PCR(JCPyV_PCR_VP1_PML_7_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_VP1_PML_8_list <- aggregate_PCR(JCPyV_PCR_VP1_PML_8_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)

JCPyV_PCR_VP23_PML_1_list <- aggregate_PCR(JCPyV_PCR_VP23_PML_1_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_VP23_PML_2_list <- aggregate_PCR(JCPyV_PCR_VP23_PML_2_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_VP23_PML_3_list <- aggregate_PCR(JCPyV_PCR_VP23_PML_3_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_VP23_PML_4_list <- aggregate_PCR(JCPyV_PCR_VP23_PML_4_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_VP23_PML_5_list <- aggregate_PCR(JCPyV_PCR_VP23_PML_5_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_VP23_PML_6_list <- aggregate_PCR(JCPyV_PCR_VP23_PML_6_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_VP23_PML_7_list <- aggregate_PCR(JCPyV_PCR_VP23_PML_7_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_VP23_PML_8_list <- aggregate_PCR(JCPyV_PCR_VP23_PML_8_raw, 5130, "late", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)

JCPyV_PCR_SuperT_PML_1_list <- aggregate_PCR(JCPyV_PCR_SuperT_PML_1_raw, 5130, "early", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_SuperT_PML_2_list <- aggregate_PCR(JCPyV_PCR_SuperT_PML_2_raw, 5130, "early", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_SuperT_PML_3_list <- aggregate_PCR(JCPyV_PCR_SuperT_PML_3_raw, 5130, "early", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_SuperT_PML_4_list <- aggregate_PCR(JCPyV_PCR_SuperT_PML_4_raw, 5130, "early", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_SuperT_PML_5_list <- aggregate_PCR(JCPyV_PCR_SuperT_PML_5_raw, 5130, "early", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_SuperT_PML_6_list <- aggregate_PCR(JCPyV_PCR_SuperT_PML_6_raw, 5130, "early", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_SuperT_PML_7_list <- aggregate_PCR(JCPyV_PCR_SuperT_PML_7_raw, 5130, "early", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)
JCPyV_PCR_SuperT_PML_8_list <- aggregate_PCR(JCPyV_PCR_SuperT_PML_8_raw, 5130, "early", all_list_IMR32, all_list_293, all_WA_list_late, all_WA_list_early)


########## filtering data ##########

JCPyV_PCR_VP1_PML_1 <- JCPyV_PCR_VP1_PML_1_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP1_PML_1_list$tx_class)
JCPyV_PCR_VP1_PML_2 <- JCPyV_PCR_VP1_PML_2_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP1_PML_2_list$tx_class)
JCPyV_PCR_VP1_PML_3 <- JCPyV_PCR_VP1_PML_3_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP1_PML_3_list$tx_class)
JCPyV_PCR_VP1_PML_4 <- JCPyV_PCR_VP1_PML_4_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP1_PML_4_list$tx_class)
JCPyV_PCR_VP1_PML_5 <- JCPyV_PCR_VP1_PML_5_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP1_PML_5_list$tx_class)
JCPyV_PCR_VP1_PML_6 <- JCPyV_PCR_VP1_PML_6_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP1_PML_6_list$tx_class)
JCPyV_PCR_VP1_PML_7 <- JCPyV_PCR_VP1_PML_7_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP1_PML_7_list$tx_class)
JCPyV_PCR_VP1_PML_8 <- JCPyV_PCR_VP1_PML_8_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP1_PML_8_list$tx_class)

JCPyV_PCR_VP23_PML_1 <- JCPyV_PCR_VP23_PML_1_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP23_PML_1_list$tx_class)
JCPyV_PCR_VP23_PML_2 <- JCPyV_PCR_VP23_PML_2_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP23_PML_2_list$tx_class)
JCPyV_PCR_VP23_PML_3 <- JCPyV_PCR_VP23_PML_3_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP23_PML_3_list$tx_class)
JCPyV_PCR_VP23_PML_4 <- JCPyV_PCR_VP23_PML_4_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP23_PML_4_list$tx_class)
JCPyV_PCR_VP23_PML_5 <- JCPyV_PCR_VP23_PML_5_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP23_PML_5_list$tx_class)
JCPyV_PCR_VP23_PML_6 <- JCPyV_PCR_VP23_PML_6_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP23_PML_6_list$tx_class)
JCPyV_PCR_VP23_PML_7 <- JCPyV_PCR_VP23_PML_7_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP23_PML_7_list$tx_class)
JCPyV_PCR_VP23_PML_8 <- JCPyV_PCR_VP23_PML_8_raw %>%
 filter(tx_class %in% JCPyV_PCR_VP23_PML_8_list$tx_class)

JCPyV_PCR_SuperT_PML_1 <- JCPyV_PCR_SuperT_PML_1_raw %>%
 filter(tx_class %in% JCPyV_PCR_SuperT_PML_1_list$tx_class)
JCPyV_PCR_SuperT_PML_2 <- JCPyV_PCR_SuperT_PML_2_raw %>%
 filter(tx_class %in% JCPyV_PCR_SuperT_PML_2_list$tx_class)
JCPyV_PCR_SuperT_PML_3 <- JCPyV_PCR_SuperT_PML_3_raw %>%
 filter(tx_class %in% JCPyV_PCR_SuperT_PML_3_list$tx_class)
JCPyV_PCR_SuperT_PML_4 <- JCPyV_PCR_SuperT_PML_4_raw %>%
 filter(tx_class %in% JCPyV_PCR_SuperT_PML_4_list$tx_class)
JCPyV_PCR_SuperT_PML_5 <- JCPyV_PCR_SuperT_PML_5_raw %>%
 filter(tx_class %in% JCPyV_PCR_SuperT_PML_5_list$tx_class)
JCPyV_PCR_SuperT_PML_6 <- JCPyV_PCR_SuperT_PML_6_raw %>%
 filter(tx_class %in% JCPyV_PCR_SuperT_PML_6_list$tx_class)
JCPyV_PCR_SuperT_PML_7 <- JCPyV_PCR_SuperT_PML_7_raw %>%
 filter(tx_class %in% JCPyV_PCR_SuperT_PML_7_list$tx_class)
JCPyV_PCR_SuperT_PML_8 <- JCPyV_PCR_SuperT_PML_8_raw %>%
 filter(tx_class %in% JCPyV_PCR_SuperT_PML_8_list$tx_class)


########## calculation of PCR product length ##########

JCPyV_PCR_VP1_PML_1_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP1_PML_1, JCPyV_PCR_VP1_PML_1_list))
JCPyV_PCR_VP1_PML_2_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP1_PML_2, JCPyV_PCR_VP1_PML_2_list))
JCPyV_PCR_VP1_PML_3_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP1_PML_3, JCPyV_PCR_VP1_PML_3_list))
JCPyV_PCR_VP1_PML_4_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP1_PML_4, JCPyV_PCR_VP1_PML_4_list))
JCPyV_PCR_VP1_PML_5_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP1_PML_5, JCPyV_PCR_VP1_PML_5_list))
JCPyV_PCR_VP1_PML_6_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP1_PML_6, JCPyV_PCR_VP1_PML_6_list))
JCPyV_PCR_VP1_PML_7_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP1_PML_7, JCPyV_PCR_VP1_PML_7_list))
JCPyV_PCR_VP1_PML_8_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP1_PML_8, JCPyV_PCR_VP1_PML_8_list))

JCPyV_PCR_VP23_PML_1_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP23_PML_1, JCPyV_PCR_VP23_PML_1_list))
JCPyV_PCR_VP23_PML_2_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP23_PML_2, JCPyV_PCR_VP23_PML_2_list))
JCPyV_PCR_VP23_PML_3_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP23_PML_3, JCPyV_PCR_VP23_PML_3_list))
JCPyV_PCR_VP23_PML_4_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP23_PML_4, JCPyV_PCR_VP23_PML_4_list))
JCPyV_PCR_VP23_PML_5_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP23_PML_5, JCPyV_PCR_VP23_PML_5_list))
JCPyV_PCR_VP23_PML_6_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP23_PML_6, JCPyV_PCR_VP23_PML_6_list))
JCPyV_PCR_VP23_PML_7_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP23_PML_7, JCPyV_PCR_VP23_PML_7_list))
JCPyV_PCR_VP23_PML_8_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_VP23_PML_8, JCPyV_PCR_VP23_PML_8_list))

JCPyV_PCR_SuperT_PML_1_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_SuperT_PML_1, JCPyV_PCR_SuperT_PML_1_list))
JCPyV_PCR_SuperT_PML_2_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_SuperT_PML_2, JCPyV_PCR_SuperT_PML_2_list))
JCPyV_PCR_SuperT_PML_3_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_SuperT_PML_3, JCPyV_PCR_SuperT_PML_3_list))
JCPyV_PCR_SuperT_PML_4_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_SuperT_PML_4, JCPyV_PCR_SuperT_PML_4_list))
JCPyV_PCR_SuperT_PML_5_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_SuperT_PML_5, JCPyV_PCR_SuperT_PML_5_list))
JCPyV_PCR_SuperT_PML_6_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_SuperT_PML_6, JCPyV_PCR_SuperT_PML_6_list))
JCPyV_PCR_SuperT_PML_7_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_SuperT_PML_7, JCPyV_PCR_SuperT_PML_7_list))
JCPyV_PCR_SuperT_PML_8_exon_length_mean <- calc_exon_length_mean(calc_exon_length(JCPyV_PCR_SuperT_PML_8, JCPyV_PCR_SuperT_PML_8_list))


JCPyV_PCR_VP1_PML_1_list_length <- combine_list_length(JCPyV_PCR_VP1_PML_1_list, WA_PCR_predicted_length, JCPyV_PCR_VP1_PML_1_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP1_IMR) | (ID_293 %in% class_WA_VP1_293))
JCPyV_PCR_VP1_PML_2_list_length <- combine_list_length(JCPyV_PCR_VP1_PML_2_list, WA_PCR_predicted_length, JCPyV_PCR_VP1_PML_2_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP1_IMR) | (ID_293 %in% class_WA_VP1_293))
JCPyV_PCR_VP1_PML_3_list_length <- combine_list_length(JCPyV_PCR_VP1_PML_3_list, WA_PCR_predicted_length, JCPyV_PCR_VP1_PML_3_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP1_IMR) | (ID_293 %in% class_WA_VP1_293))
JCPyV_PCR_VP1_PML_4_list_length <- combine_list_length(JCPyV_PCR_VP1_PML_4_list, WA_PCR_predicted_length, JCPyV_PCR_VP1_PML_4_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP1_IMR) | (ID_293 %in% class_WA_VP1_293))
JCPyV_PCR_VP1_PML_5_list_length <- combine_list_length(JCPyV_PCR_VP1_PML_5_list, WA_PCR_predicted_length, JCPyV_PCR_VP1_PML_5_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP1_IMR) | (ID_293 %in% class_WA_VP1_293))
JCPyV_PCR_VP1_PML_6_list_length <- combine_list_length(JCPyV_PCR_VP1_PML_6_list, WA_PCR_predicted_length, JCPyV_PCR_VP1_PML_6_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP1_IMR) | (ID_293 %in% class_WA_VP1_293))
JCPyV_PCR_VP1_PML_7_list_length <- combine_list_length(JCPyV_PCR_VP1_PML_7_list, WA_PCR_predicted_length, JCPyV_PCR_VP1_PML_7_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP1_IMR) | (ID_293 %in% class_WA_VP1_293))
JCPyV_PCR_VP1_PML_8_list_length <- combine_list_length(JCPyV_PCR_VP1_PML_8_list, WA_PCR_predicted_length, JCPyV_PCR_VP1_PML_8_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP1_IMR) | (ID_293 %in% class_WA_VP1_293))

JCPyV_PCR_VP23_PML_1_list_length <- combine_list_length(JCPyV_PCR_VP23_PML_1_list, WA_PCR_predicted_length, JCPyV_PCR_VP23_PML_1_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP23_IMR) | (ID_293 %in% class_WA_VP23_293))
JCPyV_PCR_VP23_PML_2_list_length <- combine_list_length(JCPyV_PCR_VP23_PML_2_list, WA_PCR_predicted_length, JCPyV_PCR_VP23_PML_2_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP23_IMR) | (ID_293 %in% class_WA_VP23_293))
JCPyV_PCR_VP23_PML_3_list_length <- combine_list_length(JCPyV_PCR_VP23_PML_3_list, WA_PCR_predicted_length, JCPyV_PCR_VP23_PML_3_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP23_IMR) | (ID_293 %in% class_WA_VP23_293))
JCPyV_PCR_VP23_PML_4_list_length <- combine_list_length(JCPyV_PCR_VP23_PML_4_list, WA_PCR_predicted_length, JCPyV_PCR_VP23_PML_4_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP23_IMR) | (ID_293 %in% class_WA_VP23_293))
JCPyV_PCR_VP23_PML_5_list_length <- combine_list_length(JCPyV_PCR_VP23_PML_5_list, WA_PCR_predicted_length, JCPyV_PCR_VP23_PML_5_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP23_IMR) | (ID_293 %in% class_WA_VP23_293))
JCPyV_PCR_VP23_PML_6_list_length <- combine_list_length(JCPyV_PCR_VP23_PML_6_list, WA_PCR_predicted_length, JCPyV_PCR_VP23_PML_6_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP23_IMR) | (ID_293 %in% class_WA_VP23_293))
JCPyV_PCR_VP23_PML_7_list_length <- combine_list_length(JCPyV_PCR_VP23_PML_7_list, WA_PCR_predicted_length, JCPyV_PCR_VP23_PML_7_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP23_IMR) | (ID_293 %in% class_WA_VP23_293))
JCPyV_PCR_VP23_PML_8_list_length <- combine_list_length(JCPyV_PCR_VP23_PML_8_list, WA_PCR_predicted_length, JCPyV_PCR_VP23_PML_8_exon_length_mean) %>%
 filter((ID_IMR %in% class_WA_VP23_IMR) | (ID_293 %in% class_WA_VP23_293))

JCPyV_PCR_SuperT_PML_1_list_length <- combine_list_length(JCPyV_PCR_SuperT_PML_1_list, WA_PCR_predicted_length, JCPyV_PCR_SuperT_PML_1_exon_length_mean) %>%
 filter((ID_IMR %in% class_SuperT_IMR) | (ID_293 %in% class_SuperT_293) | (introns %in% appendix_SuperT$introns)) %>%
 dplyr::left_join(appendix_SuperT, by = "introns") %>%
 mutate(predicted_length = ifelse(is.na(predicted_length) & !is.na(PCR_predicted_length), PCR_predicted_length, predicted_length)) %>%
 select(-PCR_predicted_length)
JCPyV_PCR_SuperT_PML_2_list_length <- combine_list_length(JCPyV_PCR_SuperT_PML_2_list, WA_PCR_predicted_length, JCPyV_PCR_SuperT_PML_2_exon_length_mean) %>%
 filter((ID_IMR %in% class_SuperT_IMR) | (ID_293 %in% class_SuperT_293) | (introns %in% appendix_SuperT$introns)) %>%
 dplyr::left_join(appendix_SuperT, by = "introns") %>%
 mutate(predicted_length = ifelse(is.na(predicted_length) & !is.na(PCR_predicted_length), PCR_predicted_length, predicted_length)) %>%
 select(-PCR_predicted_length)
JCPyV_PCR_SuperT_PML_3_list_length <- combine_list_length(JCPyV_PCR_SuperT_PML_3_list, WA_PCR_predicted_length, JCPyV_PCR_SuperT_PML_3_exon_length_mean) %>%
 filter((ID_IMR %in% class_SuperT_IMR) | (ID_293 %in% class_SuperT_293) | (introns %in% appendix_SuperT$introns)) %>%
 dplyr::left_join(appendix_SuperT, by = "introns") %>%
 mutate(predicted_length = ifelse(is.na(predicted_length) & !is.na(PCR_predicted_length), PCR_predicted_length, predicted_length)) %>%
 select(-PCR_predicted_length)
JCPyV_PCR_SuperT_PML_4_list_length <- combine_list_length(JCPyV_PCR_SuperT_PML_4_list, WA_PCR_predicted_length, JCPyV_PCR_SuperT_PML_4_exon_length_mean) %>%
 filter((ID_IMR %in% class_SuperT_IMR) | (ID_293 %in% class_SuperT_293) | (introns %in% appendix_SuperT$introns)) %>%
 dplyr::left_join(appendix_SuperT, by = "introns") %>%
 mutate(predicted_length = ifelse(is.na(predicted_length) & !is.na(PCR_predicted_length), PCR_predicted_length, predicted_length)) %>%
 select(-PCR_predicted_length)
JCPyV_PCR_SuperT_PML_5_list_length <- combine_list_length(JCPyV_PCR_SuperT_PML_5_list, WA_PCR_predicted_length, JCPyV_PCR_SuperT_PML_5_exon_length_mean) %>%
 filter((ID_IMR %in% class_SuperT_IMR) | (ID_293 %in% class_SuperT_293) | (introns %in% appendix_SuperT$introns)) %>%
 dplyr::left_join(appendix_SuperT, by = "introns") %>%
 mutate(predicted_length = ifelse(is.na(predicted_length) & !is.na(PCR_predicted_length), PCR_predicted_length, predicted_length)) %>%
 select(-PCR_predicted_length)
JCPyV_PCR_SuperT_PML_6_list_length <- combine_list_length(JCPyV_PCR_SuperT_PML_6_list, WA_PCR_predicted_length, JCPyV_PCR_SuperT_PML_6_exon_length_mean) %>%
 filter((ID_IMR %in% class_SuperT_IMR) | (ID_293 %in% class_SuperT_293) | (introns %in% appendix_SuperT$introns)) %>%
 dplyr::left_join(appendix_SuperT, by = "introns") %>%
 mutate(predicted_length = ifelse(is.na(predicted_length) & !is.na(PCR_predicted_length), PCR_predicted_length, predicted_length)) %>%
 select(-PCR_predicted_length)
JCPyV_PCR_SuperT_PML_7_list_length <- combine_list_length(JCPyV_PCR_SuperT_PML_7_list, WA_PCR_predicted_length, JCPyV_PCR_SuperT_PML_7_exon_length_mean) %>%
 filter((ID_IMR %in% class_SuperT_IMR) | (ID_293 %in% class_SuperT_293) | (introns %in% appendix_SuperT$introns)) %>%
 dplyr::left_join(appendix_SuperT, by = "introns") %>%
 mutate(predicted_length = ifelse(is.na(predicted_length) & !is.na(PCR_predicted_length), PCR_predicted_length, predicted_length)) %>%
 select(-PCR_predicted_length)
JCPyV_PCR_SuperT_PML_8_list_length <- combine_list_length(JCPyV_PCR_SuperT_PML_8_list, WA_PCR_predicted_length, JCPyV_PCR_SuperT_PML_8_exon_length_mean) %>%
 filter((ID_IMR %in% class_SuperT_IMR) | (ID_293 %in% class_SuperT_293) | (introns %in% appendix_SuperT$introns)) %>%
 dplyr::left_join(appendix_SuperT, by = "introns") %>%
 mutate(predicted_length = ifelse(is.na(predicted_length) & !is.na(PCR_predicted_length), PCR_predicted_length, predicted_length)) %>%
 select(-PCR_predicted_length)


PML_VP1_list <- list(
  JCPyV_PCR_VP1_PML_1_list_length,
  JCPyV_PCR_VP1_PML_2_list_length,
  JCPyV_PCR_VP1_PML_3_list_length,
  JCPyV_PCR_VP1_PML_4_list_length,
  JCPyV_PCR_VP1_PML_5_list_length,
  JCPyV_PCR_VP1_PML_6_list_length,
  JCPyV_PCR_VP1_PML_7_list_length,
  JCPyV_PCR_VP1_PML_8_list_length
)

JCPyV_PCR_VP1_PML_allcase <- PML_VP1_list[[1]]

for (i in 2:length(PML_VP1_list)) {
  caseNo_left <- ifelse(i == 2, 1, 0)
  JCPyV_PCR_VP1_PML_allcase <- integrate_PML_samples(JCPyV_PCR_VP1_PML_allcase, caseNo_left, PML_VP1_list[[i]], i, read_count_PML_VP1)
}


PML_VP23_list <- list(
  JCPyV_PCR_VP23_PML_1_list_length,
  JCPyV_PCR_VP23_PML_2_list_length,
  JCPyV_PCR_VP23_PML_3_list_length,
  JCPyV_PCR_VP23_PML_4_list_length,
  JCPyV_PCR_VP23_PML_5_list_length,
  JCPyV_PCR_VP23_PML_6_list_length,
  JCPyV_PCR_VP23_PML_7_list_length,
  JCPyV_PCR_VP23_PML_8_list_length
)

JCPyV_PCR_VP23_PML_allcase <- PML_VP23_list[[1]]

for (i in 2:length(PML_VP23_list)) {
  caseNo_left <- ifelse(i == 2, 1, 0)
  JCPyV_PCR_VP23_PML_allcase <- integrate_PML_samples(JCPyV_PCR_VP23_PML_allcase, caseNo_left, PML_VP23_list[[i]], i, read_count_PML_VP23)
}


PML_SuperT_list <- list(
  JCPyV_PCR_SuperT_PML_1_list_length,
  JCPyV_PCR_SuperT_PML_2_list_length,
  JCPyV_PCR_SuperT_PML_3_list_length,
  JCPyV_PCR_SuperT_PML_4_list_length,
  JCPyV_PCR_SuperT_PML_5_list_length,
  JCPyV_PCR_SuperT_PML_6_list_length,
  JCPyV_PCR_SuperT_PML_7_list_length,
  JCPyV_PCR_SuperT_PML_8_list_length
)

JCPyV_PCR_SuperT_PML_allcase <- PML_SuperT_list[[1]]

for (i in 2:length(PML_SuperT_list)) {
  caseNo_left <- ifelse(i == 2, 1, 0)
  JCPyV_PCR_SuperT_PML_allcase <- integrate_PML_samples(JCPyV_PCR_SuperT_PML_allcase, caseNo_left, PML_SuperT_list[[i]], i, read_count_PML_SuperT)
}


########## bubble chart ##########

# VP1
# Data processing
JCPyV_PCR_VP1_PML_bubble_allcase <- PML_VP1_list[[1]]
for (i in 2:length(PML_VP1_list)) {
  caseNo_left <- ifelse(i == 2, 1, 0)
  JCPyV_PCR_VP1_PML_bubble_allcase <- integrate_PML_bubble_samples(JCPyV_PCR_VP1_PML_bubble_allcase, caseNo_left, PML_VP1_list[[i]], i)
}
JCPyV_PCR_VP1_PML_bubble_data <- JCPyV_PCR_VP1_PML_bubble_allcase %>%
 select(-introns, -ID_IMR, -ID_293) %>%
 mutate(
  read_count_Case1 = as.numeric(read_count_Case1),
  read_count_Case2 = as.numeric(read_count_Case2),
  read_count_Case3 = as.numeric(read_count_Case3),
  read_count_Case4 = as.numeric(read_count_Case4),
  read_count_Case5 = as.numeric(read_count_Case5),
  read_count_Case6 = as.numeric(read_count_Case6),
  read_count_Case7 = as.numeric(read_count_Case7),
  read_count_Case8 = as.numeric(read_count_Case8)
 )
colnames(JCPyV_PCR_VP1_PML_bubble_data) <- c("Length", "Case 1", "Case 2","Case 3","Case 4","Case 5","Case 6","Case 7","Case 8")
JCPyV_PCR_VP1_PML_bubble_data_melt <- melt(JCPyV_PCR_VP1_PML_bubble_data, id.vars = c("Length"))
JCPyV_PCR_VP1_PML_bubble_data_melt$Length <- factor(JCPyV_PCR_VP1_PML_bubble_data_melt$Length, levels = rev(unique(JCPyV_PCR_VP1_PML_bubble_data_melt$Length)))

# Plot
JCPyV_PCR_VP1_PML_bubble_plot = ggplot(JCPyV_PCR_VP1_PML_bubble_data_melt, aes(x = variable, y = Length)) +
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) +
  scale_size_continuous(limits = c(1, 400), range = c(1, 15), breaks = c(1, 10, 100, 400)) +
  labs(x= "", y = expression(bold("Product Length [bp]")), size = "VP1\nRead Count", fill = "")  +
  theme(legend.key = element_blank(),
        axis.text.x = element_text(color = "#000000", size = 14, angle = 45, vjust = 0, hjust = 0),
        axis.text.y = element_text(color = "#000000", size = 14),
        legend.text = element_text(size = 12, color = "#000000"),
        legend.title = element_text(size = 12, face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(color = "#000000", fill = NA, size = 1.2),
        legend.position = "right") +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  scale_x_discrete(position = "top", limits = levels(JCPyV_PCR_VP1_PML_bubble_data_melt$variable))

ggsave("outputs/figs/fig4/JCPyV_PCR_VP1_PML_bubble_plot.tiff", plot = JCPyV_PCR_VP1_PML_bubble_plot, width = 12, height = 7, units = "cm", dpi = 300)


# VP23
# Data processing
JCPyV_PCR_VP23_PML_bubble_allcase <- PML_VP23_list[[1]]
for (i in 2:length(PML_VP23_list)) {
  caseNo_left <- ifelse(i == 2, 1, 0)
  JCPyV_PCR_VP23_PML_bubble_allcase <- integrate_PML_bubble_samples(JCPyV_PCR_VP23_PML_bubble_allcase, caseNo_left, PML_VP23_list[[i]], i)
}
JCPyV_PCR_VP23_PML_bubble_data <- JCPyV_PCR_VP23_PML_bubble_allcase %>%
 select(-introns, -ID_IMR, -ID_293) %>%
 mutate(
  read_count_Case1 = as.numeric(read_count_Case1),
  read_count_Case2 = as.numeric(read_count_Case2),
  read_count_Case3 = as.numeric(read_count_Case3),
  read_count_Case4 = as.numeric(read_count_Case4),
  read_count_Case5 = as.numeric(read_count_Case5),
  read_count_Case6 = as.numeric(read_count_Case6),
  read_count_Case7 = as.numeric(read_count_Case7),
  read_count_Case8 = as.numeric(read_count_Case8)
 )
colnames(JCPyV_PCR_VP23_PML_bubble_data) <- c("Length", "Case 1", "Case 2","Case 3","Case 4","Case 5","Case 6","Case 7","Case 8")
JCPyV_PCR_VP23_PML_bubble_data_melt <- melt(JCPyV_PCR_VP23_PML_bubble_data, id.vars = c("Length"))
JCPyV_PCR_VP23_PML_bubble_data_melt$Length <- factor(JCPyV_PCR_VP23_PML_bubble_data_melt$Length, levels = rev(unique(JCPyV_PCR_VP23_PML_bubble_data_melt$Length)))

# Plot
JCPyV_PCR_VP23_PML_bubble_plot = ggplot(JCPyV_PCR_VP23_PML_bubble_data_melt, aes(x = variable, y = Length)) +
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) +
  scale_size_continuous(limits = c(1, 200), range = c(1, 15), breaks = c(1, 10, 100, 200)) +
  labs(x= "", y = expression(bold("Product Length [bp]")), size = "VP2/3\nRead Count", fill = "")  +
  theme(legend.key = element_blank(),
        axis.text.x = element_text(color = "#000000", size = 14, angle = 45, vjust = 0, hjust = 0),
        axis.text.y = element_text(color = "#000000", size = 14),
        legend.text = element_text(size = 12, color = "#000000"),
        legend.title = element_text(size = 12, face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(color = "#000000", fill = NA, size = 1.2),
        legend.position = "right") +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  scale_x_discrete(position = "top", limits = levels(JCPyV_PCR_VP23_PML_bubble_data_melt$variable))

ggsave("outputs/figs/fig4/JCPyV_PCR_VP23_PML_bubble_plot.tiff", plot = JCPyV_PCR_VP23_PML_bubble_plot, width = 12, height = 7, units = "cm", dpi = 300)


# SuperT
# Data processing
JCPyV_PCR_SuperT_PML_bubble_allcase <- PML_SuperT_list[[1]]
for (i in 2:length(PML_SuperT_list)) {
  caseNo_left <- ifelse(i == 2, 1, 0)
  JCPyV_PCR_SuperT_PML_bubble_allcase <- integrate_PML_bubble_samples(JCPyV_PCR_SuperT_PML_bubble_allcase, caseNo_left, PML_SuperT_list[[i]], i)
}
JCPyV_PCR_SuperT_PML_bubble_data <- JCPyV_PCR_SuperT_PML_bubble_allcase %>%
 select(-introns, -ID_IMR, -ID_293) %>%
 mutate(
  read_count_Case1 = as.numeric(read_count_Case1),
  read_count_Case2 = as.numeric(read_count_Case2),
  read_count_Case3 = as.numeric(read_count_Case3),
  read_count_Case4 = as.numeric(read_count_Case4),
  read_count_Case5 = as.numeric(read_count_Case5),
  read_count_Case6 = as.numeric(read_count_Case6),
  read_count_Case7 = as.numeric(read_count_Case7),
  read_count_Case8 = as.numeric(read_count_Case8)
 )
colnames(JCPyV_PCR_SuperT_PML_bubble_data) <- c("Length", "Case 1", "Case 2","Case 3","Case 4","Case 5","Case 6","Case 7","Case 8")
JCPyV_PCR_SuperT_PML_bubble_data <- JCPyV_PCR_SuperT_PML_bubble_data %>%
  mutate(Length = case_when(
    Length == 1359 ~ "1,359",
    Length == 1206 ~ "1,206",
    Length == 1053 ~ "1,053",
    TRUE ~ as.character(Length)
  ))
JCPyV_PCR_SuperT_PML_bubble_data_melt <- melt(JCPyV_PCR_SuperT_PML_bubble_data, id.vars = c("Length"))
JCPyV_PCR_SuperT_PML_bubble_data_melt$Length <- factor(JCPyV_PCR_SuperT_PML_bubble_data_melt$Length, levels = rev(unique(JCPyV_PCR_SuperT_PML_bubble_data_melt$Length)))

# Plot
JCPyV_PCR_SuperT_PML_bubble_plot = ggplot(JCPyV_PCR_SuperT_PML_bubble_data_melt, aes(x = variable, y = Length)) +
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) +
  scale_size_continuous(limits = c(1, 50000), range = c(1, 15), breaks = c(1, 100, 10000, 50000), labels = scales::comma) +
  labs(x= "", y = expression(bold("Product Length [bp]")), size = "SuperT\nRead Count", fill = "")  +
  theme(legend.key = element_blank(),
        axis.text.x = element_text(color = "#000000", size = 14, angle = 45, vjust = 0, hjust = 0),
        axis.text.y = element_text(color = "#000000", size = 14),
        legend.text = element_text(size = 12, color = "#000000"),
        legend.title = element_text(size = 12, face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(color = "#000000", fill = NA, size = 1.2),
        legend.position = "right") +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  scale_x_discrete(position = "top", limits = levels(JCPyV_PCR_SuperT_PML_bubble_data_melt$variable))

ggsave("outputs/figs/fig5/JCPyV_PCR_SuperT_PML_bubble_plot.tiff", plot = JCPyV_PCR_SuperT_PML_bubble_plot, width = 12, height = 7, units = "cm", dpi = 300)

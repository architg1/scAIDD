#Function to get haplotype counts for each barcode

getHaplotypeCounts <- function (sample_name, chrom_names, output_folder, haplotype_alleles) {
  print(paste0("Reading in allele frequencies on ",sample_name))
  allele_freq_all <- list()
  for (i in 1:length(chrom_names)) {
    cat(chrom_names[i])
    allele_freq_all[[i]] = read.delim(paste0(output.dir,"Allele_Freq/",sample_name,"_alleleFrequencies_chr", chrom_names[i], ".txt"), sep = "\t")
  }
  allele_freq_all <- do.call(rbind, allele_freq_all)
  # #Only keep counts from accepted cells (barcodes variable)
  # allele_freq_all$Barcode <- gsub(paste0("(................)-1"), paste0(samples[s], "_\\1"),allele_freq_all$Barcode)
  # allele_freq_all <- allele_freq_all %>% filter(Barcode %in% unlist(barcodes[which(samples == samples[s])]))
  
  print(paste0("Getting counts of haplotype 1 and haplotype 2 on ",sample_name))
  all_chr_coverage <- lapply(1:length(chrom_names),function(h) {
    allele_freq_chr <- allele_freq_all %>% filter(X.CHR == chrom_names[h])
    chr_coverage <- allele_freq_chr %>%
      rename_at(vars(starts_with("Count_")), function (r) {gsub("Count_", "", r)}) %>%
      pivot_longer(cols = 4:7, names_to = "bases", values_to = "Count") %>%
      group_by(POS) %>%
      left_join(allele_freq_chr %>%
                  group_by(POS) %>%
                  summarise(Coverage = sum(Good_depth)), by = "POS")
  })
  # saveRDS(all_chr_coverage,paste0(output_folder, "/all_chr_coverage.rds"))

  haplo_1_counts_all <- lapply(1:length(chrom_names), function(x) {
    all_chr_coverage[[x]] %>%
      inner_join(haplotype_alleles[[x]] %>% dplyr::select("POS","haplo_1"),
                 by = c("POS" = "POS", "bases" = "haplo_1"))
  })

  haplo_2_counts_all <- lapply(1:length(chrom_names), function(x) {
    all_chr_coverage[[x]] %>%
      inner_join(haplotype_alleles[[x]] %>% dplyr::select("POS","haplo_2"),
                 by = c("POS" = "POS", "bases" = "haplo_2"))
  })
  
  cell_both_haplotype_counts <- do.call(rbind, haplo_1_counts_all) %>%
    dplyr::group_by(Barcode) %>% summarise(Haplo_1_Count = sum(Count)) %>%
    left_join(do.call(rbind, haplo_2_counts_all) %>%
                dplyr::group_by(Barcode) %>% summarise(Haplo_2_Count = sum(Count)), by = "Barcode")
  return(cell_both_haplotype_counts)
}

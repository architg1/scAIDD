#Function to get haplotypes
getHaplotypes <- function (sample_name, chrom_names, vcf_prefix, output_folder) {
  vcf_all_chr <- list()
  for (i in 1:length(chrom_names)) {
    vcf_all_chr[[i]] <- VariantAnnotation::readVcf(paste0(vcf_prefix,chrom_names[i],".vcf"), "grch38")
  }
  # return(vcf_all_chr)
  
  haplotype_alleles <- lapply(1:length(vcf_all_chr), function(x) {
    genotype <- as.data.frame(VariantAnnotation::geno(vcf_all_chr[[x]])$GT, rownames = "name") %>% rownames_to_column(var = "name")
    position <- as.data.frame(GenomicRanges::granges(rowRanges(vcf_all_chr[[x]]))) %>% add_column("name" = names(vcf_all_chr[[x]]))
    haplotypes <- left_join(genotype, position, by = "name") %>%
      mutate(genotype = gsub(".*_","",name)) %>% 
      mutate(POS = start) %>%
      mutate(haplo_1 = ifelse(.[,2] == "0|1", gsub("/.", "", genotype),
                              ifelse(.[,2] == "1|0", gsub("./", "", genotype), NA))) %>%
      mutate(haplo_2 = ifelse(.[,2] == "1|0", gsub("/.", "", genotype),
                              ifelse(.[,2] == "0|1", gsub("./", "", genotype), NA)))
    return(haplotypes)
  }) 
  
  saveRDS(haplotype_alleles, paste0(output_folder, "/", sample_name, "_SNP_haplotype_alleles.rds"))
  return(haplotype_alleles)
}
  

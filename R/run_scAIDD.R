#Function to run scAIDD
require(dplyr)
require(tidyr)
require(tibble)
require(GenomicRanges)
require(VariantAnnotation)
run_scAIDD_10XscRNA <- function(sample_name,
                                bam_file,
                                output_folder,
                                loci_folder,
                                vcf_prefix,
                                chrom_names = paste0("chr", 1:22),
                                min_base_qual = 20,
                                min_map_qual = 0,
                                req_flag = 0,
                                filtered_flag = 0,
                                tenX = T,
                                allelecounter.exe="alleleCounter",
                                nbcores = 8) {
  
  system(paste0("mkdir -p ",output_folder,"Allele_Freq/"))

  
  if (!all(file.exists(paste0(output_folder, "/Allele_Freq/", sample_name, "_alleleFrequencies_chr", 1:length(chrom_names), ".txt")))) {
    chr_to_run <- which(!file.exists(paste0(output.dir, "/Allele_Freq/", "test", "_alleleFrequencies_chr", 1:length(chrom.names), ".txt")))
    print(paste0("Running allelecounter on ", paste0(chrom_names[chr_to_run], collapse = " "), " for ", sample_name))
    parallel::mclapply(chr_to_run,function(i) {
      alleleCounter(bam_file = bam_file,
                    loci_file = paste0(loci_folder, chrom_names[i], ".txt"),
                    output_file = paste0(output_folder, "/Allele_Freq/", sample_name, "_alleleFrequencies_chr", chrom_names[i], ".txt"),
                    min_base_qual = min_base_qual,
                    min_map_qual = min_map_qual,
                    req_flag = req_flag,
                    filtered_flag = filtered_flag,
                    tenX = tenX,
                    allelecounter.exe = allelecounter.exe)
    },mc.cores=nbcores)
  } else {
    print(paste0("Allelecounter already run on ", paste0(chrom_names, collapse = " "), " for ", sample_name))
  }
  
  if (!file.exists(paste0(output_folder, "/", sample_name, "_SNP_haplotype_alleles.rds"))) {
    print(paste0("Getting genotypes of ", sample_name))
    haplotype_alleles <- getHaplotypes(sample_name,
                                       chrom_names, 
                                       vcf_prefix, 
                                       output_folder)
  } else {
    print(paste0("Loading genotypes of ", sample_name))
    haplotype_alleles <- readRDS(paste0(output_folder, "/", sample_name, "_SNP_haplotype_alleles.rds"))
  }
  
  cell_both_haplotype_counts <- getHaplotypeCounts(sample_name,
                                                   chrom_names, 
                                                   output_folder,
                                                   haplotype_alleles = haplotype_alleles)
  
  saveRDS(cell_both_haplotype_counts, paste0(output_folder, "/", sample_name, "_cell_both_haplotype_counts.rds"))
}

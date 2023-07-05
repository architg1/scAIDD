#Function to run scAIDD

run_scAIDD_10XscRNA <- function(sample_name,
                                bam_file,
                                output_folder,
                                loci_folder,
                                chrom_names = paste0("chr", 1:22),
                                min_base_qual = 20,
                                min_map_qual = 0,
                                req_flag = 0,
                                filtered_flag = 0,
                                tenX = T,
                                allelecounter.exe="alleleCounter",
                                nbcores = 8) {
  #system(paste0("mkdir -p ",output_folder,"Allele_Freq/",s))
  #setwd(paste0(output.dir,"Allele_Freq/",s))

  print(paste0("Running allele counter on ", bam_file))

  parallel::mclapply(1:length(chrom_names),function(i) {
    alleleCounter(bam_file = bam_file,
                  loci_file = paste0(loci_folder, chrom_names[i], ".txt"),
                  output_file = paste0(output_folder, "/", sample_name, "_alleleFrequencies_chr", i, ".txt"),
                  min_base_qual = min_base_qual,
                  min_map_qual = min_map_qual,
                  req_flag = req_flag,
                  filtered_flag = filtered_flag,
                  tenX = tenX,
                  allelecounter.exe = allelecounter.exe)
  },mc.cores=nbcores)
}

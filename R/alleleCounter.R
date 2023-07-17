#Function to run allelecounter

alleleCounter <- function (bam_file, output_file, loci_file, min_base_qual, min_map_qual,
                           req_flag = 0, filtered_flag = 0, tenX, allelecounter_path) {
  if (!file.exists(bam_file) || file.info(bam_file)$size==0) {warning('seq.file does not seem to exist or is empty'); return()}

  cmd = paste(allelecounter_path,
              "-b", bam_file,
              "-l", loci_file,
              "-o", output_file,
              "-m", min_base_qual,
              "-q", min_map_qual,
              "-f ", req_flag,
              "-F ", filtered_flag)
  # alleleCounter >= v4.0.0 is sped up considerably when run in dense-snp mode
  counter_version = system(paste0(allelecounter_path, " --version"), intern = T)
  if (as.integer(substr(x = counter_version, start = 1, stop = 1)) >= 4) {
    cmd = paste0(cmd, " --dense-snps")
  }
  if (tenX) {
    cmd = paste0(cmd, " -x")
  }
  system(cmd, wait = T)
}


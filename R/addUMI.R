#Function to add UMI to 10X scDNA bam which is needed for allelecounter

addUMI <- function (bam_file, output_folder, samtools_path) {
  bam_name = paste0("/", gsub(".*/", "", bam_file))
  print(paste0("Converting bam to sam"))
  system(paste0(samtools_path, " view -h ", bam_file, " > ", output_folder, gsub("\\.bam", ".sam", bam_name)))
  print(paste0("Adding umi"))
  system(paste0("awk '{gsub(\"RG:Z:\",\"UB:Z:\" NR \"\tRG:Z:\",$0);print}' ", output_folder, gsub("\\.bam", ".sam", bam_name), " > ", output_folder, gsub("\\.bam", "_umi.sam", bam_name))) #Adds line number as "UMI" as 10X CNVkit does not have UMIs
  print(paste0("Converting umi sam to bam"))
  system(paste0(samtools_path, " view -b ", output_folder, gsub("\\.bam", "_umi.sam", bam_name), " > ", output_folder, gsub("\\.bam", "_umi.bam", bam_name)))
  print(paste0("Indexing umi bam"))
  system(paste0(samtools_path, " index ",  output_folder, gsub("\\.bam", "_umi.bam", bam_name), " ", output_folder, gsub("\\.bam", "_umi.bai", bam_name)))
  print(paste0("Cleaning up sam"))
  system(paste0("rm ", output_folder,  gsub("\\.bam", ".sam", bam_name), " ", output_folder, gsub("\\.bam", "_umi.sam", bam_name)))
}



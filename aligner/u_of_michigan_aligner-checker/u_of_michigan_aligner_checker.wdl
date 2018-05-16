import "https://raw.githubusercontent.com/DataBiosphere/topmed-workflows/feature/umich-align-checker-wdl/aligner/umich-aligner/u_of_michigan_aligner.wdl" as TopMed_aligner
import "https://raw.githubusercontent.com/DataBiosphere/topmed-workflows/1.7.0/aligner/functional-equivalence-checker/topmed-alignment-checker.wdl" as checker

workflow checkerWorkflow {
  Int expectedNumofReads
  String docker_image

  File input_crai_file
  File input_cram_file

  File ref_alt
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac

  File ref_fasta
  File ref_fasta_index

  File dbSNP_vcf
  File dbSNP_vcf_index

  Int? increase_disk_size
  Int additional_disk = select_first([increase_disk_size, 20])
  Float bwa_disk_multiplier = 2.5

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Float ref_extra_size = size(ref_alt, "GB") + size(ref_bwt, "GB") + size(ref_pac, "GB") + size(ref_ann, "GB") + size(ref_amb, "GB") + size(ref_sa, "GB")
  Float dbsnp_size = size(dbSNP_vcf, "GB") + size(dbSNP_vcf_index, "GB")
  Float cram_size = size(input_cram_file, "GB") + size(input_crai_file, "GB")


 call TopMed_aligner.TopMedAligner as aligner { 
   input:

     input_crai_file = input_crai_file,
     input_cram_file = input_cram_file,
     
     docker_image = docker_image,

     ref_fasta = ref_fasta,
     ref_fasta_index = ref_fasta_index,
     ref_alt = ref_alt,
     ref_bwt = ref_bwt,
     ref_sa = ref_sa,
     ref_amb = ref_amb,
     ref_ann = ref_ann,
     ref_pac = ref_pac,

     dbSNP_vcf = dbSNP_vcf,
     dbSNP_vcf_index = dbSNP_vcf_index

 }


 call checker.checkerTask { input: inputCRAMFile=aligner.aligner_output, expectedNumofReads=expectedNumofReads, docker_image=docker_image }
}

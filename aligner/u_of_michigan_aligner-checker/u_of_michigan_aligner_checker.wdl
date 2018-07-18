import "https://raw.githubusercontent.com/DataBiosphere/topmed-workflows/feature/rm-dos-demo/aligner/u_of_michigan_aligner/u_of_michigan_aligner.wdl" as TopMed_aligner
import "https://raw.githubusercontent.com/DataBiosphere/topmed-workflows/feature/rm-dos-demo/aligner/u_of_michigan_aligner-checker/u_of_michigan_aligner_checker_calculation.wdl" as checker

workflow checkerWorkflow {
  String docker_image

  File input_crai_file
  File input_cram_file

  File inputTruthCRAMFile

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

  # Optional inputs for aligner checker task
  Float? est_size_inputCRAMFile
  Float? est_size_inputTruthCRAMFile
  Float? est_ref_size

  # Optional inputs for aligner task
  Float? est_ref_extra_size
  Float? est_dbsnp_size

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
     dbSNP_vcf_index = dbSNP_vcf_index,

     est_ref_size = est_ref_size,
     est_ref_extra_size = est_ref_extra_size,
     est_dbsnp_size = est_dbsnp_size,
     est_cram_size = est_size_inputCRAMFile
 }


 call checker.checkerTask { 
    input: 
        inputCRAMFile = aligner.aligner_output, 
        inputTruthCRAMFile = inputTruthCRAMFile,
        referenceFile = ref_fasta,
        docker_image = docker_image,

        est_size_inputCRAMFile = est_size_inputCRAMFile,
        est_size_inputTruthCRAMFile = est_size_inputTruthCRAMFile,
        est_ref_size = est_ref_size
    }
}

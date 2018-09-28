import "https://raw.githubusercontent.com/DataBiosphere/topmed-workflows/1.30.0/aligner/functional-equivalence-wdl/FunctionalEquivalence.wdl" as TopMed_aligner
import "https://raw.githubusercontent.com/DataBiosphere/topmed-workflows/1.30.0/aligner/functional-equivalence-checker/topmed-alignment-checker.wdl" as checker

workflow checkerWorkflow {
  Int expectedNumofReads
  String docker_image

  File wgs_evaluation_interval_list
  File wgs_coverage_interval_list

  String sample_name
  String base_file_name
  Array[File] flowcell_unmapped_bams
  String unmapped_bam_suffix

  File wgs_calling_interval_list
  Int haplotype_scatter_count
  Int break_bands_at_multiples_of
  Int? read_length

  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_alt
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac

  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  Int preemptible_tries
  Int agg_preemptible_tries

 call TopMed_aligner.PairedEndSingleSampleWorkflow as aligner { 
   input: 
     docker_image = docker_image,

     wgs_evaluation_interval_list = wgs_evaluation_interval_list,
     wgs_coverage_interval_list = wgs_coverage_interval_list,

     sample_name = sample_name,
     base_file_name = base_file_name,
     flowcell_unmapped_bams = flowcell_unmapped_bams,
     unmapped_bam_suffix = unmapped_bam_suffix,

     wgs_calling_interval_list = wgs_calling_interval_list,
     haplotype_scatter_count = haplotype_scatter_count,
     break_bands_at_multiples_of = break_bands_at_multiples_of,
     read_length = read_length,

     ref_fasta = ref_fasta,
     ref_fasta_index = ref_fasta_index,
     ref_dict = ref_dict,
     ref_alt = ref_alt,
     ref_bwt = ref_bwt,
     ref_sa = ref_sa,
     ref_amb = ref_amb,
     ref_ann = ref_ann,
     ref_pac = ref_pac,

     dbSNP_vcf = dbSNP_vcf,
     dbSNP_vcf_index = dbSNP_vcf_index,
     known_indels_sites_VCFs = known_indels_sites_VCFs,
     known_indels_sites_indices = known_indels_sites_indices,
     
     preemptible_tries = preemptible_tries,
     agg_preemptible_tries = agg_preemptible_tries

  }


 call checker.checkerTask { input: inputCRAMFile=aligner.output_cram, expectedNumofReads=expectedNumofReads, docker_image=docker_image}
}

import "https://raw.githubusercontent.com/DataBiosphere/topmed-workflows/1.14.0/variant-caller/variant-caller-wdl/topmed_freeze3_calling.wdl" as TopMed_variantcaller
import "https://raw.githubusercontent.com/DataBiosphere/topmed-workflows/1.14.0/variant-caller/variant-caller-wdl-checker/topmed-variantcaller-checker.wdl" as checker

workflow checkerWorkflow {
  File inputTruthVCFFile

  String docker_image

  Array[File] input_crai_files
  Array[File] input_cram_files

  File ref_1000G_omni2_5_b38_sites_PASS_vcf_gz
  File ref_1000G_omni2_5_b38_sites_PASS_vcf_gz_tbi
  File chr10_vcf
  File chr11_KI270927v1_alt_vcf
  File chr11_vcf
  File chr12_vcf
  File chr13_vcf
  File chr14_GL000009v2_random_vcf
  File chr14_KI270846v1_alt_vcf
  File chr14_vcf
  File chr15_vcf
  File chr16_vcf
  File chr17_KI270857v1_alt_vcf
  File chr17_KI270862v1_alt_vcf
  File chr17_KI270909v1_alt_vcf
  File chr17_vcf
  File chr18_vcf
  File chr19_KI270938v1_alt_vcf
  File chr19_vcf
  File chr1_KI270706v1_random_vcf
  File chr1_KI270766v1_alt_vcf
  File chr1_vcf
  File chr20_vcf
  File chr21_vcf
  File chr22_KI270879v1_alt_vcf
  File chr22_KI270928v1_alt_vcf
  File chr22_vcf
  File chr2_KI270773v1_alt_vcf
  File chr2_KI270894v1_alt_vcf
  File chr2_vcf
  File chr3_vcf
  File chr4_GL000008v2_random_vcf
  File chr4_vcf
  File chr5_vcf
  File chr6_vcf
  File chr7_KI270803v1_alt_vcf
  File chr7_vcf
  File chr8_KI270821v1_alt_vcf
  File chr8_vcf
  File chr9_vcf
  File chrUn_KI270742v1_vcf
  File chrX_vcf
  File ref_dbsnp_142_b38_vcf_gz
  File ref_dbsnp_142_b38_vcf_gz_tbi
  File ref_dbsnp_All_vcf_gz
  File ref_dbsnp_All_vcf_gz_tbi
  File ref_hapmap_3_3_b38_sites_vcf_gz
  File ref_hapmap_3_3_b38_sites_vcf_gz_tbi
  File ref_hs38DH_bs_umfa
  File ref_hs38DH_dict
  File ref_hs38DH_fa
  File ref_hs38DH_fa_alt
  File ref_hs38DH_fa_amb
  File ref_hs38DH_fa_ann
  File ref_hs38DH_fa_bwt
  File ref_hs38DH_fa_fai
  File ref_hs38DH_fa_pac
  File ref_hs38DH_fa_sa
  File ref_hs38DH_winsize100_gc


  call TopMed_variantcaller.TopMedVariantCaller as variantcaller {
    input:
      input_crai_files = input_crai_files,
      input_cram_files = input_cram_files,

      docker_image = docker_image,

      ref_1000G_omni2_5_b38_sites_PASS_vcf_gz = ref_1000G_omni2_5_b38_sites_PASS_vcf_gz,
      ref_1000G_omni2_5_b38_sites_PASS_vcf_gz_tbi = ref_1000G_omni2_5_b38_sites_PASS_vcf_gz_tbi,
      chr10_vcf = chr10_vcf,
      chr11_KI270927v1_alt_vcf = chr11_KI270927v1_alt_vcf,
      chr11_vcf = chr11_vcf,
      chr12_vcf = chr12_vcf,
      chr13_vcf = chr13_vcf,
      chr14_GL000009v2_random_vcf = chr14_GL000009v2_random_vcf,
      chr14_KI270846v1_alt_vcf = chr14_KI270846v1_alt_vcf,
      chr14_vcf = chr14_vcf,
      chr15_vcf = chr15_vcf,
      chr16_vcf = chr16_vcf,
      chr17_KI270857v1_alt_vcf = chr17_KI270857v1_alt_vcf,
      chr17_KI270862v1_alt_vcf = chr17_KI270862v1_alt_vcf,
      chr17_KI270909v1_alt_vcf = chr17_KI270909v1_alt_vcf,
      chr17_vcf = chr17_vcf,
      chr18_vcf = chr18_vcf,
      chr19_KI270938v1_alt_vcf = chr19_KI270938v1_alt_vcf,
      chr19_vcf = chr19_vcf,
      chr1_KI270706v1_random_vcf = chr1_KI270706v1_random_vcf,
      chr1_KI270766v1_alt_vcf = chr1_KI270766v1_alt_vcf,
      chr1_vcf = chr1_vcf,
      chr20_vcf = chr20_vcf,
      chr21_vcf = chr21_vcf,
      chr22_KI270879v1_alt_vcf = chr22_KI270879v1_alt_vcf,
      chr22_KI270928v1_alt_vcf = chr22_KI270928v1_alt_vcf,
      chr22_vcf = chr22_vcf,
      chr2_KI270773v1_alt_vcf = chr2_KI270773v1_alt_vcf,
      chr2_KI270894v1_alt_vcf = chr2_KI270894v1_alt_vcf,
      chr2_vcf = chr2_vcf,
      chr3_vcf = chr3_vcf,
      chr4_GL000008v2_random_vcf = chr4_GL000008v2_random_vcf,
      chr4_vcf = chr4_vcf,
      chr5_vcf = chr5_vcf,
      chr6_vcf = chr6_vcf,
      chr7_KI270803v1_alt_vcf = chr7_KI270803v1_alt_vcf,
      chr7_vcf = chr7_vcf,
      chr8_KI270821v1_alt_vcf = chr8_KI270821v1_alt_vcf,
      chr8_vcf = chr8_vcf,
      chr9_vcf = chr9_vcf,
      chrUn_KI270742v1_vcf = chrUn_KI270742v1_vcf,
      chrX_vcf = chrX_vcf,
      ref_dbsnp_142_b38_vcf_gz = ref_dbsnp_142_b38_vcf_gz,
      ref_dbsnp_142_b38_vcf_gz_tbi = ref_dbsnp_142_b38_vcf_gz_tbi,
      ref_dbsnp_All_vcf_gz = ref_dbsnp_All_vcf_gz,
      ref_dbsnp_All_vcf_gz_tbi = ref_dbsnp_All_vcf_gz_tbi,
      ref_hapmap_3_3_b38_sites_vcf_gz = ref_hapmap_3_3_b38_sites_vcf_gz,
      ref_hapmap_3_3_b38_sites_vcf_gz_tbi = ref_hapmap_3_3_b38_sites_vcf_gz_tbi,
      ref_hs38DH_bs_umfa = ref_hs38DH_bs_umfa,
      ref_hs38DH_dict = ref_hs38DH_dict,
      ref_hs38DH_fa = ref_hs38DH_fa,
      ref_hs38DH_fa_alt = ref_hs38DH_fa_alt,
      ref_hs38DH_fa_amb = ref_hs38DH_fa_amb,
      ref_hs38DH_fa_ann = ref_hs38DH_fa_ann,
      ref_hs38DH_fa_bwt = ref_hs38DH_fa_bwt,
      ref_hs38DH_fa_fai = ref_hs38DH_fa_fai,
      ref_hs38DH_fa_pac = ref_hs38DH_fa_pac,
      ref_hs38DH_fa_sa = ref_hs38DH_fa_sa,
      ref_hs38DH_winsize100_gc = ref_hs38DH_winsize100_gc
  }

  call checker.checkerTask { 
      input: 
          inputTruthVCFFile = inputTruthVCFFile,
          inputTestVCFFile = variantcaller.topmed_variant_caller_output, 
          docker_image = docker_image
  }
}


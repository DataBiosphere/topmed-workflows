## This is the TopMed variant caller workflow WDL for the workflow code located here:
## https://github.com/BD2KGenomics/topmed_freeze3_calling/tree/feature/dockerize-workflow
## (forked from https://github.com/statgen/topmed_freeze3_calling)
##
## It uses a Docker image built with software tools that can reproduce 
## variant calls compatible to TOPMed Freeze 3a
##
## NOTE: This workflow assumes that input CRAM files have been built with the b38
## human reference genome. In particular for the TopMed CRAM files the 
## reference genome files to use are located here:
## ftp://share.sph.umich.edu/gotcloud/ref/hs38DH-db142-v1.tgz
##
##

import "https://raw.githubusercontent.com/DataBiosphere/topmed-workflows/feature/checker-vcf/variant-caller/variant-caller-wdl/topmed_freeze3_calling.wdl" as TopMed_variantcaller

import ""https://raw.githubusercontent.com/DataBiosphere/topmed-workflows/feature/checker-vcf/variant-caller/variant-caller-wdl-checker/topmed-variantcaller-checker.wdk" as checker

workflow checkerWorkflow_vcf {

  Int expectedNumofReads
  String docker_image

  Array[File] input_crai_files
  Array[File] input_cram_files

  Float reference_files_size
  String docker_image

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
      input_crai_files = input_crai_files
      input_cram_files = input_cram_files

      reference_files_size = reference_files_size
      docker_image = docker_image

      ref_1000G_omni2_5_b38_sites_PASS_vcf_gz = ref_1000G_omni2_5_b38_sites_PASS_vcf_gz
      ref_1000G_omni2_5_b38_sites_PASS_vcf_gz_tbi = ref_1000G_omni2_5_b38_sites_PASS_vcf_gz_tbi
      chr10_vcf = chr10_vcf
      chr11_KI270927v1_alt_vcf = chr11_KI270927v1_alt_vcf
      chr11_vcf = chr11_vcf
      chr12_vcf = chr12_vcf
      chr13_vcf = chr13_vcf
      chr14_GL000009v2_random_vcf = chr14_GL000009v2_random_vcf
      chr14_KI270846v1_alt_vcf = chr14_KI270846v1_alt_vcf
      chr14_vcf = chr14_vcf
      chr15_vcf = chr15_vcf
      chr16_vcf = chr16_vcf
      chr17_KI270857v1_alt_vcf = chr17_KI270857v1_alt_vcf
      chr17_KI270862v1_alt_vcf = chr17_KI270862v1_alt_vcf
      chr17_KI270909v1_alt_vcf = chr17_KI270909v1_alt_vcf
      chr17_vcf = chr17_vcf
      chr18_vcf = chr18_vcf
      chr19_KI270938v1_alt_vcf = chr19_KI270938v1_alt_vcf
      chr19_vcf = chr19_vcf
      chr1_KI270706v1_random_vcf
      chr1_KI270766v1_alt_vcf = chr1_KI270766v1_alt_vcf
      chr1_vcf = chr1_vcf
      chr20_vcf = chr20_vcf
      chr21_vcf = chr21_vcf
      chr22_KI270879v1_alt_vcf = chr22_KI270879v1_alt_vcf
      chr22_KI270928v1_alt_vcf = chr22_KI270928v1_alt_vcf
      chr22_vcf = chr22_vcf
      File chr2_KI270773v1_alt_vcf
      chr2_KI270894v1_alt_vcf = chr2_KI270894v1_alt_vcf
      chr2_vcf = chr2_vcf
      chr3_vcf = chr3_vcf
      chr4_GL000008v2_random_vcf = chr4_GL000008v2_random_vcf
      chr4_vcf = chr4_vcf
      chr5_vcf = chr5_vcf
      chr6_vcf = chr6_vcf
      chr7_KI270803v1_alt_vcf = chr7_KI270803v1_alt_vcf
      chr7_vcf = chr7_vcf
      chr8_KI270821v1_alt_vcf = chr8_KI270821v1_alt_vcf
      chr8_vcf = chr8_vcf
      chr9_vcf = chr9_vcf
      chrUn_KI270742v1_vcf = chrUn_KI270742v1_vcf
      chrX_vcf = chrX_vcf
      ref_dbsnp_142_b38_vcf_gz
      ref_dbsnp_142_b38_vcf_gz_tbi = ref_dbsnp_142_b38_vcf_gz_tbi
      ref_dbsnp_All_vcf_gz = ref_dbsnp_All_vcf_gz
      ref_dbsnp_All_vcf_gz_tbi = ref_dbsnp_All_vcf_gz_tbi
      ref_hapmap_3_3_b38_sites_vcf_gz = ref_hapmap_3_3_b38_sites_vcf_gz
      ref_hapmap_3_3_b38_sites_vcf_gz_tbi = ref_hapmap_3_3_b38_sites_vcf_gz_tbi
      ref_hs38DH_bs_umfa = ref_hs38DH_bs_umfa
      ref_hs38DH_dict = ref_hs38DH_dict
      ref_hs38DH_fa = ref_hs38DH_fa
      ref_hs38DH_fa_alt = ref_hs38DH_fa_alt
      ref_hs38DH_fa_amb = ref_hs38DH_fa_amb
      ref_hs38DH_fa_ann = ref_hs38DH_fa_ann
      ref_hs38DH_fa_bwt = ref_hs38DH_fa_bwt
      ref_hs38DH_fa_fai = ref_hs38DH_fa_fai
      ref_hs38DH_fa_pac = ref_hs38DH_fa_pac
      ref_hs38DH_fa_sa = ref_hs38DH_fa_sa
      ref_hs38DH_winsize100_gc = ref_hs38DH_winsize100_gc


  call sumCRAMSizes {
    input:
      input_crams = input_cram_files,
      input_crais = input_crai_files,
      disk_size = reference_files_size,
      docker_image = docker_image
  }


  call variantCalling {

     input:
      input_crais = input_crai_files,
      input_crams = input_cram_files,
      disk_size = sumCRAMSizes.total_size + reference_files_size,
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
  
  output {
      File topmed_variant_caller_output = variantCalling.topmed_variant_caller_output_file
  }
}


  task sumCRAMSizes {
    Array[File] input_crams
    Array[File] input_crais
    Float disk_size
    String docker_image
  
    command {
      python <<CODE
      import os
      import functools

      cram_string = "${ sep=',' input_crams }"
      cram_list = cram_string.split(',')
      if len(cram_list) > 1:
          crams_size = functools.reduce((lambda x, y: os.stat(x).st_size + os.stat(y).st_size), cram_list)
      else:
          crams_size = os.stat(cram_list[0]).st_size

      crai_string = "${ sep=',' input_crais }"
      crai_list = crai_string.split(',')
      if len(crai_list) > 1:      
          crais_size = functools.reduce((lambda x, y: os.stat(x).st_size + os.stat(y).st_size), crai_list)
      else:
          crais_size = os.stat(crai_list[0]).st_size
         
      total_size = crams_size + crais_size
      # Shift right by 30 bits to get Gigabyte size of files
      total_size = (total_size >> 30)

      # Bump the size up 1 GB in case the total size is less than 1 GB
      print total_size + 1

      CODE
    }
    runtime {
      memory: "10 GB"
      cpu: "16"
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      zones: "us-central1-a us-central1-b us-east1-d us-central1-c us-central1-f us-east1-c"
      docker: docker_image
    }
    output {
      Float total_size = read_float(stdout())
    }
  }

  

  task variantCalling {
     # The CRAM index files are listed as an input because they are required
     # by various tools, e.g. Samtools. They should be in the same location
     # as the CRAM files when specified in the input JSON
     Array[File] input_crais
     Array[File] input_crams

     Float disk_size
     String docker_image

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


     command {
      python <<CODE

      from __future__ import print_function
      import csv
      import os
      from shutil import copyfile     

      if not os.path.exists('/root/topmed_freeze3_calling/data/'):
          os.makedirs('/root/topmed_freeze3_calling/data/')

      # Create the index file that lists the CRAM files and their location
      # We do not need to do this for the CRAM index files because the 
      # tools assume they are located in the same place as the CRAM files
      cram_string = "${ sep=',' input_crams }"
      crams_list = cram_string.split(',')

      tsv_crams_rows = []
      for cram_file in crams_list:
          # Get the Cromwell location of the CRAM file
          # The worklow will be able to access them 
          # since the Cromwell path is mounted in the
          # docker run commmand that Cromwell sets up
          base_name = os.path.basename(cram_file)
          base_name_wo_extension = base_name.split('.')[0]
          tsv_crams_rows.append([base_name_wo_extension, cram_file, '0.000'])

      # Remove the old PED file; we will not use a PED file?
      open('/root/topmed_freeze3_calling/data/trio_data.ped', 'w+').close()

      with open('/root/topmed_freeze3_calling/data/trio_data.index', 'w+') as tsv_index_file:
          writer = csv.writer(tsv_index_file, delimiter = '\t')
          for cram_info in tsv_crams_rows:
              writer.writerow(cram_info)

      CODE


      set -o pipefail
      set -e

      #echo each line of the script to stdout so we can see what is happening
      set -o xtrace
      #to turn of echo do 'set +o xtrace'

      # Make sure the directory where the reference files are supposed to be
      # located exists in the container
      mkdir -p /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38

      # Create a symlink from the where the workflow expects the reference files
      # to the Cromwell location of the reference files 
      ln -s ${ref_1000G_omni2_5_b38_sites_PASS_vcf_gz}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1000G_omni2.5.b38.sites.PASS.vcf.gz
      ln -s ${ref_1000G_omni2_5_b38_sites_PASS_vcf_gz_tbi}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1000G_omni2.5.b38.sites.PASS.vcf.gz.tbi
      ln -s ${chr10_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr10.vcf
      ln -s ${chr11_KI270927v1_alt_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr11_KI270927v1_alt.vcf
      ln -s ${chr11_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr11.vcf
      ln -s ${chr12_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr12.vcf
      ln -s ${chr13_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr13.vcf
      ln -s ${chr14_GL000009v2_random_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr14_GL000009v2_random.vcf
      ln -s ${chr14_KI270846v1_alt_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr14_KI270846v1_alt.vcf
      ln -s ${chr14_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr14.vcf
      ln -s ${chr15_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr15.vcf
      ln -s ${chr16_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr16.vcf
      ln -s ${chr17_KI270857v1_alt_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr17_KI270857v1_alt.vcf
      ln -s ${chr17_KI270862v1_alt_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr17_KI270862v1_alt.vcf
      ln -s ${chr17_KI270909v1_alt_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr17_KI270909v1_alt.vcf
      ln -s ${chr17_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr17.vcf
      ln -s ${chr18_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr18.vcf
      ln -s ${chr19_KI270938v1_alt_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr19_KI270938v1_alt.vcf
      ln -s ${chr19_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr19.vcf
      ln -s ${chr1_KI270706v1_random_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr1_KI270706v1_random.vcf
      ln -s ${chr1_KI270766v1_alt_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr1_KI270766v1_alt.vcf
      ln -s ${chr1_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr1.vcf
      ln -s ${chr20_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr20.vcf
      ln -s ${chr21_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr21.vcf
      ln -s ${chr22_KI270879v1_alt_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr22_KI270879v1_alt.vcf
      ln -s ${chr22_KI270928v1_alt_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr22_KI270928v1_alt.vcf
      ln -s ${chr22_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr22.vcf
      ln -s ${chr2_KI270773v1_alt_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr2_KI270773v1_alt.vcf
      ln -s ${chr2_KI270894v1_alt_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr2_KI270894v1_alt.vcf
      ln -s ${chr2_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr2.vcf
      ln -s ${chr3_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr3.vcf
      ln -s ${chr4_GL000008v2_random_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr4_GL000008v2_random.vcf
      ln -s ${chr4_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr4.vcf
      ln -s ${chr5_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr5.vcf
      ln -s ${chr6_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr6.vcf
      ln -s ${chr7_KI270803v1_alt_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr7_KI270803v1_alt.vcf
      ln -s ${chr7_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr7.vcf
      ln -s ${chr8_KI270821v1_alt_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr8_KI270821v1_alt.vcf
      ln -s ${chr8_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr8.vcf
      ln -s ${chr9_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr9.vcf
      ln -s ${chrUn_KI270742v1_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chrUn_KI270742v1.vcf
      ln -s ${chrX_vcf}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/1kg.pilot_release.merged.indels.sites.hg38.chrX.vcf
      ln -s ${ref_dbsnp_142_b38_vcf_gz}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/dbsnp_142.b38.vcf.gz
      ln -s ${ref_dbsnp_142_b38_vcf_gz_tbi}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/dbsnp_142.b38.vcf.gz.tbi
      ln -s ${ref_dbsnp_All_vcf_gz}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/dbsnp.All.vcf.gz
      ln -s ${ref_dbsnp_All_vcf_gz_tbi}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/dbsnp.All.vcf.gz.tbi
      ln -s ${ref_hapmap_3_3_b38_sites_vcf_gz}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hapmap_3.3.b38.sites.vcf.gz
      ln -s ${ref_hapmap_3_3_b38_sites_vcf_gz_tbi}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hapmap_3.3.b38.sites.vcf.gz.tbi
      ln -s ${ref_hs38DH_bs_umfa}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hs38DH-bs.umfa
      ln -s ${ref_hs38DH_dict}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hs38DH.dict
      ln -s ${ref_hs38DH_fa}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hs38DH.fa
      ln -s ${ref_hs38DH_fa_alt}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hs38DH.fa.alt
      ln -s ${ref_hs38DH_fa_amb}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hs38DH.fa.amb
      ln -s ${ref_hs38DH_fa_ann}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hs38DH.fa.ann
      ln -s ${ref_hs38DH_fa_bwt}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hs38DH.fa.bwt
      ln -s ${ref_hs38DH_fa_fai}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hs38DH.fa.fai
      ln -s ${ref_hs38DH_fa_pac}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hs38DH.fa.pac
      ln -s ${ref_hs38DH_fa_sa}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hs38DH.fa.sa
      ln -s ${ref_hs38DH_winsize100_gc}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hs38DH.winsize100.gc

      # Change to the directory in the container where the pipeline code is located
      pushd /root/topmed_freeze3_calling

      WORKING_DIR='./' 

      # Put the correct location of the index file into the global config file
      sed -i '/.*our $refDir.*/ c\our $refDir = "$FindBin::Bin\/..\/data\/local.org\/ref\/gotcloud.ref\/hg38";' "$WORKING_DIR"/scripts/gcconfig.pm
      sed -i '/.*our $ref = "$refDir.*/ c\our $ref = "$refDir\/hs38DH.fa";' "$WORKING_DIR"/scripts/gcconfig.pm
      sed -i '/.*our $dbsnp.*/ c\our $dbsnp = "$refDir\/dbsnp_142.b38.vcf.gz";' "$WORKING_DIR"/scripts/gcconfig.pm
      sed -i '/.*our $hapmapvcf.*/ c\our $hapmapvcf = "$refDir\/hapmap_3.3.b38.sites.vcf.gz";' "$WORKING_DIR"/scripts/gcconfig.pm
      sed -i '/.*our $omnivcf.*/ c\our $omnivcf = "$refDir\/1000G_omni2.5.b38.sites.PASS.vcf.gz";' "$WORKING_DIR"/scripts/gcconfig.pm


      # Put the correct location of the output directory into the global config file
      #sed -i '/.*our $out =.*/ c\our $out = "/root/topmed_freeze3_calling/out";' "$WORKING_DIR"/scripts/gcconfig.pm
      # Put the correct location of references into the config file
      sed -i '/.*our $md5 =.*/ c\our $md5 = "\/data\/local.org\/ref\/gotcloud.ref\/md5\/%2s\/%s\/%s";' "$WORKING_DIR"/scripts/config.pm
      sed -i '/.*our $ref =.*/ c\our $ref = "\/data\/local.org\/ref\/gotcloud.ref\/hg38\/hs38DH.fa";' "$WORKING_DIR"/scripts/config.pm


      echo "Running step1 - detect and merge variants"
      echo "Running step1 - detect and merge variants - removing old output dir if it exists"
      if [ -d "$WORKING_DIR"/out ]; then rm -Rf "$WORKING_DIR"/out; fi
      echo "Running step1 - detect and merge variants - generating Makefile"
      perl "$WORKING_DIR"/scripts/step1-detect-and-merge-variants.pl $(seq 1 22 | xargs -n 1 -I% echo chr%) chrX
      echo "Running step1 - detect and merge variants - running Makefile"
      make SHELL='/bin/bash' -f "$WORKING_DIR"/out/aux/Makefile -j 23
      

      echo "Running step2 - joing genotyping"
      echo "Running step2 - joing genotyping - removing old output dir if it exists"
      if [ -d "$WORKING_DIR"/paste ]; then rm -Rf "$WORKING_DIR"/paste; fi
      echo "Running step2 - joing genotyping - generating Makefile"
      perl "$WORKING_DIR"/scripts/step2-joint-genotyping.pl $(seq 1 22 | xargs -n 1 -I% echo chr%) chrX
      echo "Running step2 - joing genotyping - running Makefile"
      MAKEFILE_NAME="chrchr"$(seq -s '_chr' 1 22)_chrX".Makefile"
      make SHELL='/bin/bash' -f "$WORKING_DIR"/out/paste/"$MAKEFILE_NAME" -j 23

      # Pop back to the original working directory; on FireCloud this will be a
      # special directory 
      popd

      # Tar up the output directories into the output file provided in the input JSON
      tar -zcvf topmed_variant_caller_output.tar.gz /root/topmed_freeze3_calling/out/paste/ /root/topmed_freeze3_calling/out/aux/individual/

    }
     output {
      File topmed_variant_caller_output_file = "topmed_variant_caller_output.tar.gz"
    }
   runtime {
      memory: "10 GB"
      cpu: "16"
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      zones: "us-central1-a us-central1-b us-east1-d us-central1-c us-central1-f us-east1-c"
      docker: docker_image
    }
  }


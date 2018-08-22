#import "https://raw.githubusercontent.com/DataBiosphere/topmed-workflows/1.23.0/variant-caller/variant-caller-wdl/calculate_contamination.wdl" as getDNAContamination
import "/Users/waltershands/Documents/UCSC/gitroot/topmed-workflows/variant-caller/variant-caller-wdl/calculate_contamination.wdl" as getDNAContamination

## This is the U of Michigan variant caller workflow WDL for the workflow code located here:
## https://github.com/statgen/topmed_freeze3_calling
##
## It uses a Docker image built with software tools that can reproduce 
## variant calls compatible to TopMed Freeze 3a
##
## NOTE: This workflow assumes that input CRAM files have been built with the b38
## human reference genome. In particular for the TopMed CRAM files the 
## reference genome files to use are located here:
## ftp://share.sph.umich.edu/gotcloud/ref/hs38DH-db142-v1.tgz
##
##

workflow TopMedVariantCaller {
  Boolean? calculate_DNA_contamination
  Boolean calculate_contamination = select_first([calculate_DNA_contamination, true])

  Int? SumCRAMs_CPUs
  Int SumCRAMs_CPUs_default = select_first([SumCRAMs_CPUs, 2])

  Int? CalcContamination_CPUs
  Int CalcContamination_CPUs_default = select_first([CalcContamination_CPUs, 2])

  Int? VariantCaller_CPUs
  Int VariantCaller_CPUs_default = select_first([VariantCaller_CPUs, 32])

  Array[File]? input_crai_files
  Array[File] input_cram_files

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

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  Int additional_disk = select_first([increase_disk_size, 20])

  Float reference_size = ( 
  size(ref_1000G_omni2_5_b38_sites_PASS_vcf_gz, "GB") +
  size(ref_1000G_omni2_5_b38_sites_PASS_vcf_gz_tbi, "GB") +
  size(chr10_vcf, "GB") +
  size(chr11_KI270927v1_alt_vcf, "GB") +
  size(chr11_vcf, "GB") +
  size(chr12_vcf, "GB") +
  size(chr13_vcf, "GB") +
  size(chr14_GL000009v2_random_vcf, "GB") +
  size(chr14_KI270846v1_alt_vcf, "GB") +
  size(chr14_vcf, "GB") +
  size(chr15_vcf, "GB") +
  size(chr16_vcf, "GB") +
  size(chr17_KI270857v1_alt_vcf, "GB") +
  size(chr17_KI270862v1_alt_vcf, "GB") +
  size(chr17_KI270909v1_alt_vcf, "GB") +
  size(chr17_vcf, "GB") +
  size(chr18_vcf, "GB") +
  size(chr19_KI270938v1_alt_vcf, "GB") +
  size(chr19_vcf, "GB") +
  size(chr1_KI270706v1_random_vcf, "GB") +
  size(chr1_KI270766v1_alt_vcf, "GB") +
  size(chr1_vcf, "GB") +
  size(chr20_vcf, "GB") +
  size(chr21_vcf, "GB") +
  size(chr22_KI270879v1_alt_vcf, "GB") +
  size(chr22_KI270928v1_alt_vcf, "GB") +
  size(chr22_vcf, "GB") +
  size(chr2_KI270773v1_alt_vcf, "GB") +
  size(chr2_KI270894v1_alt_vcf, "GB") +
  size(chr2_vcf, "GB") +
  size(chr3_vcf, "GB") +
  size(chr4_GL000008v2_random_vcf, "GB") +
  size(chr4_vcf, "GB") +
  size(chr5_vcf, "GB") +
  size(chr6_vcf, "GB") +
  size(chr7_KI270803v1_alt_vcf, "GB") +
  size(chr7_vcf, "GB") +
  size(chr8_KI270821v1_alt_vcf, "GB") +
  size(chr8_vcf, "GB") +
  size(chr9_vcf, "GB") +
  size(chrUn_KI270742v1_vcf, "GB") +
  size(chrX_vcf, "GB") +
  size(ref_dbsnp_142_b38_vcf_gz, "GB") +
  size(ref_dbsnp_142_b38_vcf_gz_tbi, "GB") +
  size(ref_dbsnp_All_vcf_gz, "GB") +
  size(ref_dbsnp_All_vcf_gz_tbi, "GB") +
  size(ref_hapmap_3_3_b38_sites_vcf_gz, "GB") +
  size(ref_hapmap_3_3_b38_sites_vcf_gz_tbi, "GB") +
  size(ref_hs38DH_bs_umfa, "GB") +
  size(ref_hs38DH_dict, "GB") +
  size(ref_hs38DH_fa, "GB") +
  size(ref_hs38DH_fa_alt, "GB") +
  size(ref_hs38DH_fa_amb, "GB") +
  size(ref_hs38DH_fa_ann, "GB") +
  size(ref_hs38DH_fa_bwt, "GB") +
  size(ref_hs38DH_fa_fai, "GB") +
  size(ref_hs38DH_fa_pac, "GB") +
  size(ref_hs38DH_fa_sa, "GB") +
  size(ref_hs38DH_winsize100_gc, "GB")
  )

  call sumCRAMSizes {
    input:
      input_crams = input_cram_files,
      input_crais = input_crai_files,
      disk_size = reference_size + additional_disk,
      SumCRAMs_CPUs_default = SumCRAMs_CPUs_default,
      docker_image = docker_image
  }

  if (calculate_contamination) {

#      Array[Array[File]] optional_contamination_scatter_output_files
      if (defined(input_crai_files)) {
          # Create an empty array to use to get around zip's requirement that
          # the input arrays not be optional
          Array[File] no_crai_files = []
          Array[File] crai_files = select_first([input_crai_files, no_crai_files])
          Array[Pair[File, File]] cram_and_crai_files = zip(input_cram_files, crai_files)
  
          scatter(cram_or_crai_file in cram_and_crai_files) {
              call getDNAContamination.calulateDNAContamination as scatter_getContamination {
                input:
                    input_cram_file = cram_or_crai_file.left,
                    input_crai_file = cram_or_crai_file.right,
        
                    ref_fasta = ref_hs38DH_fa,
                    ref_fasta_index = ref_hs38DH_fa_fai,
              
                    CalcContamination_CPUs = CalcContamination_CPUs_default 
              }
          }
#          Array[Array[File]] optional_contamination_scatter_output_files = scatter_getContamination.calculate_DNA_contamination_output
      } 

      if (!defined(input_crai_files)) {
          scatter(cram_file in input_cram_files) {
              call getDNAContamination.calulateDNAContamination as scatter_getContamination_no_crai {
                input:
                    input_cram_file = cram_file,
        
                    ref_fasta = ref_hs38DH_fa,
                    ref_fasta_index = ref_hs38DH_fa_fai,
              
                    CalcContamination_CPUs = CalcContamination_CPUs_default 
              }
          }
#          Array[Array[File]] optional_contamination_no_crai_scatter_output_files = scatter_getContamination_no_crai.calculate_DNA_contamination_output
      }

     Array[Array[File]] optional_contamination_scatter_output_files = select_first([scatter_getContamination.calculate_DNA_contamination_output, scatter_getContamination_no_crai.calculate_DNA_contamination_output])
#     if (defined(input_crai_files)) {   
#         optional_contamination_scatter_output_files = scatter_getContamination.calculate_DNA_contamination_output
#     }
#     if (!defined(input_crai_files)) {   
#         optional_contamination_scatter_output_files = scatter_getContamination_no_crai.calculate_DNA_contamination_output
#     }
     Array[File] contamination_output_files = flatten(optional_contamination_scatter_output_files)     
  }

  Array[File]? optional_contamination_output_files = contamination_output_files

  call variantCalling {

     input:
      contamination_output_files = optional_contamination_output_files,

      input_crais = input_crai_files,
      input_crams = input_cram_files,
      disk_size = sumCRAMSizes.total_size + reference_size + additional_disk,
      VariantCaller_CPUs_default = VariantCaller_CPUs_default,
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
    Array[File]? input_crais
    Int SumCRAMs_CPUs_default
    Float disk_size
    String docker_image

    # Create empty array to get around requirement for non optional array in length function
    Array[File] no_input_crais = []
    Array[File] optional_input_crais = select_first([input_crais, no_input_crais])
#    Int crai_array_len = if (defined(input_crais)) then length(optional_input_crais) else 0
    Int crai_array_len = length(optional_input_crais)

    command <<<
      python <<CODE
      import os
      import functools

      cram_string = "${ sep=',' input_crams }"
      cram_list = cram_string.split(',')

      crams_size = 0
      for cram_file in cram_list:
          crams_size = crams_size + os.stat(cram_file).st_size


      crais_size = 0
      if ${crai_array_len} > 0:
          crai_string = "${ sep=',' input_crais }"
          crai_list = crai_string.split(',')

          for crai_file in crai_list:
              crais_size = crais_size + os.stat(crai_file).st_size

      total_size = crams_size + crais_size
      # Shift right by 30 bits to get Gigabyte size of files
      total_size = (total_size >> 30)
      # Bump the size up 1 GB in case the total size is less than 1 GB
      print total_size + 1
      CODE
    >>>
    runtime {
      memory: "10 GB"
      cpu: sub(SumCRAMs_CPUs_default, "\\..*", "")
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      zones: "us-central1-a us-central1-b us-east1-d us-central1-c us-central1-f us-east1-c"
      docker: docker_image
    }
    output {
      Float total_size = read_float(stdout())
    }
  }

  task variantCalling {
     String? chromosomes
     String chromosomes_to_process = select_first([chromosomes, "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X" ])

     Int? num_of_jobs
     Int num_of_jobs_to_run = select_first([num_of_jobs, 23 ])

     Int? discoverUnit
     Int? genotypeUnit 
     File? PED_file

     Array[File]? contamination_output_files

     # The CRAM index files are listed as an input because they are required
     # by various tools, e.g. Samtools. They should be in the same location
     # as the CRAM files when specified in the input JSON
     Array[File]? input_crais
     Array[File] input_crams

     Float disk_size
     Int VariantCaller_CPUs_default
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

     String indexFileName = "trio_data.index"

     # We have to use a trick to make Cromwell
     # skip substitution when using the bash ${<variable} syntax
     # This is necessary to get the <var>=$(<command>) sub shell 
     # syntax to work and assign the value to a variable when 
     # running in Cromwell
     # See https://gatkforums.broadinstitute.org/wdl/discussion/comment/44570#Comment_44570 
     String dollar = "$"

     command <<<
      python3.5 <<CODE

      import csv
      import os
      from shutil import copy 

      # Erase the existing PED file; if no PED file is provided as input
      # this will sidestep the pedigree operations
      open('/root/topmed_freeze3_calling/data/trio_data.ped', 'w+').close()

      # If there is a PED file input copy the contents to the PED file
      # in the location where the program expects it to be 
      if len("${PED_file}") > 0:
         copy("${PED_file}", "/root/topmed_freeze3_calling/data/trio_data.ped")

      # Convert the WDL array of strings to a python list
      # The resulting string will be empty if the contamination values
      # are not calculated
      contamination_output_file_names_string = "${ sep=',' contamination_output_files }"
      # If DNA contaminiation was calculated for input files (CRAMs)
      if len(contamination_output_file_names_string) > 0:
          contamination_output_file_names_list = contamination_output_file_names_string.split(',')
          print("variantCalling: Contamination output files list is {}".format(contamination_output_file_names_list))

   
          # The user did not provide an index file so create it 
          # DNA contamination values should have already been calculated in
          # a previous task
          tsv_crams_rows = []
          # Create list of tuples from the list of cram paths and contamination
          # Input is [/path/to/cram, contamination, /path/to/cram, contamination...]
          # see https://stackoverflow.com/questions/23286254/convert-list-to-a-list-of-tuples-python
          file_pairs_it = iter(contamination_output_file_names_list)
          file_pairs_tuples = zip(file_pairs_it, file_pairs_it)
          for file_tuple in file_pairs_tuples:
              cram_file = file_tuple[0]
              contamination_file = file_tuple[1]
              print("variantCalling: CRAM file is {} contamination file is {}".format(cram_file, contamination_file))
          
              with open(contamination_file, 'rt') as in_file:
                 contents = in_file.read()
                 broj = contents.split("Alpha:",4)
                 contamination = broj[1].rstrip()
                 print("variantCalling: Contamination is {}".format(contamination))
    
              # Get the Cromwell location of the CRAM file
              # The worklow will be able to access them
              # since the Cromwell path is mounted in the
              # docker run commmand that Cromwell sets up
              base_name = os.path.basename(cram_file)
              base_name_wo_extension = base_name.split('.')[0]
     
              # Get the Cromwell path to the input CRAM file using the filename. The
              # filename at this time consists of the TopMed DNA sample
              # unique identifier of the form NWD123456.  
              tsv_crams_rows.append([base_name_wo_extension, cram_file, contamination])
    
      else:
          tsv_crams_rows = []
          # Convert the WDL array of strings to a python list
          input_crams_file_names_string = "${ sep=',' input_crams }"
          input_crams_file_names_list = input_crams_file_names_string.split(',')
          print("variantCalling: Input CRAM files names list is {}".format(input_crams_file_names_list))
          for cram_file in input_crams_file_names_list:
              # Get the Cromwell location of the CRAM file
              # The worklow will be able to access them
              # since the Cromwell path is mounted in the
              # docker run commmand that Cromwell sets up
              base_name = os.path.basename(cram_file)
              base_name_wo_extension = base_name.split('.')[0]
    
              # Get the Cromwell path to the input CRAM file using the filename. The
              # filename at this time consists of the TopMed DNA sample
              # unique identifier of the form NWD123456.  
              tsv_crams_rows.append([base_name_wo_extension, cram_file, "0.0"])

      print("variantCalling:  Writing index file {} with contents {}".format("${indexFileName}", tsv_crams_rows))
      with open("/root/topmed_freeze3_calling/data/${indexFileName}", 'w+') as tsv_index_file:
          writer = csv.writer(tsv_index_file, delimiter = '\t')
          for cram_info in tsv_crams_rows:
              writer.writerow(cram_info)

      # Print the index file to stdout for debugging purposes
      with open("/root/topmed_freeze3_calling/data/${indexFileName}", 'r') as tsv_index_file:
          print("variantCalling: Index file is:\n")
          print(tsv_index_file.read())
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

      # Check if the variable is set
      #https://unix.stackexchange.com/questions/212183/how-do-i-check-if-a-variable-exists-in-an-if-statement
      if [[ -n "${discoverUnit}" ]]; then
         printf "Setting discoverUnit to %s in gcconfig.pm\n" ${discoverUnit}
         sed -i '/.*our $discoverUnit.*/ c\our $discoverUnit = ${discoverUnit};' "$WORKING_DIR"/scripts/gcconfig.pm
      fi

      if [[ -n "${genotypeUnit}" ]]; then
         printf "Setting genotypeUnit to %s in gcconfig.pm\n" ${genotypeUnit}
         sed -i '/.*our $genotypeUnit.*/ c\our $genotypeUnit = ${genotypeUnit};' "$WORKING_DIR"/scripts/gcconfig.pm
      fi

      # Put the correct location of the output directory into the local config file
      #sed -i '/.*our $out =.*/ c\our $out = "/root/topmed_freeze3_calling/out";' "$WORKING_DIR"/scripts/gcconfig.pm
      # Put the correct location of references into the config file
      sed -i '/.*our $md5 =.*/ c\our $md5 = "\/data\/local.org\/ref\/gotcloud.ref\/md5\/%2s\/%s\/%s";' "$WORKING_DIR"/scripts/config.pm
      sed -i '/.*our $ref =.*/ c\our $ref = "\/data\/local.org\/ref\/gotcloud.ref\/hg38\/hs38DH.fa";' "$WORKING_DIR"/scripts/config.pm

      # Format the list of chromosomes to be e.g. "chr2 chr5 chrX"
      total=$(echo ${chromosomes_to_process} | wc -w)
      formatted_chromosomes_string=$(j=0; for i in ${chromosomes_to_process}; do printf "chr""$i"; let "j=j+1"; if [ "$j" -lt "$total" ]; then printf " "; fi done)

      echo "Running step1 - detect and merge variants"
      echo "Running step1 - detect and merge variants - removing old output dir if it exists"
      if [ -d "$WORKING_DIR"/out ]; then rm -Rf "$WORKING_DIR"/out; fi
      echo "Running step1 - detect and merge variants - generating Makefile"
      perl "$WORKING_DIR"/scripts/step1-detect-and-merge-variants.pl ${dollar}{formatted_chromosomes_string} 
      echo "Running step1 - detect and merge variants - running Makefile"
      make SHELL='/bin/bash' -f "$WORKING_DIR"/out/aux/Makefile -j ${num_of_jobs_to_run}
      

      echo "Running step2 - joint genotyping"
      echo "Running step2 - joint genotyping - removing old output dir if it exists"
      if [ -d "$WORKING_DIR"/paste ]; then rm -Rf "$WORKING_DIR"/paste; fi
      echo "Running step2 - joint genotyping - generating Makefile"
      perl "$WORKING_DIR"/scripts/step2-joint-genotyping.pl ${dollar}{formatted_chromosomes_string}
      echo "Running step2 - joint genotyping - running Makefile"
      # Format makefile name to be e.g. "chrchr2_chr15_chrX.Makefile"
      MAKEFILE_NAME="chr"$(j=0; for i in ${chromosomes_to_process}; do printf "chr""$i"; let "j=j+1"; if [ "$j" -lt "$total" ]; then printf "_"; fi done)".Makefile"
      make SHELL='/bin/bash' -f "$WORKING_DIR"/out/paste/"$MAKEFILE_NAME" -j ${num_of_jobs_to_run}

      if [[ -n "${PED_file}" ]]; then
         printf "variantCalling: Performing variant filtering using pedigree information\n"
         perl "$WORKING_DIR"/scripts/step3a-compute-milk-score.pl ${dollar}{formatted_chromosomes_string}
         make SHELL='/bin/bash' -f "$WORKING_DIR"/out/aux/milk/*.Makefile -j ${num_of_jobs_to_run}
         perl "$WORKING_DIR"/scripts/step3b-run-svm-milk-filter.pl ${dollar}{formatted_chromosomes_string}
      fi

      # Pop back to the original working directory; on FireCloud this will be a
      # special directory 
      popd

      if [[ -n "${PED_file}" ]]; then
          # Tar up the output directories into the output file provided in the input JSON
          # Pedigree information is in the svm directory
          tar -zcvf topmed_variant_caller_output.tar.gz /root/topmed_freeze3_calling/out/paste/ /root/topmed_freeze3_calling/out/aux/individual/ /root/topmed_freeze3_calling/out/svm/
      else
           # Tar up the output directories into the output file provided in the input JSON
          tar -zcvf topmed_variant_caller_output.tar.gz /root/topmed_freeze3_calling/out/paste/ /root/topmed_freeze3_calling/out/aux/individual/
      fi  
    >>>
     output {
      File topmed_variant_caller_output_file = "topmed_variant_caller_output.tar.gz"
    }
   runtime {
      memory: "10 GB"
      cpu: sub(VariantCaller_CPUs_default, "\\..*", "")
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      zones: "us-central1-a us-central1-b us-east1-d us-central1-c us-central1-f us-east1-c"
      docker: docker_image
    }
  }




## This is the U of Michigan variant caller workflow WDL for the workflow code located here:
## https://github.com/statgen/topmed_variant_calling
##
## It uses a Docker image built with software tools that can reproduce
## variant calls compatible to TopMed Freeze 8
##
## NOTE: This workflow assumes that input CRAM files have been built with the b38
## human reference genome. In particular for the TopMed CRAM files the
## reference genome files to use are located here:
## ftp://share.sph.umich.edu/1000genomes/fullProject/hg38_resources
## wget ftp://share.sph.umich.edu/1000genomes/fullProject/hg38_resources/topmed_variant_calling_example_resources.tar.gz
##
##

workflow TopMedVariantCaller {
  Boolean? dynamically_calculate_file_size
  Boolean dynamically_calculate_disk_requirement = select_first([dynamically_calculate_file_size, true])

  Float? All_CRAMs_disk_size_override
  Float All_CRAMs_disk_size_override_default = select_first([All_CRAMs_disk_size_override, 1000.0])

  Float? All_CRAIs_disk_size_override
  Float All_CRAIs_disk_size_override_default = select_first([All_CRAIs_disk_size_override, 100.0])

  Float? CRAM_file_max_disk_size_override
  Float CRAM_file_max_disk_size_override_default = select_first([CRAM_file_max_disk_size_override, 200.0])

  Float? ReferenceGenome_disk_size_override
  Float ReferenceGenome_disk_size_override_default = select_first([ReferenceGenome_disk_size_override, 30.0])

  Int? SumFileSizes_preemptible_tries
  Int  SumFileSizes_preemptible_tries_default = select_first([SumFileSizes_preemptible_tries, 3])
  Int? SumFileSizes_maxretries_tries
  Int SumFileSizes_maxretries_tries_default = select_first([SumFileSizes_maxretries_tries, 3])
  Int? SumFileSizes_memory
  Int SumFileSizes_memory_default = select_first([SumFileSizes_memory, 7])
  Int? SumFileSizes_CPUs
  Int SumFileSizes_CPUs_default = select_first([SumFileSizes_CPUs, 1])
  Int? SumFileSizes_disk_size
  Int SumFileSizes_disk_size_default = select_first([SumFileSizes_disk_size, 1])

  Int? CreateCRAMIndex_preemptible_tries
  Int CreateCRAMIndex_preemptible_tries_default = select_first([CreateCRAMIndex_preemptible_tries, 3])
  Int? CreateCRAMIndex_maxretries_tries
  Int CreateCRAMIndex_maxretries_tries_default = select_first([CreateCRAMIndex_maxretries_tries, 3])
  Int? CreateCRAMIndex_memory
  Int CreateCRAMIndex_memory_default = select_first([CreateCRAMIndex_memory, 7])
  Int? CreateCRAMIndex_CPUs
  Int CreateCRAMIndex_CPUs_default = select_first([CreateCRAMIndex_CPUs, 1])

  # The variant caller typically takes more than 24 hours to run. GCP terminates
  #  preemptible tasks after 24 hours. So by using 0 for preemptible tries the 
  #  task is non preemtible
  #  if preemptible is set to 0 -- then its set to false
  #  if preemptible is set to a positive integer -- its automatically true
  Int? VariantCaller_preemptible_tries
  Int VariantCaller_preemptible_tries_default = select_first([VariantCaller_preemptible_tries, 0])
  #if preemptible is 0 and maxRetries is 3 -- then that task can be retried upto 3 times
  #if preemptible is 3 and maxRetries is 3 for a task -- that can be retried upto 6 times
  #https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#maxretries
  Int? VariantCaller_maxretries_tries
  Int VariantCaller_maxretries_tries_default = select_first([VariantCaller_maxretries_tries, 3])
  Int? VariantCaller_memory
  # Select memory and CPUs to choose a GCP n1-highmem-64 machine
  Int VariantCaller_memory_default = select_first([VariantCaller_memory, 400])
  Int? VariantCaller_CPUs
  Int VariantCaller_CPUs_default = select_first([VariantCaller_CPUs, 64])
  # For adding more disk space for the variant caller from an input file
  Int? VariantCaller_additional_disk
  Int VariantCaller_additional_disk_default = select_first([VariantCaller_additional_disk, 1])

  String? VariantCallerHomePath
  String variantCallerHomePath_default = select_first([VariantCallerHomePath, "/topmed_variant_calling" ])


  Int? ExpandRefBlob_preemptible_tries
  Int ExpandRefBlob_preemptible_tries_default = select_first([ExpandRefBlob_preemptible_tries, 3])
  Int? ExpandRefBlob_maxretries_tries
  Int ExpandRefBlob_maxretries_tries_default = select_first([ExpandRefBlob_maxretries_tries, 3])
  Int? ExpandRefBlob_memory
  Int ExpandRefBlob_memory_default = select_first([ExpandRefBlob_memory, 100 ])
  Int? ExpandRefBlob_CPUs
  Int ExpandRefBlob_CPUs_default = select_first([ExpandRefBlob_CPUs, 1])
  String? ExpandRefBlob_expected_path
  String ExpandRefBlob_expected_path_default = select_first([ExpandRefBlob_expected_path, "resources/ref"])
  String ExpandRefBlob_glob_path = ExpandRefBlob_expected_path_default + "/*"

  Array[File]? input_crai_files
  Array[File] input_cram_files
  Array[String] input_cram_files_names = input_cram_files

  String? docker_image
  String  docker_image_default = select_first([docker_image, "statgen/topmed-variant-calling:v8.0.2"])

  File? referenceFilesBlob

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size
  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  Int additional_disk = select_first([increase_disk_size, 20])

  # The number of threads to use in a particular step of the pipeline
  Int? num_of_jobs
  Int num_of_jobs_to_run = select_first([num_of_jobs, 4 ])

  Float reference_size = if (dynamically_calculate_disk_requirement) then size(referenceFilesBlob, "GB") * 3
  else ReferenceGenome_disk_size_override_default

  call expandReferenceFileBlob {
     input:
      referenceFileBlob = referenceFilesBlob,
      referenceFilesExpectedPath = ExpandRefBlob_expected_path_default,
      referenceFilesGlobPath = ExpandRefBlob_glob_path,
      preemptible_tries = ExpandRefBlob_preemptible_tries_default,
      max_retries = ExpandRefBlob_maxretries_tries_default,
      docker_image = docker_image_default,
      CPUs = ExpandRefBlob_CPUs_default,
      disk_size = reference_size,
      memory = ExpandRefBlob_memory_default
  }


  if (dynamically_calculate_disk_requirement) {
      # Use scatter to get the size of each CRAM file:
      # Add 1 GB to size in case size is less than 1 GB
      # Use an array of String instead of File so Cromwell doesn't try to download them
      scatter(cram_file in input_cram_files_names ) { Float cram_file_size = round(size(cram_file, "GB")) + 1 }
      # Gather the sizes of the CRAM files:
      Array[Float] cram_file_sizes = cram_file_size
      # Use a task to sum the array:
      call sum_file_sizes as sum_cram_file_sizes {
        input:
          file_sizes = cram_file_sizes,
          preemptible_tries = SumFileSizes_preemptible_tries_default,
          max_retries = SumFileSizes_maxretries_tries_default,
          CPUs = SumFileSizes_CPUs_default,
          disk_size = SumFileSizes_disk_size_default,
          memory = SumFileSizes_memory_default
      }
  }

  #Float cram_files_size = if (dynamically_calculate_disk_requirement) then sum_cram_file_sizes.total_size else All_CRAMs_disk_size_override_default
  Float cram_files_size = select_first([sum_cram_file_sizes.total_size, All_CRAMs_disk_size_override_default])

  # If no CRAM index files were input then
  # create the index files
  if (!defined(input_crai_files)) {
      scatter(cram_file in input_cram_files) {
          Float cram_size = if (dynamically_calculate_disk_requirement) then  size(cram_file, "GB") else  CRAM_file_max_disk_size_override_default
          call createCRAMIndex as scatter_createCRAMIndex {
            input:
              input_cram = cram_file,
              variantCallerHomePath = variantCallerHomePath_default,
              disk_size = cram_size + additional_disk,
              docker_image = docker_image_default,
              CPUs = CreateCRAMIndex_CPUs_default,
              preemptible_tries = CreateCRAMIndex_preemptible_tries_default,
              memory = CreateCRAMIndex_memory_default,
              max_retries = CreateCRAMIndex_maxretries_tries_default
           }
      }
      Array[File] generated_crai_files = scatter_createCRAMIndex.output_crai_file
  }

  # if the CRAM index files were input then capture them otherwise they must have
  # been created so save those
  Array[File]? crai_files = if (defined(input_crai_files)) then input_crai_files else generated_crai_files

  # If there is an array of input CRAM index files then select those
  # otherwise they were generated so save those as an array
  Array[File] crai_files_array = select_first([input_crai_files, generated_crai_files])
  # Get the path and names of the CRAM index files for use when getting the
  # sizes of the files so Cromwell doesn't try to download them
  Array[String] crai_files_names_array = crai_files_array


  if (dynamically_calculate_disk_requirement) {
      # Use scatter to get the size of each CRAI file:
      # Add 1 GB to size in case size is less than 1 GB
      scatter(crai_file in crai_files_names_array ) { Float crai_file_size = round(size(crai_file, "GB")) + 1 }
      # Gather the sizes of the CRAI files:
      Array[Float] crai_file_sizes_array = crai_file_size
      # Use a task to sum the array:
      call sum_file_sizes as sum_crai_file_sizes {
        input:
          file_sizes = crai_file_sizes_array,
          preemptible_tries = SumFileSizes_preemptible_tries_default,
          CPUs = SumFileSizes_CPUs_default,
          max_retries = SumFileSizes_maxretries_tries_default,
          disk_size = SumFileSizes_disk_size_default,
          memory = SumFileSizes_memory_default
      }
  }

#  Float crai_files_size = if (dynamically_calculate_disk_requirement) then sum_crai_file_sizes.total_size else  All_CRAIs_disk_size_override_default
#  !!!!Why do I need to do this:
  Float crai_files_size = select_first([sum_crai_file_sizes.total_size, All_CRAIs_disk_size_override_default])

  Array[String] discoveryCommandsToRun = [
    "cd examples/",
    "../apigenome/bin/cloudify --cmd ../scripts/run-discovery-local.cmd",
    "make -f log/discover/example-discovery.mk -k -j ${num_of_jobs_to_run}"
    ]

  String discoveryGlobPath = "examples/out/sm/*/*"

  scatter(cram_file in input_cram_files) {
      Array[File] batchCRAMFiles = [cram_file]
      call variantCalling as scatter_runVariantCallingDiscovery {
         input:
          input_crais = crai_files,
          input_crams = batchCRAMFiles,
          input_cram_files_names = batchCRAMFiles,

          variantCallerHomePath = variantCallerHomePath_default,
          commandsToRun = discoveryCommandsToRun,
          globPath = discoveryGlobPath,

          disk_size = cram_files_size + crai_files_size + reference_size + additional_disk + VariantCaller_additional_disk_default,

          CPUs = VariantCaller_CPUs_default,
          preemptible_tries = VariantCaller_preemptible_tries_default,
          max_retries = VariantCaller_maxretries_tries_default,
          memory = VariantCaller_memory_default,
          docker_image = docker_image_default,

          referenceFiles = expandReferenceFileBlob.outputReferenceFiles,
          referenceFileExpectedPath = ExpandRefBlob_expected_path_default
      }
  }

  Array[File] individualCRAMVariants = flatten(scatter_runVariantCallingDiscovery.topmed_variant_caller_output_files)


  # We have to use a trick to make Cromwell
  # skip substitution when using the bash ${<variable} syntax
  # This is necessary to get the <var>=$(<command>) sub shell 
  # syntax to work and assign the value to a variable when 
  # running in Cromwell
  # See https://gatkforums.broadinstitute.org/wdl/discussion/comment/44570#Comment_44570 
  String dollar = "$"
  String doubledollar = "${dollar}{dollar}"

  Array[String] mergeCommandsToRun = [
    "cd examples/",
    "mkdir -p out/index",
    "../apigenome/bin/cram-vb-xy-index --index index/list.107.local.crams.index --dir out/sm/ --out out/index/list.107.local.crams.vb_xy.index",
    "../apigenome/bin/cloudify --cmd ../scripts/run-merge-sites-local.cmd",
    "make -f log/merge/example-merge.mk -k -j ${num_of_jobs_to_run}",
    "../apigenome/bin/cloudify --cmd ../scripts/run-batch-genotype-local.cmd",
    "make -f log/batch-geno/example-batch-genotype.mk -k -j ${num_of_jobs_to_run}",
    "../apigenome/bin/cloudify --cmd ../scripts/run-paste-genotype-local.cmd",
    "make -f log/paste-geno/example-paste-genotype.mk -k -j ${num_of_jobs_to_run}",
    'cut -f 1,4,5 index/intervals/b38.intervals.X.10Mb.1Mb.txt | grep -v ^chrX | awk \'{print "out/genotypes/hgdp/"${dollar}1"/merged."${dollar}1"_"${dollar}2"_"${dollar}3".gtonly.minDP0.hgdp.bcf"}\' > out/index/hgdp.auto.bcflist.txt',
    "../bcftools/bcftools concat -n -f out/index/hgdp.auto.bcflist.txt -Ob -o out/genotypes/hgdp/merged.autosomes.gtonly.minDP0.hgdp.bcf",
    "plink-1.9 --bcf out/genotypes/hgdp/merged.autosomes.gtonly.minDP0.hgdp.bcf --make-bed --out out/genotypes/hgdp/merged.autosomes.gtonly.minDP0.hgdp.plink --allow-extra-chr",
    "../king/king -b out/genotypes/hgdp/merged.autosomes.gtonly.minDP0.hgdp.plink.bed --degree 4 --kinship --prefix out/genotypes/hgdp/merged.autosomes.gtonly.minDP0.hgdp.king",
    "../apigenome/bin/vcf-infer-ped --kin0 out/genotypes/hgdp/merged.autosomes.gtonly.minDP0.hgdp.king.kin0 --sex out/genotypes/merged/chr1/merged.chr1_1_1000000.sex_map.txt --out out/genotypes/hgdp/merged.autosomes.gtonly.minDP0.hgdp.king.inferred.ped",
    "../apigenome/bin/cloudify --cmd ../scripts/run-milk-local.cmd",
    'cut -f 1,4,5 index/intervals/b38.intervals.X.10Mb.1Mb.txt | awk \'{print "out/milk/"${dollar}1"/milk."${dollar}1"_"${dollar}2"_"${dollar}3".sites.vcf.gz"}\' > out/index/milk.autoX.bcflist.txt',
    "(seq 1 22; echo X;) | xargs -I {} -P 10 bash -c \"grep chr{}_ out/index/milk.autoX.bcflist.txt | ../bcftools/bcftools concat -f /dev/stdin -Oz -o out/milk/milk.chr{}.sites.vcf.gz\"",
    "(seq 1 22; echo X;) | xargs -I {} -P 10 ../htslib/tabix -f -pvcf out/milk/milk.chr{}.sites.vcf.gz",
    "mkdir out/svm",
    "../apigenome/bin/vcf-svm-milk-filter --in-vcf out/milk/milk.chr2.sites.vcf.gz --out out/svm/milk_svm.chr2 --ref resources/ref/hs38DH.fa --dbsnp resources/ref/dbsnp_142.b38.vcf.gz --posvcf resources/ref/hapmap_3.3.b38.sites.vcf.gz --posvcf resources/ref/1000G_omni2.5.b38.sites.PASS.vcf.gz --train --centromere resources/ref/hg38.centromere.bed.gz --bgzip ../htslib/bgzip --tabix ../htslib/tabix --invNorm ${doubledollar}GC/bin/invNorm --svm-train ${doubledollar}GC/bin/svm-train --svm-predict ${doubledollar}GC/bin/svm-predict",
    "(seq 1 22; echo X;) | grep -v -w 2 | xargs -I {} -P 10 ../apigenome/bin/vcf-svm-milk-filter --in-vcf out/milk/milk.chr{}.sites.vcf.gz --out out/svm/milk_svm.chr{} --ref resources/ref/hs38DH.fa --dbsnp resources/ref/dbsnp_142.b38.vcf.gz --posvcf resources/ref/hapmap_3.3.b38.sites.vcf.gz --posvcf resources/ref/1000G_omni2.5.b38.sites.PASS.vcf.gz --model out/svm/milk_svm.chr2.svm.model --centromere resources/ref/hg38.centromere.bed.gz --bgzip ../htslib/bgzip --tabix ../htslib/tabix --invNorm ${doubledollar}GC/bin/invNorm --svm-train ../libsvm/svm-train --svm-predict ../libsvm/svm-predict"
    ]

  String mergeGlobPath = "examples/out/*/* examples/out/*/*/* examples/out/*/*/*/*"

  call  variantCalling as runVariantCallingMerge {
      input:
          input_cram_files_names = input_cram_files_names,
          cramVariants = individualCRAMVariants,

          variantCallerHomePath = variantCallerHomePath_default,
          commandsToRun = mergeCommandsToRun,
          globPath = mergeGlobPath,

          disk_size = cram_files_size + crai_files_size + reference_size + additional_disk + VariantCaller_additional_disk_default,

          CPUs = VariantCaller_CPUs_default,
          preemptible_tries = VariantCaller_preemptible_tries_default,
          max_retries = VariantCaller_maxretries_tries_default,
          memory = VariantCaller_memory_default,
          docker_image = docker_image_default,

          referenceFiles = expandReferenceFileBlob.outputReferenceFiles,
          referenceFileExpectedPath = ExpandRefBlob_expected_path_default
  }

  output {
      Array[File] topmed_variant_caller_output = runVariantCallingMerge.topmed_variant_caller_output_files
  }
  meta {
      author : "Walt Shands"
      email : "jshands@ucsc.edu"
      description: "This is the workflow WDL for U of Michigan's [TOPMed Freeze 8 Variant Calling Pipeline](https://github.com/statgen/topmed_variant_calling)"
   }
}

  task expandReferenceFileBlob {
     File? referenceFileBlob
     String? referenceFilesExpectedPath
     String referenceFilesGlobPath

     Float memory
     Float disk_size
     Int CPUs
     Int preemptible_tries
     Int max_retries
     String docker_image

     # We have to use a trick to make Cromwell
     # skip substitution when using the bash ${<variable} syntax
     # See https://gatkforums.broadinstitute.org/wdl/discussion/comment/44570#Comment_44570
     String dollar = "$"

     command <<<
      # Set the exit code of a pipeline to that of the rightmost command
      # to exit with a non-zero status, or zero if all commands of the pipeline exit
      set -o pipefail
      # cause a bash script to exit immediately when a command fails
      set -e
      # cause the bash shell to treat unset variables as an error and exit immediately
      set -u
      # echo each line of the script to stdout so we can see what is happening
      set -o xtrace
      #to turn off echo do 'set +o xtrace'

      if [[ -n "${referenceFileBlob}" ]]; then
          echo "Expanding reference files blob"
          printf "Untarring reference blob ${referenceFileBlob}"
          tar xvzf ${referenceFileBlob}
      else
          mkdir -p resources
          ln -s ${referenceFilesExpectedPath} "resources/ref"
      fi

      printf "Globbing reference files at resources/ref/*"

      >>>
        output {
          Array[File] outputReferenceFiles = glob("resources/ref/*")
       }
      runtime {
         memory: sub(memory, "\\..*", "") + " GB"
         cpu: sub(CPUs, "\\..*", "")
         maxRetries: max_retries
         preemptible: preemptible_tries
         disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
         zones: "us-central1-a us-central1-b us-east1-d us-central1-c us-central1-f us-east1-c"
         docker: docker_image
       }
  }

  task createCRAMIndex {
     File input_cram
     String variantCallerHomePath

     Float memory
     Float disk_size
     Int CPUs
     Int preemptible_tries
     Int max_retries
     String docker_image

     String CRAM_basename = basename(input_cram)
     String output_crai_file_name = "${CRAM_basename}.crai"

     # We have to use a trick to make Cromwell
     # skip substitution when using the bash ${<variable} syntax
     # See https://gatkforums.broadinstitute.org/wdl/discussion/comment/44570#Comment_44570 
     String dollar = "$"

     command <<<
      # Set the exit code of a pipeline to that of the rightmost command
      # to exit with a non-zero status, or zero if all commands of the pipeline exit 
      set -o pipefail
      # cause a bash script to exit immediately when a command fails
      set -e
      # cause the bash shell to treat unset variables as an error and exit immediately
      set -u
      # echo each line of the script to stdout so we can see what is happening
      set -o xtrace
      #to turn off echo do 'set +o xtrace'

      echo "Running create CRAM index"

      printf "Creating index ${input_cram}.crai for ${input_cram}"
      ${variantCallerHomePath}/samtools/samtools index ${input_cram} ${output_crai_file_name}

      >>>
        output {
          File output_crai_file = "${output_crai_file_name}"
       }
      runtime {
         memory: sub(memory, "\\..*", "") + " GB"
         cpu: sub(CPUs, "\\..*", "")
         maxRetries: max_retries
         preemptible: preemptible_tries
         disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
         zones: "us-central1-a us-central1-b us-east1-d us-central1-c us-central1-f us-east1-c"
         docker: docker_image
       }
  }


  # Calculates sum of a list of floats
  task sum_file_sizes {
    Array[Float] file_sizes

    Float memory
    Float disk_size
    Int CPUs
    Int preemptible_tries
    #String docker_image
    Int max_retries

    command <<<
    python -c "print ${sep="+" file_sizes}"
    >>>
    output {
      Float total_size = read_float(stdout())
    }
    runtime {
      docker: "python:2.7"
      preemptible: preemptible_tries
      maxRetries: max_retries
      memory: sub(memory, "\\..*", "") + " GB"
      cpu: sub(CPUs, "\\..*", "")
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"

    }
  }


  task variantCalling {
     # The number of jobs to run is the number of cores to use
     # Typically we use n1-highmem-64 but with 32 processes (ie, -j 32)
     # These are hyperthreaded cores, so we hope to get a slight performance boost by over-allocating cpus

     # The CRAM index files are listed as an input because they are required
     # by various tools, e.g. Samtools. They should be in the same location
     # as the CRAM files when specified in the input JSON
     Array[File]? input_crais
     Array[File]? input_crams
     Array[String] input_cram_files_names
     Array[File]? cramVariants

     Array[String] commandsToRun
     String globPath

     Float memory
     Float disk_size
     Int CPUs
     Int preemptible_tries
     String docker_image
     Int max_retries

     Array[File] referenceFiles
     String referenceFileExpectedPath
     String variantCallerHomePath

     String indexFileName = "examples/index/list.107.local.crams.index"

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
      import sys
      import pathlib

      cwd = os.getcwd()
      # Create the path to where the pipeline code will be executed
      pathlib.Path("examples/index").mkdir(parents=True, exist_ok=True)
      #os.symlink("examples", "${variantCallerHomePath}/examples")

      # Erase the existing sequence batch file
      #open('/topmed_variant_calling/examples/index/seq.batches.by.20.txt', 'w+').close()
      with open("examples/index/seq.batches.by.20.txt", 'w') as the_file:
          the_file.write('1')

      # Convert the WDL array of strings to a python list
      # The resulting string will be empty if the contamination values
      # are not calculated

      tsv_crams_rows = []
      # Convert the WDL array of strings to a python list
      input_crams_file_names_string = "${ sep=' ' input_cram_files_names }"
      input_crams_file_names_list = input_crams_file_names_string.split()
      print("variantCalling: Input CRAM files names list is {}".format(input_crams_file_names_list))
      for cram_file in input_crams_file_names_list:
          # Get the Cromwell basename  of the CRAM file
          # The worklow will be able to access them
          # since the Cromwell path is mounted in the
          # docker run commmand that Cromwell sets up
          base_name = os.path.basename(cram_file)
          base_name_wo_extension = base_name.split('.')[0]

          # The ID must be unique; and this depends on the input CRAM file names
          # being unique. Test to make sure the IDs are unique and fail the
          # workflow if they are not
          if(any(tsv_entry[0] == base_name_wo_extension for tsv_entry in tsv_crams_rows)):
              error_string = "variantCalling: ERROR: Duplicate ID {}. Input CRAM file names are probably not unique".format(base_name_wo_extension)
              print(error_string)
              sys.exit(error_string)

          # Use the basename of the CRAM file without suffix as an ID
          # The filename at this time consists of the TopMed DNA sample
          # unique identifier of the form NWD123456 followed by a suffix like .realigned.cram
          tsv_crams_rows.append([base_name_wo_extension, cwd + "/examples/crams/" + base_name])

      # Symlink the CRAM index files to the Cromwell working dir so the variant
      # caller can find them
      # Create the path to where the symlinks to the reference files will be
      pathlib.Path("examples/crams").mkdir(parents=True, exist_ok=True)

      input_crais_file_names_string = "${ sep=' ' input_crais }"
      input_crais_file_names_list = input_crais_file_names_string.split()
      print("variantCalling: Input CRAM index files names list is {}".format(input_crais_file_names_list))
      for crai_file in input_crais_file_names_list:
            crai_file_basename = os.path.basename(crai_file)
            print("variantCalling: Creating symlink {} for CRAM index file {}".format(cwd + "/examples/crams/" + crai_file_basename, crai_file))
            os.symlink(crai_file, cwd + "/examples/crams/" + crai_file_basename)

      # Symlink the CRAM files to the Cromwell working dir so the variant
      # caller can find them
      input_crams_file_names_string = "${ sep=' ' input_crams }"
      input_crams_file_names_list = input_crams_file_names_string.split()
      print("variantCalling: Input CRAM files names list is {}".format(input_crams_file_names_list))
      for cram_file in input_crams_file_names_list:
            cram_file_basename = os.path.basename(cram_file)
            print("variantCalling: Creating symlink {} for CRAM file {}".format(cwd + "/examples/crams/" + cram_file_basename, cram_file))
            os.symlink(cram_file, cwd + "/examples/crams/" + cram_file_basename)

      print("variantCalling:  Writing index file {} with contents {}".format("${indexFileName}", tsv_crams_rows))
      with open("${indexFileName}", 'w+') as tsv_index_file:
          writer = csv.writer(tsv_index_file, delimiter = '\t')
          for cram_info in tsv_crams_rows:
              writer.writerow(cram_info)

      # Print the index file to stdout for debugging purposes
      with open("${indexFileName}", 'r') as tsv_index_file:
          print("variantCalling: Index file is:\n")
          print(tsv_index_file.read())


      # Create the path to where the symlinks to the reference files will be
      pathlib.Path("examples/resources/ref").mkdir(parents=True, exist_ok=True)
      # Symlink the reference files to the pipeline expected path so the variant
      # caller can find them
      referenceFiles_file_names_string = "${ sep=' ' referenceFiles }"
      referenceFiles_file_names_list = referenceFiles_file_names_string.split()
      print("variantCalling: reference files names list is {}".format(referenceFiles_file_names_list))
      for reference_file in referenceFiles_file_names_list:
            reference_file_basename = os.path.basename(reference_file)
            print("variantCalling: Creating symlink {} for reference file {}".format(reference_file_basename, reference_file))
            os.symlink(reference_file, cwd + "/examples/resources/ref/" + reference_file_basename)

      # Create the path to where the symlinks to the reference files will be
      pathlib.Path("examples/out/sm").mkdir(parents=True, exist_ok=True)
      # Symlink the CRAM variant files to the Cromwell working dir so the variant
      # caller can find them
      input_cram_variants_file_names_string = "${ sep=' ' cramVariants }"
      input_cram_variants_file_names_list = input_cram_variants_file_names_string.split()
      print("variantCalling: Input CRAM variants files names list is {}".format(input_cram_variants_file_names_list))
      for cram_variant_file in input_cram_variants_file_names_list:
            cram_variant_file_basename = os.path.basename(cram_variant_file)
            cram_variant_base_name_wo_extension = cram_variant_file_basename.split('.')[0]
            # Create the path to where the pipeline code will be executed
            pathlib.Path(cwd + "/examples/out/sm/" + cram_variant_base_name_wo_extension).mkdir(parents=True, exist_ok=True)
            print("variantCalling: Creating symlink {} for CRAM variant file {}".format(cwd + "/examples/out/sm/" + cram_variant_base_name_wo_extension + "/" + cram_variant_file_basename, cram_variant_file))
            os.symlink(cram_variant_file, cwd + "/examples/out/sm/" + cram_variant_base_name_wo_extension + "/" + cram_variant_file_basename)

      CODE

      # Set the exit code of a pipeline to that of the rightmost command
      # to exit with a non-zero status, or zero if all commands of the pipeline exit 
      set -o pipefail
      # cause a bash script to exit immediately when a command fails
      set -e
      # cause the bash shell to treat unset variables as an error and exit immediately
      set -u
      # echo each line of the script to stdout so we can see what is happening
      set -o xtrace
      #to turn off echo do 'set +o xtrace'

      ln -s ${variantCallerHomePath}/build-all.sh build-all.sh
      ln -s ${variantCallerHomePath}/invNorm ivNorm
      ln -s ${variantCallerHomePath}/libStatGen libStatGen
      ln -s ${variantCallerHomePath}/samtools samtools
      ln -s ${variantCallerHomePath}/vt-topmed vt-topmed
      ln -s ${variantCallerHomePath}/apigenome apigenome
      ln -s ${variantCallerHomePath}/bcftools bcftools
      ln -s ${variantCallerHomePath}/bamUtil bamUtil
      ln -s ${variantCallerHomePath}/cramore cramore
      ln -s ${variantCallerHomePath}/htslib htslib
      ln -s ${variantCallerHomePath}/king king
      ln -s ${variantCallerHomePath}/libsvm libsvm

      ln -s ${variantCallerHomePath}/examples/index/intervals examples/index/intervals
      ln -s ${variantCallerHomePath}/scripts scripts

      # concatenate the array of commands to run in to a bash command line
      #commandsToRunStr='
      ${ sep=' && ' commandsToRun }

      # Tar up the output directories into the output file provided in the input JSON
      #tar -zcvf topmed_variant_caller_output.tar.gz "$CROMWELL_WORKING_DIR"/out/

    >>>
     output {
      Array[File] topmed_variant_caller_output_files = glob("${globPath}")
    }
   runtime {
      preemptible: preemptible_tries
      maxRetries: max_retries
      memory: sub(memory, "\\..*", "") + " GB"
      cpu: sub(CPUs, "\\..*", "")
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      zones: "us-central1-a us-central1-b us-east1-d us-central1-c us-central1-f us-east1-c"
      docker: docker_image
    }
  }


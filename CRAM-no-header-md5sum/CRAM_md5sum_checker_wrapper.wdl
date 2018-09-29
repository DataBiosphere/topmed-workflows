import "https://raw.githubusercontent.com/DataBiosphere/topmed-workflows/1.30.0/CRAM-no-header-md5sum/md5sum/CRAM_md5sum.wdl" as f1
import "https://raw.githubusercontent.com/DataBiosphere/topmed-workflows/1.30.0/CRAM-no-header-md5sum/checker/CRAM_md5sum_checker.wdl" as f2

workflow CRAMMd5sumChecker {
  File inputCRAMFile
  File? inputCRAMIndexFile
  File referenceFile
  File referenceIndexFile
  File truthMd5sumFile

  call f1.CRAM_to_md5sum as CRAMToMd5sum {
      input: inputCRAMFile = inputCRAMFile,
                    inputCRAMIndexFile = inputCRAMIndexFile,
                    referenceFile = referenceFile,
                    referenceIndexFile = referenceIndexFile}


  call f2.CRAM_to_md5sum_checker  {
      input: inputMd5sumFile = CRAMToMd5sum.md5sum_value ,truthMd5sumFile = truthMd5sumFile
  }
}

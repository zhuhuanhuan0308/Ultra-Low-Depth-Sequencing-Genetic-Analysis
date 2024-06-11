#! R
args=commandArgs(T)
outputdir=args[1]
bamlist=args[2]
ref=args[3]
human_posfile=args[4]
human_K=as.numeric(args[5])
human_nGen=as.numeric(args[6])
nCores=as.numeric(args[7])
niterations=as.numeric(args[8])
regionStart=as.numeric(args[9])
regionEnd=as.numeric(args[10])
chr=as.character(args[11])
buffer=as.numeric(args[12])
human_reference_sample_file=as.character(args[13])
human_reference_legend_file=as.character(args[14])
human_reference_haplotype_file=as.character(args[15])
sampleName=as.character(args[16])

tempdir=outputdir
setwd(outputdir)
library("STITCH",lib.loc="./software/Rpackages-3.5.1/")
STITCH(
  bamlist = bamlist,
  reference = ref,
  outputdir = outputdir,
  method = "diploid",
  regenerateInput = TRUE,
  regionStart = regionStart,
  regionEnd = regionEnd,
  buffer = buffer,
  niterations = niterations,
  chr = chr,
  sampleNames_file = sampleName,
  inputBundleBlockSize = 100,
  reference_populations = c("CHB", "CHS", "CDX"),
  reference_haplotype_file = human_reference_haplotype_file,
  reference_sample_file = human_reference_sample_file,
  reference_legend_file = human_reference_legend_file,
  posfile = human_posfile, K = human_K, tempdir = tempdir, nCores = nCores, nGen = human_nGen)

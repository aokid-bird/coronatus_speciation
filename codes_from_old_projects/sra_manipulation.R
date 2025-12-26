#=========================#
#### SRA manipulation ####
#=========================#
# This script will prepare sequences downloaded from sequence read archive (SRA)
# This script has been modified to work in the conda environment through rstudio-server and create some scripts to be run
# in the cluster PC.

# The purpose of this script is to
# 0) download and prepare FASTQ files
# 1) trim and quality filter SRA reads
# 2) quality check
# 3) trim, quality filter and quality check for the long-read (PacBio)
# 4) map raw reads against reference genome

source("code/dataprep.R")

# Directories
datadir <- "../data/NGS/migseq_200906/"
refdir <- "/lfs/user/data/WGS/Laterallus_jamaicensis_coturniculus/JAKCOX01.fasta"
sradir <- "/lfs/user/data/sra/sra"
projdir <- getwd()

#====================#
#### 0) Data prep ####
#====================#
#===========================#
##### 0-1) downlaod sra #####

# Laterallus jamaicensis
system("~/sratoolkit.3.0.5-ubuntu64/bin/prefetch SRR18439750 --max-size 100g")
# Atlantisia rogersi
system("~/sratoolkit.3.0.5-ubuntu64/bin/prefetch SRR9853853 --max-size 100g")

#==================================#
##### 0-2) unpack sra as fastq #####
# check the available space
system("df -h .")
# conduct the following inline of the linux shell under "sradir"
system("~/sratoolkit.3.0.5-ubuntu64/bin/fasterq-dump SRR18439750 SRR9853853")

#===============================#
#### 1) Trimmomatic trimming ####
#===============================#
ref.df

mkdir("res/trimmomatic_sra")

# to run on the cluster PC
for(i in 1:nrow(ref.df)){
  logname <- ref.df[i,]
  # parameters
  resdir_trim.clust <- "/lfs/user/project_RailPhylogeography_revise/res/trimmomatic_sra"
  pe = logname[2]
  threads = 72
  phred = 33
  sw = "4:15"
  leading = ifelse(pe == "PE", 30, 20)
  trailing = ifelse(pe == "PE", 30, 20)
  minlen = logname[4]
  
  if(pe == "PE"){
    trimmomatic <- str_c("trimmomatic ", pe,
                         " -threads ", threads,
                         " -phred", phred,
                         " -trimlog ", resdir_trim.clust, "/", logname[1], "_log.txt",
                         " ", sradir, "/", logname[1], "_1.fastq",
                         " ", sradir, "/", logname[1], "_2.fastq",
                         " ", resdir_trim.clust, "/", logname[1], "_pair_R1.fastq.gz",
                         " ", resdir_trim.clust, "/", logname[1], "_unpair_R1.fastq.gz",
                         " ", resdir_trim.clust, "/", logname[1], "_pair_R2.fastq.gz",
                         " ", resdir_trim.clust, "/", logname[1], "_unpair_R2.fastq.gz",
                         " ILLUMINACLIP:$adapt/TruSeq3-PE-2.fa:2:30:10",
                         " SLIDINGWINDOW:", sw,
                         " LEADING:", leading,
                         " TRAILING:", trailing,
                         " MINLEN:", minlen)
  }else{
    trimmomatic <- str_c("trimmomatic -Xms200g ", pe,
                         " -threads ", threads,
                         " -phred", phred,
                         " -trimlog ", resdir_trim.clust, "/", logname[1], "_log.txt",
                         " ", sradir, "/", logname[1], ".fastq",
                         " ", resdir_trim.clust, "/", logname[1], ".fastq.gz",
                         " SLIDINGWINDOW:", sw,
                         " LEADING:", leading,
                         " TRAILING:", trailing,
                         " MINLEN:", minlen)
  }

  adapter <- "adapt=\"/home/user/.conda/envs/bioinfo/share/trimmomatic-0.39-2/adapters\""

  # scripts for analysis on the cluster
  cluster.scripter(filename = str_interp("trim_${logname[1]}"),
                  local.directory = "res/trimmomatic_sra",
                  argument = list(adapter, trimmomatic),
                  name = str_interp("tr${logname[1]}"),
                  memsz = 250, cpunum = 72,
                  cluster.directory = resdir_trim.clust
                  )
}

#=====================================================#
#### 2) Quality and data structure check by fastqc ####
#=====================================================#

mkdir("res/fastqc_sra")

for(i in 2:nrow(ref.df)){
  logname <- ref.df[i,]
  values <- str_c(logname[1], c("_pair_R1.fastq.gz", "_pair_R2.fastq.gz", "_unpair_R1.fastq.gz")) %>%
    str_c(collapse = " ")

  fastqcloop <- str_c(
    str_interp("for eachfastq in ${values}; do"),
    "/home/user/.conda/envs/bioinfo/bin/fastqc -o ./ -t 72 /lfs/user/project_RailPhylogeography_coturnicops/res/trimmomatic/${eachfastq}",
    "done",
    sep = "\n"
  )

  # scripts for analysis on the cluster
  cluster.scripter(filename = str_interp("fastqc_${logname[1]}"),
                  local.directory = "res/fastqc_sra",
                  argument = list(fastqcloop),
                  name = str_interp("qc${logname[1]}"),
                  memsz = 250, cpunum = 72,
                  cluster.directory = "/lfs/user/project_RailPhylogeography_coturnicops/res/fastqc")

}

#======================================================#
#### 3) Trimming & QC for PACBIO Long read sequence ####
#======================================================#
# for PacBio sequence, first check the raw reads before trimming

# fastqc
i=1
logname <- logname.df[i,]
arg <- str_interp("/home/user/.conda/envs/bioinfo/bin/fastqc -o ./ -t 72 /lfs/user/data/sra/${logname[1]}.fastq")
# scripts for analysis on the cluster
cluster.scripter(filename = str_interp("fastqc_${logname[1]}_raw"),
                local.directory = "res/fastqc_sra",
                argument = list(arg),
                name = str_interp("qc${logname[1]}"),
                memsz = 250, cpunum = 72,
                cluster.directory = "/lfs/user/project_RailPhylogeography_coturnicops/res/fastqc")

# trimmomatic
threads = 36
phred = 33
sw = "4:15"
leading = 20
trailing = 20

i=1
trimmomatic <- str_c("/home/user/.conda/envs/bioinfo/bin/trimmomatic -Xms200g ", logname[2],
                     " -threads ", threads,
                     " -phred", phred,
                     " -trimlog ", resdir_trim, logname[1], "_log.txt",
                     " ", sradir, "/", logname[1], ".fastq",
                     " ", resdir_trim, logname[1], ".fastq.gz",
                     #" ILLUMINACLIP:$adapt/TruSeq3-PE-2.fa:2:30:10",
                     " SLIDINGWINDOW:", sw,
                     " MINLEN:", logname[4],
                     " LEADING:", leading,
                     " TRAILING:", trailing)

adapter <- "adapt=\"/home/user/.conda/envs/bioinfo/share/trimmomatic-0.39-2/adapters\""

cluster.scripter(filename = str_interp("trim_${logname[1]}"),
                local.directory = "res/trimmomatic_sra",
                argument = list(adapter, trimmomatic),
                name = str_interp("tr${logname[1]}"),
                memsz = 250, cpunum = 72,
                cluster.directory = "/lfs/user/project_RailPhylogeography_coturnicops/res/trimmomatic")

i=1

fastqc.arg <- str_interp("/home/user/.conda/envs/bioinfo/bin/fastqc -o ./ -t 72 ${resdir_trim}${logname[1]}.fastq.gz")

cluster.scripter(filename = str_interp("trim_${logname[1]}"),
                local.directory = "res/fastqc_sra",
                argument = list(fastqc.arg),
                name = str_interp("qc${logname[1]}"),
                memsz = 250, cpunum = 72,
                cluster.directory = "/lfs/user/project_RailPhylogeography_coturnicops/res/fastqc")


#====================================================#
#### 4) Reference mapping and bam merge, sam sort ####
#====================================================#
##### 4-1) Bwa-mem2/bwa-mem for reference mapping #####
# paired seqs

for(i in 1:nrow(ref.df)){
  logname <- ref.df[i,]

  # parameters
  algorithm = "mem"
  thread = 72
  resdir_bwa <- "/lfs/user/project_RailPhylogeography_revise/res/bwa_sra"

  if(logname$type == "PE"){

    bwap <-
      str_c("bwa ", algorithm, " -t ", thread, " ", refdir,
            " /lfs/user/project_RailPhylogeography_revise/res/trimmomatic_sra/", logname[1], "_pair_R1.fastq.gz",
            " /lfs/user/project_RailPhylogeography_revise/res/trimmomatic_sra/", logname[1], "_pair_R2.fastq.gz",
            " | samtools sort -@ ", thread, " -o ", resdir_bwa, "/", logname[1], "_p.bam")

    bwa1 <-
      str_c("bwa ", algorithm, " -t ", thread, " ", refdir,
            " /lfs/user/project_RailPhylogeography_revise/res/trimmomatic_sra/", logname[1], "_unpair_R1.fastq.gz",
            " | samtools sort -@ ", thread, " -o ", resdir_bwa, "/", logname[1], "_R1.bam")

    bwa2 <-
      str_c("bwa ", algorithm, " -t ", thread, " ", refdir,
            " /lfs/user/project_RailPhylogeography_revise/res/trimmomatic_sra/", logname[1], "_unpair_R2.fastq.gz",
            " | samtools sort -@ ", thread, " -o ", resdir_bwa, "/", logname[1], "_R2.bam")

  }else{
    bwa1 <- bwa2 <- NULL
    bwap <-
      str_c("bwa ", algorithm, " -t ", 50, " ", refdir,
            " /lfs/user/project_RailPhylogeography_revise/res/trimmomatic_sra/", logname[1], ".fastq.gz",
            " | samtools sort -@ ", thread, " -o ", resdir_bwa, "/", logname[1], ".bam")

  }

  cluster.scripter(filename = str_interp("bwa_${logname[1]}"),
                   local.directory = "res/bwa_sra",
                   argument = list(bwap, bwa1, bwa2),
                   name = str_interp("bw${logname[1]}"),
                   memsz = 250, cpunum = 72,
                   cluster.directory = resdir_bwa)

}


##### 4-2) indexing #####
cd.arg <- "cd /lfs/user/data/WGS/Laterallus_jamaicensis_coturniculus/"
index.arg <- "bwa-mem2 index JAKCOX01.fasta"

cluster.scripter(filename = str_interp("bwa_index2"),
                 local.directory = "res/bwa_sra",
                 argument = list(cd.arg, index.arg),
                 name = str_interp("bwindex"),
                 memsz = 250, cpunum = 72,
                 cluster.directory = resdir_bwa)

##### 4-3) combine bams to one #####
for(i in 1:nrow(ref.df)){
  logname <- ref.df[i,]
  resdir_bwa.clust <- "/lfs/user/project_RailPhylogeography_revise/res/bwa_sra"

  if(logname$type == "PE"){
    out <- str_interp("${resdir_bwa.clust}/${logname[1]}.bam")
    pbam <- str_interp("${resdir_bwa.clust}/${logname[1]}_p.bam")
    unpbam1 <- str_interp("${resdir_bwa.clust}/${logname[1]}_R1.bam")
    unpbam2 <- str_interp("${resdir_bwa.clust}/${logname[1]}_R2.bam")

    run <- str_interp("samtools merge ${out} ${pbam} ${unpbam1} ${unpbam2}")
    #rename <- str_interp("rename -v s/_p.bam/.bam/ ${pbam}")
    sort <- str_interp("samtools sort -@ ${thread} -o ${resdir_bwa.clust}/${logname[1]}.bam ${resdir_bwa.clust}/${logname[1]}.bam")

  }else{
    run <- NULL
    #rename <- NULL
    sort <- str_interp("samtools sort -@ ${thread} -o ${resdir_bwa.clust}/${logname[1]}.bam ${resdir_bwa.clust}/${logname[1]}.bam")
  }

  cluster.scripter(filename = str_interp("bammerge_${logname[1]}"),
                   local.directory = "res/bwa_sra",
                   argument = list(run, sort),
                   name = str_interp("mg${logname[1]}"),
                   memsz = 250, cpunum = 72,
                   cluster.directory = resdir_bwa.clust)
}


#!/usr/bin/env Rscript

######
# LOAD PACKAGES
######

# Load R packages
packages <- c("openxlsx", "stringr", "dplyr", "colorspace", "argparse")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# Make argparse object and arguments
parser <- ArgumentParser(description="Annotate SNPs with functional information via AnnoVar.")
parser$add_argument("-i", "--input_file", type="character",
                    help="Path to single-column file containing unique SNP rsIDs")
parser$add_argument("-a", "--annovar_folder", type="character",
                    help="Path to AnnoVar folder")
parser$add_argument("-b", "--build", type="character", default="hg19",
                    help="Genome build to use [default %(default)s]")
parser$add_argument("-o", "--output_folder", type="character", default="output",
                    help="Path to output folder [default %(default)s]")
args <- parser$parse_args()

######
# PARAMETER SETTING
######

# Sets user parameters
input_file <- args$input_file
annovar_folder <- args$annovar_folder
build <- args$build
output_folder <- args$output_folder

## NOTE for testing
input_file <- "input/snps_unique_hc.txt"
annovar_folder <- "~/Software/annovar"
build <- "hg19"
output_folder <- "output"
###

# Check if required files and folders exist
if (!file.exists(input_file)) {
  stop("Input SNP file does not exist")
}
if (!file.exists(annovar_folder)) {
  stop("AnnoVar folder does not exist (download from web)")
}

# Define other important inputs
humandb_folder <- file.path(annovar_folder, "humandb")
CONVERT2ANNO <- file.path(annovar_folder, "convert2annovar.pl")
TAB_ANNO <- file.path(annovar_folder, "table_annovar.pl")
ANNO_VAR <- file.path(annovar_folder, "annotate_variation.pl")

# Generate parent output folder if it doesn't already exist
if (!dir.exists(output_folder)) { dir.create(output_folder) }

# Print user parameters to console
cat(sprintf("Hostname: %s\n", Sys.info()["nodename"]))
cat(sprintf("Start time: %s\n", format(Sys.time(), "%a %b %d %X %Y")))
cat(sprintf("Working directory: %s\nOutput directory: %s\n", getwd(), output_folder))
cat(sprintf("Input SNP rsID file: %s\n", input_file))

######
# DOWNLOAD ANNOTATION DATABASES
######

cat("Downloading AnnoVar annotation databases\n")

# refGene FASTA sequences for all annotated transcripts in RefSeq Gene
system(sprintf("perl %s -build %s -downdb -webfrom annovar refGene %s", ANNO_VAR, build, humandb_folder))

# Cytogenetic bands
system(sprintf("perl %s -build %s -downdb cytoBand %s", ANNO_VAR, build, humandb_folder))

# SNP138
system(sprintf("perl %s -build %s -downdb snp138 %s", ANNO_VAR, build, humandb_folder))

# gnomAD exome collection
system(sprintf("perl %s -build %s -downdb -webfrom annovar gnomad211_exome %s", ANNO_VAR, build, humandb_folder))

# PhastCons
system(sprintf("perl %s -build %s -downdb phastConsElements100way %s", ANNO_VAR, build, humandb_folder))

######
# PREPARE ANNOVAR INPUT
######

# Check if SNPs are in proper rsID format
if (!unique(grepl("^rs[0-9]", readLines(input_file)))) {
  stop("SNPs are not in proper rsID format")
}

# Re-formatted input file name
input_name <- paste0("output/", sub(".txt", "", basename(input_file)), ".avinput")

# Convert SNP list to AnnoVar input format using SNP138 reference
cat("Formatting SNPs for AnnoVar input\n")
str1 <- sprintf("perl %s -format rsid %s", CONVERT2ANNO, input_file)
str2 <- sprintf("-dbsnpfile %s --outfile hg19_snp138.txt", dbsnp_file, input_name)
cmd1 <- paste(str1, str2)
system(cmd1)

######
# ANNOTATE SNPS
######

# Annotate variants in avinput file using table_annovar.pl
str1 <- sprintf("perl %s %s %s", TAB_ANNO, input_name, humandb_folder)
str2 <- sprintf("-buildver %s -out %s/snp_anno -remove", build, output_folder)
str3 <- sprintf("-protocol refGeneWithVer,cytoBand,gnomad211_exome -operation g,r,f")
str4 <- sprintf("-nastring . -polish")
cmd2 <- paste(str1, str2, str3, str4)
system(cmd2)

# Run gene-based annotations (ie. classify them as intergenic, intronic,
# non-synonymous SNP, frameshift deletion, large-scale duplication, etc.)
cmd2 <- sprintf("perl %s -geneanno -buildver hg19 %s %s -outfile %s", ANNO_VAR, input_name, humandb_folder, input_name)
system(cmd2)

# Run region-based annotations (ie. annotate variants that fall within conserved regions)
str1 <- sprintf("perl %s -regionanno -build %s -outfile %s", ANNO_VAR, build, input_name)
str2 <- sprintf("-dbtype phastConsElements100way %s %s", input_name, humandb_folder)
cmd3 <- paste(str1, str2)
system(cmd3)

# Run filter-based annotations (ie. identify a subset of variants input file
# that are not observed in 1000G version 2014 Oct and those that are observed
# with allele frequencies)
#cat("Running filter-based annotations...")
#str1 <- sprintf("perl %s -filter", ANNO_VAR)
#str2 <- sprintf("-dbtype 1000g2014oct_all -buildver %s %s/snplist.avinput %s",
#                annoBuild, outDir, annodb)
#cmd <- sprintf("%s %s", str1, str2)
#system(cmd)
#cat(" done.\n")

######
# ANNOTATION SUMMARY
######

# Print summary of SNP annotations
dat <- read.delim(sprintf("%s/snp_anno.hg19_multianno.txt", outDir), h=TRUE, as.is=TRUE)

# Write output to excel format
write.xlsx(dat, file=sprintf("%s/snp_anno.hg19_multianno.xlsx", outDir), col=TRUE, row=FALSE)

 # Pie chart of SNP functional annotations (all annotated SNPs)
 annoTable <- table(dat$Func.refGene)
 annoTabledf <- as.data.frame(annoTable)
 pred_anno <- annoTabledf$Var1
 freq_anno <- annoTabledf$Freq
   freq_anno_pct <- round(freq_anno/sum(freq_anno) * 100)
   freq_anno_pct <- paste(freq_anno_pct, "%", sep="")
   freq_anno_pct <- sprintf("(%s)", freq_anno_pct)
 plot_anno_lbls <- paste(pred_anno, freq_anno_pct, sep=" ")

 png(file=sprintf("%s/%s_function_piechart.png",
     outDir, snpGroup[i]))
 pie(annoTable, labels=plot_anno_lbls,
     col=rainbow_hcl(length(annoTable)),
     main=paste0(sprintf("Proportion of %s variant function", snpGroup[i])))
 dev.off()

# Rearrange data to get columns in the following order: "avsnp144 (rsID)",
# "Func.refGene", "ExonicFunc.refGene", "SIFT_score", "SIFT_pred",
# "Polyphen2_HDIV_score", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score",
# "Polyphen2_HVAR_pred", etc..
dat2 <- dat[c(13,6,9,14:ncol(dat))]
cat(sprintf("Total number of annotated SNPs: %i\n", nrow(dat2)))
cat("Function of all annotated SNPs in refGene database:")
print(table(dat2[2]))

# Determining predictions
cat("Functional predictions per annotation database:\n")
# PolyPhen2 HDIV prediction
type_of_change <- data.frame()    # Initialize empty data frame
exonic <- dat2[which(dat2$Func.refGene == "exonic"),]

polyphenPred <- exonic$Polyphen2_HDIV_pred
for (i in 1:nrow(exonic)) {
  if(polyphenPred[i] == ".") {
    type_of_change[i, "Polyphen2_HDIV_pred"] <- "missing"
  }
  else if(polyphenPred[i] %in% "D") {
    type_of_change[i, "Polyphen2_HDIV_pred"] <- "probably_damaging"
  }
  else if(polyphenPred[i] %in% "P") {
    type_of_change[i, "Polyphen2_HDIV_pred"] <- "possibly_damaging"
  }
  else if(polyphenPred[i] %in% "B") {
    type_of_change[i, "Polyphen2_HDIV_pred"] <- "benign"
  }
  else {
    type_of_change[i, "Polyphen2_HDIV_pred"] <- "other"
  }
}

# PolyPhen2 HVAR prediction
polyphenPred2 <- exonic$Polyphen2_HVAR_pred
for (i in 1:nrow(exonic)) {
  if(polyphenPred2[i] == ".") {
    type_of_change[i, "Polyphen2_HVAR_pred"] <- "missing"
  }
  else if(polyphenPred2[i] %in% "D") {
    type_of_change[i, "Polyphen2_HVAR_pred"] <- "probably_damaging"
  }
  else if(polyphenPred2[i] %in% "P") {
    type_of_change[i, "Polyphen2_HVAR_pred"] <- "possibly_damaging"
  }
  else if(polyphenPred2[i] %in% "B") {
    type_of_change[i, "Polyphen2_HVAR_pred"] <- "benign"
  }
  else {
    type_of_change[i, "Polyphen2_HVAR_pred"] <- "other"
  }
}

# SIFT prediction
SIFTpred <- exonic$SIFT_pred
for (i in 1:nrow(exonic)) {
 if(SIFTpred[i] == ".") {
   type_of_change[i, "SIFT_pred"] <- "missing"
 }
 else if(SIFTpred[i] == "D") {
   type_of_change[i, "SIFT_pred"] <- "damaging"
 }
 else if(SIFTpred[i] == "T") {
   type_of_change[i, "SIFT_pred"] <- "tolerated"
 }
 else {
   type_of_change[i, "SIFT_pred"] <- "other"
 }
}

# LRT prediction
LRTpred <- exonic$LRT_pred
for (i in 1:nrow(exonic)) {
  if(LRTpred[i] == ".") {
    type_of_change[i, "LRT_pred"] <- "missing"
  }
  else if(LRTpred[i] == "D") {
    type_of_change[i, "LRT_pred"] <- "deleterious"
  }
  else if(LRTpred[i] == "N") {
    type_of_change[i, "LRT_pred"] <- "neutral"
  }
  else if(LRTpred[i] == "U") {
    type_of_change[i, "LRT_pred"] <- "unknown"
  }
}

 # MutationTaster prediction
 MTpred <- exonic$MutationTaster_pred
 for (i in 1:nrow(exonic)) {
   if(MTpred[i] == ".") {
     type_of_change[i, "MutationTaster_pred"] <- "missing"
   }
   else if(MTpred[i] == "A") {
     type_of_change[i, "MutationTaster_pred"] <- "disease_causing_automatic"
   }
   else if(MTpred[i] == "D") {
     type_of_change[i, "MutationTaster_pred"] <- "disease_causing"
   }
   else if(MTpred[i] == "N") {
     type_of_change[i, "MutationTaster_pred"] <- "polymorphism"
   }
   else if(MTpred[i] == "P") {
     type_of_change[i, "MutationTaster_pred"] <- "polymorphism_automatic"
   }
   else {
     type_of_change[i, "MutationTaster_pred"] <- "other"
   }
 }

 # MutationAssessor prediction
 MApred <- exonic$MutationAssessor_pred
 for (i in 1:nrow(exonic)) {
   if(MApred[i] == ".") {
     type_of_change[i, "MutationAssessor_pred"] <- "missing"
   }
   else if(MApred[i] == "H") {
     type_of_change[i, "MutationAssessor_pred"] <- "high_functional"
   }
   else if(MApred[i] == "M") {
     type_of_change[i, "MutationAssessor_pred"] <- "medium_functional"
   }
   else if(MApred[i] == "L") {
     type_of_change[i, "MutationAssessor_pred"] <- "low_non-functional"
   }
   else if(MApred[i] == "N") {
     type_of_change[i, "MutationAssessor_pred"] <- "neutral_non-functional"
   }
   else {
     type_of_change[i, "MutationAssessor_pred"] <- "other"
   }
 }

 # FATHMM prediction
 FATHMMpred <- exonic$FATHMM_pred
 for (i in 1:nrow(exonic)) {
   if(FATHMMpred[i] == ".") {
     type_of_change[i, "FATHMM_pred"] <- "missing"
    }
   else if(FATHMMpred[i] == "D") {
      type_of_change[i, "FATHMM_pred"] <- "damaging"
    }
   else if(FATHMMpred[i] == "T") {
     type_of_change[i, "FATHMM_pred"] <- "tolerated"
   }
   else {
     type_of_change[i, "FATHMM_pred"] <- "other"
   }
 }

 # GERP++_RS prediction
 GERPpred <- as.numeric(exonic$GERP)
 for (i in 1:nrow(exonic)) {
   if(GERPpred[i] %in% NA) {
     type_of_change[i, "GERP"] <- "missing"
    }
   else if(GERPpred[i] > 0) {
      type_of_change[i, "GERP"] <- "conserved"
     }
   else if(GERPpred[i] < 0) {
      type_of_change[i, "GERP"] <- "not_conserved"
     }
   else {
      type_of_change[i, "GERP"] <- "other"
     }
  }

# Print table of AnnoVar predictions for each exonic SNV to log
exonic_SNVfunc <- exonic[c(1,3)]
colnames(exonic_SNVfunc) <- c("Exonic_SNV", "Function")
exonic_predTable <- cbind(exonic_SNVfunc, type_of_change)
print(exonic_predTable)

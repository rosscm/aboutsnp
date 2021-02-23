#' Run annovar (taken from http://annovar.openbioinformatics.org/en/latest/user-guide/startup/)

# Directory containing annovar software
cd /media/catherine/DATAPART1/Software/annovar

# Directory containing SNP input files + location of annovar output files
dataDir=/media/catherine/DATAPART1/Data/PopulationPathways/results/annovar/HM3/top13pathways

# Download required databases (uncomment if needed)
# perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
# perl annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
# perl annotate_variation.pl -buildver hg19 -downdb genomicSuperDups humandb/
# perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all humandb/
# perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2014oct humandb/
# perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp138 humandb/
# perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ljb26_all humandb/

# Take dbSNP rs identifiers and annotate functionality of these SNPs
perl convert2annovar.pl -format rsid ${dataDir}/topSnpList.txt \
			-dbsnpfile humandb/hg19_snp138.txt > ${dataDir}/snplist.avinput

# Run the table_annovar.pl program to annotate the variants in the .avinput file
perl table_annovar.pl ${dataDir}/snplist.avinput humandb/ \
		 	-buildver hg19 \
			-out ${dataDir}/topSnp_anno \
			-remove \
			-protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all \
			-operation g,r,r,f,f,f,f,f,f,f \
			-nastring . \
			-csvout

# Run the following three commands that correspond to gene-based, region-based
# and filter-based annotations
perl annotate_variation.pl -geneanno \
			   -buildver hg19 ${dataDir}/snplist.avinput humandb/
perl annotate_variation.pl -regionanno \
			   -dbtype cytoBand \
			   -buildver hg19 ${dataDir}/snplist.avinput humandb/
perl annotate_variation.pl -filter \
			   -dbtype 1000g2014oct_all \
			   -buildver hg19 ${dataDir}/snplist.avinput humandb/

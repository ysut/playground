ANNOVAR="/usr/local/bio/src/annovar202203/table_annovar.pl"
INPUT_VCF="/path/to/input.vcf"
ASSEMBLY="hg19"
ANNOVAR_DB="/betelgeuse04/analysis/utsu/resources/annovar/hg19"

"${ANNOVAR}" "${INPUT_VCF}" \
    --buildver ${ASSEMBLY} \
	--out "${INPUT_VCF}".anv.tsv \
	# --protocol refGene,esp6500siv2_all,dbnsfp30a,exac03,snp138NonFlagged,vcf,vcf,vcf,vcf,snp20171005_tommo3.5k_passed,gerp++gt2,cosmic70,clinvar_20221231,bed,bed,bed,bed,dbscsnv11,spidex \
	--protocol refGene,clinvar20221231,bed,bed,revel,mcap,dbnsfp42a,gnomad211_genome,esp6500siv2_all,vcf,vcf \
	--operation g,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r,f,f \
	--vcfdbfile hg19_snpHGVD20170303.txt,hg19_ABRaOM_fix.txt \
	--vcfinput ${ANNOVAR_DB} \
	--bedfile hg19_hgmdAllMut.bed,hg19_hgmdAllMut-collapse.bed \
	-arg '--splicing_threshold $(SPLICING_THRESHOLD),,-colsWanted 4,-colsWanted 4,,,,-infoasscore,-infoasscore,,,,,,' \
	--remove \
	--nastring . \
	--otherinfo


<< COMMENTOUT
Gene-based annotation
1. refGene
2. ensGene41

Filter-based annotation
1. clinvar_20221231
2. dbnsfp42a
3. gnomad211_exome
4. gnomda211_genome
5. esp6500siv2_all
6. gme
7. revel
8. mcap

VCF-based annotation
1. hg19_snpHGVD20170303.txt		-infoasscore
2. hg19_ABRaOM_fix.txt			-infoasscore

BED-based annotation
1. hg19_hgmdAllMut				-colsWanted 4
2. hg19_hgmdAllMut-collapse		-colsWanted 4



COMMENTOUT
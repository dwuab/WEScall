# create g1k reference panel used for WEScall

set -euo pipefail

export g1k_path
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )"/../ >/dev/null 2>&1 && pwd )"
export DIR
echo $DIR

# first check existence of g1k reference panel

if [ -d $DIR/resources/1000G_ref_panel ]; then
	echo "It looks like you already have the reference in place."
	printf "In case you want to create the reference panel again, please delete %s first and try again.\n" $DIR/resources/1000G_ref_panel
	exit
fi

url_g1k=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
printf "WEScall needs 1000G phase 3 data as reference panel. If you don't have 1000G phase 3 data, you can download it from %s\n" $url_g1k
read -p "Please enter the path to 1000G phase 3 data: " g1k_path
while [ ! -d ${g1k_path} ]; do
	printf "directory %s does not exists!\n" ${g1k_path}
	read -p "Please enter the path to 1000G phase 3 data: " g1k_path
done

function filter_autosome {
	chr=$1
	bcftools annotate ${g1k_path}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
		-h <(printf "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">") -Ou | \
		bcftools norm --rm-dup none -Ou | \
		bcftools view -m 2 -M 2 -Ou | \
		bcftools view -i 'MAC>=5 && INFO/VT!="SV"' -Oz --threads 10 > $DIR/resources/1000G_ref_panel/ALL.chr${chr}.phase3.20130502.SNP.indel.biallelic.mac5.vcf.gz
}

function filter_X_chr {
	bcftools annotate ${g1k_path}/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1a.20130502.genotypes.vcf.gz \
		-h <(printf "%s\n%s" "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" "##INFO=<ID=OLD_VARIANT,Number=1,Type=String,Description=\"unknown\">") -Ou | \
		bcftools norm --rm-dup none -Ou | \
		bcftools view -m 2 -M 2 -Ou | \
		bcftools view -i 'MAC>=5 && INFO/VT!="SV"' -Oz --threads 10 > $DIR/resources/1000G_ref_panel/ALL.chrX.phase3.20130502.SNP.indel.biallelic.mac5.vcf.gz

	bcftools convert --haploid2diploid -h $DIR/resources/1000G_ref_panel/chrX_tmp $DIR/resources/1000G_ref_panel/ALL.chrX.phase3.20130502.SNP.indel.biallelic.mac5.vcf.gz
	bcftools convert -H $DIR/resources/1000G_ref_panel/chrX_tmp -Oz > $DIR/resources/1000G_ref_panel/ALL.chrX.phase3.20130502.SNP.indel.biallelic.mac5.vcf.gz
	rm $DIR/resources/1000G_ref_panel/chrX_tmp*

}

export -f filter_autosome
export -f filter_X_chr

mkdir $DIR/resources/1000G_ref_panel

if hash parallel 2>/dev/null; then
	printf "Great! parallel detected!\nrunning!\n"
	parallel -j8 --env g1k_path,DIR,filter_autosome 'echo processing chr {} && filter_autosome {}' ::: {1..22}
	echo processing chr X && filter_X_chr

	parallel -j8 --env DIR 'bcftools view -v snps -i "MAF>=0.01" -G resources/1000G_ref_panel/ALL.chr{}.phase3.20130502.SNP.indel.biallelic.mac5.vcf.gz -Oz --threads 2 > $DIR/resources/1000G_ref_panel/ALL.chr{}.phase3.20130502.SNP.biallelic.MAF0.01.sites.vcf.gz' ::: {X,{1..22}}
	parallel -j8 --env DIR 'bcftools index $DIR/resources/1000G_ref_panel/ALL.chr{}.phase3.20130502.SNP.biallelic.MAF0.01.sites.vcf.gz' ::: {X,{1..22}} 
else
	for chr in {1..22}; do echo processing chr${chr}; filter_autosome ${chr}; done
	echo processing chr X && filter_X_chr

	for chr in {X,{1..22}}; do bcftools view -v snps -i "MAF>=0.01" -G resources/1000G_ref_panel/ALL.chr${chr}.phase3.20130502.SNP.indel.biallelic.mac5.vcf.gz -Oz --threads 2 > $DIR/resources/1000G_ref_panel/ALL.chr${chr}.phase3.20130502.SNP.biallelic.MAF0.01.sites.vcf.gz; done
	for chr in {X,{1..22}}; do bcftools index $DIR/resources/1000G_ref_panel/ALL.chr${chr}.phase3.20130502.SNP.biallelic.MAF0.01.sites.vcf.gz; done
fi

printf "All done!"
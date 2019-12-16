# create g1k reference panel used for WEScall

set -euo pipefail

# path to g1k vcf files. Assume vcf files from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
g1k_path=/home/clwang/resource/1000G_phase3/data_v5a/
export g1k_path
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export DIR

function filter_autosome {
	chr=$1
	bcftools annotate ${g1k_path}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
		-h <(printf "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">") -Ou | \
		bcftools norm --rm-dup none -Ou | \
		bcftools view -m 2 -M 2 -Ou | \
		bcftools view -i 'MAC>=5 && INFO/VT!="SV"' -Oz --threads 10 > $DIR/resources/1000G_v5a_filtered/ALL.chr${chr}.phase3.20130502.SNP.indel.biallelic.mac5.vcf.gz
}

function filter_X_chr {
	bcftools annotate ${g1k_path}/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1a.20130502.genotypes.vcf.gz \
		-h <(printf "%s\n%s" "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" "##INFO=<ID=OLD_VARIANT,Number=1,Type=String,Description=\"unknown\">") -Ou | \
		bcftools norm --rm-dup none -Ou | \
		bcftools view -m 2 -M 2 -Ou | \
		bcftools view -i 'MAC>=5 && INFO/VT!="SV"' -Oz --threads 10 > $DIR/resources/1000G_v5a_filtered/ALL.chrX.phase3.20130502.SNP.indel.biallelic.mac5.vcf.gz

	bcftools convert --haploid2diploid -h $DIR/resources/1000G_v5a_filtered/chrX_tmp $DIR/resources/1000G_v5a_filtered/ALL.chrX.phase3.20130502.SNP.indel.biallelic.mac5.vcf.gz
	bcftools convert -H $DIR/resources/1000G_v5a_filtered/chrX_tmp -Oz > $DIR/resources/1000G_v5a_filtered/ALL.chrX.phase3.20130502.SNP.indel.biallelic.mac5.vcf.gz
	rm $DIR/resources/1000G_v5a_filtered/chrX_tmp*

}

export -f filter_autosome
export -f filter_X_chr

if [ ! -d $DIR/resources/1000G_v5a_filtered ]; then
	mkdir $DIR/resources/1000G_v5a_filtered
fi

echo {1..22} | tr " " "\n" | parallel --env g1k_path,DIR,filter_autosome 'filter_autosome {}'
filter_X_chr

echo {X,{1..22}} | tr " " "\n" | parallel --env DIR 'bcftools view -v snps -i "MAF>=0.01" -G resources/1000G_v5a_filtered/ALL.chr{}.phase3.20130502.SNP.indel.biallelic.mac5.vcf.gz -Oz --threads 2 > $DIR/resources/1000G_v5a_filtered/ALL.chr{}.phase3.20130502.SNP.biallelic.MAF0.01.sites.vcf.gz'
echo {X,{1..22}} | tr " " "\n" | parallel --env DIR 'bcftools index $DIR/resources/1000G_v5a_filtered/ALL.chr{}.phase3.20130502.SNP.biallelic.MAF0.01.sites.vcf.gz'
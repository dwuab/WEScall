# create g1k reference panel used for WEScall

set -euo pipefail

# path to g1k vcf files. Assume vcf files from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
g1k_path=/home/clwang/resource/1000G_phase3/data_v5a/
export g1k_path

if [ ! -d resources/data_v5a_filtered ]; then
	mkdir resources/data_v5a_filtered
fi

echo {1..22} | tr " " "\n" | parallel --env g1k_path 'bcftools reheader ${g1k_path}/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -h resources/hdr.txt | bcftools view -m 2 -M 2 -Ou | bcftools view -i "MAC>=5" -Oz --threads 10 > resources/data_v5a_filtered/ALL.chr{}.phase3.20130502.SNP.indel.biallelic.mac5.vcf.gz'
bcftools reheader ${g1k_path}/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1a.20130502.genotypes.vcf.gz -h resources/hdr_X.txt | bcftools view -m 2 -M 2 -Ou | bcftools view -i "MAC>=5" -Oz --threads 10 > resources/data_v5a_filtered/ALL.chrX.phase3.20130502.SNP.indel.biallelic.mac5.vcf.gz

echo {X,{1..22}} | tr " " "\n" | parallel 'bcftools view -v snps -i "MAF>=0.01" -G resources/data_v5a_filtered/ALL.chr{}.phase3.20130502.SNP.indel.biallelic.mac5.vcf.gz -Oz --threads 2 > resources/data_v5a_filtered/ALL.chr{}.phase3.20130502.SNP.biallelic.MAF0.01.sites.vcf.gz'
echo {X,{1..22}} | tr " " "\n" | parallel 'bcftools index resources/data_v5a_filtered/ALL.chr{}.phase3.20130502.SNP.biallelic.MAF0.01.sites.vcf.gz'
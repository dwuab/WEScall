# this script shows how we generate the indel interval list from g1k file
# first we find the longest allele (in terms of number of bases), let it be, N_max
# then the interval will be [POS,POS+N_max]

set -euo pipefail

bcftools reheader /home/clwang/resource/1000G_phase3/data_v5a/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz -h hdr_X_no_sample.txt | \
	bcftools view --type indels -Ou | \
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' | \
	awk -F'[\t,]' '{x=length($3);for(i=4;i<=NF;i++){if(length($i)>x)x=length($i)}print $1"\t"$2"\t"$2+x}' | bgzip -c > g1k.indel.interval.bed.gz
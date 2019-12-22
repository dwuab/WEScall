# a script to help you check whether necessary resources are located at the right place
# this script will not attempt to download ther resrouces from the Internet as they
# are so common that they are likely already on your computer/server, and you just 
# need to make links to resource files

set -euo pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" >/dev/null 2>&1 && pwd )"

printf "Checking whether required resources existed in the proper place.\n"
printf "You should see \"Checking completed!\" if everything is fine.\n"

if [ ! -e ${DIR}/resources/SeqCap_EZ_Exome_v3_primary.bed ]; then
	printf "unzipping resources/SeqCap_EZ_Exome_v3_primary.bed.gz\n"
	gunzip ${DIR}/resources/SeqCap_EZ_Exome_v3_primary.bed.gz -c > ${DIR}/resources/SeqCap_EZ_Exome_v3_primary.bed
fi

genetic_map_location=${DIR}/resources/geneticMap_GRCh37
printf "Expecting genetic map at %s\n" ${genetic_map_location}
if [ ! -d ${genetic_map_location} ]; then
	printf "Genetic maps not found at %s! Genetic maps can be downloaded from %s\n" \
		${genetic_map_location} \
		"http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip"
	exit
else
	echo OK!
fi

g1k_ref_location=${DIR}/resources/1000G_ref_panel
printf "Expecting 1KG reference panel at %s\n" ${g1k_ref_location}
if [ ! -d ${g1k_ref_location} ]; then
	printf "1KG reference panel not found at %s!\n" ${g1k_ref_location}
	printf "You can run the script, create_g1k_ref.sh, accompanying this package\n to generate the reference panel"
	exit
else
	echo OK!
fi

gotcloud_url=ftp://anonymous@share.sph.umich.edu/gotcloud/ref/hs37d5-db142-v1.tgz
g1k_omni25_location=${DIR}/pipelines/varCall/gotcloud.ref/1000G_omni2.5.b37.sites.PASS.vcf.gz
printf "Expecting 1000G_omni2.5.b37.sites.PASS.vcf.gz at %s\n" $g1k_omni25_location
if [ ! -e ${g1k_omni25_location} ]; then
	printf "%s not found!\n" ${g1k_omni25_location}
	printf "You can get it from gotcloud resource bundle at %s\n" ${gotcloud_url}
	exit
else
	echo OK!
fi

g1k_omni25_index_location=${DIR}/pipelines/varCall/gotcloud.ref/1000G_omni2.5.b37.sites.PASS.vcf.gz.tbi
printf "Expecting 1000G_omni2.5.b37.sites.PASS.vcf.gz.tbi at %s\n" $g1k_omni25_index_location
if [ ! -e ${g1k_omni25_index_location} ]; then
	printf "%s not found!\n" ${g1k_omni25_index_location}
	printf "You can index %s by tabix" ${g1k_omni25_index_location}
	exit
else
	echo OK!
fi

hapmap_sites_location=${DIR}/pipelines/varCall/gotcloud.ref/hapmap_3.3.b37.sites.vcf.gz
printf "Expecting hapmap_3.3.b37.sites.vcf.gz at %s\n" $hapmap_sites_location
if [ ! -e ${hapmap_sites_location} ]; then
	printf "%s not found!\n" ${hapmap_sites_location}
	printf "You can get it from gotcloud resource bundle at %s\n" ${gotcloud_url}
	exit
else
	echo OK!
fi

hapmap_sites_index_location=${DIR}/pipelines/varCall/gotcloud.ref/hapmap_3.3.b37.sites.vcf.gz.tbi
printf "Expecting hapmap_3.3.b37.sites.vcf.gz.tbi at %s\n" $hapmap_sites_index_location
if [ ! -e ${hapmap_sites_index_location} ]; then
	printf "%s not found!\n" ${hapmap_sites_index_location}
	printf "You can index %s by tabix" hapmap_3.3.b37.sites.vcf.gz
	exit
else
	echo OK!
fi

dbsnp_location=${DIR}/pipelines/varCall/gotcloud.ref/dbsnp_142.b37.vcf.gz
printf "Expecting dbsnp_142.b37.vcf.gz at %s\n" $dbsnp_location
if [ ! -e $dbsnp_location ]; then
	printf "%s not found!\n" $dbsnp_location
	printf "You can get it from gotcloud resource bundle at %s\n" $gotcloud_url
	exit
else
	echo OK!
fi

dbsnp_index_location=${DIR}/pipelines/varCall/gotcloud.ref/dbsnp_142.b37.vcf.gz.tbi
printf "Expecting dbsnp_142.b37.vcf.gz.tbi at %s\n" $dbsnp_index_location
if [ ! -e $dbsnp_index_location ]; then
	printf "%s not found!\n" $dbsnp_index_location
	printf "You can index %s by tabix" dbsnp_142.b37.vcf.gz.tb
	exit
else
	echo OK!
fi

refGenome_location=${DIR}/pipelines/varCall/gotcloud.ref/hs37d5.fa
printf "Expecting hs37d5.fa at %s\n" $refGenome_location
if [ ! -e ${refGenome_location} ]; then
	printf "%s not found!\n" ${refGenome_location}
	printf "You can get it from gotcloud resource bundle at %s\n" ${gotcloud_url}
	exit
else
	echo OK!
fi

refGenome_index_location=${DIR}/pipelines/varCall/gotcloud.ref/hs37d5.fa.fai
printf "Expecting hs37d5.fa.fai at %s\n" $refGenome_index_location
if [ ! -e ${refGenome_index_location} ]; then
	printf "%s not found!\n" ${refGenome_index_location}
	printf "Please run samtools index to get the reference genome index.\n"
	exit
else
	echo OK!
fi

printf "Checking completed!\n"

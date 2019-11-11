# a script to help you check whether necessary resources are located at the right place
# this script will not attempt to download ther resrouces from the Internet as they
# are so common that they are likely already on your computer/server, and you just 
# need to make links to resource files

set -e

printf "Checking whether required resources existed in the proper place.\n"
printf "You should see \"Checking completed!\" if everything is fine.\n"

if [ ! -e resources/SeqCap_EZ_Exome_v3_primary.bed ]; then
	printf "unzipping resources/SeqCap_EZ_Exome_v3_primary.bed.gz\n"
	gunzip resources/SeqCap_EZ_Exome_v3_primary.bed.gz -c > resources/SeqCap_EZ_Exome_v3_primary.bed
fi

genetic_map_location=resources/geneticMap_GRCh37
printf "Expecting genetic map at %s\n" ${genetic_map_location}
if [ ! -d ${genetic_map_location} ]; then
	printf "Genetic maps not found at %s! Genetic maps can be downloaded from %s\n" \
		${genetic_map_location} \
		"http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip"
	exit
else
	echo OK!
fi

g1k_ref_location=resources/data_v5a_filtered
printf "Expecting 1KG reference panel at %s\n" ${g1k_ref_location}
if [ ! -d ${g1k_ref_location} ]; then
	printf "1KG reference panel not found at %s!\n" ${g1k_ref_location}
	printf "You can run the script, create_g1k_ref.sh, accompanying this package\n to generate the reference panel"
	exit
else
	echo OK!
fi

url_gotcloud_bundle=ftp://share.sph.umich.edu/vt/grch37/
g1k_omni25_location=topMed/gotcloud.ref/1000G_omni2.5.b37.sites.PASS.vcf.gz
printf "Expecting %s at topMed/gotcloud.ref\n" $g1k_omni25_location
if [ ! -e ${g1k_omni25_location} ]; then
	printf "%s not found!\n" ${g1k_omni25_location}
	printf "You can download it from %s\n" "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_omni2.5.b37.vcf.gz"
	exit
else
	echo OK!
fi

refGenome_location=topMed/gotcloud.ref/hs37d5.fa
printf "Expecting %s at topMed/gotcloud.ref\n" $refGenome_location
if [ ! -e ${refGenome_location} ]; then
	printf "%s not found!\n" ${refGenome_location}
	printf "You can download it from gotcloud bundle at %s\n" ftp://share.sph.umich.edu/vt/grch37/hs37d5.fa
	exit
else
	echo OK!
fi

refGenome_index_location=topMed/gotcloud.ref/hs37d5.fa.fai
printf "Expecting %s at topMed/gotcloud.ref\n" $refGenome_index_location
if [ ! -e ${refGenome_index_location} ]; then
	printf "%s not found!\n" ${refGenome_index_location}
	printf "Please run samtools index to get the reference genome index.\n"
	exit
else
	echo OK!
fi

printf "Checking completed!\n"

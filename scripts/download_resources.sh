set -euo pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" >/dev/null 2>&1 && pwd )"

printf "downloading beagle genetic maps ...\n"

if [ ! -d $DIR/resources/geneticMap_GRCh37 ]; then
        mkdir $DIR/resources/geneticMap_GRCh37
fi

pushd $DIR/resources/geneticMap_GRCh37

wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip

unzip plink.GRCh37.map.zip

rm plink.GRCh37.map.zip

popd

if [ ! -d $DIR/pipelines/varCall/gotcloud.ref/ ]; then
	mkdir $DIR/pipelines/varCall/gotcloud.ref/
fi

printf "downloading gotcloud bundles ...\n"

pushd $DIR/pipelines/varCall/

wget ftp://anonymous@share.sph.umich.edu/gotcloud/ref/hs37d5-db142-v1.tgz

printf "uncompressing hs37d5-db142-v1.tgz...\n"

tar zxvf hs37d5-db142-v1.tgz

rm hs37d5-db142-v1.tgz


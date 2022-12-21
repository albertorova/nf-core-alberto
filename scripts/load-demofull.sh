#!/bin/sh

ADMIN_PASSWORD=$1
OPENCGA_HOME=${2:-/opt/opencga}
USER_PASSWORD='Demo2022_'

OPENCGA_BIN=$OPENCGA_HOME/bin
if [ ! -d $OPENCGA_BIN ]; then
    echo "Opencga executable folder $OPENCGA_BIN not found!"
    exit 1
fi

set -x

echo "Removing existent information of demofull user and demo_id study in MongoDB"
mongo $OPENCGA_HOME/misc/demo/mongo_remove_data.js <<< $ADMIN_PASSWORD

echo "Removing existent core in Solr"
/opt/solr/bin/solr delete -c opencga_demofull_demo_project

echo "Removing OpenCGA session information for demofull"
sudo rm -rf /app/data/sessions/users/demofull
sudo rm -rf /app/data/sessions/jobs/JOBS/demofull

echo "Creating user for OpenCGA Catalog ....."
echo "$ADMIN_PASSWORD" | $OPENCGA_BIN/opencga-admin.sh users create -u demofull --email demo@t-systems.com --name "Demo User" --user-password "$USER_PASSWORD" --organization "Demo health service"
echo "Login user demofull ...."
echo $USER_PASSWORD | $OPENCGA_BIN/opencga.sh users login -u demofull -p

    
echo "Creating demofull@demo_project:demo_study ...."
$OPENCGA_BIN/opencga.sh projects create --id 'demo_project' --name 'Demo Family Studies GRCh38'  \
    --organism-scientific-name 'homo sapiens' \
    --organism-assembly 'GRCh38'
$OPENCGA_BIN/opencga.sh studies create --project 'demofull@demo_project' --name 'Demo Family' --id 'demo_study' \
    --description 'This study simulates two disorders and some phenotypes in the Corpas family for training purposes'

echo "Moving corpasome yml data to /tmp (job won't work if placed in $OPENCGA_HOME/misc/demo/corpasome)"
rm -rf /tmp/corpasome
cp -r $OPENCGA_HOME/misc/demo/corpasome /tmp/
TEMPLATE=`$OPENCGA_BIN/opencga.sh studies templates-upload -i /tmp/corpasome/ --study 'demofull@demo_project:demo_study'`
echo "Waiting for OpenCGA to update template info"
sleep 10
$OPENCGA_BIN/opencga.sh studies templates-run --id "$TEMPLATE" --study 'demofull@demo_project:demo_study' --overwrite

$OPENCGA_BIN/opencga.sh files create --study 'demofull@demo_project:demo_study' --path 'data' --type 'DIRECTORY'

# $OPENCGA_BIN/opencga.sh files fetch --study 'demofull@demo_project:demo_study' --path 'data' --url 'http://resources.opencb.org/datasets/corpasome/data/quartet.variants.annotated.vcf.gz' \
#     --job-id 'download_quartet.variants.annotated.vcf.gz'
$OPENCGA_BIN/opencga.sh files link -i /var/lib/mongo/tmp/corpasome/quartet.variants.annotated.hg38.vcf.gz --study 'demofull@demo_project:demo_study' --path 'data'

$OPENCGA_BIN/opencga.sh operations variant-index --file 'quartet.variants.annotated.hg38.vcf.gz' --family \
    --job-id 'variant_index'  # --job-depends-on 'download_quartet.variants.annotated.vcf.gz'
$OPENCGA_BIN/opencga.sh operations variant-annotation-index --project 'demofull@demo_project' \
    --job-id 'variant_annotation' --job-depends-on 'variant_index'
$OPENCGA_BIN/opencga.sh operations variant-stats-index --study 'demofull@demo_project:demo_study' --cohort 'ALL' \
    --job-id 'variant_stats' --job-depends-on 'variant_annotation'
$OPENCGA_BIN/opencga.sh operations variant-secondary-index --project 'demofull@demo_project' \
    --job-id 'variant_secondary_index' --job-depends-on 'variant_stats,variant_annotation'

echo "Creating disease panels..."
for panel_json in $OPENCGA_HOME/misc/demo/disease_panels/*.json; do
    id=`grep '"id":' $panel_json | head -1 | awk '{print substr($2, 1, length($2)-1)}'`
    echo $OPENCGA_BIN/opencga.sh panels create --json-file $panel_json --id $id
    $OPENCGA_BIN/opencga.sh panels create --json-file $panel_json --id $id
done

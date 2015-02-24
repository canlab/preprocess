#!/bin/bash
today=`date +%Y_%m_%d`

tools='SCN_Core_Support 3DheadUtility_lite Atlas_Localization_Tools densityUtility diagnostics hewma_utility Logit_HRF_Toolbox mediation_toolbox OptimizeDesign11 preprocess robust_toolbox spmUtility'

for t in $tools
do
if [ -f ${t}_${today}.zip ]; then
    rm ${t}_${today}.zip
fi

echo ${t}_${today}.zip

zip -r ${t}_${today}.zip ${t}/ -x \*~ -x \*.svn\*

cp ${t}_${today}.zip ${t}_latest.zip

done

ls -lt *zip

currentdir=`pwd` 
echo "Copying *latest.zip to tor@canlab.colorado.edu:~/canlab_tools_zips"
echo "Enter CANlab password for user tor."
echo "Then on CANlab use: sudo cp -R ~/canlab_tools_zips/*latest.zip /Library/WebServer/External/files/tools"

echo "Copying today's files to canlab_tools_zips/archive"
rsync -avu ${currentdir}/*${today}.zip tor@canlab.colorado.edu:~/canlab_tools_zips/archive



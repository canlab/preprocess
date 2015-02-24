#! /bin/tcsh 
# created 10-17-91 dcn ftp raw files from signa2
# modified by Tor Wager 
# USE THIS on a remote computer to zip and PUT files to Kosslyn_G5
# ENTER wildcard as argument

##***To run one subject*****

#if ($#argv < 1) then
# echo "$0 [-h hostname] subjcode"
#  echo "please type in: sftp_put_to_kosslyn wildcard"
#  echo "subjidcode is subject id, e.g. sftp_put_to_kosslyn 0*"
#  exit 1
#endif

#if ( "$argv[1]" == "-h") then
#  set hostname = $argv[2]
#   shift
#  shift
#else
    set hostname = "dyn-233-79.dyn.columbia.edu"
#endif
*******************************

set sids = $argv[1]
echo $sids
set sids = `ls | grep $sids`

echo " "
echo "working on subjects:"
echo "_________________________"
echo $sids
echo "_________________________"
echo " "
set pid = $$
echo "script number is $pid"


echo cd /Users/scnlab/Documents/Data_and_Tools/IMAGING_DATA/dropbox >> script.$pid

foreach sid = ( $sids )

	# create archive
	echo "\!gtar cfz $sid.tar.gz $sid" >> script.$pid

	echo "put $sid*.gz" >> script.$pid

	echo "\!rm $sid.tar.gz" >> script.$pid
end

echo quit >> script.$pid

sftp scnlab@$hostname < script.$pid
\rm -f script.$$


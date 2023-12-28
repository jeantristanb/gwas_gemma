scl_enabled java-11

checklog() {
if [ ! -f $1 ]
then
echo $1" file doesn't exist"
return 0
fi
ErrorMessage=`grep "Error executing process" $1|wc -l`
ErrorMessage2=`grep "Completed at: " $1|wc -l`
if [ "$ErrorMessage" != "0" ]
then
 echo $1" error "$ErrorMessage2
 return 0
else 
 echo $1" ok "$ErrorMessage2" "
 return 1
fi
}

runsim() {
echo "-------------begin test "$1" -------------------------"
checklog $1".log";
if [ $? == "0" ]
then
echo "rerun "$1
echo " ./run_example.bash $1"
rm .nextflow.log
bash ./run_example.bash $1 &> $1".log"
if [ ! -f .nextflow.log ]
then
echo $1" : nextflow did not run "
exit
fi
cp .nextflow.log nf_$1".log"
echo $?
checklog $1".log"
if [ $? == "0" ]
then
echo "after rerun "$1" doesn't work"
echo "----------------------------- stop -----------------"
exit
fi
else 
echo " no rerun "$1
fi
echo "-------------end test "$1" -------------------------"
}
#runsim sr1_loco0_dosage0_vcfi_addpc
#runsim sr1_loco0_dosage0_vcfi
#runsim sr1_loco0_dosage0_vcfi_addpc_res_rin
#runsim sr1_loco0_dosage0_vcfi_addpc_nores
#runsim sr1_loco0_dosage1_vcfi_addpc_res_rin
#runsim sr1_loco1_dosage1_vcfi
#runsim sr1_loco0_dosage0_vcfmulti
#runsim sr1_loco1_dosage0_vcfi
#runsim sr1_loco0_dosage1_vcfi
#runsim sr1_loco0_dosage1_vcfi_covres
#runsim sr1_loco0_dosage1_vcfi_covnores
runsim sr1_loco1_dosage1_bimbammulti


for filelog in `ls *.log`
do
checklog $filelog
done

#!/bin/bash
#$ -l h_rt=48:00:00

#How many qsubs
nprocesses=400
#Limit of sim processes
proc_lim=400
dir=/fastdata-sharc/sm1epk/SDGVM_glob_run
#Input file for run.Several values should be set to ARGUMENT
inputfile=/home/sm1epk/SDGVM/SDGVM_ant/input_f/glob_run_new.dat
#exe directory
exedir=/home/sm1epk/SDGVM/SDGVM_ant/source/f03_ant/bin
#Output final directory
outputdir=$dir/fin_outs
#File with sites at resolution
sites=/data/sm1epk/SDGVM_drivers_new/land_use/sage_isam/land_sites_1d.dat

#Make temporary output folder
rm -fr $dir
tmpdir=$dir/tempoutput
mkdir -p $tmpdir

#Count sites and adjust number of processes
count=`wc -l < $sites`
if [ $nprocesses -gt $count ]
then
  nprocesses=$count
fi

#How many sites per process
chunk=$[($count-1)/$nprocesses+1]
echo sites=$count processes=$nprocesses chunk=$chunk

#Create batch file and folder for each chunk
#Creates the strings needed for the prompt
for i in $(seq 1 $nprocesses)
do
  outd=$tmpdir/r$i/
  batchfile=$tmpdir/batch-$i
  mkdir $outd
  let start=($i-1)*$chunk+1
  let end=$i*$chunk
  if [ $end -gt $count ]
  then
    end=$count
  fi
  echo "#!/bin/bash">> $batchfile
  echo "#$ -l h_rt=22:59:00">> $batchfile
  echo "module load compilers/gcc/5.2">> $batchfile
  echo "$exedir/sdgvm.exe $inputfile $outd $sites $start $end">> $batchfile
done


#Launch the jobs
for i in $(seq 1 $nprocesses)
do
  prunning=`Qstat | grep 'batch-'|grep 'sm1epk'| wc -l`
  while [ $prunning -ge $proc_lim ]
  do
      prunning=`Qstat | grep 'batch-'|grep 'sm1epk'| wc -l`
      echo $prunning
      sleep 60
  done

  outtest=$tmpdir/r$i/
  while [ ! -d "$outtest" ] 
  do
      echo directory $outtest does not exist
      sleep 10
  done   

  echo submitting job $i
  qsub $tmpdir/batch-$i
  sleep 10
done

#Wait for final job to finish
prunning=`Qstat | grep 'batch-'|grep 'sm1epk'| wc -l`
while [ $prunning -gt 0 ]
do

  for i in $(seq 1 $nprocesses)
  do
    if grep -Fq "PROGRAM TERMINATED" batch-$i.o*
    then
      rm batch-$i.o*
      rm batch-$i.e* 
      qsub $tmpdir/batch-$i   
    fi
  done

  sleep 300
  prunning=`Qstat | grep 'batch-'|grep 'sm1epk'| wc -l`
done  

echo Finished jobs starting merging files
  
#Create output folder
mkdir -p $outputdir    

#Removes first line from site_info.dat file
#It just keeps the header from the 1st one
for i in $(seq 2 $nprocesses)
do
  sed 1d $tmpdir/r$i/site_info.dat > $tmpdir/r$i/temp
  cp $tmpdir/r$i/temp $tmpdir/r$i/site_info.dat
  rm $tmpdir/r$i/temp
done
 
#Gets the names of the output files
files=`ls $tmpdir/r1`

#Merges output files
for i in $(seq 1 $nprocesses)
do
  echo merging batch $i
  for file in $files
  do
    cat $tmpdir/r$i/$file >> $outputdir/$file
  done
done
  





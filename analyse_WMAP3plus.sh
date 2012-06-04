#!/bin/bash
# for all the xfrac3d files in the directory given as 
# input do the following steps:
# - convert into single precision                        "0"
# - list the ion  redshifts of z in file "ion_redshift" "A1"
# - list the dens redshfits of z in file "dens_redshift""A2"
# - call the program fof_bubbles and write output to file
#                    fof_HII_z.out  and   fof_HI_z.out   "B"
# - call the program zahn_bubbles and write output to file
#                    zahn_z.out                          "C"
# - call the program minkowski and collect outputs in file
#                    mink_history                        "D"
# - call the program cc_spherical and write ouput to file 
#                    cc_z.out                            "E"

# as second input give the directory where output should be 
# saved
# to do the cross correltation.

# Some parameters-------------------------------------------
N=256  # number of cells in one direction
h=0.7  #WMAP3+ values as on the group site
# ----------------------------------------------------------

while [ -n "$1" ]; do    # as input parameters give a list 
                         # like this: "file" boxsize 
                         # example:
                         # "100Mpc_f250C_203_new" 100

# Name some directories-------------------------------------
dirneu=/home/martina/WMAP3plus_results/"$1" # this is the dir where the outputs are saved
olddir2="/home/martina/WMAP3plus"/"$1"/results # dir where the xfrac files are
if [ $2 = "114" ] 
   then olddir="/home/martina/WMAP3plus/114Mpc_densities" # dir where density files are
   elif [ $2 = "37" ]
   then olddir="/home/martina/WMAP3plus/37Mpc_densities"
fi
 mkdir $dirneu
# ----------------------------------------------------------

#---------------------------------------------------------A1
echo $olddir2
for g in "$olddir2"/xfrac3d_*.bin 
do
   base=`basename "$g" .bin `
   z="${base#xfrac3d_}" 
   echo $z>>$dirneu/redshiftlist_inter 
done  
mv $dirneu/redshiftlist_inter $dirneu/redshiftlist2
sort -n -r $dirneu/redshiftlist2 > $dirneu/ion_redshift
rm $dirneu/redshiftlist2
#-----------------------------------------------------------

#---------------------------------------------------------A2
for g in "$olddir/"*n_all.dat  
do
   base=`basename "$g" .dat`
   base2="${base%n_all}"
   echo $base2 >> $dirneu/redshift_inter
done
mv $dirneu/redshift_inter $dirneu/redshift2
sort -n -r $dirneu/redshift2 > $dirneu/dens_redshift
rm $dirneu/redshift2
#-----------------------------------------------------------
#-----------------------------------------------------------

# read in the whole redshift list from the density files----
count=0
old_IFS=$IFS     #to set field seperator variable
IFS=$'\n'        # to end of line and set it back later at (!)
allz=($(cat $dirneu/dens_redshift)) 
#-----------------------------------------------------------
# line by line, assign variables and call the correlation function
while read z1
do
read z2
z3=${allz[$count]}
let count=$count+1

./double2single256 "$olddir2/xfrac3d_"$z1".bin" intermed 
#--------------------------------------------------------"E"
cc_spherical "$olddir/${z3}n_all.dat"  "intermed" "$dirneu/cc_$z1.out"  $N $2 $h
#-----------------------------------------------------------
#--------------------------------------------------------"B"
  fof_bubbles "intermed" $dirneu/fof_HII_$z1.out 0.5 0 $2 $h    #for HII bubbles
  fof_bubbles "intermed" $dirneu/fof_HI_$z1.out 0.5 1  $2 $h    #for HI  bubbles
#-----------------------------------------------------------
#--------------------------------------------------------"C"
 zahn_bubbles "intermed" $dirneu/zahn_$z1.out $N $2 $h
#-----------------------------------------------------------

#--------------------------------------------------------"D"
  tail -c 67108868 intermed > "junk"   # to cut away the byte-size info from 
  head -c 67108864 "junk" > "$base"    # unformatted fortran output
  rm junk
  minkowski -x $N -y $N -z $N -b1 -l0.5 -h0.501 -s0 -o $base.mink -i"$base"
  echo $z1 >> $base.mink
  cat $base.mink >> $dirneu/mink_history_inter
  rm $base.mink
  rm $base

 ./redshiftlist "intermed" "$olddir/${z3}n_all.dat"  zout  #this calc the mass av ion frac
 while read a; do
 echo -n $a >> zout2
 done < zout
 echo -n -e ' ' $z1 '\n' >> zout2
 cat zout2 >> $dirneu/ionfracs_inter
 rm zout 
 rm zout2
echo end
#-----------------------------------------------------------
###ONLY FOR NEW SIMULATION THAT HAS ONE ION STEP PER DENS STEP
#########################3if [$z3 -gt 9.
#z3=${allz[$count]}   # that is for f0.4
#count=$count+1
############################3fi
###

# Now, do the same things again for the second xfrac file that is read for the 
# same density file:

sleep 2       # the option with read does not work because I "jump" over two redshifts then 
./double2single256 "$olddir2/xfrac3d_"$z2".bin" intermed           
cc_spherical "$olddir/${z3}n_all.dat"  "intermed" "$dirneu/cc_$z2.out" $N $2 $h
#read -t 2   
sleep 1

#--------------------------------------------------------"B"
  fof_bubbles "intermed" $dirneu/fof_HII_$z2.out 0.5 0 $2 $h    #for HII bubbles
  fof_bubbles "intermed" $dirneu/fof_HI_$z2.out 0.5 1  $2 $h    #for HI  bubbles
#-----------------------------------------------------------

#--------------------------------------------------------"C"
  zahn_bubbles "intermed" $dirneu/zahn_$z2.out $N $2 $h
#-----------------------------------------------------------
#--------------------------------------------------------"D"
  tail -c 67108868 intermed > "junk"
  head -c 67108864 "junk" > "$base"
  minkowski -x $N -y $N -z $N -b1 -l0.5 -h0.501 -s0 -o $base.mink -i"$base"
  echo $z2 >> $base.mink
  cat $base.mink >> $dirneu/mink_history_inter
  rm $base.mink
  rm junk
  rm $base
 ./redshiftlist "intermed" "$olddir/${z3}n_all.dat"  zout
 while read a; do
 echo -n $a >> zout2
 done < zout
 echo -n -e ' ' $z2 '\n' >> zout2
 cat zout2 >> $dirneu/ionfracs_inter
 rm zout 
 rm zout2
echo end 
#-----------------------------------------------------------

done < "$dirneu/ion_redshift"

IFS=$old_IFS   # (!)

mv $dirneu/mink_history_inter $dirneu/mink_history
mv $dirneu/ionfracs_inter $dirneu/ion_history

shift 2
echo NEW SIMULATION
done
rm intermed
    

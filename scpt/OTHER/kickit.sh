#!/bin/sh
# Scripts to run ocean spinup (phase3 and phase4)
# Written for Richard Yablonsky's research
# Author : Biju Thomas on 07/20/05
# Modified to run on cyclone with archived data and multiple climatologies
# Author : Richard Yablonsky on 07/27/05

#################################################################################
# SET PARAMETERS FOR EACH NEW RUN                                               #
#################################################################################
##for climo in gd50 gd25 lv25  # e.g. (gd50 gd25 lv25) OR (gd50m gd25m lv25m)
##do
##echo $climo

stormid=18L                  # e.g. 05L
pdy=20050922                 # e.g. 20050718
cyc=00                       # e.g. 00
name=rita                    # e.g. emily,rita,dennis
region=united                # e.g. united
sharp=gfdlm                  # i.e. gfdl,gfdlm
sharp2=sharpmcs              # i.e. sharpop,sharpmcs
loc=ncep                     # i.e. ncep,cyc
climo=gd50m                  # i.e. gd50,gd25,lv25,mdas,
                             #      gd50m,gd25m,lv25m
loop=rl2                     # i.e. NA,276,264,rl2
wind=msg                     # i.e. msg,hrd
trunc=newcd2_92pct           # i.e. orig,tcd90,hwnd,newcd,newcd2,newcd2_92pct
suffx=std                    # i.e. std,long,vlong
#################################################################################
# CREATE AND MODIFY 'gfdl_$region.parm.4.$pdy$cyc-$stormid-$suffx' IN '$parm_d' #
#################################################################################

# Define directory setup
export stormid
export pdy
export cyc
export region
export climo
home_d=/hurricane5/richard/oceanspin
scpt_d=$home_d/scpt
sstd_d=$home_d/sstd
exec_d=$home_d/exec
fixd_d=$home_d/fixd
clim_d=$fixd_d/$sharp2$loop$loc
parm_d=$home_d/parm
vitl_d=$home_d/vitl
hwnd_d=$home_d/hwnd
read_d=$home_d/read
work_d=$home_d/work/$pdy$cyc-$stormid-$climo-$loop-$wind-$trunc-$suffx
pha3_d=$work_d/phase3
pha4_d=$work_d/phase4
mm=`echo $pdy | cut -c5-6`
yymmdd=`echo $pdy | cut -c3-8`
yyyy=`echo $pdy | cut -c1-4`
export mm
export yymmdd
export yyyy

# Clean up work dir:!!!
if [   -d $work_d ];   then rm -rf   $work_d; fi
if [ ! -d $work_d ];   then mkdir -p $work_d; fi
if [ ! -d $pha3_d ];   then mkdir -p $pha3_d; fi
if [ ! -d $pha4_d ];   then mkdir -p $pha4_d; fi

# Use correct executable depending on CD truncation choice
cd $exec_d
cp gfdl_ocean_united_$trunc gfdl_ocean_united

# Extract gfs realtime sst and mask
cd $work_d
cp $sstd_d/lonlat.gfs.$yyyy lonlat.gfs
cp $sstd_d/sst.gfs.$pdy.dat sst.gfs.dat
cp $sstd_d/mask.gfs.$pdy.dat mask.gfs.dat

# Define all run instructions for phase3
cd $pha3_d
cp $parm_d/gfdl_$region.parm.3          parameters.inp.shell
sed -e "s/_date_/${yymmdd}${cyc}/g"  \
      -e "s/_name1_/${name}/g"     \
      -e "s/_name2_/${name}/g"     \
      parameters.inp.shell >parameters.inp
/bin/rm parameters.inp.shell

# Link all files for phase3 and run
ln -s -f parameters.inp                               fort.10
ln -s -f $clim_d/$climo"_initdata."$region.$mm        fort.13
ln -s -f nullfile                                     fort.15
ln -s -f $work_d/sst.gfs.dat                          fort.21
ln -s -f $work_d/mask.gfs.dat                         fort.22
ln -s -f $work_d/lonlat.gfs                           fort.23
ln -s -f $fixd_d/$sharp"_Hdeepgsu."$region            fort.66

$exec_d/gfdl_ocean_$region >gfdl_ocean.$region.phase3.out
mv RST.* RST.phase3.$region
/bin/rm fort.*
echo "Phase3 is finished"

# Do phase3 postprocessing
$read_d/readp1_$name
echo "Phase3 postprocessing done... check for .dat files"

# Prepare files for phase4 
cd $pha4_d
echo "$yymmdd$cyc" >date2day.inp
ln -s -f date2day.dat fort.61
$exec_d/gfdl_date2day <date2day.inp
tmp=`cat date2day.dat`
julday=`echo $tmp | awk -F, '{print $1}'`
julday=`echo "2k $julday 3 - p" | dc`
year=`echo $tmp | awk '{print $2}'`
echo "$year $julday" >day2date.inp
ln -s -f day2date.dat fort.61
$exec_d/gfdl_day2date <day2date.inp
start_date=`cat day2date.dat`
grep $stormid $vitl_d/syndat_tcvitals.$yyyy | $scpt_d/gfdl_pre_sortvit.sh > track
echo \
"000  000 000000    000000 0000 000  000 000  00  000 0000  000 00  00 0000 0000 0000 0000 0000 0000 0000 0000"\
>> track

cp $parm_d/gfdl_$region.parm.4.$pdy$cyc-$stormid-$suffx parameters.inp

# Copy all necessary HRD wind files for phase4 (if applicable)
##cp $hwnd_d/parameter1.inp.$name parameter1.inp
##cp $hwnd_d/$name/*.* .

sed "s/_date_/$start_date/g" parameters.inp     > tmp
sed "s/_name1_/track/g"      tmp                > parameters.inp
sed "s/_name2_/track/g"      parameters.inp     > tmp
sed "s/_rstph3file_/RST.phase3.$region/g" tmp > parameters.inp

# Link all files for phase4 and run
ln -s -f parameters.inp                               fort.10
ln -s -f $clim_d/$climo"_initdata."$region.$mm        fort.13
ln -s -f $work_d/sst.gfs.dat                          fort.21
ln -s -f $work_d/mask.gfs.dat                         fort.22
ln -s -f $work_d/lonlat.gfs                           fort.23
ln -s -f $fixd_d/$sharp"_Hdeepgsu."$region            fort.66
ln -s -f $pha3_d/RST.phase3.$region                   fort.14
ln -s -f track                                        fort.15

/bin/rm RST.*
$exec_d/gfdl_ocean_$region >gfdl_ocean.$region.phase4.out
/bin/rm fort.*
mv RST.* RST.final
echo " phase4 finished"
if [ -s RST.final ]
then
  echo "Probable phase 4 ocean spin-up success."
  # Do phase4 postprocessing
  $read_d/readp2_$name
  echo "Phase4 postprocessing done... check for .dat files"
else
  echo "ERROR: Phase 4 ocean spin-up failure.  Not my fault ..Good Luck"
fi

cd $scpt_d

##done

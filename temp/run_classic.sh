# Purpose of script:

# Copy classic_submit.sh from generalTools to output directory
# Edit classic_submit.sh
# Copy initialization file to output directory and rename it to rsFile.nc
# Copy job_options.txt from source code directory to output directory
# Copy run_parameters.txt from source code directory to output directory
# Edit job_options.txt
# Submit run

#-------------------------------------------------------
# (I) Provide inputs below
#-------------------------------------------------------
# classic_submit.sh
#-------------------------------------------------------

# simulation ID
simulationID=daisyRun

# Meteorological forcing
metForcing=CRUJRAv2.4.5

# spinup or transient?
simulationType=transient # current options: spinup, transient, ssp585

# Is this the last 100 years of the spinup with spinfast = 1?
lastPart=TRUE # TRUE or FALSE

# Restart file
# init_file='/home/aadams/projects/def-cseiler-ab/aadams/classic_inputs/initfiles/global/T63/rsFilev02.nc'
# init_file='/home/aadams/projects/def-cseiler-ab/aadams/data-assimilation/simulations/spinup_CRUJRAv2.4.5/netcdf_files/rsFile_modified.nc'
# init_file='/home/aadams/projects/def-cseiler-ab/aadams/data-assimilation/simulations/spinup_CRUJRAv2.4.5_lastPart/netcdf_files/rsFile_modified.nc'
# init_file='/home/aadams/projects/def-cseiler-ab/aadams/data-assimilation/simulations/transient_CRUJRAv2.4.5_2000/netcdf_files/rsFile_modified.nc'
init_file='/home/aadams/projects/def-cseiler-ab/aadams/data-assimilation/rsFile_daisy.nc'
# init_file='/home/aadams/projects/def-cseiler-ab/aadams/data-assimilation/simulations/spinup_CRUJRAv2.4.5_opt_lastPart/netcdf_files/rsFile_modified.nc'
# For DRA cluster: Make sure that there are no NaN's in FLND, set them to zero.

# Run CLASSIC globally (global), a subset of grid cells (multipleGridCells) or a single grid cell only (local)?
spatialCoverage=multipleGridCells # global # multipleGridCells

# Switches
PFTCompetition=.false.
dofire=.true.
Ncycle_on=.false. 

# Path of CLASSIC output directory 
outputDir="/home/aadams/projects/def-cseiler-ab/aadams/data-assimilation/simulations/$simulationID"

# Path of CLASSIC executable
sourceCodeDir="/home/aadams/classic/classic_code/classic"

# email
email=19aja3@queensu.ca

xmlFile='/home/aadams/classic/classic_code/classic/configurationFiles/outputVariableDescriptors_daisy.xml'
# xmlFile='/home/aadams/classic/classic_code/classic/configurationFiles/outputVariableDescriptors.xml'

# met data
metDataDir="/home/aadams/projects/def-cseiler-ab/aadams/classic_inputs/meteorology/global/T63"

if [ $spatialCoverage == local ]
then
lon=120
lat=65
fi

#-------------------------------------------------------
# job_options.txt
#-------------------------------------------------------
Comment='develop' # e.g. branch name
#-------------------------------------------------------

# PFT competition
start_bare=.false.
inibioclim=.true.

if [ $PFTCompetition == .true. ] && [ $simulationType == spinup ] && [ $lastPart == FALSE ]
then
start_bare=.true.
inibioclim =.false.
fi

if [ $simulationType == spinup ]
then
fertilizeron=.true.
transientFER=.false.
depositionon=.true.
transientDEP=.false.
fi

if [ $simulationType != spinup ]
then
fertilizeron=.true.
transientFER=.true.
depositionon=.true.
transientDEP=.true.
fi


if [ $metForcing == CRUJRAv2.4.5 ]
then
# Met forcings
metFileFss=$metDataDir'/CRUJRAv2.4.5/tswrf_crujra_v2.4.5d_T63_1700_2022.nc'
metFileFdl=$metDataDir'/CRUJRAv2.4.5/dlwrf_crujra_v2.4.5d_T63_1700_2022.nc '
metFilePre=$metDataDir'/CRUJRAv2.4.5/pre_crujra_v2.4.5d_T63_1700_2022.nc'
metFileTa=$metDataDir'/CRUJRAv2.4.5/tmp_crujra_v2.4.5d_T63_1700_2022.nc'
metFileQa=$metDataDir'/CRUJRAv2.4.5/spfh_crujra_v2.4.5d_T63_1700_2022.nc'
metFileUv=$metDataDir'/CRUJRAv2.4.5/wind_crujra_v2.4.5d_T63_1700_2022.nc'
metFilePres=$metDataDir'/CRUJRAv2.4.5/pres_crujra_v2.4.5d_T63_1700_2022.nc'

# Other forcings
FERFile='/home/aadams/projects/def-cseiler-ab/aadams/classic_inputs/nitrogen/global/T63/nfer_TRENDYv12_t63.nc'
DEPFile='/home/aadams/projects/def-cseiler-ab/aadams/classic_inputs/nitrogen/global/T63/ndep_TRENDYv12_t63.nc'
CO2File='/home/aadams/projects/def-cseiler-ab/aadams/classic_inputs/CO2/TRENDY_v12_CO2_1700-2022_GCP2023_fixed.nc'
# POPDFile='/home/aadams/projects/def-cseiler-ab/aadams/classic_inputs/population/global/T63/popd_trendy_v12_1700_2023_t63.nc'
POPDFile='/home/aadams/projects/def-cseiler-ab/aadams/classic_inputs/population/global/T63/popd_trendy_v12_1700_2023_t63_Test03.nc'
LUCFile='/home/aadams/projects/def-cseiler-ab/aadams/classic_inputs/landcover/global/T63/ESA_CCI_land_cover_with_LUH_TRENDY_v12_crops_9_PFTs_T63_1700_2022.nc'
LGHTFile='/home/aadams/projects/def-cseiler-ab/aadams/classic_inputs/lightning/global/T63/lisotd_1995_2014_climtlgl_lghtng_as_ts_1700_2050_chunked.nc'
fi

# Additional inputs
CH4File='/home/aadams/projects/def-cseiler-ab/aadams/classic_inputs/CH4/CH4_1700_2017_for_GCP_2018.nc'
OBSWETFFile='/home/aadams/projects/def-cseiler-ab/aadams/classic_inputs/wetland/global/T63/gcp-ch4_wetlands_1838-2017_t63_final_daily.nc'
alb4BandParamsFile='/home/aadams/projects/def-cseiler-ab/aadams/classic_inputs/albedo_lookup/classic4BandAlbedoLookUpFile.nc'

# Model spin up
if [ $metForcing == CRUJRAv2.4.5 ] && [ $simulationType == spinup ]
then
fixedYear=1700
transientCO2=.false.
lnduseon=.false.
transientLGHT=.false.
readMetStartYear=1901
readMetEndYear=1920
spinfast=10
metLoop=20 # 20 for 400 years spin up

if [ $Ncycle_on == .true. ]
    then
    spinfast=10
    metLoop=100
fi

if [ $lastPart == TRUE ] && [ $Ncycle_on == .false. ]
    then
    spinfast=1
    metLoop=5
fi

if [ $lastPart == TRUE ] && [ $Ncycle_on == .true. ]
    then
    spinfast=1
    metLoop=25
fi

leap=.false.
allLocalTime=.false.
transientPOPD=.false.

domonthoutput=.false.
doperpftoutput=.false.
doAnnualOutput=.true.
JMOSTY=$readMetStartYear
fi

if [ $metForcing == CRUJRAv2.4.5 ] && [ $simulationType == transient ]
then
fixedYear=1700
transientCO2=.true.
lnduseon=.true.
transientLGHT=.true.
readMetStartYear=2001 # 1701 # 2001-2010 for tuning
readMetEndYear=2010 # 2022 # 2010 # 2001-2010 for tuning
spinfast=1
metLoop=1
leap=.false.
allLocalTime=.false.
transientPOPD=.true.

domonthoutput=.true.
doperpftoutput=.false.
doAnnualOutput=.false.
JMOSTY=1970
fi

#-------------------------------------------------------
# Info on experimental protocol (TRENDY, v9, S3)
#-------------------------------------------------------

# 1700 Model spin up:
# CO2: 1700 CO2 concentration (276.59ppm)
# Climate: recycling climate mean and variability from the early decades of the 20th century (i.e. 1901-1920)
# LULCC: constant 1700 LUC (crops and pasture distribution)

# 1701-1900 transient simulation:
# CO2: varying CO2
# Climate: continue recycling spin up climate (i.e. 1901-1920)
# LULCC: varying LUC

# 1901-2019 transient simulation:
# CO2: varying CO2
# Climate: varying climate
# LULCC: varying LUC

# Transient runs for other simulations:
# CanESM5 and CanESM5-ISIMIP3b: transient from 1850-2014
# GSWP3-W5E5: transient from 1901-2019

#-------------------------------------------------------
# (II) Code below should work as is, no changes required
#-------------------------------------------------------
# Create output directory
rm -rf $outputDir
mkdir $outputDir
mkdir $outputDir/outputFiles

executableDir=$sourceCodeDir"/bin/CLASSIC_parallel_intel"

echo
echo "CLASSIC outputs will be stored here:" 
echo
echo $outputDir
echo
echo "CLASSIC source code is located here:"
echo
echo $sourceCodeDir
echo

#-------------------------------------------------------
# Copy classic_submit.sh from generalTools to output directory
#-------------------------------------------------------
cp /home/aadams/classic/classic_production/templates/classic_submit_dra_ntask20.sh $outputDir/classic_submit.sh
if [ $spatialCoverage == multipleGridCells ]
then
cp /home/aadams/classic/classic_production/templates/classic_submit_dra_daisy.sh $outputDir/classic_submit.sh
fi
echo "Copy classic_submit.sh to output directory completed"
echo

#-------------------------------------------------------
# Copy the initialization file and rename it as restart file
#-------------------------------------------------------
cp $init_file $outputDir/rsFile.nc
echo "Copy initialization file and rename it as restart file completed"
echo
#-------------------------------------------------------

#-------------------------------------------------------
# Edit classic_submit.sh
#-------------------------------------------------------

# Set your output directory in classic_submit.sh
# sed -i "s|old|new|g" classic_submit.sh

sed -i "s|output_directory=.*|output_directory=$outputDir|g" $outputDir/classic_submit.sh
echo "Output directory updated in classic_submit.sh"
echo
# set your executable directory in classic_submit.sh
sed -i "s|executable=.*|executable=$executableDir|g" $outputDir/classic_submit.sh
echo "Executable directory updated in classic_submit.sh"
echo

# set your email in classic_submit.sh
sed -i "s|email=.*|email=$email|g" $outputDir/classic_submit.sh
echo "Email updated in classic_submit.sh"
echo

# Edit name of job options file:
# sed -i "s|template_job_options_file.txt|job_options.txt|g" $outputDir/classic_submit.sh

#-------------------------------------------------------
# Copy run_parameters.txt
#-------------------------------------------------------
# Copy template_run_parameters.txt to CLASSIC output folder
# cp $sourceCodeDir/configurationFiles/default_run_parameters.txt $outputDir/run_parameters.txt
# cp /home/aadams/projects/def-cseiler-ab/aadams/data-assimilation/template_run_parameters.txt $outputDir/run_parameters.txt
# cp /home/aadams/projects/def-cseiler-ab/aadams/data-assimilation/run_parameters.txt $outputDir/run_parameters.txt
cp run_parameters.txt $outputDir/run_parameters.txt
sleep 1
rm run_parameters.txt

runparams_file=$outputDir/run_parameters.txt
echo "Copy run_parameters.txt to output directory completed"
echo

#-------------------------------------------------------
# Copy and edit job_options.txt
#-------------------------------------------------------

# Copy job_options.txt and ouputVariablesDescriptior to CLASSIC output folder
cp $sourceCodeDir/configurationFiles/template_job_options_file.txt $outputDir/job_options.txt
cp $sourceCodeDir/configurationFiles/outputVariableDescriptors.xml $outputDir/outputVariableDescriptors.xml

echo "Copy job_options.txt from source code directory to output directory completed"
echo

#---------------------------------------------------
# CRUJRAv2.4.5
#---------------------------------------------------

echo "CLASSIC will use the following met forcing files:"
echo $metFileFss
echo $metFileFdl
echo $metFilePre
echo $metFileTa
echo $metFileQa
echo $metFileUv
echo $metFilePres

# Update Met path:
sed -i "s|metFileFss =.*|metFileFss ='$metFileFss'|g" $outputDir/job_options.txt
sed -i "s|metFileFdl =.*|metFileFdl ='$metFileFdl'|g" $outputDir/job_options.txt
sed -i "s|metFilePre =.*|metFilePre ='$metFilePre'|g" $outputDir/job_options.txt
sed -i "s|metFileTa =.*|metFileTa ='$metFileTa'|g" $outputDir/job_options.txt
sed -i "s|metFileQa =.*|metFileQa ='$metFileQa'|g" $outputDir/job_options.txt
sed -i "s|metFileUv =.*|metFileUv ='$metFileUv'|g" $outputDir/job_options.txt
sed -i "s|metFilePres =.*|metFilePres ='$metFilePres'|g" $outputDir/job_options.txt

echo
echo "Met file path updated in job_options.txt"
echo

sed -i "s|Ncycle_on =.*|Ncycle_on =$Ncycle_on|g" $outputDir/job_options.txt
echo
echo "Ncycle_on is set to" $Ncycle_on
echo

# Update 
rs_file_to_overwrite=$outputDir/rsFile.nc
sed -i "s|rs_file_to_overwrite =.*|rs_file_to_overwrite ='$rs_file_to_overwrite'|g" $outputDir/job_options.txt
sed -i "s|init_file =.*|init_file ='$init_file'|g" $outputDir/job_options.txt
sed -i "s|CH4File =.*|CH4File ='$CH4File'|g" $outputDir/job_options.txt
sed -i "s|OBSWETFFile =.*|OBSWETFFile ='$OBSWETFFile'|g" $outputDir/job_options.txt
sed -i "s|FERFile =.*|FERFile ='$FERFile'|g" $outputDir/job_options.txt
sed -i "s|DEPFile =.*|DEPFile ='$DEPFile'|g" $outputDir/job_options.txt
sed -i "s|fertilizeron =.*|fertilizeron =$fertilizeron|g" $outputDir/job_options.txt
sed -i "s|transientFER =.*|transientFER =$transientFER|g" $outputDir/job_options.txt
sed -i "s|depositionon =.*|depositionon =$depositionon|g" $outputDir/job_options.txt
sed -i "s|transientDEP =.*|transientDEP =$transientDEP|g" $outputDir/job_options.txt
sed -i "s|alb4BandParamsFile =.*|alb4BandParamsFile ='$alb4BandParamsFile'|g" $outputDir/job_options.txt
sed -i "s|fixedYearCO2 =.*|fixedYearCO2 =$fixedYear|g" $outputDir/job_options.txt
sed -i "s|fixedYearPOPD =.*|fixedYearPOPD =$fixedYear|g" $outputDir/job_options.txt
sed -i "s|POPDFile =.*|POPDFile ='$POPDFile'|g" $outputDir/job_options.txt
sed -i "s|fixedYearLGHT =.*|fixedYearLGHT =$fixedYear|g" $outputDir/job_options.txt
sed -i "s|fixedYearLUC =.*|fixedYearLUC =$fixedYear|g" $outputDir/job_options.txt
sed -i "s|fixedYearFER =.*|fixedYearFER =$fixedYear|g" $outputDir/job_options.txt
sed -i "s|fixedYearDEP =.*|fixedYearDEP =$fixedYear|g" $outputDir/job_options.txt
sed -i "s|transientCO2 =.*|transientCO2 =$transientCO2|g" $outputDir/job_options.txt
sed -i "s|tracerCO2file =.*|tracerCO2file ='$metFileFss'|g" $outputDir/job_options.txt # dummy file
sed -i "s|CO2File =.*|CO2File ='$CO2File'|g" $outputDir/job_options.txt
sed -i "s|PFTCompetition =.*|PFTCompetition =$PFTCompetition|g" $outputDir/job_options.txt
sed -i "s|start_bare =.*|start_bare =$start_bare|g" $outputDir/job_options.txt
sed -i "s|lnduseon =.*|lnduseon =$lnduseon|g" $outputDir/job_options.txt
sed -i "s|LUCFile =.*|LUCFile ='$LUCFile'|g" $outputDir/job_options.txt
sed -i "s|LGHTFile =.*|LGHTFile ='$LGHTFile'|g" $outputDir/job_options.txt
sed -i "s|transientLGHT=.*|transientLGHT =$transientLGHT|g" $outputDir/job_options.txt
sed -i "s|inibioclim =.*|inibioclim =$inibioclim|g" $outputDir/job_options.txt
sed -i "s|dofire =.*|dofire =$dofire|g" $outputDir/job_options.txt
sed -i "s|readMetStartYear =.*|readMetStartYear =$readMetStartYear|g" $outputDir/job_options.txt
sed -i "s|readMetEndYear =.*|readMetEndYear =$readMetEndYear|g" $outputDir/job_options.txt
sed -i "s|metLoop =.*|metLoop =$metLoop|g" $outputDir/job_options.txt
sed -i "s|leap =.*|leap =$leap|g" $outputDir/job_options.txt
sed -i "s|spinfast =.*|spinfast =$spinfast|g" $outputDir/job_options.txt
sed -i "s|transientPOPD =.*|transientPOPD =$transientPOPD|g" $outputDir/job_options.txt
sed -i "s|domonthoutput =.*|domonthoutput =$domonthoutput|g" $outputDir/job_options.txt
sed -i "s|doperpftoutput =.*|doperpftoutput =$doperpftoutput|g" $outputDir/job_options.txt
sed -i "s|doAnnualOutput =.*|doAnnualOutput =$doAnnualOutput|g" $outputDir/job_options.txt
sed -i "s|JMOSTY =.*|JMOSTY =$JMOSTY|g" $outputDir/job_options.txt
sed -i "s|Comment =.*|Comment ='$Comment'|g" $outputDir/job_options.txt
sed -i "s|runparams_file =.*|runparams_file ='$runparams_file'|g" $outputDir/job_options.txt
sed -i "s|xmlFile =.*|xmlFile ='$xmlFile'|g" $outputDir/job_options.txt 
sed -i "s|allLocalTime =.*|allLocalTime =$allLocalTime|g" $outputDir/job_options.txt

# Replace default path with no path for files that you don't need
sed -i "s|blackCarbonFile =.*|blackCarbonFile ="",|g" $outputDir/job_options.txt
sed -i "s|prescribedFireFile =.*|prescribedFireFile ="",|g" $outputDir/job_options.txt
sed -i "s|timberHarvestFile =.*|timberHarvestFile ="",|g" $outputDir/job_options.txt


# Change output directory so that output files are written to virtual temporary file system
# sed -i "s|output_directory =.*|output_directory = '/tmp/$simulationID'|g" $outputDir/job_options.txt

echo
echo "The joboptions file is now set up for a" $simulationType
echo

#---------------------------------------------------
# Submit run
#---------------------------------------------------

if [ $spatialCoverage == global ]
then
# Copy this script into output folder
cp run_classic.sh $outputDir/run_classic.sh
sleep 1
# Remove the script from the current location so that there are no duplicates.
rm run_classic.sh
# Submit classic run
cd $outputDir
chmod +x classic_submit.sh
./classic_submit.sh
fi

if [ $spatialCoverage == multipleGridCells ]
then
# Copy this script into output folder
cp run_classic.sh $outputDir/run_classic.sh
sleep 1
# Remove the script from the current location so that there are no duplicates.
rm run_classic.sh
# Submit classic run
cd $outputDir
chmod +x classic_submit.sh
./classic_submit.sh
fi

if [ $spatialCoverage == local ]
then
# Copy this script into output folder
cp run_classic.sh $outputDir/run_classic.sh
sleep 1
# Remove the script from the current location so that there are no duplicates.
rm run_classic.sh
# Create the submission script
cd $outputDir

touch classic_submit_gridCell.sh

echo "#!/bin/sh
#SBATCH --account=def-cseiler-ab 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# #SBATCH --exclusive
#SBATCH --mem=10G
#SBATCH --time=00:30:00
#SBATCH --job-name=$simulationID
#SBATCH --output=/home/aadams/projects/def-cseiler-ab/aadams/data-assimilation/simulations/classic.out
#SBATCH --error=errors_CLASSIC
# #SBATCH --mail-user=$email
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

cd $outputDir
/home/aadams/classic/classic_code/classic/bin/CLASSIC_parallel_intel job_options.txt $lon/$lat
touch /home/aadams/projects/def-cseiler-ab/aadams/data-assimilation/simulations/classicFinished" >> classic_submit_gridCell.sh

chmod +x $outputDir/classic_submit_gridCell.sh

# Submit classic run
sbatch  $outputDir/classic_submit_gridCell.sh
fi 

#!/bin/bash

LOG_TO_IFDH=0

LOGYLOG () {
  echo ${1}
  if [ ${LOG_TO_IFDH} == "1" ]; then
    ifdh log ${1}
  fi
}

PNFS_PATH_APPEND=${1}
if [ -z ${1} ]; then
  LOGYLOG "[ERROR]: Failed to find PNFS_PATH_APPEND passed on command line."
  exit 2
fi

ANAVERSION="4"
if [ ! -z ${2} ]; then
  ANAVERSION="${2}"
fi

#if [ ! -e ${INPUT_TAR_DIR_LOCAL}/CAFAna/InputStateFiles.list ]; then
#  LOGYLOG "[ERROR]: Expected to recieve a CAF file list @ CAFAna/InputStateFiles.list but didn't."
#  ls ${INPUT_TAR_DIR_LOCAL}/CAFAna
#  exit 2
#fi

printenv

set -x #start bash debugging at this point
LOGYLOG "Start $(date)"
LOGYLOG "Site:${GLIDEIN_ResourceName}"
LOGYLOG "The worker node is " `hostname` "OS: " `uname -a`
LOGYLOG "The user id is $(whoami)"
LOGYLOG "The output of id is: $(id)"
set +x #stop bash debugging at this point

if [ -z ${GRID_USER} ]; then
  GRID_USER=$(basename $X509_USER_PROXY | cut -d "_" -f 2)
fi

if [ -z ${GRID_USER} ]; then
  LOGYLOG "Failed to get GRID_USER."
  exit 2
fi

echo "Start to move files, but will see warning mv: cannot remove ${INPUT_TAR_DIR_LOCAL}/CAFAna/*': Read-only file system (ignore for now)"
mv ${INPUT_TAR_DIR_LOCAL}/CAFAna $_CONDOR_SCRATCH_DIR/

cd $_CONDOR_SCRATCH_DIR

export CAFANA=$(readlink -f CAFAna)
source ${CAFANA}/CAFAnaEnv.sh

voms-proxy-info --all
setup ifdhc
source ${CAFANA}/CAFAnaEnv.sh

export IFDH_CP_UNLINK_ON_ERROR=1;
export IFDH_CP_MAXRETRIES=2;

ifdh ls ${PNFS_OUTDIR}
PNFS_OUTDIR=/pnfs/dune/scratch/users/${GRID_USER}/${PNFS_PATH_APPEND}
LOGYLOG "Output dir is ${PNFS_OUTDIR}"

myinfilehadron=""
# PROCESS starts from 0, 1, ... N-1
(( LINE_N = ${PROCESS} + 1 ))

for ifile in $(cat ${CAFANA}/myOutputFDGeoEffHadron.txt | head -${LINE_N} | tail -1); do
  myinfilehadron=${ifile}
done

myinfilemu=""
#loop over ntuples with muon combined
for ifile2 in $(cat ${CAFANA}/myOutputFDGeoEffMu.txt | head -${LINE_N} | tail -1); do
  myinfilemu=${ifile2}
done

echo "Got input file hadron: $myinfilehadron"
echo "Got input file mu: $myinfilemu"

if [ $? -ne 0 ]; then
  LOGYLOG "Unable to read ${PNFS_OUTDIR}. Make sure that you have created this directory and given it group write permission (chmod g+w ${PNFS_OUTDIR})."
  exit 10
fi

# LOGYLOG "Building interps @ $(date)"
#
# OUTFILENAME=""
# OUTFILENAME=Hadded_State.${CLUSTER}.${PROCESS}.root
#
# LOGYLOG "Output file name: ${OUTFILENAME}"

export CAFANA_ANALYSIS_VERSION=${ANAVERSION}
echo "CAFANA_ANALYSIS_VERSION=${CAFANA_ANALYSIS_VERSION}"

source ${CAFANA}/CAFAnaEnv.sh

cp ${CAFANA}/FileWithCoeffsNuMu.root .
echo "copying root files to node : ifdh cp -D $myinfilehadron ."
ifdh cp -D $myinfilehadron .
echo "copying root files to node : ifdh cp -D $myinfilemu ."
ifdh cp -D $myinfilemu .
ls -lhrt

myFileHadronCut="${myinfilehadron##*/}"
echo "this file should be list in ls $myFileHadronCut"
myFileMuonCut="${myinfilemu##*/}"
echo "this file should be listed in ls $myFileMuonCut"

##for NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA
LOGYLOG "Start of code at  $(date)"
LOGYLOG "NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA $myFileHadronCut $myFileMuonCut"
NtupleOutVetoAndTrimE_AssumeEqualEffAtAllOA $myFileHadronCut $myFileMuonCut

# Unique name in case we send multiple jobs.
OUTFILE=EtrimHistosFile_CoeffsNoOsc_OAPosSameAsCAFs_MuonCombined_VisEtrim_${CLUSTER}_${PROCESS}.root

if [ -f FileWithHistEtrim_CoeffsAndOAPoswithSameEff_FirstMuResultsIncludedWithVisEtrim_NuOscProc.root ]; then

  echo "mv FileWithHistEtrim_CoeffsAndOAPoswithSameEff_FirstMuResultsIncludedWithVisEtrim_NuOscProc.root $OUTFILE"
  mv FileWithHistEtrim_CoeffsAndOAPoswithSameEff_FirstMuResultsIncludedWithVisEtrim_NuOscProc.root $OUTFILE
LOGYLOG "End of code at  $(date)"
  # and copy our output file back
  echo "copy command: ifdh cp -D $OUTFILE $PNFS_OUTDIR"
  ifdh cp -D $OUTFILE $PNFS_OUTDIR
  ls $PNFS_OUTDIR
  # check the exit status to see if the copyback actually worked. Print a message if it did not.
  IFDH_RESULT=$?
  #if [ $IFDH_RESULT -ne 0 ]; then
   # echo "Error during output copyback. See output logs."
   # exit $IFDH_RESULT
  #fi
fi

echo " remove root file from the batch cluster: rm -r $myFileHadronCut rm -r $myFileMuonCut"
rm -r $myFileHadronCut 
rm -r $myFileMuonCut
pwd
ls -lhrt 
#If we got this far, we succeeded.
echo "Completed successfully."
exit 0

#LOGYLOG "hadd_cafana ${OUTFILENAME} ${INPFILE}"

#hadd_cafana ${OUTFILENAME} ${INPFILE}


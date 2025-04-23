# Instruction for running translation and rotations on FD n-tuples

Prerequisite: [Produce Ntuple from DUNE FD MC files](https://github.com/weishi10141993/myntuples#produce-ntuple-from-dune-fd-mc-files). The produced FD n-tuples will be used as input files for the following program to run.

[First time only]
```
cd /dune/app/users/<your_username>
mkdir NDEff (first time only)
cd NDEff
git clone --recurse-submodules -b FD_Wei https://github.com/weishi10141993/DUNE_ND_GeoEff.git      # Get geoEff library
# Note for git version (git --version) before 2.13, use: git clone --recursive -b FD_Wei https://github.com/weishi10141993/DUNE_ND_GeoEff.git
cd DUNE_ND_GeoEff
source setup.sh                                                                                    # Necessary setups for build
cmake -DCMAKE_CXX_STANDARD=17 -DPYTHON_EXECUTABLE:FILEPATH=`which python` .                        # -DCMAKE_CXX_STANDARD=17 for running with root 6.28.10
make -j geoEff                                                                                     # Build geoEff (can also use: make -j pyGeoEff)
```

To (re)compile
```
cd /dune/app/users/<your_username>/NDEff/DUNE_ND_GeoEff/
#
# In case you log out, need to source setup.sh to setup ROOT
#
source setup.sh         

# Compile program
cd app
make runGeoEffFDEvtSim                                                                       
```

To (re)run program,
```
cd ../bin
# Usage: ./runGeoEffFDEvtSim inputFDntuple
./runGeoEffFDEvtSim /dune/app/users/weishi/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/myntuple.root
```
this will produce a root file containing throws and the hadron throw result.

If the source files in src are changed, recompile:

```
source setup.sh
cmake -DPYTHON_EXECUTABLE:FILEPATH=`which python` .
make -j geoEff    
```

Tips on DUNE FNAL machines: if want to run interactively for longer time even when terminal connection is lost, use screen option:

```
screen
# do the enviroment setup, in this case: source setup.sh
nohup ./runGeoEffFDEvtSim >& out_throws_nohup.log &                                                # Check status: jobs -l
# To detach from the screen session, press Ctrl+a (release) and then d to detach the process/screen.
# To resume detached process, use: screen -r
# 10k evts: 6.20pm start, end second day 4:52am, 10hrs32mins
```

## Instruction for calculate FD event efficiency

The output root file from running ```runGeoEffFDEvtSim``` can be used to calculate FD event hadron containment efficiency by running:

```
cd /dune/app/users/weishi/NDEff/DUNE_ND_GeoEff
source setup.sh
cd app
root -l -b -q FDEffCalc.C
# 5k evts: 10mins
```

The output root file from running ```runGeoEffFDEvtSim``` can also be used for lepton NN training.

## Event displays

The hadronic hits can be plotted in 2D event displays at FD via:
```
echo 'gROOT->ProcessLine(".L ReadHadronHitNtuple.cpp"); ReadHadronHitNtuple_FD()'| root -l -b
```

## Run on Grid

First get the work env setup:
```
cd /dune/app/users/weishi
wget https://raw.githubusercontent.com/weishi10141993/NeutrinoPhysics/main/setupNDEff-grid.sh --no-check-certificate
```

Suppose the input FD ntuples are in this directory,
```
/pnfs/dune/scratch/users/weishi/myFDntuples
```
write the list to txt file,
```
ls -d "/pnfs/dune/scratch/users/weishi/myFDntuples"/* | sed "s\/pnfs\root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr\g" > myFDntuples.txt
# it also changes pnfs to xrootd so that worker node can access
```
which is going to be sent to the grid in a tarball below.

Now make the tarball,
```
tar -czvf NDEff.tar.gz setupNDEff-grid.sh myFDntuples.txt
# Check the tarball *.tar.gz is indeed created and open with: tar -xf *.tar.gz
```

Get the running script,
```
wget https://raw.githubusercontent.com/weishi10141993/NeutrinoPhysics/main/run_NDEff_autogrid.sh --no-check-certificate

# this submits N jobs (N = number of input files, so each job runs 1 file)
jobsub_submit -G dune -N 1 --memory=500MB --disk=1GB --expected-lifetime=20m --cpu=1 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE --tar_file_name=dropbox:///dune/app/users/weishi/NDEff.tar.gz --use-cvmfs-dropbox -l '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\"' --append_condor_requirements='(TARGET.HAS_Singularity==true&&TARGET.HAS_CVMFS_dune_opensciencegrid_org==true&&TARGET.HAS_CVMFS_larsoft_opensciencegrid_org==true&&TARGET.CVMFS_dune_opensciencegrid_org_REVISION>=1105&&TARGET.HAS_CVMFS_fifeuser1_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser2_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser3_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser4_opensciencegrid_org==true)' file:///dune/app/users/weishi/run_NDEff_autogrid.sh

```
Reference:
3 files each 100 events:``` -N 3 --memory=500MB --disk=1GB --expected-lifetime=20m --cpu=1```
# DUNE_ND_GeoEff

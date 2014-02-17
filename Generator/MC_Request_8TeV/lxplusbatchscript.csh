#!/bin/tcsh

##TO RUN: bsub -q 8nm -J job1 < lxplusbatchscript.csh

## Lxplus Batch Job Script
set PROJECTDIR="/afs/cern.ch/work/s/santanas/Releases/CMSSW_5_3_14_Generator/src/Test_DijetScouting8TeV_GENSIM/"
set MASS=1000
#Qstar
#set CFGFILE="QstarToJJ_M_"$MASS"_Tune4C_8TeV_pythia8_cff_py_GEN_SIM".py
#set ROOTFILE="QstarToJJ_M_"$MASS"_Tune4C_8TeV_pythia8_cff_py_GEN_SIM".root
#RSToQQbar
#set CFGFILE="RSGravitonToQQbar_kMpl01_M_"$MASS"_Tune4C_8TeV_pythia8_cff_py_GEN_SIM".py
#set ROOTFILE="RSGravitonToQQbar_kMpl01_M_"$MASS"_Tune4C_8TeV_pythia8_cff_py_GEN_SIM".root
#RStoGG
set CFGFILE="RSGravitonToGG_kMpl01_M_"$MASS"_Tune4C_8TeV_pythia8_cff_py_GEN_SIM".py
set ROOTFILE="RSGravitonToGG_kMpl01_M_"$MASS"_Tune4C_8TeV_pythia8_cff_py_GEN_SIM".root

## Working directory on lxbatch machine
set TOP="$PWD"
echo "we are in" $TOP
echo "ls"
ls

## Setup CMSSW area
echo "cd" $PROJECTDIR
cd $PROJECTDIR
echo "ls"
ls
echo "Setting CMSSW..."
eval `scramv1 runtime -csh`
echo "CMSSW area set!"

## Run job from working directory
echo "cd" $TOP
cd $TOP
echo "ls"
ls
echo "time cmsRun" $PROJECTDIR$CFGFILE
time cmsRun $PROJECTDIR$CFGFILE

## Copy output 
echo "ls"
ls
echo "cp $ROOTFILE $PROJECTDIR"
cp $ROOTFILE $PROJECTDIR


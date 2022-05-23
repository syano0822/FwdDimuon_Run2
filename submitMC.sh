#!/bin/bash

MODE=$1
JDL=$2
DATA=$3
MIX=$4
YEAR=$5

function main_func () {

    period=$1
    mode=$2
    jdl=$3
    data=$4
    mix=$5

    if [ ${data} == LHC20f10a ]; then    
	DIR_OUT=/Users/syano_mbp2021/analysis/run2/PWGLF/FwdDimuon_Run2/LHC20f10/$period
	PER_DIR=/Users/syano_mbp2021/analysis/run2/PWGLF/mass_mc/LHC20f10/input/Default
    elif [ ${data} == LHC20f10b ]; then    
	DIR_OUT=/Users/syano_mbp2021/analysis/run2/PWGLF/FwdDimuon_Run2/LHC20f10/$period
	PER_DIR=/Users/syano_mbp2021/analysis/run2/PWGLF/mass_mc/LHC20f10/input/Default
    elif [ ${data} == LHC20f10c ]; then    
	DIR_OUT=/Users/syano_mbp2021/analysis/run2/PWGLF/FwdDimuon_Run2/LHC20f10/$period
	PER_DIR=/Users/syano_mbp2021/analysis/run2/PWGLF/mass_mc/LHC20f10/input/Default
    else
	DIR_OUT=/Users/syano_mbp2021/analysis/run2/PWGLF/FwdDimuon_Run2/${DATA}/$period
	PER_DIR=/Users/syano_mbp2021/analysis/run2/PWGLF/mass_mc/${data}/input/Default
    fi
    
    mkdir -pv $DIR_OUT

    cd $DIR_OUT
    
    ln -s /Users/syano_mbp2021/analysis/run2/PWGLF/FwdDimuon_Run2/* ${DIR_OUT}
    
    aliroot -l -b -q runAnalysisMC.C\(\"${period}\",\"${mode}\",${jdl},\"${data}\",0,${mix}\)        

    if [ ${MODE} == "terminate" ] && [ $JDL == 0 ]; then
	mv Dimuon.root ${period}.root
    fi

    
    
    

    chmod 755 ./${period}.root
    mv ./${period}.root $PER_DIR

    rm -rf ${DIR_OUT}
}


if [ ${YEAR} == 2016 ]; then
    PERIOD=(LHC16h LHC16j LHC16k LHC16l LHC16o LHC16p)
elif [ ${YEAR} == 2017 ]; then
    PERIOD=(LHC17h LHC17i LHC17k LHC17l LHC17m LHC17o LHC17r)
else 
    PERIOD=(LHC18c LHC18d LHC18e LHC18f LHC18l LHC18m LHC18o LHC18p)
fi

for PERIOD in ${PERIOD[@]}
do
    main_func $PERIOD $MODE $JDL $DATA $MIX $YEAR &
done

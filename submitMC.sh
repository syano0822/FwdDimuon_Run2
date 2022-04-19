#!/bin/bash

MODE=$1
JDL=$2
DATA=$3
MIX=$4
YEAR=$5

if [ ${YEAR} == 2016 ]; then
    PERIOD=(LHC16h LHC16j LHC16k LHC16l LHC16o LHC16p)
elif [ ${YEAR} == 2017 ]; then
    PERIOD=(LHC17h LHC17i LHC17k LHC17l LHC17m LHC17o LHC17r)
else 
    PERIOD=(LHC18c LHC18d LHC18e LHC18f LHC18l LHC18m LHC18o LHC18p)
fi

for i in ${PERIOD[@]}
do
    echo ${i}
    aliroot -l -b -q runAnalysisMC.C\(\"${i}\",\"${MODE}\",${JDL},\"${DATA}\",0,${MIX}\)

    if [ ${MODE} == "terminate" ] && [ $JDL == 0 ]; then
	mv Dimuon.root ${i}.root
    fi

done

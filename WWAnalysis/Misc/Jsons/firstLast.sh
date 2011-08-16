#!/usr/bin/env bash

triggers="
HLT_DoubleMu7_v
HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v
HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v
HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v
HLT_Ele52_CaloIdVT_TrkIdT_v
HLT_IsoMu17_eta2p1_v
HLT_IsoMu17_v
HLT_IsoMu20_eta2p1_v
HLT_IsoMu20_v
HLT_Mu13_Mu8_v
HLT_Mu15_v
HLT_Mu17_Ele8_CaloIdL_v
HLT_Mu24_v
HLT_Mu30_v
HLT_Mu30_v
HLT_IsoMu17_v
HLT_Mu8_Ele17_CaloIdL_v
"

#cat Mu15 | tee temp.out | grep Mu15 | awk -F'|' -v t=Mu15 '
for x in $triggers; do 
    echo "beginning to calculate $x"
    lumiCalc2.py -hltpath $x* --nowarning recorded | tee out/$x.out | grep HLT | awk -F'|' -v t=$x '
        BEGIN {i=0;} 
        { 
            if($4 == $7 && !($2 in arr)) {
                arr[i]=$2;++i
            } 
        } 
        END { 
            n=asort(arr); 
            printf("%s*:%d-%d\n",t,arr[1],arr[length(arr)]);
            for(j=1;j<=n;++j) print arr[j];
        }
    '  >> fullSummary 2>&1 &
done    


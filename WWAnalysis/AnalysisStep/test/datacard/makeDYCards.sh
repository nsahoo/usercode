for C in ee mm; do for J in 0 1; do
    echo " ---- DYCard_${C}_${J}j.txt ----"
    echo "#Yield  Err(Data)   Err(MC) - both to be gmM" | tee DYCard_${C}_${J}j.txt
    grep "N.Z->$C, all cuts, ${J}j" ../step3/dyEstimate.cuts*txt | \
        perl -npe 's/.*cuts(\d+).*=\s+(\S+) .*-\s+(\S+).*-\s+(\S+).*/$1  $2 $3 $4/' | \
        tee -a DYCard_${C}_${J}j.txt
    done
done

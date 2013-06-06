#!/bin/bash    

createJob() {
if [ -d $4 ]; then 
    echo "$4 exists"
else
    mkdir $4
fi
cd $4
echo "#!/bin/bash" > $5
echo "cd $PWD" >> $5
echo "export SCRAM_ARCH=slc5_amd64_gcc462" >> $5
echo "eval \`scramv1 runtime -sh\`" >> $5

for X in $(seq $6 $7 $8)
    do
        cp ../`echo "card*m"$X"*"` ./
        if [ $1 == '7TeV' -o $1 == '8TeV' ]; then
            combineCards.py `echo "card_"$2"_m"$X"_"$1"_4mu.txt"` `echo "card_"$2"_m"$X"_"$1"_4e.txt"` `echo "card_"$2"_m"$X"_"$1"_2e2mu.txt"` > `echo "card_"$2"_m"$X"_"$1".txt"`  
        elif [ $1 == '78TeV' ]; then
            combineCards.py `echo "card_"$2"_m"$X"_7TeV_4mu.txt"` `echo "card_"$2"_m"$X"_7TeV_4e.txt"` `echo "card_"$2"_m"$X"_7TeV_2e2mu.txt"` `echo "card_"$2"_m"$X"_8TeV_4mu.txt"` `echo "card_"$2"_m"$X"_8TeV_4e.txt"` `echo "card_"$2"_m"$X"_8TeV_2e2mu.txt"`> `echo "card_"$2"_m"$X"_"$1".txt"`
        fi

        if [ $3 == 'limits' ]; then
            echo "combine -M Asymptotic -m "$X" card_"$2"_m"$X"_"$1".txt > limits_"$2"_m"$X"_"$1".txt" >> $5
        elif [ $3 == 'signif'  ]; then echo "combine -M ProfileLikelihood -m "$X" --signif card_"$2"_m"$X"_"$1".txt > signif_"$2"_m"$X"_"$1".txt" >> $5
        else echo "combine -M ProfileLikelihood --significance -m "$X" card_"$2"_m"$X"_"$1".txt -t -1 --expectSignal=1 --toysFreq > expsignif_"$2"_m"$X"_"$1".txt" >> $5
        fi
        chmod 755 $5
    done

bsub -q 1nh $5
cd ../
}

createJob $1 $2 $3 jobwreg job1.sh   114 1 118
createJob $1 $2 $3 jobwreg job2.sh   119 1 123
createJob $1 $2 $3 jobwreg job3.sh   124 1 128
createJob $1 $2 $3 jobwreg job4.sh   129 1 133
createJob $1 $2 $3 jobwreg job5.sh   134 1 138
createJob $1 $2 $3 jobwreg job6.sh   139 1 143
createJob $1 $2 $3 jobwreg job7.sh   144 1 148
createJob $1 $2 $3 jobwreg job8.sh   149 1 153
createJob $1 $2 $3 jobwreg job9.sh   154 1 160
createJob $1 $2 $3 jobwreg job10.sh  162 2 164


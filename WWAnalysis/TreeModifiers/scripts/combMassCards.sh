#!/bin/bash

dirs='cards/c1DNoMassErr cards/c1DMassErr cards/c2DMassErr'

echo "Combining m(4l) cards..."
for dir in $dirs 
  do
  cd $dir
  pwd=$PWD
  echo $pwd
  suffix=''
  if [[ $dir =~ 'cards/c1DNoMassErr' ]]
      then
      suffix='1D'
  elif [[ $dir =~ 'cards/c1DMassErr' ]]
      then
      suffix='1Dmerr'
  elif [[ $dir =~ 'cards/c2DMassErr' ]]
      then
      suffix='2Dmerr'
  fi
  sed -i 's|card|'$pwd'/card|g' card*txt
  combineCards.py -S hzz4l_7TeV_4mu=card_$suffix\_DCBApprox_m126_7TeV_4mu.txt     hzz4l_8TeV_4mu=card_$suffix\_DCBApprox_m126_8TeV_4mu.txt     > card_$suffix\_DCBApprox_m126_4mu.txt
  combineCards.py -S hzz4l_7TeV_4e=card_$suffix\_DCBApprox_m126_7TeV_4e.txt       hzz4l_8TeV_4e=card_$suffix\_DCBApprox_m126_8TeV_4e.txt       > card_$suffix\_DCBApprox_m126_4e.txt
  combineCards.py -S hzz4l_7TeV_2e2mu=card_$suffix\_DCBApprox_m126_7TeV_2e2mu.txt hzz4l_8TeV_2e2mu=card_$suffix\_DCBApprox_m126_8TeV_2e2mu.txt > card_$suffix\_DCBApprox_m126_2e2mu.txt
  
  combineCards.py -S hzz4l_7TeV_4mu=card_$suffix\_DCBApprox_m126_7TeV_4mu.txt hzz4l_7TeV_4e=card_$suffix\_DCBApprox_m126_7TeV_4e.txt hzz4l_7TeV_2e2mu=card_$suffix\_DCBApprox_m126_7TeV_2e2mu.txt \
      hzz4l_8TeV_4mu=card_$suffix\_DCBApprox_m126_8TeV_4mu.txt hzz4l_8TeV_4e=card_$suffix\_DCBApprox_m126_8TeV_4e.txt hzz4l_8TeV_2e2mu=card_$suffix\_DCBApprox_m126_8TeV_2e2mu.txt \
      > card_$suffix\_DCBApprox_m126_comb.txt
  
  text2workspace.py card_$suffix\_DCBApprox_m126_4mu.txt    -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass --PO higgsMassRange=122,132 -o FloatMass_4mu_hzz.root
  text2workspace.py card_$suffix\_DCBApprox_m126_4e.txt     -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass --PO higgsMassRange=122,132 -o FloatMass_4e_hzz.root
  text2workspace.py card_$suffix\_DCBApprox_m126_2e2mu.txt  -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass --PO higgsMassRange=122,132 -o FloatMass_2e2mu_hzz.root
  
  text2workspace.py card_$suffix\_DCBApprox_m126_comb.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass --PO higgsMassRange=122,132 -o FloatMass_comb_hzz.root
  cd -
done;

echo "done"

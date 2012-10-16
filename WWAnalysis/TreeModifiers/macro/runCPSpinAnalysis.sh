

############################
#Make Some Plots
############################
root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/MakeMELAPlots.C+



############################
#Make MELA 2D shape
############################
root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/create2DHistogram.C+

root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/createPSMelaVsMassHistogram.C+
root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/createSpinTwoMinimalMelaVsMassHistogram.C+



############################
#Calculate Bkg Model parameters
############################
root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/getParams.C+



############################
#Make Data Cards
############################
root -l WWAnalysis/TreeModifiers/macro/rootlogon.C WWAnalysis/TreeModifiers/macro/makeSpinCPMeasurementDatacards.C+'(1)'
combineCards.py card_pseudoscalar_2D_m125_8TeV_4mu.txt card_pseudoscalar_2D_m125_8TeV_4e.txt card_pseudoscalar_2D_m125_8TeV_2e2mu.txt >! card_pseudoscalar_2D_m125_8TeV_combined.txt

############################
#Do Hypothesis Tests
############################

################
#Make Workspace
################
text2workspace.py -m 125 /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/card_pseudoscalar_2D_m125_8TeV_combined.txt -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs -o workspace_HCP_pseudoscalar_fixedMu.root
text2workspace.py -m 125 /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/card_pseudoscalar_2D_m125_8TeV_4mu.txt -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs -o workspace_4mu_HCP_pseudoscalar_fixedMu.root
text2workspace.py -m 125 /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/card_pseudoscalar_2D_m125_8TeV_comb.txt -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs -o workspace_comb_HCP_pseudoscalar_fixedMu.root

################
#Do Asimov Toy
################
combine -m 125 -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/workspace_HCP_pseudoscalar_fixedMu.root --singlePoint 1 --saveHybridResult --fork 4 -T 4000 -i 1 --clsAcc 0 --fullBToys

combine -m 125 -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 /afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/workspace_HCP_pseudoscalar_fixedMu.root --singlePoint 1 --saveHybridResult -T 4000 -i 1 --clsAcc 0 --fullBToys

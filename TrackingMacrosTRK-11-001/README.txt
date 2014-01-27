1) Initialization:
lunch root and compile doubleCrystaBall function:
- root.exe
- .L rooDoubleCB.cc++
The latter creates the library rooDoubleCB_cc.so


2) Running the macros:
- root.exe
- gSystem->Load("./rooDoubleCB_cc.so");
- gROOT->LoadMacro("./tdrstyle.C");
- setTDRStyle();
- compile and run the macro as: .x macroName.C++(argument)

Purpose of the different macros:

- makePlotsSingleParticles.C:
Used for efficiency and fake rate plots as a function of eta and pt. Used
for muons, pions and electrons.

Used also for resolution plots (as a function of eta and pt) for muons and pions.


- makeElectronResolutionPlots.C:
Used for resolution plots (as a function of eta) for electrons. It produces
resolution plots separately for the following configuration: 
tail-side, gaussian-side and inclusive.

It also produces the bias plots.

- makePlotsTTbar.C:
Used for efficiency and fake rate plots of ttbar, both generalTracks vs highPurity and
highPurityWithPU vs highPurityWithoutPU

- makePlotsTTbarSplit.C:
Used to make resolution vs pt plots for ttbar sample. It splits the resolution in 
barrel,endcap and barrel-endcap transition regions.


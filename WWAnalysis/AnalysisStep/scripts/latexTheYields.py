import sys, os, math, re, copy

from WWAnalysis.AnalysisStep.setupSamples_cff import *

def printFinalYields(filePrefix,masses,suffix,description,cut=-1,lumi=1000.,addData=False):
    for mass in masses:
        f = TFile("data/{0}{1}.root".format(filePrefix,mass))
        outFile = open("tex/yield{0}_{1}_{2:.0f}.tex".format(filePrefix.replace('.',''),mass,lumi),'w')
#         outFile = sys.stdout
        print >> outFile,  "\\begin{table}[h!]\\begin{center}"
        print >> outFile,  "\\caption{{Yields for $m_{{H}}$ = {0} with {3:.2f} $\mathrm{{pb}}^{{-1}}$, {1} \\label{{tab:yield_{2}_{0} }}}}".format(mass,description,filePrefix.replace('.',''),lumi)
        print >> outFile,  "\\begin{tabular}{|c|c|c|c|c|c|} \\hline"

        print >> outFile,  " final state ",
        for label in latexLabels[:5]:
            print >> outFile,  " & ",label,
        print >> outFile,  " \\\\ \\hline"
        for chan in channels:
            print >> outFile,  channelLabels[chan],
            for label in latexLabels[:5]:
                yieldAndError = getYieldAndError(f,latexMassSamples[mass][label],[chan],cut,suffix,lumi)
                print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format(yieldAndError[0],yieldAndError[1]),
            print >> outFile,  "\\\\ " 
        print >> outFile,  "\\hline",
        print >> outFile,  "all     ",
        for label in latexLabels[:5]:
            yieldAndError = getYieldAndError(f,latexMassSamples[mass][label],channels,cut,suffix,lumi)
            print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format( yieldAndError[0],yieldAndError[1] ),
        print >> outFile,  "\\\\ \\hline" 

        print >> outFile,  " final state ",
        for label in latexLabels[5:]:
            print >> outFile,  " & ",label,
        print >> outFile,  " & all bkg & data \\\\ \\hline"
        for chan in channels:
            print >> outFile,  channelLabels[chan],
            for label in latexLabels[5:]:
                yieldAndError = getYieldAndError(f,latexMassSamples[mass][label],[chan],cut,suffix,lumi)
                print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format(yieldAndError[0],yieldAndError[1]),
            yieldAndError = getYieldAndError(f,allBackgroundSamples,[chan],cut,suffix,lumi)
            print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format(yieldAndError[0],yieldAndError[1]),
            if addData:
                yieldAndError = getYieldAndError(f,dataSamples,[chan],cut,suffix,-1)
                print >> outFile,  " & ${0:6.0f}$".format(yieldAndError[0]),
            else:
                print >> outFile, " & N/A ",
            print >> outFile,  "\\\\ " 
        print >> outFile,  "\\hline",
        print >> outFile,  "all     ",
        for label in latexLabels[5:]:
            yieldAndError = getYieldAndError(f,latexMassSamples[mass][label],channels,cut,suffix,lumi)
            print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format( yieldAndError[0],yieldAndError[1] ),
        yieldAndError = getYieldAndError(f,allBackgroundSamples,channels,cut,suffix,lumi)
        print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format(yieldAndError[0],yieldAndError[1]),
        if addData:
            yieldAndError = getYieldAndError(f,dataSamples,channels,cut,suffix,-1)
            print >> outFile,  " & ${0:6.0f}$".format(yieldAndError[0]),
        else:
            print >> outFile, " & N/A ",
        print >> outFile,  "\\\\ \\hline" 
        print >> outFile,  "\\end{tabular}\\end{center}\\end{table}"
        print >> outFile
        print >> outFile
        outFile.close()

def buildNewDict(sampleIds,mass,lumi):
    newDict = {}
    for id in sampleIds:
        modifiedID = id if int(id) < 101 else id+mass
        for sample in latexMassSamples[mass]:
            if lumi < 0:
                if modifiedID in dataSamples:
                    newDict[modifiedID] = copy.deepcopy(dataSamples[modifiedID])
            else:
                if modifiedID in latexMassSamples[mass][sample]:
                    newDict[modifiedID] = copy.deepcopy(latexMassSamples[mass][sample][modifiedID])
    return newDict

def printCutProgression(filePrefix,samplePrefix,channels,masses,suffix,sampleIds,description,lumi = 1000):
    for mass in masses:
        fullSet = cloneVPSet(defaultWW)
        addMassDependentCuts(fullSet,getattr(sys.modules[__name__],"h{0}{1}".format(samplePrefix,mass)))
#         injectTightIsolationCut(fullSet)
        injectIPCut(fullSet)

        newDict = buildNewDict(sampleIds,mass,lumi)
    
        f = TFile("data/{0}{1}.root".format(filePrefix,mass))
        outFile = open("tex/break_{0}_{1}_{2}_{3:.0f}.tex".format(filePrefix.replace('.',''),sampleIds[0],mass,lumi if lumi > 0 else 0),'w')
    
        print >> outFile,  "\\begin{table}[h!]\\begin{center}"
        print >> outFile,  "\\caption{{Breakdown for {0} for $m_{{H}}$ = {1} \\label{{tab:break_{2}_{3}_{4} }}}}".format(description,mass,filePrefix.replace('.',''),sampleIds[0],mass)
        print >> outFile,  "\\begin{tabular}{|c|c|c|c|c|c|} \\hline"
        print >> outFile,  " cut level ",
        for chan in channels:
            print >> outFile,  " & ",channelLabels[chan],
        print >> outFile,  "& all \\\\ \\hline"
        for cut in range(1,len(fullSet)+1):
            print >> outFile,  convertToRealLatex(fullSet[cut-1].label.value()),
            for chan in channels:
                yieldAndError = getYieldAndError(f,newDict,[chan],cut,suffix,lumi) 
                if lumi < 0:
                    print >> outFile,  " & ${0:6.0f}$".format( yieldAndError[0]),
                else:
                    print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format( yieldAndError[0],yieldAndError[1] ),
            yieldAndError = getYieldAndError(f,newDict,channels,cut,suffix,lumi) 
            if lumi < 0:
                print >> outFile,  " & ${0:6.0f}$ \\\\".format( yieldAndError[0]),
            else:
                print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$ \\\\".format( yieldAndError[0],yieldAndError[1] )
        print >> outFile,  "\\hline" 
        print >> outFile,  "\\end{tabular}\\end{center}\\end{table}"
        print >> outFile
        print >> outFile
        outFile.close()
    
def printAllProgression(filePrefix,samplePrefix,channels,masses,suffix,description,lumi = 1000,cutMin=1,cutMax=-1):
    for mass in masses:
        fullSet = cloneVPSet(defaultWW)
        addMassDependentCuts(fullSet,getattr(sys.modules[__name__],"h{0}{1}".format(samplePrefix,mass)))
#         injectTightIsolationCut(fullSet)
        injectIPCut(fullSet)

        f = TFile("data/{0}{1}.root".format(filePrefix,mass))
        for chan in channels:
            outFile = open("tex/breakAll_{0}_{1}_{2}_{3:.0f}.tex".format(filePrefix.replace('.',''),mass,chan,lumi if lumi > 0 else 0),'w')
#             outFile = sys.stdout
    
            print >> outFile,  "\\begin{table}[h!]\\begin{center}"
            print >> outFile,  "\\caption{{Breakdown for {0} for $m_{{H}}$ = {1} in the {3} channel with {5}/pb \\label{{tab:breakAll_{2}_{5}_{4} }}}}".format(
                description,mass,filePrefix.replace('.',''),channelLabels[chan],mass,lumi,chan)
            print >> outFile,  "\\begin{tabular}{|c|c|c|c|c|c|c|} \\hline"
            print >> outFile,  " cut level ",
            for samp in latexLabels[:5]:
                print >> outFile, " & {0}".format(samp),
            print >> outFile, " \\\\ \\hline"
            for cut in range(cutMin,cutMax if cutMax!=-1 else len(fullSet)+1):
                print >> outFile,  convertToRealLatex(fullSet[cut-1].label.value()),
                for samp in latexLabels[:5]:
                    yieldAndError = getYieldAndError(f,latexMassSamples[mass][samp],[chan],cut,suffix,lumi) 
                    print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format( yieldAndError[0],yieldAndError[1] ),
                print >> outFile, " \\\\"
            print >> outFile,  " \\hline \\hline",

            print >> outFile,  " cut level ",
            for samp in latexLabels[5:]:
                print >> outFile, " & {0}".format(samp),
            print >> outFile, " & all bkg & data \\\\ \\hline"
            for cut in range(cutMin,cutMax if cutMax!=-1 else len(fullSet)+1):
                print >> outFile,  convertToRealLatex(fullSet[cut-1].label.value()),
                for samp in latexLabels[5:]:
                    yieldAndError = getYieldAndError(f,latexMassSamples[mass][samp],[chan],cut,suffix,lumi) 
                    print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format( yieldAndError[0],yieldAndError[1] ),
                yieldAndError = getYieldAndError(f,allBackgroundSamples,[chan],cut,suffix,lumi) 
                print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format( yieldAndError[0],yieldAndError[1] ),
                yieldAndError = getYieldAndError(f,dataSamples,[chan],cut,suffix,-1) 
                print >> outFile,  " & ${0:6.0f}$".format( yieldAndError[0]),
                print >> outFile, " \\\\"
            print >> outFile,  "\\hline" 
            print >> outFile,  "\\end{tabular}\\end{center}\\end{table}"
            print >> outFile
            print >> outFile
            outFile.close()

        outFile = open("tex/breakAll_{0}_{1}_{2}_{3:.0f}.tex".format(filePrefix.replace('.',''),mass,"all",lumi if lumi > 0 else 0),'w')
#         outFile = sys.stdout
    
        print >> outFile,  "\\begin{table}[h!]\\begin{center}"
        print >> outFile,  "\\caption{{Breakdown for {0} for $m_{{H}}$ = {1} in all channels with {4}/pb \\label{{tab:breakAll_{2}_{4}_{3} }}}}".format(description,mass,filePrefix.replace('.',''),mass,lumi)
        print >> outFile,  "\\begin{tabular}{|c|c|c|c|c|c|c|} \\hline"
        print >> outFile,  " cut level ",
        for samp in latexLabels[:5]:
            print >> outFile, " & {0}".format(samp),
        print >> outFile, " \\\\ \\hline"
        for cut in range(cutMin,cutMax if cutMax!=-1 else len(fullSet)+1):
            print >> outFile,  convertToRealLatex(fullSet[cut-1].label.value()),
            for samp in latexLabels[:5]:
                yieldAndError = getYieldAndError(f,latexMassSamples[mass][samp],channels,cut,suffix,lumi) 
                if samp == "HWW" and cut < 6:
                    print >> outFile, " & ",
                else:
                    print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format( yieldAndError[0],yieldAndError[1] ),
            print >> outFile, " \\\\"
        print >> outFile,  " \\hline \\hline",

        print >> outFile,  " cut level ",
        for samp in latexLabels[5:]:
            print >> outFile, " & {0}".format(samp),
        print >> outFile, " & all bkg & data \\\\ \\hline"
        for cut in range(cutMin,cutMax if cutMax!=-1 else len(fullSet)+1):
            print >> outFile,  convertToRealLatex(fullSet[cut-1].label.value()),
            for samp in latexLabels[5:]:
                yieldAndError = getYieldAndError(f,latexMassSamples[mass][samp],channels,cut,suffix,lumi) 
                print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format( yieldAndError[0],yieldAndError[1] ),
            yieldAndError = getYieldAndError(f,allBackgroundSamples,channels,cut,suffix,lumi) 
            print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format( yieldAndError[0],yieldAndError[1] ),
            yieldAndError = getYieldAndError(f,dataSamples,channels,cut,suffix,-1) 
            print >> outFile,  " & ${0:6.0f}$".format( yieldAndError[0]),
            print >> outFile, " \\\\"
        print >> outFile,  "\\hline" 
        print >> outFile,  "\\end{tabular}\\end{center}\\end{table}"
        print >> outFile
        print >> outFile
        outFile.close()
    
    

# ____             _         __  __       _       
#|  _ \           (_)       |  \/  |     (_)      
#| |_) | ___  __ _ _ _ __   | \  / | __ _ _ _ __  
#|  _ < / _ \/ _` | | '_ \  | |\/| |/ _` | | '_ \ 
#| |_) |  __/ (_| | | | | | | |  | | (_| | | | | |
#|____/ \___|\__, |_|_| |_| |_|  |_|\__,_|_|_| |_|
#             __/ |                               
#            |___/                                


                     #fn prefix #masses                                 # branch    # description                            cut lumi addData
# printFinalYields(   'note.125.',  [ m for m in allMasses if int(m) < 700 ] , "CONVLHT",  "using $2 \cdot \gamma \cdot m_{R}^{*}$", -1, 125.6, True)
# printFinalYields(   'note.146.',  [ m for m in allMasses if int(m) < 700 ] , "CONVLHT",  "using $2 \cdot \gamma \cdot m_{R}^{*}$", -1, 1000, False)
# printFinalYields(   'note.146.',  [ m for m in allMasses if int(m) < 700 ] , "CONVLHT",  "using $2 \cdot \gamma \cdot m_{R}^{*}$", -1, 146.1, True)

# printAllProgression('note.125.','GammaMR',channels,[ m for m in allMasses if int(m) < 700 ],"CONVLHT","all samples and data", 125.6, 6, -1)
# printAllProgression('note.146.','GammaMR',channels,[ m for m in allMasses if int(m) < 700 ],"CONVLHT","all samples and data", 146.1, 6, -1)

# printCutProgression('note', 'GammaMR', channels, [ m for m in allMasses if int(m) < 700 ] , "CONVLHT", ['101'], "gg to Higgs to $2\\ell 2\\nu$",1000)
# printCutProgression('note', 'GammaMR', channels, [ m for m in allMasses if int(m) < 700 ] , "CONVLHT", ['014'], "WW",1000)
# printCutProgression('note', 'GammaMR', channels, [ m for m in allMasses if int(m) < 700 ] , "CONVLHT", ['023'], "$t\\bar{t}$",1000)
# printCutProgression('note', 'GammaMR', channels, [ m for m in allMasses if int(m) < 700 ] , "CONVLHT", ['026'], "W+Jets",1000)
# printCutProgression('note', 'GammaMR', ['wwelel'], [ m for m in allMasses if int(m) < 700 ] , "CONVLHT", ['003'], "Drell-Yan to ee",1000)
# printCutProgression('note', 'GammaMR', ['wwmumu'], [ m for m in allMasses if int(m) < 700 ] , "CONVLHT", ['004'], "Drell-Yan to $\\mu\\mu$",1000)
# 
# printCutProgression('note', 'GammaMR', channels, [ m for m in allMasses if int(m) < 700 ] , "CONVLHT", ['101'], "gg to Higgs to $2\\ell 2\\nu$",146.1)
# printCutProgression('note', 'GammaMR', channels, [ m for m in allMasses if int(m) < 700 ] , "CONVLHT", ['014'], "WW",146.1)
# printCutProgression('note', 'GammaMR', channels, [ m for m in allMasses if int(m) < 700 ] , "CONVLHT", ['023'], "$t\\bar{t}$",146.1)
# printCutProgression('note', 'GammaMR', channels, [ m for m in allMasses if int(m) < 700 ] , "CONVLHT", ['026'], "W+Jets",146.1)
# printCutProgression('note', 'GammaMR', ['wwelel'], [ m for m in allMasses if int(m) < 700 ] , "CONVLHT", ['003'], "Drell-Yan to ee",146.1)
# printCutProgression('note', 'GammaMR', ['wwmumu'], [ m for m in allMasses if int(m) < 700 ] , "CONVLHT", ['004'], "Drell-Yan to $\\mu\\mu$",146.1)
# 
# printCutProgression('note', 'GammaMR', channels, [ m for m in allMasses if int(m) < 700 ] , "CONVLHT", [k for k in dataSamples.keys()], "Data", -1)


# printFinalYields(   'backNew',  [ m for m in allMasses if int(m) == 150 ] , "CONVLHT",  "using $2 \cdot \gamma \cdot m_{R}^{*}$", -1, 125.6, True)
# printAllProgression('backNew','GammaMR',channels,[ m for m in allMasses if int(m) == 150 ],"CONVLHT","all samples and data", 125.6, 1, -1)

printFinalYields(   'triggerFirst.125.',  [ m for m in allMasses if int(m) == 150 ] , "CONVLHT",  "using $2 \cdot \gamma \cdot m_{R}^{*}$", -1, 125.6, True)
printAllProgression('triggerFirst.125.','GammaMR',channels,[ m for m in allMasses if int(m) == 150 ],"CONVLHT","all samples and data", 125.6, 1, -1)

printFinalYields(   'triggerFirst.NotTriggerFirst.125.',  [ m for m in allMasses if int(m) == 150 ] , "CONVLHT",  "using $2 \cdot \gamma \cdot m_{R}^{*}$", -1, 125.6, True)
printAllProgression('triggerFirst.NotTriggerFirst.125.','GammaMR',channels,[ m for m in allMasses if int(m) == 150 ],"CONVLHT","all samples and data", 125.6, 1, -1)



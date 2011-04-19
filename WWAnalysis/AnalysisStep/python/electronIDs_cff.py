
import FWCore.ParameterSet.Config as cms

ELE_NOCONV =("( (abs(userFloat('convValueMapProd:dist')) > 0.02   ||" +
             "   abs(userFloat('convValueMapProd:dcot')) > 0.02 ) &&"
             "   gsfTrack.trackerExpectedHitsInner.numberOfLostHits == 0 )")


ELE_IP =("( abs(userFloat('dxyPV')) < 0.02 && " +
         "  abs(userFloat('dzPV'))  < 1.0  )" )



# _      _    _  __      ______     ___  _____ 
#| |    | |  | | \ \    / /___ \   / _ \| ____|
#| |    | |__| |  \ \  / /  __) | | (_) | |__  
#| |    |  __  |   \ \/ /  |__ <   \__, |___ \ 
#| |____| |  | |    \  /   ___) |    / / ___) |
#|______|_|  |_|     \/   |____/    /_/ |____/ 
#                                              

ELE_ID_LH_95_2011=("( (  isEB  && numberOfBrems == 0 && electronID('egammaIDLikelihood') > -4.274 ) ||" +           
                   "  (  isEB  && numberOfBrems  > 0 && electronID('egammaIDLikelihood') > -3.773 ) ||" + 
                   "  ((!isEB) && numberOfBrems == 0 && electronID('egammaIDLikelihood') > -5.092 ) ||" + 
                   "  ((!isEB) && numberOfBrems  > 0 && electronID('egammaIDLikelihood') > -2.796 ) )" )

ELE_ISO_LH_95_2011=("(( isEB && (dr03TkSumPt +" +
                    "  max(0,dr03EcalRecHitSumEt - 1.0) + hcalOverEcal()*superCluster.energy() + " +
                    "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.113) ||"+
                    "( isEE && (dr03TkSumPt + hcalOverEcal*superCluster.energy() + " +
                    "  dr03EcalRecHitSumEt +" +
                    "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.109) ) ")

# _      _    _  __      ______     ___   ___  
#| |    | |  | | \ \    / /___ \   / _ \ / _ \ 
#| |    | |__| |  \ \  / /  __) | | (_) | | | |
#| |    |  __  |   \ \/ /  |__ <   \__, | | | |
#| |____| |  | |    \  /   ___) |    / /| |_| |
#|______|_|  |_|     \/   |____/    /_/  \___/ 
#                                              
#

ELE_ID_LH_90_2011=("( (  isEB  && numberOfBrems == 0 && electronID('egammaIDLikelihood') > -1.497 ) ||" +           
                   "  (  isEB  && numberOfBrems  > 0 && electronID('egammaIDLikelihood') > -1.521 ) ||" + 
                   "  ((!isEB) && numberOfBrems == 0 && electronID('egammaIDLikelihood') > -2.571 ) ||" + 
                   "  ((!isEB) && numberOfBrems  > 0 && electronID('egammaIDLikelihood') > -0.657 ) )" )

ELE_ISO_LH_90_2011=("(( isEB && (dr03TkSumPt +" +
                    "  max(0,dr03EcalRecHitSumEt - 1.0) + hcalOverEcal()*superCluster.energy() + " +
                    "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.070) ||"+
                    "( isEE && (dr03TkSumPt + hcalOverEcal*superCluster.energy() + " +
                    "  dr03EcalRecHitSumEt +" +
                    "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.069) ) ")



# _      _    _  __      ______     ___  _____ 
#| |    | |  | | \ \    / /___ \   / _ \| ____|
#| |    | |__| |  \ \  / /  __) | | (_) | |__  
#| |    |  __  |   \ \/ /  |__ <   > _ <|___ \ 
#| |____| |  | |    \  /   ___) | | (_) |___) |
#|______|_|  |_|     \/   |____/   \___/|____/ 
#                                              

ELE_ID_LH_85_2011=("( (  isEB  && numberOfBrems == 0 && electronID('egammaIDLikelihood') >  0.163 ) ||" +           
                   "  (  isEB  && numberOfBrems  > 0 && electronID('egammaIDLikelihood') >  0.065 ) ||" + 
                   "  ((!isEB) && numberOfBrems == 0 && electronID('egammaIDLikelihood') > -0.683 ) ||" + 
                   "  ((!isEB) && numberOfBrems  > 0 && electronID('egammaIDLikelihood') >  1.564 ) )" )

ELE_ISO_LH_85_2011=("(( isEB && (dr03TkSumPt +" +
                    "  max(0,dr03EcalRecHitSumEt - 1.0) + hcalOverEcal()*superCluster.energy() + " +
                    "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.046) ||"+
                    "( isEE && (dr03TkSumPt + hcalOverEcal*superCluster.energy() + " +
                    "  dr03EcalRecHitSumEt +" +
                    "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.046) ) ")


#__      ______ _______ ______ ___  _____   ___   ___  __ __ 
#\ \    / /  _ \__   __|  ____/ _ \| ____| |__ \ / _ \/_ /_ |
# \ \  / /| |_) | | |  | |__ | (_) | |__      ) | | | || || |
#  \ \/ / |  _ <  | |  |  __| \__, |___ \    / /| | | || || |
#   \  /  | |_) | | |  | |      / / ___) |  / /_| |_| || || |
#    \/   |____/  |_|  |_|     /_/ |____/  |____|\___/ |_||_|
#                                                            
                                                            
ELE_ID_CB_95_2011=("((isEB && sigmaIetaIeta < 0.012 &&" +           
                   " deltaPhiSuperClusterTrackAtVtx > -0.800 && deltaPhiSuperClusterTrackAtVtx < 0.800 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.007 && deltaEtaSuperClusterTrackAtVtx < 0.007) || " + 
                   " ( (!isEB) && sigmaIetaIeta < 0.031  &&  " +
                   " deltaPhiSuperClusterTrackAtVtx > -0.700 && deltaPhiSuperClusterTrackAtVtx < 0.700 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.011 && deltaEtaSuperClusterTrackAtVtx < 0.011 )) ")

ELE_ISO_CB_95_2011=("(( isEB && (dr03TkSumPt +" +
                    "  max(0,dr03EcalRecHitSumEt - 1.0) + hcalOverEcal()*superCluster.energy() + " +
                    "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.150) ||"+
                    "( isEE && (dr03TkSumPt + hcalOverEcal*superCluster.energy() + " +
                    "  dr03EcalRecHitSumEt +" +
                    "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.100) ) ")


#__      ______ _______ ______ ___   ___    ___   ___  __ __ 
#\ \    / /  _ \__   __|  ____/ _ \ / _ \  |__ \ / _ \/_ /_ |
# \ \  / /| |_) | | |  | |__ | (_) | | | |    ) | | | || || |
#  \ \/ / |  _ <  | |  |  __| \__, | | | |   / /| | | || || |
#   \  /  | |_) | | |  | |      / /| |_| |  / /_| |_| || || |
#    \/   |____/  |_|  |_|     /_/  \___/  |____|\___/ |_||_|
#

ELE_ID_CB_90_2011=("((isEB && sigmaIetaIeta < 0.01 &&" +           
                   " deltaPhiSuperClusterTrackAtVtx > -0.071 && deltaPhiSuperClusterTrackAtVtx < 0.071 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.007 && deltaEtaSuperClusterTrackAtVtx < 0.007) || " + 
                   " ( (!isEB) && sigmaIetaIeta < 0.031  &&  " +
                   " deltaPhiSuperClusterTrackAtVtx > -0.047 && deltaPhiSuperClusterTrackAtVtx < 0.047 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.011 && deltaEtaSuperClusterTrackAtVtx < 0.011 )) ")

ELE_ISO_CB_90_2011=("(( isEB && (dr03TkSumPt +" +
                    "  max(0,dr03EcalRecHitSumEt - 1.0) + hcalOverEcal()*superCluster.energy() + " +
                    "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.085) ||"+
                    "( isEE && (dr03TkSumPt + hcalOverEcal*superCluster.energy() + " +
                    "  dr03EcalRecHitSumEt +" +
                    "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.051) ) ")


#__      ______ _______ ______ ___  _____   ___   ___  __ __ 
#\ \    / /  _ \__   __|  ____/ _ \| ____| |__ \ / _ \/_ /_ |
# \ \  / /| |_) | | |  | |__ | (_) | |__      ) | | | || || |
#  \ \/ / |  _ <  | |  |  __| > _ <|___ \    / /| | | || || |
#   \  /  | |_) | | |  | |   | (_) |___) |  / /_| |_| || || |
#    \/   |____/  |_|  |_|    \___/|____/  |____|\___/ |_||_|
#                                                            

ELE_ID_CB_85_2011=("((isEB && sigmaIetaIeta < 0.01 &&" +           
                   " deltaPhiSuperClusterTrackAtVtx > -0.039 && deltaPhiSuperClusterTrackAtVtx < 0.039 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.005 && deltaEtaSuperClusterTrackAtVtx < 0.005) || " + 
                   " ( (!isEB) && sigmaIetaIeta < 0.031  &&  " +
                   " deltaPhiSuperClusterTrackAtVtx > -0.028 && deltaPhiSuperClusterTrackAtVtx < 0.028 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.007 && deltaEtaSuperClusterTrackAtVtx < 0.007 )) ")

ELE_ISO_CB_85_2011=("(( isEB && (dr03TkSumPt +" +
                 "  max(0,dr03EcalRecHitSumEt - 1.0) + hcalOverEcal()*superCluster.energy() + " +
                 "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.053) ||"+
                 "( isEE && (dr03TkSumPt + hcalOverEcal*superCluster.energy() + " +
                 "  dr03EcalRecHitSumEt +" +
                 "  dr03HcalTowerSumEt - userFloat('rhoEl')*3.14159265*0.3*0.3)/pt < 0.042) ) ")


#__      ______ _______ ______ ___   ___    ___   ___  __ __ 
#\ \    / /  _ \__   __|  ____/ _ \ / _ \  |__ \ / _ \/_ /_ |
# \ \  / /| |_) | | |  | |__ | (_) | | | |    ) | | | || || |
#  \ \/ / |  _ <  | |  |  __| > _ <| | | |   / /| | | || || |
#   \  /  | |_) | | |  | |   | (_) | |_| |  / /_| |_| || || |
#    \/   |____/  |_|  |_|    \___/ \___/  |____|\___/ |_||_|
#

ELE_ID_CB_80_2011=("((isEB && sigmaIetaIeta < 0.01 &&" +           
                   " deltaPhiSuperClusterTrackAtVtx > -0.027 && deltaPhiSuperClusterTrackAtVtx < 0.027 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.005 && deltaEtaSuperClusterTrackAtVtx < 0.005) || " + 
                   " ( (!isEB) && sigmaIetaIeta < 0.031  &&  " +
                   " deltaPhiSuperClusterTrackAtVtx > -0.021 && deltaPhiSuperClusterTrackAtVtx < 0.021 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.006 && deltaEtaSuperClusterTrackAtVtx < 0.006 )) ")


#  ____  _      _____     _____ _    _ _____ _______ 
# / __ \| |    |  __ \   / ____| |  | |_   _|__   __|
#| |  | | |    | |  | | | (___ | |__| | | |    | |   
#| |  | | |    | |  | |  \___ \|  __  | | |    | |   
#| |__| | |____| |__| |  ____) | |  | |_| |_   | |   
# \____/|______|_____/  |_____/|_|  |_|_____|  |_|   
#                                                    
                                                    

ELE_ID_CB_80_2010=("((isEB && sigmaIetaIeta < 0.01 &&" +           
                   " deltaPhiSuperClusterTrackAtVtx > -0.06 && deltaPhiSuperClusterTrackAtVtx < 0.06 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.004 && deltaEtaSuperClusterTrackAtVtx < 0.004 &&" +
                   " hadronicOverEm < 0.04) ||" +
                   " ( (!isEB) && sigmaIetaIeta < 0.03  &&  " +
                   " deltaPhiSuperClusterTrackAtVtx > -0.03 && deltaPhiSuperClusterTrackAtVtx < 0.03 &&" +
                   " deltaEtaSuperClusterTrackAtVtx > -0.007 && deltaEtaSuperClusterTrackAtVtx < 0.007 " + 
                   " hadronicOverEm < 0.025 )) ")

ELE_ISO_CB_80_2010=("(( isEB && (dr03TkSumPt +" +
                 "  max(0,dr03EcalRecHitSumEt - 1.0) +" +
                 "  dr03HcalTowerSumEt)/pt < 0.1) ||"+
                 "( isEE && (dr03TkSumPt +" +
                 "  dr03EcalRecHitSumEt +" +
                 "  dr03HcalTowerSumEt)/pt < 0.1) ) ")


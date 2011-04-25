#include "PhysicsTools/UtilAlgos/interface/BasicAnalyzer.h"
#include "WWAnalysis/AnalysisStep/interface/EventBitsFiller.h"
#include "PhysicsTools/SelectorUtils/interface/RunLumiSelector.h"

#include<vector>

class CreateEventHists : public edm::BasicAnalyzer {

    public:
        enum cutBits {
            SKIMMED         =  0,
            FIDUCIAL        ,//=  1,
//             ID              ,//=  2,
//             ISO             ,//=  3,
//             CONVERSION      ,//=  4,
//             IP              ,//=  5,
            LEPTONVETO      ,//=  2,
            MET             ,//=  3,
            MLL             ,//=  4,
            MZ              ,//=  5,
            PROJMET         ,//=  6,
            JETVETO         ,//=  7,
            SOFTMU          ,//=  8,
            BJETS           ,//=  9,
            MLL2            ,//= 10,
            PTMAX           ,//= 11,
            PTMIN           ,//= 12,
            DPHI             //= 13
        };                    

        CreateEventHists(const edm::ParameterSet& cfg, TFileDirectory& fs);
        virtual ~CreateEventHists(){};

        void beginJob(){};
        void endJob();
        void analyze(const edm::EventBase& event);

    private:
        EventBitsFiller eventFiller_;
        edm::ParameterSet selectionParams_;
        RunLumiSelector myLumiSel_;
        std::vector<std::string> hypoTypes_;

};







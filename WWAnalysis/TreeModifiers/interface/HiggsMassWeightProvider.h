#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;


class HiggsMassWeightProvider {

    private :

        std::vector<std::vector<float> > weights;
        bool intOnly;
        bool applyWeights;

        const std::string trim(const std::string& pString, const std::string& pWhitespace = " \t") {
            const size_t beginStr = pString.find_first_not_of(pWhitespace);
            if (beginStr == std::string::npos) return "";
        
            const size_t endStr = pString.find_last_not_of(pWhitespace);
            const size_t range = endStr - beginStr + 1;
        
            return pString.substr(beginStr, range);
        }
        
        const std::string reduce(const std::string& pString, const std::string& pFill = " ", const std::string& pWhitespace = " \t") {
            std::string result(trim(pString, pWhitespace));
            size_t beginSpace = result.find_first_of(pWhitespace);
            while (beginSpace != std::string::npos) {
                const size_t endSpace = result.find_first_not_of(pWhitespace, beginSpace);
                const size_t range = endSpace - beginSpace;
                result.replace(beginSpace, range, pFill);
                const size_t newStart = beginSpace + pFill.length();
                beginSpace = result.find_first_of(pWhitespace, newStart);
            }
            return result;
        }


    public : 
        
        HiggsMassWeightProvider(std::string filename, bool iOnly) : intOnly(iOnly) {
            string line;
            string word;
            if (filename != "") {
                ifstream ifs(filename.c_str());
                while (!ifs.eof()) {
                    getline(ifs, line);
                    line = trim(line);
                    line = reduce(line);
                    istringstream iss(line);
                    std::vector<float> weightrow;
                    while(!iss.eof() && line.size() > 0) {
                        getline(iss, word, ' ');
                        float f;
                        istringstream wss(word);
                        wss >> f;
                        weightrow.push_back(f);
                    }
                    if (weightrow.size() != 0) weights.push_back(weightrow);
                }
                if (weights.size() != 0) {
                    for (std::size_t i = 1; i < weights[0].size(); i++) {
                        float norm = 0.0;
                        for (std::size_t j = 0; j < weights.size(); j++) {
                            norm += weights[j][i];
                        }
                        for (std::size_t j = 0; j < weights.size(); j++) {
                            weights[j][i] /= norm;
                        }
                    }
                }
            }


        }

        float getWeight(float mass) {
        
            if (mass >= 2395) return 0.0;
            if (weights.size() == 0) return 0.0;
            
            int bin = -1;
            
            for (std::size_t i = 0; i < weights.size()-1; i++) {
                if (weights[i].size() == 8 && mass >= weights[i][0] && mass < weights[i+1][0]) {
                    if (mass*2.0 >= (weights[i][0] + weights[i+1][0])) bin = i+1;
                    else bin = i;
                }
            }
            
            if (bin == -1) return 0.0;
            else {
                if (weights[bin][1] == 0) return 0.0;
                else {
                    if      ( intOnly && weights[bin][2] != 0) return weights[bin][3] / weights[bin][2];
                    else if (!intOnly && weights[bin][1] != 0) return weights[bin][3] / weights[bin][1];
                    else return 0;
                }
            }
        }


        float getWeightUp(float mass) {

            if (mass >= 2395) return 0.0;
            if (weights.size() == 0) return 0.0;

            int bin = -1;

            for (std::size_t i = 0; i < weights.size()-1; i++) {
                if (weights[i].size() == 8 && mass >= weights[i][0] && mass < weights[i+1][0]) {
                    if (mass*2.0 >= (weights[i][0] + weights[i+1][0])) bin = i+1;
                    else bin = i;
                }
            }

            if (bin == -1) return 0.0;
            else {
                if (weights[bin][1] == 0) return 0.0;
                else {
                    if      ( intOnly && weights[bin][2] != 0) return weights[bin][6] / weights[bin][2];
                    else if (!intOnly && weights[bin][1] != 0) return weights[bin][6] / weights[bin][1];
                    else return 0;
                }
            }
        }

        float getWeightDown(float mass) {

            if (mass >= 2395) return 0.0;
            if (weights.size() == 0) return 0.0;

            int bin = -1;

            for (std::size_t i = 0; i < weights.size()-1; i++) {
                if (weights[i].size() == 8 && mass >= weights[i][0] && mass < weights[i+1][0]) {
                    if (mass*2.0 >= (weights[i][0] + weights[i+1][0])) bin = i+1;
                    else bin = i;
                }
            }

            if (bin == -1) return 0.0;
            else {
                if (weights[bin][1] == 0) return 0.0;
                else {
                    if      ( intOnly && weights[bin][2] != 0) return weights[bin][7] / weights[bin][2];
                    else if (!intOnly && weights[bin][1] != 0) return weights[bin][7] / weights[bin][1];
                    else return 0;
                }
            }
        }

};

#include <iostream>
#include <vector>
#include <limits>
void TestDummy() {
    const int nDetPos = 6;
    const int nDetSlices = 5;
    double nDetLowestSlice = -2;
    double nDetHighestSlice = 2;
    double detPos[nDetPos] = {0, 1, 4, 5, 8, 13};
    double OAPos[nDetPos][nDetSlices];
    vector<double> LowBinEdge;
    vector<double> HighBinEdge;

    // Populate the OAPos array
    for(int iDetPos = 0; iDetPos < nDetPos; iDetPos++) {
        for(int iSlice = 0; iSlice < nDetSlices; iSlice++) {
            OAPos[iDetPos][iSlice] = detPos[iDetPos] + nDetLowestSlice + iSlice;
        }
    }

    // //fisrt LowBinEdge will always be the first value of the first detPos:
    // LowBinEdge.push_back(OAPos[0][0]);

    for (int iDetPos = 0; iDetPos < nDetPos; iDetPos++) {
        for (int iCountSlice = 0; iCountSlice < nDetSlices; iCountSlice++) {
           std::cout << OAPos[iDetPos][iCountSlice] << " ";
        }
        cout<<"  "<<endl;
    }

    // different algorithm for first bin (i.e first 2 detector positions)
    // Boolean array to track covered slice indices
    bool covered[nDetSlices] = {false};
    // Variables to track the indices of minimum iSlice1 and maximum iSlice2
    int minISlice1 = nDetSlices; // Initialize to an invalid index
    int maxISlice2 = -1;         // Initialize to an invalid index
    // Iterate through possible pairs to check for equal OAPos values
    for (int iSlice1 = 0; iSlice1 < nDetSlices; iSlice1++) {
        for (int iSlice2 = 0; iSlice2 < nDetSlices; iSlice2++) {
            if (OAPos[0][iSlice1] == OAPos[1][iSlice2]) {
                std::cout << "Equal overlap between bin 0 and bin 1 at bin0 det Slice = "
                          << iSlice1 << ", bin1 det Slice = " << iSlice2 << std::endl;
                // Mark the corresponding indices as covered
                covered[iSlice1] = true;
                covered[iSlice2] = true;
                // Update minISlice1 and maxISlice2
                if (iSlice1 < minISlice1) {
                    minISlice1 = iSlice1;
                }
                if (iSlice2 > maxISlice2) {
                    maxISlice2 = iSlice2;
                }
            }
        }
    }
    // Check if all indices are covered
    bool allCovered = true;
    for (int i = 0; i < nDetSlices; i++) {
        if (!covered[i]) {
            allCovered = false;
            break;
        }
    }

    int minSliceIndexNotCovered = nDetSlices; // Initialize to an invalid index
    if (allCovered) {
        std::cout << "All slice indices are covered." << std::endl;
        std::cout << "Value corresponding to minimum iSlice1 in OAPos[0][iSlice1]: "
                  << OAPos[0][minISlice1] << std::endl;
        std::cout << "Value corresponding to maximum iSlice2 in OAPos[1][iSlice2]: "
                  << OAPos[1][maxISlice2] << std::endl;
        LowBinEdge.push_back(OAPos[0][minISlice1]);
        HighBinEdge.push_back(OAPos[1][maxISlice2]);
    } else {
        std::cout << "The following slice indices are not covered:" << std::endl;
        for (int i = 0; i < nDetSlices; i++) {
            if (!covered[i]) {
                std::cout << "Slice index " << i << " (OAPos[0][" << i << "] = "
                          << OAPos[0][i] << ")" << std::endl;
                if(i < minSliceIndexNotCovered){
                  minSliceIndexNotCovered = i;
                }
                // cout<<" mis index: "<<minSliceIndexNotCovered<<" val at missing index: "<< OAPos[0][minSliceIndexNotCovered]<<endl;
            }
        }
        cout<<" mis index: "<<minSliceIndexNotCovered<<" val at missing index: "<< OAPos[0][minSliceIndexNotCovered]<<endl;
        LowBinEdge.push_back(OAPos[0][minSliceIndexNotCovered]);
        HighBinEdge.push_back(OAPos[0][nDetSlices-1]); //always end at last det slice in 0m axis if all indexes are not covered by both (in fact i think we always end at this last det slice in the case of the first 2 bins)
    }

    for(int i = 0; i< LowBinEdge.size(); i++){
      cout<< " bin for the first 2 detector postions at "<<detPos[0]<<" m and "<<detPos[1]<<" m are from: "<<LowBinEdge[i]<<" to "<<HighBinEdge[i]<<endl;
    }



    double currentThreshold = -std::numeric_limits<double>::infinity();
    bool allColumnsCovered = true;
    do {
        bool columnsCovered[nDetSlices] = {false};
        std::vector<std::pair<double, std::pair<int, int>>> minElements;
        double maxValueInIteration = -std::numeric_limits<double>::infinity();
        // Find the minimum elements covering all columns with the current threshold
        while (std::any_of(std::begin(columnsCovered), std::end(columnsCovered), [](bool covered) { return !covered; })) {
            double minValue = std::numeric_limits<double>::infinity();
            int minRow = -1;
            int minCol = -1;
            for (int i = 0; i < nDetPos; ++i) {
                for (int j = 0; j < nDetSlices; ++j) {
                    if (!columnsCovered[j] && OAPos[i][j] > currentThreshold && OAPos[i][j] < minValue) {
                        minValue = OAPos[i][j];
                        minRow = i;
                        minCol = j;
                    }
                }
            }
            if (minCol != -1) {
                columnsCovered[minCol] = true;
                minElements.push_back({minValue, {minRow, minCol}});
                if (minValue > maxValueInIteration) {
                    maxValueInIteration = minValue;
                }
            } else {
                allColumnsCovered = false;
                break;
            }
        }
        if (allColumnsCovered) {
            // Output the results for the current iteration
            std::cout << "Selected minimum elements covering all columns with threshold > " << currentThreshold << ":\n";
            for (const auto& elem : minElements) {
                std::cout << "Value: " << elem.first << " at position [" << elem.second.first << "][" << elem.second.second << "]\n";
            }
            auto minElemIt = std::min_element(minElements.begin(), minElements.end(),
            [](const auto& lhs, const auto& rhs) {
                return lhs.first < rhs.first;
            });
            // Access the minimum element
           const auto& minElem = *minElemIt;

           auto maxElemIt = std::max_element(minElements.begin(), minElements.end(),
           [](const auto& lhs, const auto& rhs) {
               return lhs.first < rhs.first;
           });
           const auto& maxElem = *maxElemIt;
           // Output the result
           if(currentThreshold!=(-std::numeric_limits<double>::infinity()) ){
             std::cout << "Minimum Value: " << minElem.first <<" max value: "<<maxElem.first
                       <<"\n";
             LowBinEdge.push_back(minElem.first);
             HighBinEdge.push_back(maxElem.first);
           }


            std::cout << std::endl;
            // Update the threshold for the next iteration
            currentThreshold = maxValueInIteration;
        }
    } while (allColumnsCovered);


    // //last HighBinEdge will always be the last value of the last detPos:
    // HighBinEdge.push_back(OAPos[nDetPos-1][nDetSlices-1]);
    for(int i = 0; i< LowBinEdge.size(); i++){
      cout<< " bins from: "<<LowBinEdge[i]<<" to "<<HighBinEdge[i]<<endl;
    }
}

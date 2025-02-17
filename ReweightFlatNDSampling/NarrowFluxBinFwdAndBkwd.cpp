#include <vector>
#include <set>
#include <iostream>
#include <limits>
#include <algorithm>
#include <TGraph.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TRootCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <cstdio>

void FunctionToNarrowBin(bool forward, int nDetPos, int nDetSlices, const std::vector<double>& detPos, std::vector<double>& LowBinEdge, std::vector<double>& HighBinEdge) {
    // Initialize the 2D array of OA positions (OAPos)
    std::vector<std::vector<double>> OAPos(nDetPos, std::vector<double>(nDetSlices));
    for (int i = 0; i < nDetPos; ++i) {
        for (int j = 0; j < nDetSlices; ++j) {
            // Assuming you have a way to calculate OAPos from detPos
            OAPos[i][j] = detPos[i] + (-2) + 0.5*j ;

        }
    }

    std::set<double> processedValues;
  // std::vector<double> HighBinEdgeForward;
  // std::vector<double> LowBinEdgeForward;
  // std::vector<double> HighBinEdgeBackward;
  // std::vector<double> LowBinEdgeBackward;

  if (forward) {
      // Forward binning logic
      for (int startDetPos = 0; startDetPos < nDetPos; ++startDetPos) {
          for (int startDetSlice = 0; startDetSlice < nDetSlices; ++startDetSlice) {
              double currentValue = OAPos[startDetPos][startDetSlice];

              // Skip if the value has already been processed
              if (processedValues.find(currentValue) != processedValues.end()) {
                  continue;
              }

              // Use std::vector<bool> to track covered columns
              std::vector<bool> columnsCovered(nDetSlices, false);
              std::vector<std::pair<double, std::pair<int, int>>> minElements;

              // Consider the starting point as covered
              columnsCovered[startDetSlice] = true;
              minElements.push_back({OAPos[startDetPos][startDetSlice], {startDetPos, startDetSlice}});

              // Find the remaining minimum elements to cover other columns
              while (std::any_of(std::begin(columnsCovered), std::end(columnsCovered), [](bool covered) { return !covered; })) {
                  double minValue = std::numeric_limits<double>::infinity();
                  int minRow = -1;
                  int minCol = -1;

                  // Only look at elements that are after the current starting point
                  for (int i = startDetPos; i < nDetPos; ++i) {
                      for (int j = 0; j < nDetSlices; ++j) {

                          // Skip columns before the starting point
                          if (i == startDetPos && j <= startDetSlice) continue;

                          // Get the minimum values that cover all j indexes
                          if (!columnsCovered[j] && OAPos[i][j] < minValue && OAPos[i][j] >= currentValue) {
                              minValue = OAPos[i][j];
                              minRow = i;
                              minCol = j;
                          }
                      }
                  }

                  if (minRow == -1 || minCol == -1) {
                      break;
                  }

                  // Mark the column as covered and store the element
                  if (minCol != -1) {
                      columnsCovered[minCol] = true;
                      minElements.push_back({minValue, {minRow, minCol}});
                  }
              }

              // Only add the result to minElements if all columns are covered
              if (std::all_of(std::begin(columnsCovered), std::end(columnsCovered), [](bool covered) { return covered; })) {
                  auto maxElemIt = std::max_element(minElements.begin(), minElements.end(),
                  [](const auto& lhs, const auto& rhs) {
                      return lhs.first < rhs.first;
                  });
                  const auto& maxElem = *maxElemIt;
                  HighBinEdge.push_back(maxElem.first);
                  LowBinEdge.push_back(OAPos[startDetPos][startDetSlice]);

                  //std::cout << "Bin from " << OAPos[startDetPos][startDetSlice] << " to " << maxElem.first << std::endl;
              }

              // Add the current value to processedValues to prevent reprocessing
              processedValues.insert(currentValue);
          }
      }
  } else {
      // Backward binning logic
      for (int startDetPos = nDetPos - 1; startDetPos >= 0; --startDetPos) {
          for (int startDetSlice = nDetSlices - 1; startDetSlice >= 0; --startDetSlice) {
              double currentValue = OAPos[startDetPos][startDetSlice];

              // Skip if the value has already been processed
              if (processedValues.find(currentValue) != processedValues.end()) {
                  continue;
              }

              // Use std::vector<bool> to track covered columns
              std::vector<bool> columnsCovered(nDetSlices, false);
              std::vector<std::pair<double, std::pair<int, int>>> minElements;

              // Consider the starting point as covered
              columnsCovered[startDetSlice] = true;
              minElements.push_back({OAPos[startDetPos][startDetSlice], {startDetPos, startDetSlice}});

              // Find the remaining minimum elements to cover other columns
              while (std::any_of(std::begin(columnsCovered), std::end(columnsCovered), [](bool covered) { return !covered; })) {
                  double maxValue = -std::numeric_limits<double>::infinity();
                  int minRow = -1;
                  int minCol = -1;

                  // Only look at elements that are before the current starting point (since we're going backward)
                  for (int i = startDetPos; i >= 0; --i) {
                      for (int j = nDetSlices - 1; j >= 0; --j) {
                          if (i == startDetPos && j >= startDetSlice) continue;  // Skip columns before the starting point

                          if (!columnsCovered[j] && OAPos[i][j] > maxValue && OAPos[i][j] <= currentValue) {
                              maxValue = OAPos[i][j];
                              minRow = i;
                              minCol = j;
                          }
                      }
                  }

                  if (minRow == -1 || minCol == -1) {
                      break;
                  }

                  // Mark the column as covered and store the element
                  if (minCol != -1) {
                      columnsCovered[minCol] = true;
                      minElements.push_back({maxValue, {minRow, minCol}});
                  }
              }

              // Only add the result to minElements if all columns are covered
              if (std::all_of(std::begin(columnsCovered), std::end(columnsCovered), [](bool covered) { return covered; })) {
                  auto minElemIt = std::min_element(minElements.begin(), minElements.end(),
                  [](const auto& lhs, const auto& rhs) {
                      return lhs.first < rhs.first;
                  });
                  const auto& minElem = *minElemIt;
                  LowBinEdge.push_back(minElem.first);
                  HighBinEdge.push_back(OAPos[startDetPos][startDetSlice]);

                  //std::cout << "BACKWARDS Bin from " << minElem.first << " to " << OAPos[startDetPos][startDetSlice] << std::endl;
              }

              // Add the current value to processedValues to prevent reprocessing
              processedValues.insert(currentValue);
          }
      }
  }
}

std::pair<std::vector<double>, std::vector<double>> ChosenBins(
    const std::vector<double>& LowBinEdgeBackward,
    const std::vector<double>& HighBinEdgeBackward,
    const std::vector<double>& LowBinEdgeForward,
    const std::vector<double>& HighBinEdgeForward
) {
    // Vectors to store the equal bin edges
    std::vector<double> EqualLowBinEdge;
    std::vector<double> EqualHighBinEdge;

    // Sort the input vectors to ensure proper comparison
    std::vector<double> SortedLowBinEdgeBackward = LowBinEdgeBackward;
    std::vector<double> SortedHighBinEdgeBackward = HighBinEdgeBackward;
    std::vector<double> SortedLowBinEdgeForward = LowBinEdgeForward;
    std::vector<double> SortedHighBinEdgeForward = HighBinEdgeForward;

    std::sort(SortedLowBinEdgeBackward.begin(), SortedLowBinEdgeBackward.end());
    std::sort(SortedHighBinEdgeBackward.begin(), SortedHighBinEdgeBackward.end());
    std::sort(SortedLowBinEdgeForward.begin(), SortedLowBinEdgeForward.end());
    std::sort(SortedHighBinEdgeForward.begin(), SortedHighBinEdgeForward.end());

    // Save equal edges (where low and high edges match)
    for (int i = 0; i < SortedLowBinEdgeBackward.size(); i++) {
        for (int j = 0; j < SortedLowBinEdgeForward.size(); j++) {
            if (SortedLowBinEdgeForward[j] == SortedLowBinEdgeBackward[i] &&
                SortedHighBinEdgeForward[j] == SortedHighBinEdgeBackward[i]) {
                EqualLowBinEdge.push_back(SortedLowBinEdgeBackward[i]);
                EqualHighBinEdge.push_back(SortedHighBinEdgeForward[j]);
            }
        }
    }

    // Vectors to store the remaining values (bin edges that are individual to forward or backward)
    std::vector<double> RemainingLowBinEge;
    std::vector<double> RemainingHighBinEge;

    // Add remaining low edges from LowBinEdgeForward that were not found in EqualLowBinEdge
    for (size_t i = 0; i < SortedLowBinEdgeForward.size(); ++i) {
        if (std::find(EqualLowBinEdge.begin(), EqualLowBinEdge.end(), SortedLowBinEdgeForward[i]) == EqualLowBinEdge.end()) {
            RemainingLowBinEge.push_back(SortedLowBinEdgeForward[i]);
            RemainingHighBinEge.push_back(SortedHighBinEdgeForward[i]);
        }
    }

    // Add remaining high edges from HighBinEdgeBackward that were not found in EqualHighBinEdge
    for (size_t i = 0; i < SortedHighBinEdgeBackward.size(); ++i) {
        if (std::find(EqualHighBinEdge.begin(), EqualHighBinEdge.end(), SortedHighBinEdgeBackward[i]) == EqualHighBinEdge.end()) {
            RemainingHighBinEge.push_back(SortedHighBinEdgeBackward[i]);
            RemainingLowBinEge.push_back(SortedLowBinEdgeBackward[i]);
        }
    }

    // Vectors to store the final chosen low and high bin edges
    std::vector<double> ChosenLowBinEdge;
    std::vector<double> ChosenHighBinEdge;

    // Insert remaining values into ChosenLowBinEdge and ChosenHighBinEdge
    ChosenLowBinEdge.insert(ChosenLowBinEdge.end(), RemainingLowBinEge.begin(), RemainingLowBinEge.end());
    ChosenLowBinEdge.insert(ChosenLowBinEdge.end(), EqualLowBinEdge.begin(), EqualLowBinEdge.end());

    ChosenHighBinEdge.insert(ChosenHighBinEdge.end(), RemainingHighBinEge.begin(), RemainingHighBinEge.end());
    ChosenHighBinEdge.insert(ChosenHighBinEdge.end(), EqualHighBinEdge.begin(), EqualHighBinEdge.end());

    // Return the two vectors as a pair
    return {ChosenLowBinEdge, ChosenHighBinEdge};
}


int main(int argc, char **argv) {

  using namespace std;

  //for interactive canvas root
    TApplication theApp("App", &argc, argv);
    int nDetPos = 2;
    int nDetSlices = 9;
    std::vector<double> detPos = {0,2};

    double OAPos[nDetPos][nDetSlices];
    double nDetLowestSlice = -2.0;
    // Populate the OAPos array
    for (int iDetPos = 0; iDetPos < nDetPos; iDetPos++) {
        for (int iSlice = 0; iSlice < nDetSlices; iSlice++) {
            OAPos[iDetPos][iSlice] = detPos[iDetPos] + nDetLowestSlice + 0.5*iSlice;
        }
    }

    for (int iDetPos = 0; iDetPos < nDetPos; iDetPos++) {
        for (int iCountSlice = 0; iCountSlice < nDetSlices; iCountSlice++) {
           std::cout << OAPos[iDetPos][iCountSlice] << " ";
        }
        cout<<"  "<<endl;
    }

    vector<double> uniqueOAPos; //need for diagonal
    for (int iDetPos = 0; iDetPos < nDetPos; iDetPos++) {
        for (int iSlice = 0; iSlice < nDetSlices; iSlice++) {
            int element = OAPos[iDetPos][iSlice];
            if(std::find(uniqueOAPos.begin(), uniqueOAPos.end(), element) == uniqueOAPos.end()){
              uniqueOAPos.push_back(element);
            }
        }
    }

    // Vectors to store the Low and High bin edges
    std::vector<double> LowBinEdgeForward;
    std::vector<double> HighBinEdgeForward;

    std::vector<double> LowBinEdgeBackward;
    std::vector<double> HighBinEdgeBackward;


    FunctionToNarrowBin(true, nDetPos, nDetSlices, detPos, LowBinEdgeForward, HighBinEdgeForward);  // For forward binning
    FunctionToNarrowBin(false, nDetPos, nDetSlices, detPos, LowBinEdgeBackward, HighBinEdgeBackward); // For backward binning

    // Optionally, print the results
    std::cout << "Forward Binning:" << std::endl;
    for (size_t i = 0; i < LowBinEdgeForward.size(); ++i) {
        std::cout << "Low: " << LowBinEdgeForward[i] << ", High: " << HighBinEdgeForward[i] << std::endl;
    }

    std::cout << "Backward Binning:" << std::endl;
    for (size_t i = 0; i < LowBinEdgeBackward.size(); ++i) {
        std::cout << "Low: " << LowBinEdgeBackward[i] << ", High: " << HighBinEdgeBackward[i] << std::endl;
    }


    // Call the ChosenBins function
     auto finalbins = ChosenBins(LowBinEdgeBackward, HighBinEdgeBackward, LowBinEdgeForward, HighBinEdgeForward);

    // Output the results
    std::cout << "Chosen Low Bin Edges: ";
    for (double val : finalbins.first) {
       std::cout << val << " ";
    }
    std::cout << std::endl;

    std::cout << "Chosen High Bin Edges: ";
    for (double val : finalbins.second) {
       std::cout << val << " ";
    }
    std::cout << std::endl;


    TCanvas* c1 = new TCanvas("c1", "Graph Example", 800, 600);

    TGraph* LowVsHighEdgeGraphForward = new TGraph(LowBinEdgeForward.size(), &LowBinEdgeForward[0], &HighBinEdgeForward[0]);
    LowVsHighEdgeGraphForward->SetMarkerStyle(20);
    LowVsHighEdgeGraphForward->SetMarkerSize(2.0);
    LowVsHighEdgeGraphForward->SetMarkerColor(kRed-8);
    LowVsHighEdgeGraphForward->SetMinimum(-3);
    LowVsHighEdgeGraphForward->SetMaximum(34);
    LowVsHighEdgeGraphForward->GetXaxis()->SetLimits(-3, 34);
    LowVsHighEdgeGraphForward->Draw("AP");
    TGraph* LowVsHighEdgeGraphBackward = new TGraph(LowBinEdgeBackward.size(), &LowBinEdgeBackward[0], &HighBinEdgeBackward[0]);
    LowVsHighEdgeGraphBackward->SetMarkerStyle(20);
    LowVsHighEdgeGraphBackward->SetMarkerColor(4);
    LowVsHighEdgeGraphBackward->SetMinimum(-3);
    LowVsHighEdgeGraphBackward->Draw("Psame");
    TGraph* ChosenLowAndHigBins = new TGraph(finalbins.first.size(), &finalbins.first[0], &finalbins.second[0]);
    ChosenLowAndHigBins->SetMarkerStyle(4);
    ChosenLowAndHigBins->SetMarkerColor(kGreen+2);
    ChosenLowAndHigBins->SetMarkerSize(3);
    ChosenLowAndHigBins->Draw("Psame");
    TGraph* DiagonalGraph = new TGraph(LowBinEdgeForward.size(), &uniqueOAPos[0], &uniqueOAPos[0]);
    DiagonalGraph->SetLineColor(1);
    DiagonalGraph->SetLineWidth(2);
    DiagonalGraph->Draw("Lsame");

    TLegend* leg = new TLegend();
    leg->AddEntry(LowVsHighEdgeGraphForward, "FORWARD", "p");
    leg->AddEntry(LowVsHighEdgeGraphBackward, "BACKWARD", "p");
    leg->AddEntry(ChosenLowAndHigBins, Form("final bins = %d" , finalbins.second.size()), "p");
    leg->Draw("same");

    TLatex latex;
    latex.SetTextSize(0.04);
    char buffer[200];

    for(int i = 0; i< nDetPos; i++){
      if(i==0)
        sprintf(buffer,  "det stops: ");
      sprintf(buffer + strlen(buffer), " %.2f m,", detPos[i]);

    }

    latex.DrawLatex(0.1, 1.0 , buffer);
    c1->Update();




    //connect canvas close to terminate TApplication
    TRootCanvas *rc = (TRootCanvas *)c1->GetCanvasImp();
    rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    theApp.Run();




    //
    // for(int i = 0; i<LowBinEdgeBackward.size(); i++){
    //   for(int j = 0; j<LowBinEdgeForward.size(); j++){
    //     if(LowBinEdgeBackward[i] == LowBinEdgeForward[j] &&  HighBinEdgeBackward[i]== HighBinEdgeForward[j]){
    //       cout<<" points with same bins , only save once=== same low bin : "<<LowBinEdgeBackward[i]<<" high bin backw: "<<HighBinEdgeBackward[i]<<" high bin fwd: "<<HighBinEdgeForward[j]<<endl;
    //     }
    //     else if(LowBinEdgeBackward[i] != LowBinEdgeForward[j] && HighBinEdgeBackward[i] != HighBinEdgeForward[j]){
    //       cout<<" low bin bkwd: "<<LowBinEdgeBackward[i]<<" high bin bwd: "<<HighBinEdgeBackward[i]<<" low bin fwd: "<<LowBinEdgeForward[j]<<" high bin fwd: "<<HighBinEdgeForward[j]<<endl;
    //     }
    //   }
    // }

    // theApp.Terminate();


    return 0;
}

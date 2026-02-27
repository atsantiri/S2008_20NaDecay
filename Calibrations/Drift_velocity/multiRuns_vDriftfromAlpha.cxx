// Drift velocity from multiple runs to make a curve

#include "ActDataManager.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActModularData.h"
#include "ActModularData.h"
#include "ActMergerData.h"
#include "ActTPCData.h"
#include "ActCluster.h"
#include "ActVoxel.h"
#include "ActCutsManager.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include "TLegend.h"

double vDriftfromAlpha(int run)
{
    bool inspect{false};
    ROOT::EnableImplicitMT();
    // Get the data run 123
    auto df{ROOT::RDataFrame("GETTree", TString::Format("../../RootFiles/Cluster/Clusters_Run_%04d.root", run))};
    auto vdrift{0.};
    // df.Describe().Print();
    std::cout << "Processing Run " << run << std::endl;

    // Define last point of cluster in x y z, as the projection of the alpha track
    auto dfLastPoint = df.Define("fLastPoint", [](ActRoot::TPCData &d)
                                 {
        if(d.fClusters.size() != 1)
        {
            ROOT::Math::XYZPointF lastPoint {-1000, -1000, -1000};
            return lastPoint;
        }
        else
        {
            auto cluster {d.fClusters[0]};
            auto line {cluster.GetRefToLine()};
            auto dir {line.GetDirection()};
            cluster.SortAlongDir(dir);
            auto lastVoxel {cluster.GetRefToVoxels().back()};
            auto projectionPointLine {line.ProjectionPointOnLine(lastVoxel.GetPosition())};
            return projectionPointLine;
        } }, {"TPCData"});

    // Define other point
    auto dfOtherPoint = dfLastPoint.Define("fOtherPoint", [](ActRoot::TPCData &d)
                                           {
        if(d.fClusters.size() != 1)
        {
            ROOT::Math::XYZPointF otherPoint {-1000, -1000, -2000};
            return otherPoint;
        }
        else
        {
            auto cluster {d.fClusters[0]};
            auto line {cluster.GetRefToLine()};
            auto otherPoint {line.MoveToX(-50)};
            return otherPoint;
        } }, {"TPCData"});

    // Define the coordinates of the points
    auto dfXY = dfOtherPoint
                    .Define("fLastX", "fLastPoint.X()")
                    .Define("fLastY", "fLastPoint.Y()")
                    .Define("fOtherX", "fOtherPoint.X()")
                    .Define("fOtherY", "fOtherPoint.Y()");

    // // Plot
    // auto c = new TCanvas("c", "Points XY", 1000, 500);
    // auto hLast = dfXY.Histo2D({"hLast", "LastPoint XY;X [pads];Y [pads]", 1000, -100, 100, 1000, -100, 100}, "fLastX", "fLastY");
    // auto hOther = dfXY.Histo2D({"hOther", "OtherPoint XY;X [pads];Y [pads]", 1000, -100, 100, 1000, -100, 100}, "fOtherX", "fOtherY");
    // hLast->DrawClone("colz");
    // hOther->DrawClone("same");

    // // Create the lines and draw them with the foreach
    // int counter = 0;
    // dfXY.Foreach(
    //     [&](float otherX, float otherY, float lastX, float lastY)
    //     {
    //         counter++;
    //         auto line = new TLine(otherX, otherY, lastX, lastY);
    //         line->SetLineColorAlpha(kBlue, 0.3); // transparente para ver cruces
    //         if (counter % 50 == 0 && lastX > 5 && lastX < 60)
    //             line->Draw("same");
    //     },
    //     {"fOtherX", "fOtherY", "fLastX", "fLastY"});

    // With the plot I guess the alpha source is at (-27, 41)
    float xSource{-30.};
    float ySource{38.};

    auto dfDrift = dfXY.Define("fDeltaZ", [&](ActRoot::TPCData &d)
                               {
        if(d.fClusters.size() != 1)
            return -1000.;
        else
        {
            auto cluster {d.fClusters[0]};
            auto line {cluster.GetRefToLine()};
            auto dir {line.GetDirection()};
            cluster.SortAlongDir(dir);
            //auto firstVoxel {cluster.GetRefToVoxels().front()};
            //auto projectionFirstPointLine {line.ProjectionPointOnLine(firstVoxel.GetPosition())};
            auto lastVoxel {cluster.GetRefToVoxels().back()};
            auto projectionLastPointLine {line.ProjectionPointOnLine(lastVoxel.GetPosition())};
            auto zSource {line.MoveToX(xSource).Z()};
            double deltaZ = projectionLastPointLine.Z() - zSource;
            return deltaZ *0.32; // Conversion factor from bin time bucket to micro seconds. 4 time buckets in 1 bin time bucket. 1 time bucket is 12.5MHz.
        } }, {"TPCData"})
                       .Define("fLxy", [&](ActRoot::TPCData &d)
                               {
        if(d.fClusters.size() != 1)
            return -1000.;
        else
        {
            auto cluster {d.fClusters[0]};
            auto line {cluster.GetRefToLine()};
            auto dir {line.GetDirection()};
            cluster.SortAlongDir(dir);
            auto lastVoxel {cluster.GetRefToVoxels().back()};
            auto projectionPointLine {line.ProjectionPointOnLine(lastVoxel.GetPosition())};
            double lxy = TMath::Sqrt(TMath::Power(projectionPointLine.X() - xSource, 2) + TMath::Power(projectionPointLine.Y() - ySource, 2));
            return (lxy * 2)/10; // Conversion factor from pads to cm
        } }, {"TPCData"})
                       .Define("fDeltaZSquare", "fDeltaZ * fDeltaZ")
                       .Define("fLxySquare", "fLxy * fLxy");

    // Plot the DeltaZ and lxy
    auto graphDrift = dfDrift.Graph("fDeltaZ", "fLxy");
    graphDrift->SetTitle("Delta Z vs Lxy;#Delta Z [#mus]; Lxy [cm]");
    graphDrift->GetXaxis()->SetRangeUser(-30, 30);
    graphDrift->GetYaxis()->SetRangeUser(-10, 30);

    // Linearize the graph
    auto graphDriftLinear = dfDrift.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLinear->SetTitle("Delta Z^2 vs Lxy^2;#Delta Z^2 [#mus^2]; Lxy^2 [cm^2]");

    // auto c1 = new TCanvas("c1", "Delta Z vs Lxy", 1400, 800);
    // c1->DivideSquare(2);
    // c1->cd(1);
    // graphDrift->DrawClone("AP");
    // c1->cd(2);
    // graphDriftLinear->DrawClone("AP");

    // Cuts for good events (no broad region) and for each line
    ActRoot::CutsManager<std::string> cuts;
    // Gas PID
    cuts.ReadCut("goodEvents", "./Inputs/cut_DriftVelocity_GoodAlphaEvents.root");
    cuts.ReadCut("first", TString::Format("./Inputs/cut_firstPeak_%d.root", run).Data());
    cuts.ReadCut("second", TString::Format("./Inputs/cut_secondPeak_%d.root", run).Data());
    cuts.ReadCut("third", TString::Format("./Inputs/cut_thirdPeak_%d.root", run).Data());

    auto dfFiltered = dfDrift.Filter([&](double lxy, double deltaZ)
                                     { return cuts.IsInside("goodEvents", deltaZ, lxy); },
                                     {"fLxy", "fDeltaZ"});

    auto graphDriftFiltered = dfFiltered.Graph("fDeltaZ", "fLxy");
    graphDriftFiltered->SetTitle("Delta Z vs Lxy ;#Delta Z [#mus]; Lxy [cm]");
    auto graphDriftFilteredLinear = dfFiltered.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftFilteredLinear->SetTitle("Delta Z^2 vs Lxy^2;#Delta Z^2 [#mus^2];Lxy^2 [cm^2]");

    // Do graphs for each peak
    auto dfFirst = dfFiltered.Filter([&](double lxy2, double deltaZ2)
                                     { return cuts.IsInside("first", deltaZ2, lxy2); },
                                     {"fLxySquare", "fDeltaZSquare"});
    auto graphDriftLineFirst = dfFirst.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLineFirst->SetTitle("Delta Z^2 vs Lxy^2 (first peak);#Delta Z^2 [#mus^2];Lxy^2 [cm^2]");
    graphDriftLineFirst->Fit("pol1", "Q");
    auto f1{graphDriftLineFirst->GetFunction("pol1")};
    f1->SetLineColor(kRed);
    auto dfSecond = dfFiltered.Filter([&](double lxy2, double deltaZ2)
                                      { return cuts.IsInside("second", deltaZ2, lxy2); },
                                      {"fLxySquare", "fDeltaZSquare"});
    auto graphDriftLineSecond = dfSecond.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLineSecond->SetTitle("Delta Z^2 vs Lxy^2 (second peak);#Delta Z^2 [#mus^2];Lxy^2 [cm^2]");
    graphDriftLineSecond->Fit("pol1", "Q");
    auto f2{graphDriftLineSecond->GetFunction("pol1")};
    auto dfThird = dfFiltered.Filter([&](double lxy2, double deltaZ2)
                                     { return cuts.IsInside("third", deltaZ2, lxy2); },
                                     {"fLxySquare", "fDeltaZSquare"});
    auto graphDriftLineThird = dfThird.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLineThird->SetTitle("Delta Z^2 vs Lxy^2 (third peak);#Delta Z^2 [#mus^2];Lxy^2 [cm^2]");
    graphDriftLineThird->Fit("pol1", "Q");
    auto f3{graphDriftLineThird->GetFunction("pol1")};

    // auto c3 = new TCanvas("c3", "Delta Z vs Lxy lines", 2100, 700);
    // c3->DivideSquare(3);
    // c3->cd(1);
    // graphDriftLineFirst->DrawClone("AP");
    // f1->Draw("same");
    // c3->cd(2);
    // graphDriftLineSecond->DrawClone("AP");
    // f2->Draw("same");
    // c3->cd(3);
    // graphDriftLineThird->DrawClone("AP");
    // f3->Draw("same");

    // Draw them also in the filtered plot
    if (inspect)
    {
        auto c2 = new TCanvas("c2", "Delta Z vs Lxy filtered", 1400, 800);
        c2->DivideSquare(2);
        c2->cd(1);
        graphDriftFiltered->DrawClone("AP");
        c2->cd(2);
        graphDriftFilteredLinear->DrawClone("AP");

        c2->cd(1);
        cuts.DrawCut("goodEvents");
        c2->cd(2);
        cuts.DrawCut("first");
        cuts.DrawCut("second");
        cuts.DrawCut("third");
        f1->DrawClone("same");
        f2->SetLineColor(kGreen);
        f2->DrawClone("same");
        f3->SetLineColor(kBlue);
        f3->DrawClone("same");
        // Text of the fit parameters
        auto t1 = new TLatex(50, 30000, TString::Format("First peak: Vdrift = %.2f#pm%.2f ", TMath::Sqrt(-f1->GetParameter(1)), TMath::Sqrt(f1->GetParError(1))));
        auto t2 = new TLatex(50, 28000, TString::Format("Second peak: Vdrift = %.2f#pm%.2f", TMath::Sqrt(-f2->GetParameter(1)), TMath::Sqrt(f2->GetParError(1))));
        auto t3 = new TLatex(50, 25000, TString::Format("Third peak: Vdrift = %.2f#pm%.2f", TMath::Sqrt(-f3->GetParameter(1)), TMath::Sqrt(f3->GetParError(1))));
        t1->DrawClone();
        t2->DrawClone();
        t3->DrawClone();

        // Sleep while waiting for user
        c2->cd();
        c2->Update();
        c2->WaitPrimitive("dummy", "more dummy");
        c2->Update();
    }
    vdrift = TMath::Mean(3, (double[]){TMath::Sqrt(-f1->GetParameter(1)), TMath::Sqrt(-f2->GetParameter(1)), TMath::Sqrt(-f3->GetParameter(1))});

    return vdrift;
}

void multiRuns_vDriftfromAlpha()
{
    /* List of runs for drift measurements
    Run# Vdt (V) Vm (V) Pad Gain (fC) + Notes
    19 	6000 	490 	120
    20 	5500 	490 	120 - 1 CoBo failed
    21 	5500 	490 	120
    22 	4500 	490 	120
    23 	4000 	490 	120
    24 	3000 	490 	120
    */

    std::map<int, double> data = {
        {19, (6000. - 490.) / 23.5},
        {21, (5500. - 490.) / 23.5},
        {22, (4500. - 490.) / 23.5},
        {23, (4000. - 490.) / 23.5},
        {24, (3000. - 490.) / 23.5}};

    std::vector<double> drifts;
    auto *hDrifts = new TGraph();
    auto nPoints{0};

    for (const auto &m : data)
    {
        auto drift = vDriftfromAlpha(m.first);
        drifts.push_back(drift);
        std::cout << m.second << " " << drift << std::endl;
        hDrifts->SetPoint(nPoints, m.second, drift);
        nPoints++;
    }

    auto c3 = new TCanvas("c3", "Drifts vs E", 1400, 800);
    c3->cd();

    hDrifts->SetMarkerStyle(20);
    hDrifts->SetMarkerSize(1.2);
    hDrifts->SetMarkerColor(kBlue);
    hDrifts->SetTitle("Drift vs E;E [V/cm];Drift [cm/us]");

    hDrifts->Draw("AP");
    c3->Update();

    // Compare with graph from Garfield
    std::ifstream fIn("./Inputs/H2-95_iC4H10-5_950mbar.print");
    if (!fIn.is_open())
    {
        std::cout << "Did not find Garfield file" << std::endl;
        return;
    }
    auto *hTheo = new TGraph();
    int n{0};
    std::string line;
    int iline{0};
    while (std::getline(fIn, line))
    {
        iline++;
        if (iline < 13)
            continue;
        if (iline > 32)
            break;
        std::istringstream iss(line);
        double E, V;
        if (!(iss >> E >> V))
            continue;
        hTheo->SetPoint(n, E, V);
        n++;
    }
    fIn.close();
    c3->cd();
    hTheo->SetLineColor(kRed);
    hTheo->Draw("L same");

    auto l = new TLegend(0.65, 0.1, 0.9, 0.3);
    l->AddEntry(hDrifts, "Drift measurements");
    l->AddEntry(hTheo, "Garfield");
    l->Draw();

    auto run129 = new TLine((6000. - 480.) / 23.5, 0., (6000. - 480.) / 23.5, 1.);
    run129->SetLineColor(kMagenta);
    run129->SetLineWidth(2);
    run129->SetLineStyle(kDashed);
    run129->Draw("same");
    c3->SaveAs("vdrift_runs19-24.png");
}
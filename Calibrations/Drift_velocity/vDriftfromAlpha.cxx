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
#include "TEllipse.h"

void vDriftfromAlpha()
{
    ROOT::EnableImplicitMT();
    // Get the data run 123
    int run = 19;
    auto df{ROOT::RDataFrame("GETTree", TString::Format("../../RootFiles/Cluster/Clusters_Run_00%d.root",run))};

    // df.Describe().Print();

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

    // Plot
    auto c = new TCanvas("c", "Points XY", 1000, 500);
    auto hLast = dfXY.Histo2D({"hLast", "LastPoint XY;X [pads];Y [pads]", 1000, -100, 150, 1000, -100, 150}, "fLastX", "fLastY");
    auto hOther = dfXY.Histo2D({"hOther", "OtherPoint XY;X [pads];Y [pads]", 1000, -100, 150, 1000, -100, 150}, "fOtherX", "fOtherY");
    hLast->DrawClone("colz");
    hOther->DrawClone("same");

    // Create the lines and draw them with the foreach
    int counter = 0;
    dfXY.Foreach(
        [&](float otherX, float otherY, float lastX, float lastY)
        {
            counter++;
            auto line = new TLine(otherX, otherY, lastX, lastY);
            line->SetLineColorAlpha(kBlue, 0.3); // transparente para ver cruces
            if (counter % 50 == 0 && lastX > 5 && lastX < 60)
                line->Draw("same");
        },
        {"fOtherX", "fOtherY", "fLastX", "fLastY"});

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
            return deltaZ *0.32; // Conversion factor from btb to micro seconds
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
            return (lxy * 2) / 10; // Conversion factor from pads to cm
        } }, {"TPCData"})
                       .Define("fDeltaZSquare", "fDeltaZ * fDeltaZ")
                       .Define("fLxySquare", "fLxy * fLxy");

    // Plot the DeltaZ and lxy
    auto graphDrift = dfDrift.Graph("fDeltaZ", "fLxy");
    graphDrift->SetTitle("Delta Z vs Lxy;#Deltat [#mus]; #Deltaxy [cm]");

    // Linearize the graph
    auto graphDriftLinear = dfDrift.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLinear->SetTitle("Delta Z^2 vs Lxy^2;(#Deltat)^{2} [#mus^{2}];(#Deltaxy)^{2} [cm^{2}]");


    auto c1 = new TCanvas("c1", "Delta Z vs Lxy", 1400, 800);
    c1->DivideSquare(2);
    c1->cd(1);
    graphDrift->GetXaxis()->SetRangeUser(-20,20);
    graphDrift->GetYaxis()->SetRangeUser(0,20);
    graphDrift->DrawClone("AP");
    c1->cd(2);
    graphDriftLinear->DrawClone("AP");

    // Cuts for good events (no broad region) and for each line
    ActRoot::CutsManager<std::string> cuts;
    // Gas PID
    cuts.ReadCut("goodEvents", "./Inputs/cut_DriftVelocity_GoodAlphaEvents.root");
    cuts.ReadCut("first", TString::Format("./Inputs/cut_firstPeak_%d.root", run).Data());
    cuts.ReadCut("second", TString::Format("./Inputs/cut_secondPeak_%d.root", run).Data());
    cuts.ReadCut("third", TString::Format("./Inputs/cut_thirdPeak_%d.root", run).Data());

    c1->cd(1);
    cuts.DrawCut("goodEvents");

    auto dfFiltered = dfDrift.Filter([&](double lxy, double deltaZ)
                                { return cuts.IsInside("goodEvents", deltaZ, lxy); },
                                {"fLxy", "fDeltaZ"});

    auto graphDriftFiltered = dfFiltered.Graph("fDeltaZ", "fLxy");
    graphDriftFiltered->SetTitle("Delta Z vs Lxy ;#Deltat [#mus]; #Deltaxy [cm]");
    auto graphDriftFilteredLinear = dfFiltered.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftFilteredLinear->SetTitle("Delta Z^2 vs Lxy^2;(#Deltat)^{2} [#mus^{2}];(#Deltaxy)^{2} [cm^{2}]");
    auto c2 = new TCanvas("c2", "Delta Z vs Lxy filtered", 1400, 800);
    c2->DivideSquare(2);
    c2->cd(1);
    graphDriftFiltered->DrawClone("AP");
    c2->cd(2);
    graphDriftFilteredLinear->DrawClone("AP");

    // Do graphs for each peak
    auto dfFirst = dfFiltered.Filter([&](double lxy2, double deltaZ2)
                                { return cuts.IsInside("first", deltaZ2, lxy2); },
                                {"fLxySquare", "fDeltaZSquare"});
    auto graphDriftLineFirst = dfFirst.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLineFirst->SetTitle("Delta Z^2 vs Lxy^2 (first peak);(#Deltat)^{2} [#mus^{2}];(#Deltaxy)^{2} [cm^{2}]");
    graphDriftLineFirst->Fit("pol1");
    auto f1 {graphDriftLineFirst->GetFunction("pol1")};
    f1->SetLineColor(kRed);
    auto dfSecond = dfFiltered.Filter([&](double lxy2, double deltaZ2)
                                { return cuts.IsInside("second", deltaZ2, lxy2); },
                                {"fLxySquare", "fDeltaZSquare"});
    auto graphDriftLineSecond = dfSecond.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLineSecond->SetTitle("Delta Z^2 vs Lxy^2 (second peak);(#Deltat)^{2} [#mus^{2}];(#Deltaxy)^{2} [cm^{2}]");
    graphDriftLineSecond->Fit("pol1");
    auto f2 {graphDriftLineSecond->GetFunction("pol1")};
    auto dfThird = dfFiltered.Filter([&](double lxy2, double deltaZ2)
                                { return cuts.IsInside("third", deltaZ2, lxy2); },
                                {"fLxySquare", "fDeltaZSquare"});
    auto graphDriftLineThird = dfThird.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLineThird->SetTitle("Delta Z^2 vs Lxy^2 (third peak);(#Deltat)^{2} [#mus^{2}];(#Deltaxy)^{2} [cm^{2}]");
    graphDriftLineThird->Fit("pol1");
    auto f3 {graphDriftLineThird->GetFunction("pol1")};
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
    c2->cd(2);
    f1->DrawClone("same");
    cuts.DrawCut("first");
    f2->SetLineColor(kGreen);
    f2->DrawClone("same");
    cuts.DrawCut("second");
    f3->SetLineColor(kBlue);
    f3->DrawClone("same");
    cuts.DrawCut("third");

    // Text of the fit parameters
    auto t1 = new TLatex(60, 200, TString::Format("First peak: Vdrift = %.2f#pm%.2f ", TMath::Sqrt(-f1->GetParameter(1)), TMath::Sqrt(f1->GetParError(1))));
    auto t2 = new TLatex(60, 180, TString::Format("Second peak: Vdrift = %.2f#pm%.2f", TMath::Sqrt(-f2->GetParameter(1)), TMath::Sqrt(f2->GetParError(1))));
    auto t3 = new TLatex(60, 160, TString::Format("Third peak: Vdrift = %.2f#pm%.2f", TMath::Sqrt(-f3->GetParameter(1)), TMath::Sqrt(f3->GetParError(1))));
    t1->DrawClone();
    t2->DrawClone();
    t3->DrawClone();
    
}

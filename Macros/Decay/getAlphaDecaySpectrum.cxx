#include "ActCluster.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSRIM.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"

#include <fstream>
// run 129 is processed without any user functions enabled

double calcTL(ROOT::Math::XYZPointF A, ROOT::Math::XYZPointF B, double drift)
{
    A += ROOT::Math::XYZVector {0.5, 0.5,0.5}; // when converting a bin point to physical units which wasnt already corrected
    A.SetX(A.X() * 2.);               // pad in mm
    A.SetY(A.Y() * 2.);
    A.SetZ(A.Z() * drift);

    B += ROOT::Math::XYZVector {0.5, 0.5,
                                0.5}; // when converting a bin point to physical units which wasnt already corrected
    B.SetX(B.X() * 2.);               // pad in mm
    B.SetY(B.Y() * 2.);
    B.SetZ(B.Z() * drift);

    return (A - B).R();
}

void getAlphaDecaySpectrum()
{
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EFilter};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain3.get());


    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d {*chain};

    ROOT::TThreadedObject<TH2D> h2d("hPad", "Pad plane;X [pad];Y [pad]", 128, 0, 128, 128, 0, 128);
    ROOT::TThreadedObject<TH1D> hEne("hEne", "Alpha Energy; E_{#alpha} [MeV]; ", 200, 0, 10);

    ROOT::Math::XYZPointF beamEndPoint {105, 65, 0}; // visually inspected output of Pipe0
    double maxDist {2.};                             // max distance from beam stop point
    double edge {2.}; // min distance from x wall to remove tracks that out of the pad plane

    ActRoot::InputParser parserDet {"../../configs/detector.conf"};
    auto bl1 {parserDet.GetBlock("Merger")};
    auto drift {bl1->GetDouble("DriftFactor")};
    std::cout << drift << std::endl;

    auto* srim {new ActPhysics::SRIM()};
    srim->ReadTable("HeInGas", "../../Calibrations/SRIM/4He_950mbar_95-5.txt");

    d.Foreach(
        [&](ActRoot::TPCData& tpc)
        {
            if(tpc.fRPs.empty())
            {
                auto& clusters = tpc.fClusters;
                for(auto& c : clusters)
                {
                    if(!c.GetIsBeamLike())
                    {
                        auto line {c.GetRefToLine()};
                        auto dir {line.GetDirection()};
                        c.SortAlongDir();
                        auto firstPoint {line.ProjectionPointOnLine(c.GetRefToVoxels().front().GetPosition())};
                        auto lastPoint {line.ProjectionPointOnLine(c.GetRefToVoxels().back().GetPosition())};
                        double lxyA = TMath::Sqrt(TMath::Power(firstPoint.X() - beamEndPoint.X(), 2) +
                                                  TMath::Power(firstPoint.Y() - beamEndPoint.Y(), 2));
                        double lxyB = TMath::Sqrt(TMath::Power(lastPoint.X() - beamEndPoint.X(), 2) +
                                                  TMath::Power(lastPoint.Y() - beamEndPoint.Y(), 2));

                        //remove tracks that don't stop in the pad plane, note this will affect relative peak intensities
                        bool isNearEdge {(128. - firstPoint.X()) <= edge || (128. - firstPoint.Y()) <= edge ||
                                         (128. - lastPoint.X()) <= edge || (128. - lastPoint.Y()) <= edge};
                        if((lxyA <= maxDist || lxyB <= maxDist) && !isNearEdge)
                        {
                            for(const auto& v : c.GetVoxels())
                            {
                                auto& pos {v.GetPosition()};
                                h2d->Fill(pos.X(), pos.Y());
                            }
                            double TL = calcTL(firstPoint, lastPoint, drift);
                            double ene {srim->EvalInitialEnergy("HeInGas", 0, TL)};
                            hEne->Fill(ene);
                        }
                    }
                }
            }
        },
        {"TPCData"});

    // Draw
    auto* c0 {new TCanvas {"c0", ""}};
    h2d.Merge()->DrawClone("colz");

    auto* c1 {new TCanvas {"c1", ""}};
    hEne.Merge()->DrawClone();

    //   alpha energies from nds to draw
    std::vector<double> energies {{2153., 2482., 3806., 4434., 4891.}};
    std::vector<double> intensities {{15.96, 0.583, 0.241, 2.88, 0.174}};
    int i = 0;
    for(auto& s : energies)
    {
        TLine* st = new TLine(s * 1e-3, 0, s * 1e-3, 7000);
        st->SetLineColor(2);
        st->SetLineStyle(2);
        st->Draw("same");
        auto t = new TLatex(s * 1e-3 - 0.02, 7000 - i * 500, Form("%.2f (%.2f%%)", s * 1e-3,intensities[i]));
        t->SetTextColor(2);
        t->SetTextSize(0.03);
        t->Draw("same");
        i++;
    }
}
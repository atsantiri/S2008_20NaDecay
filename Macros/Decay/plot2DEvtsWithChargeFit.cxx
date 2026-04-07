#include "ActCluster.h"
#include "ActDataManager.h"
#include "ActDetectorManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSRIM.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"
#include "ActTypes.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include <TF1.h>

#include <fstream>

#include "./miscFunctions.cxx"

void plot2DEvtsWithChargeFit()
{
    // ========================================================================================================================================
    // Input Stuff
    // ========================================================================================================================================
    ROOT::Math::XYZPointF beamEndPoint {105, 65, 0}; // visually inspected output of Pipe0
    double maxDist {2.};                             // max distance from beam stop point
    double edge {5.};                                // min distance from Z wall to remove tracks
                                                     // that come from the pad plane - 20Na drifted before decaying
    int kmax {10};                                   // how many events to look at
    std::string srimfile {"../../Calibrations/SRIM/4He_950mbar_95-5.txt"};
    double qval {-4.73};
    double ma {4.0026};      // amu
    double mo16 {15.994914}; // amu
    // ========================================================================================================================================
    // ========================================================================================================================================

    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EFilter};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain3.get());

    ROOT::RDataFrame d {*chain};

    ActRoot::TPCParameters TPCparams {"Actar"};
    ActRoot::InputParser parserDet {"../../configs/detector.conf"};
    auto bl1 {parserDet.GetBlock("Merger")};
    auto drift {bl1->GetDouble("DriftFactor")};
    std::cout << drift << std::endl;


    auto* srim {new ActPhysics::SRIM()};
    const std::string& particleKey = "HeInGas";
    srim->ReadTable(particleKey, srimfile);
    auto srimSpline {buildSRIMspline(srim, 100, particleKey)};

    ROOT::Math::XYZVectorF beam = {1, 0, 0};
    using SQ = std::vector<std::pair<double, double>>;
    int k = 0;
    auto* ct {new TCanvas {"ct", "Testing", 1500, 600}};
    ct->Divide(2);
    int badfits {0};

    d.Foreach(
        [&](ActRoot::TPCData& tpc)
        {
            if(tpc.fRPs.empty())
            {
                auto& clusters = tpc.fClusters;
                for(auto& c : clusters)
                {
                    if((k < kmax))
                    {
                        if(!c.GetIsBeamLike())
                        {
                            auto line {c.GetRefToLine()};
                            auto dir {line.GetDirection()};
                            c.SortAlongDir();

                            auto startVoxelPos = findStartVoxel(c, true);
                            auto endVoxelPos = findStartVoxel(c, false);

                            double lxy = calcDistXY(startVoxelPos, beamEndPoint);

                            // remove tracks that come from pad plane, note this will affect relative peak
                            // intensities
                            bool isNearEdge {
                                (startVoxelPos.Z() * drift <= edge) || (startVoxelPos.Z() * drift >= (272 - edge)) ||
                                (endVoxelPos.Z() * drift <= edge) || (endVoxelPos.Z() * drift >= (272 - edge))};

                            // remove a forward solid angle of 65 deg, for events that make it out of the pad plane
                            double cosAngle = beam.Unit().Dot(dir.Unit());
                            bool forwardCone =
                                (cosAngle > std::cos(65. * TMath::DegToRad())) && (endVoxelPos.X() > beamEndPoint.X());

                            if((lxy <= maxDist) && !forwardCone && !isNearEdge)
                            {

                                // get charge profile
                                SQ sq = getChargeProf(c, startVoxelPos, endVoxelPos, dir, drift);

                                auto hcharge =
                                    new TH1D {TString::Format("hCharge%d", k), "Charge profile; L [mm]; Q", 50, 0, 150};

                                for(size_t i = 0; i < sq.size(); ++i)
                                    hcharge->Fill(sq[i].first, sq[i].second);

                                auto fitSpline =
                                    FitSRIMtoChargeProfile(hcharge, srimSpline, TString::Format("fSRIM%d", k));
                                hcharge->Fit(fitSpline, "QR0", "", hcharge->GetXaxis()->GetXmin(),
                                             hcharge->GetXaxis()->GetXmax());


                                double chi2 = fitSpline->GetChisquare();
                                double ndf = fitSpline->GetNDF();
                                double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
                                double TLfromCharge = calcTLfromChargeFit(fitSpline, hcharge->GetXaxis()->GetXmin(),
                                                                          hcharge->GetXaxis()->GetXmax());
                                double ene {srim->EvalInitialEnergy("HeInGas", 0, TLfromCharge)};
                                double TLfromVoxel = calcTLfromVoxel(startVoxelPos, endVoxelPos, line, drift);
                                double ene2 {srim->EvalInitialEnergy("HeInGas", 0, TLfromVoxel)};

                                double diff {(TLfromCharge - TLfromVoxel) / TLfromVoxel * 100.};
                                if(std::abs(diff) > 10)
                                    badfits++;

                                auto hCharge2D = new TH2F(TString::Format("hCharge2D%d", k), ";x [pad], y [pad]",
                                                          TPCparams.GetNPADSX(), 0, TPCparams.GetNPADSX(),
                                                          TPCparams.GetNPADSY(), 0, TPCparams.GetNPADSY());
                                for(const auto& v : c.GetVoxels())
                                    hCharge2D->Fill(v.GetPosition().X(), v.GetPosition().Y(), v.GetCharge());

                                ct->cd(1);
                                hcharge->DrawCopy("hist");
                                fitSpline->SetLineWidth(3);
                                fitSpline->DrawCopy("same");

                                ct->cd(2);
                                hCharge2D->Draw("LEGO2");

                                std::cout << "\n===== Fit Results: " << fitSpline->GetName() << " =====\n";
                                int npar = fitSpline->GetNpar();
                                for(int i = 0; i < npar; ++i)
                                {
                                    std::cout << std::setw(10) << fitSpline->GetParName(i) << " = " << std::setw(12)
                                              << fitSpline->GetParameter(i) << " ± " << std::setw(12)
                                              << fitSpline->GetParError(i) << "\n";
                                }
                                std::cout << "-----------------------------------\n";
                                std::cout << "Chi2 / NDF = " << chi2 << " / " << ndf << " = " << chi2ndf << "\n";

                                std::cout << "TL from charge: " << TLfromCharge
                                          << " vs TL from last voxel: " << TLfromVoxel << " ---- ( " << diff << " %)"
                                          << std::endl;

                                ct->Update();
                                ct->WaitPrimitive();
                                ct->Clear();
                                ct->Divide(2);

                                delete hCharge2D;

                                delete hcharge;

                                k++;
                            }
                        }
                    }
                }
            }
        },
        {"TPCData"});
}
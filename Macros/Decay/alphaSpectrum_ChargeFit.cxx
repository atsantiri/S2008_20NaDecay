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
#include "TFitResult.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include <TF1.h>

#include <fstream>

#include "./miscFunctions.cxx"
// run 129 is processed without any user functions enabled

void alphaSpectrum_ChargeFit()
{
    // ========================================================================================================================================
    // Input Stuff
    // ========================================================================================================================================
    bool debug {false}; // set to true if we're testing. Otherwise we will do only few events
    int kmax {10};      // how many events to look at
    std::string srimfile {"../../Calibrations/SRIM/4He_950mbar_95-5.txt"};
    double qval {-4.73};
    double ma {4.0026};      // amu
    double mo16 {15.994914}; // amu
    // ========================================================================================================================================
    // ========================================================================================================================================

    std::string fIn {"../../PostAnalysis/Outputs/tree_20Na_decays.root"};
    ROOT::RDataFrame d {"Decay_Tree", fIn};

    ActRoot::TPCParameters TPCparams {"Actar"};
    ActRoot::InputParser parserDet {"../../configs/detector.conf"};
    auto bl1 {parserDet.GetBlock("Merger")};
    auto drift {bl1->GetDouble("DriftFactor")};
    std::cout << drift << std::endl;

    auto hPad = new TH2F("hPad", "Pad;X [pad];Y [pad]", TPCparams.GetNPADSX(), 0, TPCparams.GetNPADSX(),
                         TPCparams.GetNPADSY(), 0, TPCparams.GetNPADSY());

    auto hSide = new TH2F("hSide", "Side;X [pad];Z [tb]", TPCparams.GetNPADSX(), 0, TPCparams.GetNPADSX(),
                          TPCparams.GetNPADSZ() / 4, 0, TPCparams.GetNPADSZ());

    auto hEne = new TH1F("hEne", "Alpha Energy s; E_{#alpha} [MeV]; ", 200, 0, 7);
    auto hEne2 = new TH1F("hEneVoxel", "Alpha Energy Voxel; E_{#alpha} [MeV]; ", 200, 0, 7);
    auto hEne3 = new TH1F("hEne3", "Alpha Energy lambda; E_{#alpha} [MeV]; ", 200, 0, 7);
    auto hEx = new TH1F("hEx", "Excitation Energy; E_{20Ne*} [MeV]; ", 200, 4, 14);
    auto hEx2 = new TH1F("hExVoxel", "Excitation Energy; E_{20Ne*} [MeV]; ", 200, 4, 14);
    std::cout << TPCparams.GetNPADSZ() << " " << TPCparams.GetREBINZ() << std::endl;

    auto* srim {new ActPhysics::SRIM()};
    const std::string& particleKey = "HeInGas";
    srim->ReadTable(particleKey, srimfile);
    double splinerange = 150.;
    auto srimSpline {buildSRIMspline(srim, splinerange, particleKey)};
    auto srimSplineJ {buildSRIMSplineJ(srim, particleKey, splinerange)};

    ROOT::Math::XYZVectorF beam = {1, 0, 0};

    std::vector<std::shared_ptr<TPolyLine>> ls;
    std::vector<double> xsA, ysA, xsB, ysB;
    using SQ = std::vector<std::pair<double, double>>;
    int k = 0;
    auto* ct {new TCanvas {"ct", "Testing", 800, 600}};
    ct->Divide(2);
    int badfits {0};
    auto nTot = *d.Count();

    if(debug)
    {
        std::cout << "========= Track Lengths =============" << std::endl;
        std::cout << "Charge fit s \t charge fit λ \t Last Voxel" << std::endl;
    }

    d.Foreach(
        [&](ActRoot::TPCData& tpc, std::vector<bool> isDecay)
        {
            auto& clusters = tpc.fClusters;
            for(auto& c : clusters)
            {
                if(((k < kmax) || !debug) && isDecay[c.GetClusterID()])
                {
                    auto line {c.GetRefToLine()};
                    auto dir {line.GetDirection()};
                    c.SortAlongDir();

                    auto startVoxelPos = findStartVoxel(c, true);
                    auto endVoxelPos = findStartVoxel(c, false);

                    for(const auto& v : c.GetVoxels())
                    {
                        auto& pos {v.GetPosition()};
                        hPad->Fill(pos.X(), pos.Y());
                        hSide->Fill(pos.X(), pos.Z() * drift);
                    }

                    // get charge profile
                    SQ sq = getChargeProf(c, startVoxelPos, endVoxelPos, dir, drift);
                    auto hcharge = new TH1D {TString::Format("hCharge%d", k), "Charge profile; L [mm]; Q", 50, 0, 150};

                    for(size_t i = 0; i < sq.size(); ++i)
                        hcharge->Fill(sq[i].first, sq[i].second);

                    auto fitSpline = FitSRIMtoChargeProfile(hcharge, srimSpline, TString::Format("fSRIM%d", k));
                    hcharge->Fit(fitSpline, "QR0", "", hcharge->GetXaxis()->GetXmin(), hcharge->GetXaxis()->GetXmax());

                    TH1D* hchargeNorm =
                        ConvertToLambdaSpace(hcharge, splinerange, TString::Format("hChargeLambda%d", k));
                    auto fitSplineJ =
                        FitSRIMtoChargeProfileJ(hchargeNorm, srimSplineJ, splinerange, TString::Format("fSRIMJ%d", k));
                    auto r = hchargeNorm->Fit(fitSplineJ, "QR0S", "", hchargeNorm->GetXaxis()->GetXmin(),
                                              hchargeNorm->GetXaxis()->GetXmax());

                    // auto ndf = fitSpline->GetNDF();
                    // double chi2 = ComputeChi2(hcharge, fitSpline, ndf);
                    // auto ndfJ = fitSplineJ->GetNDF();
                    // double chi2J = ComputeChi2(hcharge, fitSplineJ, ndfJ);

                    double TLfromCharge =
                        calcTLfromChargeFit(fitSpline, hcharge->GetXaxis()->GetXmin(), hcharge->GetXaxis()->GetXmax());
                    double ene {srim->EvalInitialEnergy("HeInGas", 0, TLfromCharge)};
                    double TLfromVoxel = calcTLfromVoxel(startVoxelPos, endVoxelPos, line, drift);
                    double ene2 {srim->EvalInitialEnergy("HeInGas", 0, TLfromVoxel)};
                    double TLfromCharge2 = 0.;
                    double ene3 {0.};

                    if(r.Get() && r->Status() == 0)
                    {
                        TLfromCharge2 = splinerange * fitSplineJ->GetParameter(1);
                        ene3 = srim->EvalInitialEnergy("HeInGas", 0, TLfromCharge2);
                    }

                    double diff {(TLfromCharge - TLfromVoxel) / TLfromVoxel * 100.};
                    if(std::abs(diff) > 10)
                        badfits++;

                    hEne->Fill(ene);
                    // double Ex {ene * (1 + ma / mo16) - qval};
                    // hEx->Fill(Ex);
                    hEne2->Fill(ene2);
                    hEne3->Fill(ene3);
                    // double Ex2 {ene2 * (1 + ma / mo16) - qval};
                    // hEx2->Fill(Ex2);
                    if(debug)
                    {
                        xsA.push_back(startVoxelPos.X());
                        ysA.push_back(startVoxelPos.Y());

                        ct->cd(1);
                        hcharge->DrawCopy("hist");
                        fitSpline->SetLineWidth(3);
                        fitSpline->SetLineColor(94);
                        fitSpline->DrawCopy("same");

                        ct->cd(2);
                        hchargeNorm->Draw("hist");
                        fitSplineJ->SetLineWidth(3);
                        fitSplineJ->SetLineColor(6);
                        fitSplineJ->DrawCopy("same");

                        std::cout << TLfromCharge << "\t\t" << TLfromCharge2 << "\t\t" << TLfromVoxel << std::endl;

                        ct->Update();
                        ct->WaitPrimitive();
                        ct->Clear();
                        ct->Divide(2);
                    }
                    delete hcharge;

                    k++;
                    if(!debug)
                        printProgress(nTot, k);
                }
            }
        },
        {"TPCData", "isDecay"});
    delete ct;

    // Draw
    std::cout << "size " << k << std::endl;
    std::cout << "Bad fits: " << badfits << " out of " << k << std::endl;

    auto* c0 {new TCanvas {"c0", "", 1500, 600}};
    c0->Divide(2);
    c0->cd(1);
    hPad->DrawClone("colz");

    if(debug)
    {
        auto* gBegins = new TGraph(xsA.size(), xsA.data(), ysA.data());
        gBegins->SetMarkerStyle(20);
        gBegins->SetMarkerSize(0.8);
        gBegins->SetMarkerColor(kRed);
        gBegins->Draw("P same");
    }

    c0->cd(2);
    hSide->DrawClone("colz");

    // auto* c1 {new TCanvas {"c1", ""}};
    // c1->DivideSquare(2);
    // c1->cd(1);
    // hEne->SetLineWidth(2);
    // hEne->DrawClone();
    // hEne2->SetLineColor(2);
    // hEne2->SetLineWidth(2);
    // hEne2->DrawClone("same");
    // hEne3->SetLineColor(9);
    // hEne3->SetLineWidth(2);
    // hEne3->DrawClone("same");
    // c1->cd(2);
    // hEx->DrawClone();
    // hEx2->SetLineColor(2);
    // hEx2->DrawClone("same");
    // auto l = new TLegend(0.6, 0.8, 0.9, 0.9);
    // l->AddEntry(hEx, "from charge fit s");
    // l->AddEntry(hEne3, "from charge fit λ");
    // l->AddEntry(hEx2, "from last voxel");
    // l->Draw();
    auto* c1 {new TCanvas {"c1", "", 1200, 500}};
    c1->DivideSquare(3);
    c1->cd(1);
    hEne->SetLineWidth(2);
    hEne->DrawClone();
    c1->cd(2);
    hEne2->SetLineWidth(2);
    hEne2->DrawClone();
    c1->cd(3);
    hEne3->SetLineWidth(2);
    hEne3->DrawClone();

    //   alpha energies from nds to draw
    std::vector<double> energies {{2153., 2482., 3806., 4434., 4683., 4891.}};
    std::vector<double> intensities {{15.96, 0.583, 0.241, 2.88, 0.09, 0.174}};
    std::vector<int> colors {{6, 2, 4, 6, 2, 4}};

    for(int cv = 1; cv < 4; cv++)
    {
        int i = 0;
        c1->cd(cv);
        for(auto& s : energies)
        {
            TLine* st = new TLine(s * 1e-3, 0, s * 1e-3, 1800 - i * 200);
            st->SetLineColor(colors[i]);
            st->SetLineStyle(2);
            st->Draw("same");
            // auto t = new TLatex(s * 1e-3 - 0.02, 1800 - i * 200, Form("%.2f (%.2f%%)", s * 1e-3, intensities[i]));
            auto t = new TLatex(s * 1e-3 - 0.02, 1800 - i * 200, Form("%.2f", s * 1e-3));
            t->SetTextColor(colors[i]);
            t->SetTextSize(0.03);
            t->Draw("same");
            i++;
        }
    }


    // //  20Ne levels
    // std::vector<double> exs {{7421., 7833., 9487., 10273., 10584., 10843.}};
    // // std::vector<double> intensities {{15.96, 0.583, 0.241, 2.88, 0.09, 0.174}};
    // // std::vector<int> colors {{6,2,4,6,2,4}};
    // i = 0;
    // c1->cd(2);
    // for(auto& s : exs)
    // {
    //     TLine* st = new TLine(s * 1e-3, 0, s * 1e-3, 1900 - i * 200);
    //     st->SetLineColor(colors[i]);
    //     st->SetLineStyle(2);
    //     st->Draw("same");
    //     auto t = new TLatex(s * 1e-3 - 0.02, 1900 - i * 200, Form("%.2f (%.2f%%)", s * 1e-3, intensities[i]));
    //     t->SetTextColor(colors[i]);
    //     t->SetTextSize(0.03);
    //     t->Draw("same");
    //     i++;
    // }
}
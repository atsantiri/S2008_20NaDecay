#include "ActCluster.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSRIM.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"
#include "ActTypes.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"

#include <fstream>

#include "./miscFunctions.cxx"

void getAlphaDecaySpectrum_LastVoxel()
{
    ROOT::EnableImplicitMT();
    std::string fIn {"../../PostAnalysis/Outputs/tree_20Na_decays.root"};
    ROOT::RDataFrame d {"Decay_Tree", fIn};

    ROOT::TThreadedObject<TH2F> h2d("hPad", "Pad plane;X [pad];Y [pad]", 128, 0, 128, 128, 0, 128);
    ROOT::TThreadedObject<TH2F> h2dS("hPadS", "Pad plane (Sample events);X [pad];Y [pad]", 128, 0, 128, 128, 0, 128);
    ROOT::TThreadedObject<TH1F> hEne("hEne", "Alpha Energy; E_{#alpha} [MeV]; ", 400, 0, 7);
    ROOT::TThreadedObject<TH1F> hEx("hEx", "Excitation Energy; E_{20Ne*} [MeV]; ", 400, 4, 14);


    ActRoot::InputParser parserDet {"../../configs/detector.conf"};
    auto bl1 {parserDet.GetBlock("Merger")};
    auto drift {bl1->GetDouble("DriftFactor")};
    std::cout << drift << std::endl;

    auto* srim {new ActPhysics::SRIM()};
    srim->ReadTable("HeInGas", "../../Calibrations/SRIM/4He_950mbar_95-5.txt");
    // srim->ReadTable("HeInGas", "../../Calibrations/SRIM/He_actar_192g_1.07_38.8stoich.txt");

    double qval {-4.73};
    double ma {4.0026};      // amu
    double mo16 {15.994914}; // amu

    std::vector<double> xsA, ysA, xsB, ysB;

    d.Foreach(
        [&](ActRoot::TPCData& tpc, std::vector<bool> isDecay)
        {
            if(tpc.fRPs.empty())
            {
                auto& clusters = tpc.fClusters;
                for(auto& c : clusters)
                {
                    if(isDecay[c.GetClusterID()])
                    {
                        auto line {c.GetRefToLine()};
                        c.SortAlongDir();

                        auto firstPoint {c.GetRefToVoxels().front().GetPosition()};
                        auto lastPoint {c.GetRefToVoxels().back().GetPosition()};

                        for(const auto& v : c.GetVoxels())
                        {
                            auto& pos {v.GetPosition()};
                            h2d->Fill(pos.X(), pos.Y());
                        }
                        double TL = calcTLfromVoxel(firstPoint, lastPoint, line, drift);
                        double ene {srim->EvalInitialEnergy("HeInGas", 0, TL)};
                        hEne->Fill(ene);
                        double Ex {ene * (1 + ma / mo16) - qval};
                        hEx->Fill(Ex);
                    }
                }
            }
        },
        {"TPCData", "isDecay"});


    // Draw
    auto* c0 {new TCanvas {"c0", "", 800, 600}};
    h2d.Merge()->DrawClone("colz");

    auto* c1 {new TCanvas {"c1", "", 1200,600}};
    c1->DivideSquare(2);
    c1->cd(1);
    hEne.Merge()->DrawClone();
    c1->cd(2);
    hEx.Merge()->DrawClone();

    //   alpha energies from nds to draw
    std::vector<double> energies {{2153., 2482., 3806., 4434., 4683., 4891.}};
    std::vector<double> intensities {{15.96, 0.583, 0.241, 2.88, 0.09, 0.174}};
    std::vector<int> colors {{6, 2, 4, 6, 2, 4}};
    int i = 0;
    c1->cd(1);
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


    //  20Ne levels
    std::vector<double> exs {{7421., 7833., 9487., 10273., 10584., 10843.}};
    // std::vector<double> intensities {{15.96, 0.583, 0.241, 2.88, 0.09, 0.174}};
    // std::vector<int> colors {{6,2,4,6,2,4}};
    i = 0;
    c1->cd(2);
    for(auto& s : exs)
    {
        TLine* st = new TLine(s * 1e-3, 0, s * 1e-3, 1900 - i * 200);
        st->SetLineColor(colors[i]);
        st->SetLineStyle(2);
        st->Draw("same");
        auto t = new TLatex(s * 1e-3 - 0.02, 1900 - i * 200, Form("%.2f (%.2f%%)", s * 1e-3, intensities[i]));
        t->SetTextColor(colors[i]);
        t->SetTextSize(0.03);
        t->Draw("same");
        i++;
    }
}
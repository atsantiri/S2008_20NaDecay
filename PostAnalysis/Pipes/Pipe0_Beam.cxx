#include "ActCluster.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TROOT.h"

#include <atomic>
#include <utility>

void Pipe0_Beam()
{
    std::string dataconf {"./../configs/data.conf"};

    ROOT::EnableImplicitMT();
    // Read data
    ActRoot::DataManager datman {dataconf, ActRoot::ModeType::EReadSilMod};
    auto chain {datman.GetJoinedData()};
    auto chain2 {datman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain3 {datman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());
    ROOT::RDataFrame df {*chain};

    // Get GATCONF values
    auto defGat {df.Define("GATCONF", [](ActRoot::ModularData& mod)
                           { return static_cast<int>(mod.fLeaves["GATCONF"]); }, {"ModularData"})};

    // Book histograms
    auto hGATCONF {defGat.Histo1D("GATCONF")};
    ROOT::TThreadedObject<TH2D> h2d("hPad", "Pad plane;X [pad];Y [pad]", 128, 0, 128, 128, 0, 128);

    // And cound CFA triggers
    std::atomic<unsigned long int> cfa {};
    defGat.Foreach(
        [&](const int& gatconf, ActRoot::TPCData& tpc)
        {
            if(gatconf == 64)
                cfa++;

            for(const auto& cluster : tpc.fClusters)
            {
                for(const auto& v : cluster.GetVoxels())
                {
                    auto& pos {v.GetPosition()};
                    h2d->Fill(pos.X(), pos.Y());
                }
            }
        },
        {"GATCONF", "TPCData"});


    // Draw
    auto* c0 {new TCanvas {"c00", "Pipe 0 canvas 0"}};
    c0->DivideSquare(2);
    c0->cd(1);
    h2d.Merge()->DrawClone("colz");
    c0->cd(2);
    hGATCONF->DrawClone();

    // Print report
    std::cout << "===== GATCONF report =====" << '\n';
    std::cout << "-> CFA/div = " << cfa << '\n';
    std::cout << "==========================" << '\n';
}

#include "ActCluster.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TMath.h"

#include <fstream>

void getDecayEvts()
// This macro assumes that data was processed with user filter [CorrelateImplantDecay] enabled and
// "EnableDeleteInvalidCluster" of [FindRP] set to false!
{
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EFilter};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain3.get());

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    std::ofstream streamer {"../Outputs/correlatedDecays20Na.txt"};
    df.Foreach(
        [&](ActRoot::MergerData& mer, ActRoot::TPCData& tpc)
        {
            auto& clusters = tpc.fClusters;
            if(!clusters.empty())
                mer.Stream(streamer);
        },
        {"MergerData", "TPCData"});
}
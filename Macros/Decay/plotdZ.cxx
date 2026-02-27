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

// This macro assumes that data was processed with user filter [CorrelateImplantDecay] enabled and
// "EnableDeleteInvalidCluster" of [FindRP] set to false!
void plotdZ()
{
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EFilter};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain3.get());

    // grab drift factor, we'll need it for conversions
    ActRoot::InputParser parserDet {"../../configs/detector.conf"};
    auto bl1 {parserDet.GetBlock("Merger")};
    auto drift {bl1->GetDouble("DriftFactor")};
    std::cout << drift << std::endl;

    ActRoot::InputParser parserMA {"../../configs/multiaction.conf"};
    auto bl2 {parserMA.GetBlock("CorrelateImplantDecay")};
    auto minLxy {bl2->GetDouble("MinLength")};
    std::cout<<minLxy<<std::endl;

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d {*chain};

    auto df {d.Filter(
        [](ActRoot::TPCData& tpc)
        {
            auto& clusters = tpc.fClusters;
            if(!clusters.empty())
                return true;
            return false;
        },
        {"TPCData"})};

    auto dff =
        df.Define("beamLastPoint",
                 [](ActRoot::TPCData& tpc)
                 {
                     ROOT::Math::XYZPointF lastPoint {-1000, -1000, -1000};
                     auto& clusters = tpc.fClusters;
                     ActRoot::Cluster beamLike;

                     for(auto& c : clusters)
                     {
                         if(c.GetIsBeamLike())
                             beamLike = c;
                     }
                     auto line {beamLike.GetRefToLine()};
                     auto dir {line.GetDirection()};
                     beamLike.SortAlongDir(dir);
                     auto lastVoxel {beamLike.GetRefToVoxels().back()};
                     lastPoint = line.ProjectionPointOnLine(lastVoxel.GetPosition());

                     return lastPoint;
                 },
                 {"TPCData"})
            .Define("decayFirstPoint",
                    [&](ActRoot::TPCData& tpc, ROOT::Math::XYZPointF beamPoint)
                    {
                        ROOT::Math::XYZPointF firstPoint {-1000, -1000, -1000};
                        auto& clusters = tpc.fClusters;
                        bool hasDecay {false};

                        for(auto& c : clusters)
                        {
                            if(!c.GetIsBeamLike())
                            {
                                ActRoot::Cluster thisCluster = c;
                                auto line {thisCluster.GetRefToLine()};
                                auto dir {line.GetDirection()};
                                thisCluster.SortAlongDir(dir);
                                const auto& vs = thisCluster.GetRefToVoxels();
                                if(!vs.empty())
                                {
                                    auto firstVoxel {vs.front()};
                                    auto projectionPointLine {line.ProjectionPointOnLine(firstVoxel.GetPosition())};
                                    double lxy = TMath::Sqrt(TMath::Power(projectionPointLine.X() - beamPoint.X(), 2) +
                                                             TMath::Power(projectionPointLine.Y() - beamPoint.Y(), 2));

                                    if(lxy <= minLxy)
                                    {
                                        hasDecay = true;
                                        firstPoint = projectionPointLine;
                                    }
                                }
                            }
                        }
                        if (!hasDecay){
                            throw std::runtime_error("I'm confusion, I didn't find a decay. Is your filter working??");
                        }
                        return firstPoint;
                    },
                    {"TPCData", "beamLastPoint"})
            .Define("dZ",
                    [&](ROOT::Math::XYZPointF beamPoint, ROOT::Math::XYZPointF decayPoint)
                    {
                        double dZ = beamPoint.Z() - decayPoint.Z();
                        // std::cout << dZ << std::endl;
                        return dZ;
                    },
                    {"beamLastPoint", "decayFirstPoint"});


    auto hdZ = dff.Histo1D({"hdZ", "dZ;dZ [btb];", 1000, -500, 500}, "dZ");
    auto* c {new TCanvas("c", "c", 800, 600)};
    hdZ->DrawClone("hist");
}
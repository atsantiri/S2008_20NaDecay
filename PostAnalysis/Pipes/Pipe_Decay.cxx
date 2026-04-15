#include "ActCluster.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TROOT.h"

#include <atomic>
#include <utility>

#include "../../Macros/Decay/miscFunctions.cxx"

// For decay run 129. Assuming data processed with user action [FindDecays] enabled and
// "EnableDeleteInvalidCluster" of [FindRP] set to false


// functions commented out to import from file
// double calcDistXY(ROOT::Math::XYZPointF A, ROOT::Math::XYZPointF B)
// {
//     return TMath::Sqrt(TMath::Power(A.X() - B.X(), 2) + TMath::Power(A.Y() - B.Y(), 2));
// }

// ROOT::Math::XYZPointF
// findStartVoxel(ActRoot::Cluster& c, bool returnStart) // compute the distance to the center mass of the charge. The one
//                                                       // closer to the center mass is the end point
// {
//     float qtot = 0;
//     float xcm = 0., ycm = 0., zcm = 0.;

//     for(const auto& v : c.GetVoxels())
//     {
//         ROOT::Math::XYZPointF pos = v.GetPosition();
//         auto q {v.GetCharge()};

//         qtot += q;
//         xcm += q * pos.X();
//         ycm += q * pos.Y();
//         zcm += q * pos.Z();
//     }

//     if(qtot == 0)
//         return {};
//     xcm /= qtot;
//     ycm /= qtot;
//     zcm /= qtot;
//     ROOT::Math::XYZPointF rcm(xcm, ycm, zcm);
//     ROOT::Math::XYZPointF A = c.GetRefToVoxels().front().GetPosition();
//     ROOT::Math::XYZPointF B = c.GetRefToVoxels().back().GetPosition();
//     float dA = (A - rcm).R();
//     float dB = (B - rcm).R();

//     ROOT::Math::XYZPointF start = (dA > dB) ? A : B;
//     ROOT::Math::XYZPointF end = (dA > dB) ? B : A;

//     return returnStart ? start : end;
// }

void Pipe_Decay()
{
    std::string dataconf {"./../configs/data.conf"};

    ROOT::EnableImplicitMT();

    // Read data
    ActRoot::DataManager datman {dataconf, ActRoot::ModeType::EReadSilMod};
    auto chain {datman.GetJoinedData()};
    auto chain2 {datman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain3 {datman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());
    ROOT::RDataFrame d {*chain};

    // Parameters for event selection
    ActRoot::InputParser parserMultiAction {"./../configs/multiaction.conf"};
    auto bl {parserMultiAction.GetBlock("FindDecays")};
    // ROOT::Math::XYZPointF beamEndPoint {bl->GetDouble("beamEndPointX"), bl->GetDouble("beamEndPointY"),
    // bl->GetDouble("beamEndPointZ")}; // visually inspected output of Pipe0
    ROOT::Math::XYZPointF beamEndPoint {static_cast<float>(bl->GetDouble("beamEndPointX")),
                                        static_cast<float>(bl->GetDouble("beamEndPointY")),
                                        static_cast<float>(bl->GetDouble("beamEndPointZ"))};
    double maxDist {bl->GetDouble("maxDistBeamEnd")}; // max distance from beam stop point
    double edge {bl->GetDouble("maxDistEdge")}; // min distance from x wall to remove tracks that out of the pad plane
    double forwardAngle {bl->GetDouble("forwardAngle")}; // angle of forward angle cone to be removed
    ROOT::Math::XYZVectorF beam = {1, 0, 0};
    std::cout << "Input from filter function:" << std::endl;
    std::cout << "BeamStopPoint:" << bl->GetDouble("beamEndPointX") << ", " << bl->GetDouble("beamEndPointY") << ", "
              << bl->GetDouble("beamEndPointZ") << std::endl;
    std::cout << "Max Dist from Beam Stop: " << maxDist << std::endl;
    std::cout << "Max Dist from Edge: " << edge << std::endl;
    std::cout << "Forward angle cone: " << forwardAngle << std::endl;

    ActRoot::TPCParameters TPCparams {"Actar"};

    ActRoot::InputParser parserDet {"./../configs/detector.conf"};
    auto bl1 {parserDet.GetBlock("Merger")};
    auto drift {bl1->GetDouble("DriftFactor")};
    std::cout << "Drift velocity: " << drift << std::endl;

    // Define lambda function
    auto lambdaIsDecay {[&](const ActRoot::Cluster& c)
                        {
                            if(!c.GetIsBeamLike())
                            {
                                ActRoot::Cluster thisCluster = c; // can't sort a const object
                                auto line {thisCluster.GetRefToLine()};
                                auto dir {line.GetDirection()};
                                thisCluster.SortAlongDir(dir);

                                auto startPos = line.ProjectionPointOnLine(findStartVoxel(thisCluster, true));
                                auto endPos = line.ProjectionPointOnLine(findStartVoxel(thisCluster, false));

                                double lxy = calcDistXY(startPos, beamEndPoint);

                                // remove tracks that don't stop in the pad plane, note this will affect relative peak
                                // intensities
                                bool isNearEdge {(startPos.Z() <= edge) || (startPos.Z() >= (272 - edge)) ||
                                                 (endPos.Z() <= edge) || (endPos.Z() >= (272 - edge))};

                                // alternatively I can remove a cone for forward angles up to 65 deg
                                auto dot {beam.Unit().Dot(dir.Unit())};
                                bool forwardCone = (dot > std::cos(forwardAngle * TMath::DegToRad())) &&
                                                   (endPos.X() > beamEndPoint.X());

                                if((lxy <= maxDist) && !forwardCone && !isNearEdge)
                                    return true;
                                else
                                    return false;
                            }
                            else
                                return false;
                        }};

    // Histograms
    auto hPad = new TH2F("hPad", "Pad;X [pad];Y [pad]", TPCparams.GetNPADSX(), 0, TPCparams.GetNPADSX(),
                         TPCparams.GetNPADSY(), 0, TPCparams.GetNPADSY());
    auto hSide = new TH2F("hSide", "Side;X [pad];Z [tb]", TPCparams.GetNPADSX(), 0, TPCparams.GetNPADSX(),
                          TPCparams.GetNPADSZ() / 4, 0, TPCparams.GetNPADSZ());

    auto df = d.Define("isDecay",
                       [&](const ActRoot::TPCData& tpc)
                       {
                           auto& clusters = tpc.fClusters;
                           std::vector<bool> ret(clusters.size(), false);
                           for(auto& c : clusters)
                               ret[c.GetClusterID()] = lambdaIsDecay(c);

                           return ret;
                       },
                       {"TPCData"});
    int decaysTot {0};
    df.Foreach(
        [&](ActRoot::TPCData& tpc, std::vector<bool> isDecay)
        {
            auto& clusters = tpc.fClusters;
            for(auto& c : clusters)
            {
                if(isDecay[c.GetClusterID()])
                {
                    decaysTot++;
                    for(const auto& v : c.GetVoxels())
                    {
                        auto& pos {v.GetPosition()};
                        hPad->Fill(pos.X(), pos.Y());
                        hSide->Fill(pos.X(), pos.Z() * drift);
                    }
                }
            }
        },
        {"TPCData", "isDecay"});

    std::cout << "Found " << decaysTot << " total number of alpha decays" << std::endl;

    auto* c0 {new TCanvas {"c0", "Decay Pipe canvas 0", 1500, 600}};
    c0->Divide(2);
    c0->cd(1);
    hPad->DrawClone("colz");
    c0->cd(2);
    hSide->DrawClone("colz");

    // save decays
    auto name {"./Outputs/tree_20Na_decays.root"};
    std::cout << "Saving Decay_Tree in file : " << name << '\n';
    df.Snapshot("Decay_Tree", name, {"TPCData", "isDecay"});
}

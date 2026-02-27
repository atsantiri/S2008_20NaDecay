#include "CorrelateImplantDecay.h"

#include "ActCluster.h"
#include "ActColors.h"
#include "ActInputParser.h"
#include "ActTPCData.h"

#include "TMath.h"

#include <memory>

// AT Feb 2026
// Run 129 was 20Na decay only. We want to find events where the implant of 20Na and its decay are in the same event
// window to try and extract a halflife. If enabled, make sure that "EnableDeleteInvalidCluster" of [findRP] is set to
// false in multiaction.conf. We don't want to delete the events with invalid RP, these are the ones we're interested in
// :)

void ActAlgorithm::CorrelateImplantDecay::CorrelateImplantDecay::ReadConfiguration(
    std::shared_ptr<ActRoot::InputBlock> block)
{
    fIsEnabled = block->GetBool("IsEnabled");
    if(!fIsEnabled)
        return;
    if(block->CheckTokenExists("MinLength"))
        fMinLength = block->GetDouble("MinLength");
}

void ActAlgorithm::CorrelateImplantDecay::Run()
{
    if(!fIsEnabled)
        return;

    if(fIsVerbose)
        std::cerr << BOLDMAGENTA << "==================== CorrelateImplantDecay ====================" << RESET << '\n';

    if(!fTPCData)
    {
        if(fIsVerbose)
            std::cerr << BOLDMAGENTA << "Nope: fTPCData is null" << RESET << '\n';
        return;
    }

    // copy/ref clusters
    auto& clusters = fTPCData->fClusters;
    if(clusters.empty())
    {
        if(fIsVerbose)
            std::cout << BOLDMAGENTA << "Nope: no clusters" << RESET << '\n';
        clusters.clear();
        return;
    }

    if(clusters.size() == 1)
    {
        if(fIsVerbose)
            std::cout << BOLDMAGENTA << "Nope: only one cluster" << RESET << '\n';
        clusters.clear();
        return;
    }

    ActRoot::Cluster beamLike;
    int nOfbeamLikes = 0;
    for(const auto& c : clusters)
    {
        if(c.GetIsBeamLike())
        {
            beamLike = c;
            nOfbeamLikes++;
        }
    }

    if(nOfbeamLikes == 0)
    {
        if(fIsVerbose)
            std::cout << BOLDMAGENTA << "Nope: no BeamLike found" << RESET << '\n';
        clusters.clear();
        return;
    }
    else if(nOfbeamLikes > 1)
    {
        if(fIsVerbose)
            std::cout << BOLDMAGENTA << "Nope: multiple BeamLike tracks" << RESET << '\n';
        clusters.clear();
        return;
    }

    auto beamLine {beamLike.GetRefToLine()};
    auto beamDir {beamLine.GetDirection()};
    beamLike.SortAlongDir(beamDir);

    // guard voxels
    if(beamLike.GetRefToVoxels().empty())
    {
        if(fIsVerbose)
            std::cout << BOLDMAGENTA << "Nope: selected beam cluster has no voxels" << RESET << '\n';
        clusters.clear();
        return;
    }

    // verify there's no RP
    if(fTPCData->fRPs.size() > 0)
    {
        if(fIsVerbose)
            std::cout << BOLDMAGENTA << "Nope: There are RPs" << RESET << '\n';
        clusters.clear();
        return;
    }

    // need last point of the beam
    auto beamLastVoxel {beamLike.GetRefToVoxels().back()};
    auto beamEndPoint {beamLine.ProjectionPointOnLine(beamLastVoxel.GetPosition())};

    // find decay product
    ActRoot::Cluster decay;
    ActRoot::Line::XYZPointF decayStartPoint;

    bool hasDecay {false};
    for(const auto& c : clusters)
    {
        int idx = 0;
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
                double lxy = TMath::Sqrt(TMath::Power(projectionPointLine.X() - beamEndPoint.X(), 2) +
                                         TMath::Power(projectionPointLine.Y() - beamEndPoint.Y(), 2));
                if(fIsVerbose)
                {
                    std::cout << BOLDMAGENTA << "-------- Cluster " << idx << '\n';
                    std::cout << "Distance from BeamLike end Lxy: " << lxy << RESET << '\n';
                }
                if(lxy <= fMinLength)
                {
                    decay = thisCluster;
                    decayStartPoint = projectionPointLine;
                    hasDecay = true;
                }
            }
        }
        idx++;
    }
    if(hasDecay)
    {
        if(fIsVerbose)
        {
            std::cout << BOLDMAGENTA << " -------- I made it! ----- \n";
            std::cout << "Lxy Implant/Decay: "
                      << TMath::Sqrt(TMath::Power(decayStartPoint.X() - beamEndPoint.X(), 2) +
                                     TMath::Power(decayStartPoint.Y() - beamEndPoint.Y(), 2))
                      << " pixels \n";
            std::cout << "Beam-like Z: " << beamEndPoint.Z() << '\n';
            std::cout << "Z decay: " << decayStartPoint.Z() << '\n';
            std::cout << "Delta Z: " << beamEndPoint.Z() - decayStartPoint.Z() << RESET << '\n';
        }
    }
    else
    {
        if(fIsVerbose)
            std::cout << BOLDMAGENTA << "Nope: didn't find decay" << RESET << '\n';
        clusters.clear();
    }
    if(fIsVerbose)
        std::cerr << BOLDMAGENTA << "===============================================================" << RESET << '\n';
}

void ActAlgorithm::CorrelateImplantDecay::Print() const
{
    std::cout << BOLDCYAN << "····· " << GetActionID() << " ·····" << '\n';
    if(!fIsEnabled)
    {
        std::cout << "······························" << RESET << '\n';
        return;
    }
    std::cout << "  MinLength      : " << fMinLength << '\n';
}

// Create symbol to load class from .so
extern "C" ActAlgorithm::CorrelateImplantDecay* CreateUserAction()
{
    return new ActAlgorithm::CorrelateImplantDecay;
}

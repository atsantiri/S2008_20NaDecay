#include "FindDecays.h"

#include "ActCluster.h"
#include "ActColors.h"
#include "ActInputParser.h"
#include "ActTPCData.h"

#include "TMath.h"

#include <memory>

// AT 2026
// Run 129 was 20Na decay only. We want to find decay events that deposit all their energy inside the pad plane. If
// enabled, make sure that "EnableDeleteInvalidCluster" of [findRP] is set to false in multiaction.conf. We don't want
// to delete the events with invalid RP, these are the ones we're interested in
// :)
double calcDistXY(ROOT::Math::XYZPointF A, ROOT::Math::XYZPointF B)
{
    return TMath::Sqrt(TMath::Power(A.X() - B.X(), 2) + TMath::Power(A.Y() - B.Y(), 2));
}

ROOT::Math::XYZPointF
findStartVoxel(ActRoot::Cluster& c, bool returnStart) // compute the distance to the center mass of the charge. The one
                                                      // closer to the center mass is the end point
{
    float qtot = 0;
    float xcm = 0., ycm = 0., zcm = 0.;

    for(const auto& v : c.GetVoxels())
    {
        ROOT::Math::XYZPointF pos = v.GetPosition();
        auto q {v.GetCharge()};

        qtot += q;
        xcm += q * pos.X();
        ycm += q * pos.Y();
        zcm += q * pos.Z();
    }

    if(qtot == 0)
        return {};
    xcm /= qtot;
    ycm /= qtot;
    zcm /= qtot;
    ROOT::Math::XYZPointF rcm(xcm, ycm, zcm);
    ROOT::Math::XYZPointF A = c.GetRefToVoxels().front().GetPosition();
    ROOT::Math::XYZPointF B = c.GetRefToVoxels().back().GetPosition();
    float dA = (A - rcm).R();
    float dB = (B - rcm).R();

    ROOT::Math::XYZPointF start = (dA > dB) ? A : B;
    ROOT::Math::XYZPointF end = (dA > dB) ? B : A;

    return returnStart ? start : end;
}

void ActAlgorithm::FindDecays::FindDecays::ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block)
{
    fIsEnabled = block->GetBool("IsEnabled");
    if(!fIsEnabled)
        return;
    if(block->CheckTokenExists("beamEndPointX"))
        beamEndPointX = block->GetDouble("beamEndPointX");
    if(block->CheckTokenExists("beamEndPointY"))
        beamEndPointY = block->GetDouble("beamEndPointY");
    if(block->CheckTokenExists("beamEndPointZ"))
        beamEndPointZ = block->GetDouble("beamEndPointZ");
    if(block->CheckTokenExists("maxDistBeamEnd"))
        maxDistBeamEnd = block->GetDouble("maxDistBeamEnd");
    if(block->CheckTokenExists("maxDistEdge"))
        maxDistEdge = block->GetDouble("maxDistEdge");
    if(block->CheckTokenExists("forwardAngle"))
        forwardAngle = block->GetDouble("forwardAngle");
}

void ActAlgorithm::FindDecays::Run()
{
    if(!fIsEnabled)
        return;

    if(fIsVerbose)
        std::cerr << BOLDMAGENTA << "==================== FindDecays ====================" << RESET << '\n';

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

    // verify there's no RP
    if(fTPCData->fRPs.size() > 0)
    {
        if(fIsVerbose)
            std::cout << BOLDMAGENTA << "Nope: There are RPs" << RESET << '\n';
        clusters.clear();
        return;
    }

    // find decay product
    ActRoot::Cluster decay;
    int decayIdx {};
    ActRoot::Line::XYZPointF decayStartPoint;
    ROOT::Math::XYZPointF beamEndPoint {beamEndPointX, beamEndPointY, beamEndPointZ}; // implantation point
    ROOT::Math::XYZVectorF beam = {1, 0, 0};

    bool hasDecay {false};
    for(const auto& c : clusters)
    {
        auto idx {c.GetClusterID()};
        if(!c.GetIsBeamLike())
        {
            ActRoot::Cluster thisCluster = c;
            auto line {thisCluster.GetRefToLine()};
            auto dir {line.GetDirection()};
            thisCluster.SortAlongDir(dir);

            auto startPos = line.ProjectionPointOnLine(findStartVoxel(thisCluster, true));
            auto endPos = line.ProjectionPointOnLine(findStartVoxel(thisCluster, false));

            double lxy = calcDistXY(startPos, beamEndPoint);

            // remove tracks that come from pad plane, note this will affect relative peak
            // intensities
            bool isNearEdge {(startPos.Z() <= maxDistEdge) || (startPos.Z() >= (272 - maxDistEdge)) ||
                             (endPos.Z() <= maxDistEdge) || (endPos.Z() >= (272 - maxDistEdge))};

            // remove a forward solid angle of forwardangle deg, for events that make it out of the pad plane
            double cosAngle = beam.Unit().Dot(dir.Unit());
            bool forwardCone =
                (cosAngle > std::cos(forwardAngle * TMath::DegToRad())) && (endPos.X() > beamEndPoint.X());
            if(fIsVerbose)
            {
                std::cout << BOLDMAGENTA << "-------- Cluster " << idx << '\n';
                std::cout << "Distance from BeamLike end Lxy: " << lxy << '\n';
                std::cout << "Is in forward cone: " << forwardCone << '\n';
                std::cout << "Is near edge: " << isNearEdge << RESET << '\n';
            }
            if((lxy <= maxDistBeamEnd) && !forwardCone && !isNearEdge)
            {
                decay = thisCluster;
                decayStartPoint = startPos;
                hasDecay = true;
                decayIdx = idx;
            }
        }
    }
    if(hasDecay)
    {
        if(fIsVerbose)
        {
            std::cout << BOLDMAGENTA << " -------- I made it! ----- \n";
            std::cout << "Cluster " << decayIdx << RESET << '\n';
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

void ActAlgorithm::FindDecays::Print() const
{
    std::cout << BOLDCYAN << "····· " << GetActionID() << " ·····" << '\n';
    if(!fIsEnabled)
    {
        std::cout << "······························" << RESET << '\n';
        return;
    }
    std::cout << "  BeamEndPoint    : " << beamEndPointX << ", " << beamEndPointY << ", " << beamEndPointZ << '\n';
    std::cout << "  maxDistBeamEnd  : " << maxDistBeamEnd << '\n';
    std::cout << "  maxDistEdge     : " << maxDistEdge << '\n';
    std::cout << "  ForwardAngle    : " << forwardAngle << '\n';
}

// Create symbol to load class from .so
extern "C" ActAlgorithm::FindDecays* CreateUserAction()
{
    return new ActAlgorithm::FindDecays;
}

#include "ActVAction.h"

namespace ActAlgorithm
{
class FindDecays : public VAction
{
public:
    float beamEndPointX {};
    float beamEndPointY {};
    float beamEndPointZ {};
    double maxDistBeamEnd {};
    double maxDistEdge {};
    double forwardAngle {};
    // double drift {};

public:
    FindDecays() : VAction("FindDecays") {}

    void ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block) override;
    void Run() override;
    void Print() const override;
};
} // namespace ActAlgorithm

#include "ActVAction.h"

namespace ActAlgorithm
{
class CorrelateImplantDecay : public VAction
{
public:
    double fMinLength {}; //!< Min Lxy between implant and decay.

public:
    CorrelateImplantDecay() : VAction("CorrelateImplantDecay") {}

    void ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block) override;
    void Run() override;
    void Print() const override;
};
} // namespace ActAlgorithm

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

#include "TLine.h"
#include "TMath.h"
#include <TF1.h>

#include <fstream>

// Scale, ScalePoint and calcTl functions from ActMergerDetector and ActLine classes
void ActRoot::Line::Scale(float xy, float z)
{
    // Point
    fPoint.SetX(fPoint.X() * xy);
    fPoint.SetY(fPoint.Y() * xy);
    fPoint.SetZ(fPoint.Z() * z);
    // Direction
    fDirection.SetX(fDirection.X() * xy);
    fDirection.SetY(fDirection.Y() * xy);
    fDirection.SetZ(fDirection.Z() * z);
}

void ScalePoint(ROOT::Math::XYZPointF& point, float xy, float z)
{
    point += ROOT::Math::XYZVector {0.5, 0.5, 0.5};
    point.SetX(point.X() * xy);
    point.SetY(point.Y() * xy);
    point.SetZ(point.Z() * z);
}

double calcTLfromVoxel(ROOT::Math::XYZPointF A, ROOT::Math::XYZPointF B, ActRoot::Line line, double drift)
{
    ActRoot::TPCParameters params;
    double xy = params.GetPadSide();

    ScalePoint(A, xy, drift);
    ScalePoint(B, xy, drift);
    line.Scale(xy, drift);

    auto projBegin {line.ProjectionPointOnLine(A)};
    auto projEnd {line.ProjectionPointOnLine(B)};
    return (projBegin - projEnd).R();
}

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

std::vector<std::pair<double, double>> getChargeProf(ActRoot::Cluster& c, ROOT::Math::XYZPointF start,
                                                     ROOT::Math::XYZPointF end, ActRoot::Line::XYZVectorF dir,
                                                     double drift)
{
    ActRoot::TPCParameters params;
    double xy = params.GetPadSide();

    ScalePoint(start, xy, drift);
    ScalePoint(end, xy, drift);

    std::vector<std::pair<double, double>> sq;
    if((end - start).Dot(dir) < 0) // ensure positive dir
        dir *= -1.0;


    for(const auto& v : c.GetVoxels())
    {
        ROOT::Math::XYZPointF pos = v.GetPosition();
        double q = v.GetCharge();
        ScalePoint(pos, xy, drift);
        double s = (pos - start).Dot(dir.Unit());
        sq.emplace_back(s, q);
    }
    std::sort(sq.begin(), sq.end());
    return sq;
}

TSpline3* buildSRIMspline(
    ActPhysics::SRIM* srim, double range, const std::string& particleKey, double step = 0.5,
    double sOffset = 0.0) // function adapted from Ivan's fitTPCevent.cxx, but instead of dE here we do dE/dx to remove
                          // step size dependence
                          // Source:https://github.com/ibc12/S2384/blob/main/Simulation/Macros/fitTPCevent.cxx
{
    std::vector<double> s_pts;
    std::vector<double> y_pts;

    for(double s = 0; s < range; s += step)
    {
        double R = range - s; // remaining range

        if(R <= 0)
            break;

        double E = srim->EvalInverse(particleKey, R);
        double dEdx = srim->EvalStoppingPower(particleKey, E);

        if(dEdx < 0)
            dEdx = 0;

        s_pts.push_back(s + 0.5 * step + sOffset);
        y_pts.push_back(dEdx);
    }
    //  ensure endpoint goes to zero
    if(!y_pts.empty() && y_pts.back() > 0)
    {
        s_pts.push_back(range);
        y_pts.push_back(0.0);
    }

    if(s_pts.size() < 3)
    {
        std::cout << "Not enough points to build spline (nSteps=" << s_pts.size() << "). Need at least 3 points.\n";
        return nullptr;
    }

    auto* sp = new TSpline3(("spSRIM_" + particleKey).c_str(), s_pts.data(), y_pts.data(), s_pts.size());

    sp->SetNpx(3000);
    std::cout << "Spline max: " << sp->GetXmax() << "\n";
    std::cout << "Value near max: " << sp->Eval(range - 1) << "\n";
    return sp;
}

double FindPositionFromChargeFraction(
    TH1* h, double frac) // function from Ivan's fitTPCevent.cxx,
                         // Source:https://github.com/ibc12/S2384/blob/main/Simulation/Macros/fitTPCevent.cxx
{
    double total = h->Integral("width");
    double accum = 0.0;

    for(int i = 1; i <= h->GetNbinsX(); ++i)
    {
        accum += h->GetBinContent(i) * h->GetBinWidth(i);

        if(accum / total >= frac)
            return h->GetBinCenter(i);
    }

    return h->GetXaxis()->GetXmax();
}

TF1* FitSRIMtoChargeProfile(TH1* hcharge, TSpline3* srimSpline, TString fname) // a bit different from Ivan's version
{
    if(!srimSpline || !hcharge)
        return nullptr;

    double xmin = FindPositionFromChargeFraction(hcharge, 0.01);
    double xmax = FindPositionFromChargeFraction(hcharge, 1.);
    double xmaxSpline = srimSpline->GetXmax();

    // std::cout << "Fitting range: " << xmin << " " << xmax << std::endl;

    TF1* fit = new TF1(
        fname,
        [=](double* x, double* p)
        {
            double xs = x[0];
            double xShifted = xs + p[1]; // match endpoint of spline to charge profile

            if(xShifted < srimSpline->GetXmin() || xShifted > srimSpline->GetXmax())
                return 0.0;

            double val = srimSpline->Eval(xShifted);

            if(val < 0)
                val = 0;

            return p[0] * val;
        },
        xmin, xmax, 2);

    fit->SetParameters(hcharge->Integral("width") / hcharge->GetNbinsX(), // amplitude
                       xmaxSpline - xmax);                                // shift

    fit->SetParName(0, "amplitude");
    fit->SetParName(1, "shift");
    fit->SetParLimits(0, 0, 1e8);
    fit->SetParLimits(1, -150, 150);

    fit->SetNpx(3000);

    return fit;
}

double
calcTLfromChargeFit(TF1* f, double xMin, double xMax) // find location where fit curve drops to 1/5 of the Bragg peak
{
    if(!f)
        return -1;
    int nStepsMax = 3000;
    double dx = (xMax - xMin) / nStepsMax;

    // find maximum - Bragg peak
    double xPeak = xMin;
    double yPeak = f->Eval(xMin);
    for(int i = 1; i <= nStepsMax; ++i)
    {
        double x = xMin + i * dx;
        double y = f->Eval(x);
        if(y > yPeak)
        {
            yPeak = y;
            xPeak = x;
        }
    }

    // find point after peak where height drops below 1/5
    double threshold = yPeak / 5.0;
    int nStepsTail = 1000;
    dx = (xMax - xPeak) / nStepsTail;

    for(int i = 1; i <= nStepsTail; ++i)
    {
        double x = xPeak + i * dx;
        double y = f->Eval(x);
        if(y <= threshold)
        {
            // Linear interp
            double x0 = x - dx;
            double y0 = f->Eval(x0);
            if(y0 == y)
                return x;
            double xInterp = x0 + (threshold - y0) * (x - x0) / (y - y0);
            return xInterp;
        }
    }

    return -1;
}

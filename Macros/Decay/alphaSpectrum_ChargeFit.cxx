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

#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include <TF1.h>

#include <fstream>
// run 129 is processed without any user functions enabled

// Scale, ScalePoint and calcTl function from ActMergerDetector and ActLine classes
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

double calcTL(ROOT::Math::XYZPointF A, ROOT::Math::XYZPointF B, ActRoot::Line line, double drift)
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

TF1* FitSRIMtoChargeProfile(TH1* hcharge, TSpline3* srimSpline, TString fname)
{
    if(!srimSpline || !hcharge)
        return nullptr;

    double xmin = FindPositionFromChargeFraction(hcharge, 0.01);
    double xmax = FindPositionFromChargeFraction(hcharge, 0.999);
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
            // Linear interpolation
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

void alphaSpectrum_ChargeFit()
{
    // ========================================================================================================================================
    ROOT::Math::XYZPointF beamEndPoint {105, 65, 0}; // visually inspected output of Pipe0
    double maxDist {2.};                             // max distance from beam stop point
    double edge {5.};                                // min distance from Z wall to remove tracks
                                                     // that come from the pad plane - 20Na drifted before decaying
    bool debug {false}; // set to true if we're testing. Otherwise we will do only few events
    int kmax {15};     // how many events to look at
    std::string srimfile {"../../Calibrations/SRIM/4He_950mbar_95-5.txt"};    
    double qval {-4.73};
    double ma {4.0026};      // amu
    double mo16 {15.994914}; // amu
    // ========================================================================================================================================

    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EFilter};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain3.get());

    ROOT::RDataFrame d {*chain};

    ActRoot::TPCParameters TPCparams {"Actar"};
    ActRoot::InputParser parserDet {"../../configs/detector.conf"};
    auto bl1 {parserDet.GetBlock("Merger")};
    auto drift {bl1->GetDouble("DriftFactor")};
    std::cout << drift << std::endl;

    auto hPad = new TH2F("hPad", "Pad;X [pad];Y [pad]", TPCparams.GetNPADSX(), 0, TPCparams.GetNPADSX(),
                         TPCparams.GetNPADSY(), 0, TPCparams.GetNPADSY());

    auto hSide = new TH2F("hSide", "Side;X [pad];Z [tb]", TPCparams.GetNPADSX(), 0, TPCparams.GetNPADSX(),
                          TPCparams.GetNPADSZ() / 4, 0, TPCparams.GetNPADSZ());

    auto hEne = new TH1F("hEne", "Alpha Energy; E_{#alpha} [MeV]; ", 200, 0, 7);
    auto hEne2 = new TH1F("hEneVoxel", "Alpha Energy; E_{#alpha} [MeV]; ", 200, 0, 7);
    auto hEx = new TH1F("hEx", "Excitation Energy; E_{20Ne*} [MeV]; ", 200, 4, 14);
    auto hEx2 = new TH1F("hExVoxel", "Excitation Energy; E_{20Ne*} [MeV]; ", 200, 4, 14);
    std::cout << TPCparams.GetNPADSZ() << " " << TPCparams.GetREBINZ() << std::endl;

    auto* srim {new ActPhysics::SRIM()};
    const std::string& particleKey = "HeInGas";
    srim->ReadTable(particleKey, srimfile);
    auto srimSpline {buildSRIMspline(srim, 100, particleKey)};

    ROOT::Math::XYZVectorF beam = {1, 0, 0};

    std::vector<std::shared_ptr<TPolyLine>> ls;
    std::vector<double> xsA, ysA, xsB, ysB;
    using SQ = std::vector<std::pair<double, double>>;
    std::vector<SQ> allChargeProfiles;
    std::vector<double> TLsfromVoxels;
    int k = 0;

    d.Foreach(
        [&](ActRoot::TPCData& tpc)
        {
            if(tpc.fRPs.empty())
            {
                auto& clusters = tpc.fClusters;
                for(auto& c : clusters)
                {
                    if((k < kmax) || !debug)
                    {
                        if(!c.GetIsBeamLike())
                        {
                            auto line {c.GetRefToLine()};
                            auto dir {line.GetDirection()};
                            c.SortAlongDir();

                            auto startVoxelPos = findStartVoxel(c, true);
                            auto endVoxelPos = findStartVoxel(c, false);

                            double lxy = calcDistXY(startVoxelPos, beamEndPoint);

                            // remove tracks that come from pad plane, note this will affect relative peak
                            // intensities
                            bool isNearEdge {
                                (startVoxelPos.Z() * drift <= edge) || (startVoxelPos.Z() * drift >= (272 - edge)) ||
                                (endVoxelPos.Z() * drift <= edge) || (endVoxelPos.Z() * drift >= (272 - edge))};

                            // remove a forward solid angle of 65 deg, for events that make it out of the pad plane
                            double cosAngle = beam.Unit().Dot(dir.Unit());
                            bool forwardCone =
                                (cosAngle > std::cos(65. * TMath::DegToRad())) && (endVoxelPos.X() > beamEndPoint.X());

                            if((lxy <= maxDist) && !forwardCone && !isNearEdge)
                            {
                                for(const auto& v : c.GetVoxels())
                                {
                                    auto& pos {v.GetPosition()};
                                    hPad->Fill(pos.X(), pos.Y());
                                    hSide->Fill(pos.X(), pos.Z() * drift);
                                }

                                // get charge profile
                                SQ sq = getChargeProf(c, startVoxelPos, endVoxelPos, dir, drift);
                                allChargeProfiles.push_back(std::move(sq));
                                TLsfromVoxels.push_back(calcTL(startVoxelPos, endVoxelPos, line, drift));

                                if(debug)
                                {
                                    xsA.push_back(startVoxelPos.X());
                                    ysA.push_back(startVoxelPos.Y());
                                }
                                k++;
                            }
                        }
                    }
                    // k++;
                }
            }
        },
        {"TPCData"});


    // Draw
    std::cout << "size " << allChargeProfiles.size() << std::endl;
    int badfits {0};

    if(debug)
    {
        auto* ct {new TCanvas {"ct", "Testing", 800, 600}};
        int p {0};
        for(const auto& sq : allChargeProfiles)
        {
            ct->cd();
            auto hcharge = new TH1D {TString::Format("hCharge%d", p), "Charge profile; L [mm]; Q", 50, 0, 150};
            for(size_t i = 0; i < sq.size(); ++i)
                hcharge->Fill(sq[i].first, sq[i].second);
            hcharge->DrawCopy("hist");
            auto fitSpline = FitSRIMtoChargeProfile(hcharge, srimSpline, TString::Format("fSRIM%d", p));
            hcharge->Fit(fitSpline, "QR0", "", hcharge->GetXaxis()->GetXmin(), hcharge->GetXaxis()->GetXmax());
            fitSpline->SetLineWidth(3);
            fitSpline->DrawCopy("same");

            double chi2 = fitSpline->GetChisquare();
            double ndf = fitSpline->GetNDF();
            double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
            std::cout << "\n===== Fit Results: " << fitSpline->GetName() << " =====\n";
            int npar = fitSpline->GetNpar();
            for(int i = 0; i < npar; ++i)
            {
                std::cout << std::setw(10) << fitSpline->GetParName(i) << " = " << std::setw(12)
                          << fitSpline->GetParameter(i) << " ± " << std::setw(12) << fitSpline->GetParError(i) << "\n";
            }
            std::cout << "-----------------------------------\n";
            std::cout << "Chi2 / NDF = " << chi2 << " / " << ndf << " = " << chi2ndf << "\n";


            double TLfromVoxel = TLsfromVoxels[p];
            double TLfromCharge =
                calcTLfromChargeFit(fitSpline, hcharge->GetXaxis()->GetXmin(), hcharge->GetXaxis()->GetXmax());
            double ene {srim->EvalInitialEnergy("HeInGas", 0, TLfromCharge)};
            double diff {(TLfromCharge - TLfromVoxel) / TLfromVoxel * 100.};
            if(std::abs(diff) > 10)
                badfits++;
            std::cout << "TL from charge: " << TLfromCharge << " vs TL from last voxel: " << TLfromVoxel << " ---- ( "
                      << diff << " %)" << std::endl;

            hEne->Fill(ene);
            double Ex {ene * (1 + ma / mo16) - qval};
            hEx->Fill(Ex);


            p++;

            ct->Update();
            ct->WaitPrimitive();
            ct->Clear();
            delete hcharge;
            // delete fitSpline;
        }
        delete ct;
    }
    else
    {
        int p {0};
        for(const auto& sq : allChargeProfiles)
        {
            auto hcharge = new TH1D {TString::Format("hCharge%d", p), "Charge profile; L [mm]; Q", 50, 0, 150};
            hcharge->SetDirectory(nullptr);
            for(size_t i = 0; i < sq.size(); ++i)
                hcharge->Fill(sq[i].first, sq[i].second);
            auto fitSpline = FitSRIMtoChargeProfile(hcharge, srimSpline, TString::Format("fSRIM%d", p));
            hcharge->Fit(fitSpline, "QR0", "", hcharge->GetXaxis()->GetXmin(), hcharge->GetXaxis()->GetXmax());
            double chi2 = fitSpline->GetChisquare();
            double ndf = fitSpline->GetNDF();
            double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
            double TLfromCharge =
                calcTLfromChargeFit(fitSpline, hcharge->GetXaxis()->GetXmin(), hcharge->GetXaxis()->GetXmax());
            double ene {srim->EvalInitialEnergy("HeInGas", 0, TLfromCharge)};
            double TLfromVoxel = TLsfromVoxels[p];
            double ene2 {srim->EvalInitialEnergy("HeInGas", 0, TLfromVoxel)};

            double diff {(TLfromCharge - TLfromVoxel) / TLfromVoxel * 100.};
            if(std::abs(diff) > 10)
                badfits++;

            hEne->Fill(ene);
            double Ex {ene * (1 + ma / mo16) - qval};
            hEx->Fill(Ex);
            hEne2->Fill(ene2);
            double Ex2 {ene2 * (1 + ma / mo16) - qval};
            hEx2->Fill(Ex2);

            p++;
            delete hcharge;
            // delete fitSpline;
        }
    }
    std::cout << "Bad fits: " << badfits << " out of " << allChargeProfiles.size() << std::endl;


    auto* c0 {new TCanvas {"c0", "", 1500, 600}};
    c0->Divide(2);
    c0->cd(1);
    hPad->DrawClone("colz");

    if(debug)
    {
        auto* gBegins = new TGraph(xsA.size(), xsA.data(), ysA.data());
        gBegins->SetMarkerStyle(20);
        gBegins->SetMarkerSize(0.8);
        gBegins->SetMarkerColor(kRed);
        gBegins->Draw("P same");
    }

    c0->cd(2);
    hSide->DrawClone("colz");

    auto* c1 {new TCanvas {"c1", ""}};
    c1->DivideSquare(2);
    c1->cd(1);
    hEne->SetLineWidth(2);
    hEne->DrawClone();
    hEne2->SetLineColor(2);
    hEne2->SetLineWidth(3);
    hEne2->DrawClone("same");
    c1->cd(2);
    hEx->DrawClone();
    hEx2->SetLineColor(2);
    hEx2->DrawClone("same");
    auto l = new TLegend(0.6, 0.8, 0.9, 0.9);
    l->AddEntry(hEx, "from charge fit");
    l->AddEntry(hEx2, "from last voxel");
    l->Draw();

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
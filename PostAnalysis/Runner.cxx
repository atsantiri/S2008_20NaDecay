#include "ActColors.h"

#include "TROOT.h"
#include "TString.h"

#include <iostream>
#include <string>

void Runner(TString what = "")
{
    std::string beam {"20Na"};
    std::string target {"p"};
    std::string light {"p"};

    std::cout << BOLDGREEN << "···· Runner ····" << '\n';
    std::cout << "-> Beam   : " << beam << '\n';
    std::cout << "-> Target : " << target << '\n';
    std::cout << "-> Light  : " << light << '\n';
    std::cout << "-> What   : " << what << '\n';
    std::cout << "······························" << RESET << '\n';

    auto args {TString::Format("(\"%s\", \"%s\", \"%s\")", beam.c_str(), target.c_str(), light.c_str())};
    TString path {"./Pipes/"};
    TString func {};
    TString ext {".cxx"};

    // CFA counter
    if(what.Contains("0"))
    {
        func = "Pipe0_Beam";
        gROOT->LoadMacro(path + func + ext);
        gROOT->ProcessLine(func + "()");
    }
    // PID
    if(what.Contains("1"))
    {
        func = "Pipe1_PID";
        gROOT->LoadMacro(path + func + ext);
        gROOT->ProcessLine(func + args);
    }
    // Kin + Ex
    if(what.Contains("2"))
    {
        func = "Pipe2_Ex";
        gROOT->LoadMacro(path + func + ext);
        gROOT->ProcessLine(func + args);
    }
    if(what.Contains("3"))
    {
        func = "Pipe3_RPCuts";
        gROOT->LoadMacro(path + func + ext);
        gROOT->ProcessLine(func + args);
    }
    // 20Na Decay
        if(what.Contains("d"))
    {
        func = "Pipe_Decay";
        gROOT->LoadMacro(path + func + ext);
        gROOT->ProcessLine(func+"()");
    }
}

// 
// Generate unweighted three-body phasespace decay events
// 
#include "TComplex.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TH2D.h"
#include "TCanvas.h"

#include <iostream>
#include <cmath>

using namespace std;

const double mass_lambdac2625 = 2.62811;
const double mass_pi = 0.13957;
const double mass_lambdac = 2.28646;

const Double_t mass_sigmac0 = 2.45375;
const Double_t gamma_sigmac0 = 0.00183;

const Double_t g13 = 0.0;
const Double_t g5 = 1.0;

TComplex BreitWigner(Double_t mass, Double_t mass0, Double_t gamma)
{
    return TComplex(1) / TComplex(mass - mass0, gamma / 2);
}

Double_t T1Squared(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3)
{
    Double_t m23 = (p2 + p3).M();
    TComplex f1 = g13 * BreitWigner(m23, mass_sigmac0, gamma_sigmac0);

    // cout << "T1 squared: p1.mag = " << p1.Vect().Mag()
    //      << "f1 rho2 = " << f1.Rho2() << endl;

    return 4 / 3 * f1.Rho2() * pow(p1.Vect().Mag(), 4.) * 
        pow(p2.Vect().Mag(), 2.);
}

Double_t T3Squared(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3)
{
    Double_t m13 = (p1 + p3).M();
    TComplex f3 = g13 * BreitWigner(m13, mass_sigmac0, gamma_sigmac0);

    // cout << "T1 squared: p1.mag = " << p1.Vect().Mag()
    //      << "f1 rho2 = " << f1.Rho2() << endl;

    return 4 / 3 * f3.Rho2() * pow(p2.Vect().Mag(), 4.) * 
        pow(p1.Vect().Mag(), 2.);
}

Double_t T5Squared(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3)
{
    TVector3 vec1 = p1.Vect();
    TVector3 vec2 = p2.Vect();

    return 4 / 3 * pow(g5, 2) * (vec1.Mag2() + vec2.Mag2() + 2 * vec1.Dot(vec2));
}

Double_t T15(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3)
{
    Double_t m23 = (p2 + p3).M();
    TComplex f1 = g13 * BreitWigner(m23, mass_sigmac0, gamma_sigmac0);
    TComplex f5(g5);

    TVector3 vec1 = p1.Vect();
    TVector3 vec2 = p2.Vect();

    return 4 * (f1 * f5).Re() * 
        (2 * vec1.Mag2() * vec1.Dot(vec2) + 3 * pow(vec1.Dot(vec2), 2) - vec1.Mag2() * vec2.Mag2());
}

Double_t T35(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3)
{
    Double_t m13 = (p1 + p3).M();
    TComplex f3 = g13 * BreitWigner(m13, mass_sigmac0, gamma_sigmac0);
    TComplex f5(g5);

    TVector3 vec1 = p1.Vect();
    TVector3 vec2 = p2.Vect();

    return 4 * (f3 * f5).Re() * 
        (2 * vec2.Mag2() * vec1.Dot(vec2) + 3 * pow(vec1.Dot(vec2), 2) - vec1.Mag2() * vec2.Mag2());
}

Double_t T13(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3)
{
    Double_t m13 = (p1 + p3).M();
    Double_t m23 = (p2 + p3).M();
    
    TComplex f3 = g13 * BreitWigner(m13, mass_sigmac0, gamma_sigmac0);
    TComplex f1 = g13 * BreitWigner(m23, mass_sigmac0, gamma_sigmac0);
    
    TComplex f5(g5);

    TVector3 vec1 = p1.Vect();
    TVector3 vec2 = p2.Vect();

    return 8/3 * (f3 * f1).Re() * 
       vec1.Dot(vec2) * (3 * pow(vec1.Dot(vec2), 2) - 2 * vec1.Mag2() * vec2.Mag2());
}

Double_t AmpSquared(TLorentzVector rest_p1, TLorentzVector rest_p2, TLorentzVector rest_p3)
{
    TLorentzVector resonance = rest_p2 + rest_p3;
    TLorentzRotation resonanceBoost = -resonance.BoostVector();
    TLorentzVector p1 = resonanceBoost * rest_p1;
    TLorentzVector p2 = resonanceBoost * rest_p2;
    TLorentzVector p3 = resonanceBoost * rest_p3;
    
    return T1Squared(p1, p2, p3) + T3Squared(p1, p2, p3) + T5Squared(p1, p2, p3) + 
        T15(p1, p2, p3) + T35(p1, p2, p3) + T13(p1, p2, p3);
}

void gen_dalitz()
{
    TLorentzVector mom(0, 0, 0, mass_lambdac2625);
    Double_t daug_masses[3] = {mass_pi, mass_pi, mass_lambdac};

    TGenPhaseSpace gen;
    gen.SetDecay(mom, 3, daug_masses);
    double max_weight = gen.GetWtMax();

    

    TH2D * hist = new TH2D("hist", "", 100, 0.075, 0.12, 100, 5.87, 6.20);
    TH1D * hist_amp2 = new TH1D("hist_amp2", "", 100, 0, 10);

    //
    // Find max amplitude squared
    //
    double max_amp2 = 0.04;
    int max_amp2_search_steps = 10000;
    for(int i = 0; i < max_amp2_search_steps; i++)
    {
        gen.Generate();
        TLorentzVector p1 = *gen.GetDecay(0);
        TLorentzVector p2 = *gen.GetDecay(1);
        TLorentzVector p3 = *gen.GetDecay(2);
        double s12 = (p1 + p2).M2();
        double s23 = (p2 + p3).M2();
        double s13 = (p1 + p3).M2();
        double amp2 = AmpSquared(p1, p2, p3);
        max_amp2 = max(amp2, max_amp2);
    }
    cout << "Maximum amplitude square = " << max_amp2 << endl;

    const int maxn = 10000;
    for(int i=0; i < maxn; i++)
    {
        // This part get rid of the event weight
        if(drand48() < gen.Generate() / max_weight)
        {
            i--;
            continue;
        }

        TLorentzVector p1 = *gen.GetDecay(0);
        TLorentzVector p2 = *gen.GetDecay(1);
        TLorentzVector p3 = *gen.GetDecay(2);
        double s12 = (p1 + p2).M2();
        double s23 = (p2 + p3).M2();
        double s13 = (p1 + p3).M2();
        double amp2 = AmpSquared(p1, p2, p3);

        if(drand48() > amp2  / (max_amp2 * 1.1))
        {
            i--;
            continue;    
        }
        
        hist->Fill(s12, s13);

        // std::cout << "T1^2 = " << std::setw(10) << T1Squared(p1, p2, p3)
        //      << "\tT15 = " << T15(p1, p2, p3) << endl;
    }
    // std::cout << "max_gen_amp2 = " << max_gen_amp2 << endl;

    TCanvas * can = new TCanvas("can", "can", 800, 600);
    hist->Draw("");
    hist->SetTitle(TString::Format("g13 = %.2f g5 = %.2f", g13, g5));
    can->Draw();
}
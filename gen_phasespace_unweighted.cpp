// 
// Generate unweighted three-body phasespace decay events
// 
#include "TComplex.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TH2D.h"

#include <iostream>
#include <cmath>

using namespace std;

const double mass_lambdac2625 = 2.62811;
const double mass_pi = 0.13957;
const double mass_lambdac = 2.28646;

void gen_phasespace_unweighted()
{
    TLorentzVector mom(0, 0, 0, mass_lambdac2625);
    Double_t daug_masses[3] = {mass_pi, mass_pi, mass_lambdac};

    TGenPhaseSpace gen;
    gen.SetDecay(mom, 3, daug_masses);
    double max_weight = gen.GetWtMax();

    const int maxn = 10000;

    TH2D * hist = new TH2D("hist", "", 100, 0.075, 0.12, 100, 5.87, 6.20);

    for(int i=0; i < maxn; i++)
    {
        // This part get rid of the event weight
        if(drand48() > gen.Generate() / max_weight)
        {
            i--;
            continue;
        }

        TLorentzVector p1 = *gen.GetDecay(0);
        TLorentzVector p2 = *gen.GetDecay(1);
        TLorentzVector p3 = *gen.GetDecay(2);
        double s13 = (p1 + p3).M2();
        double s12 = (p1 + p2).M2();
        hist->Fill(s12, s13);
    }

    hist->Draw();
}
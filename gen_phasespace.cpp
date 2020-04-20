//
// Generate weighted three-body phasespace decay events
//
#include "TComplex.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

#include <iostream>
#include <cmath>

using namespace std;

const double mass_lambdac2625 = 2.62811;
const double mass_pi = 0.13957;
const double mass_lambdac = 2.28646;

void gen_phasespace()
{
    TLorentzVector mom(0, 0, 0, mass_lambdac2625);
    Double_t daug_masses[3] = {mass_pi, mass_pi, mass_lambdac};

    TGenPhaseSpace gen;
    gen.SetDecay(mom, 3, daug_masses);

    const int maxn = 10;
    for (int i = 0; i < maxn; i++)
    {
        // The weights are NOT constant for the generated events!!!
        cout << "weight = " << gen.Generate() << endl;
        TLorentzVector p1 = *gen.GetDecay(0);
        cout << p1.X() << " " << p1.Y() << " " << p1.Z() << " "
             << p1.T() << " " << p1.M() << endl;
    }
}
#include <iostream>

#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"

#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/IdealGasReactor.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/ReactorSurface.h"

#include "cantera/thermo/SurfPhase.h"

#include "cantera/base/Interface.h"

#include "fvCFD.H"

using namespace Cantera;

double simple_demo1
(
    const std::string& CanteraMechanismFile_,
    const std::string& PhaseName1_,
    const std::string& PhaseName2_,
    const std::string& PhaseName3_,
    const std::string& PhaseName4_
)
{
    // First, create mixtures, gas and surf
    //this is a test
    

    // Create a new 'Solution' object that provides access to ThermoPhase, Kinetics and
    // Transport objects.
    auto gas_Al = newSolution(CanteraMechanismFile_, PhaseName1_);
    auto bulk_Al = newSolution(CanteraMechanismFile_, PhaseName2_);
    auto bulk_Al2O3 = newSolution(CanteraMechanismFile_, PhaseName3_);

    //const std::vector<std::string>& adjacent)
    //const std::vector<shared_ptr<Solution>>& adjacent
    std::vector<shared_ptr<Cantera::Solution>> adjacent;
    adjacent.push_back(gas_Al);
    adjacent.push_back(bulk_Al);
    adjacent.push_back(bulk_Al2O3);

    auto CanteraInterface_ = Cantera::newInterface(CanteraMechanismFile_, PhaseName4_, adjacent);

    //auto CanteraReactorSurf_ = new Cantera::ReactorSurface();
    Cantera::ReactorSurface CanteraReactorSurf_;
    //! Get the surface phase Kinetics object, shared_ptr<InterfaceKinetics>
    auto surfkin = CanteraInterface_->kinetics();
    //auto *temp = surfkin.get();//convert from shared_ptr to *


    auto surf = CanteraInterface_->thermo();
    auto gas = gas_Al->thermo();
    //gas->setState_TPX(1001.0, Cantera::OneAtm, "H2:2.0, O2:1.0, N2:4.0");
    
    auto kin = gas_Al->kinetics();
    
    Info<< nl
        << "surfkin->surfacePhaseIndex() = " << surfkin->surfacePhaseIndex() << nl
        << "gas->nSpecies() = " << gas->nSpecies() << nl
        << "gas->speciesName(0) = " << gas->speciesName(0) << nl 
        << "kin->nReactions() = " << kin->nReactions() << nl
        << "kin->reactionString(0) = " << kin->reactionString(0) << nl
        << "surf->nSpecies() = " << surf->nSpecies() << nl
        << "surf->speciesName(0) = " << surf->speciesName(0) << nl
        << "surfkin->nReactions() = " << surfkin->nReactions() << nl
        << "surfkin->reactionString(0) = " << surfkin->reactionString(0) << nl
        << endl;
        


    scalar relTol_ = 1.0e-6;
    
    scalar absTol_ = 1.0e-12;

    const scalar V = 0.001;
    const scalar A = V*2.0/5.0;

    label gasSpecie_ = gas->nSpecies();
    label surfSpecie_ = surf->nSpecies();

    scalarField c0(gasSpecie_, 0.0);
    scalarField c1(gasSpecie_, 0.0);

    scalarField surfC0(surfSpecie_, 0.0);
    scalarField surfC1(surfSpecie_, 0.0);

    


    surfC0[0] = 1.0;//AL
    surf->setCoverages(surfC0.begin());
    //CanteraReactorSurf_->setCoverages(surfC0.begin());

    // create a reactor
    //Cantera::Reactor CanteraReactor_;
    Cantera::IdealGasReactor CanteraReactor_;

    // set the volume for the reactor
    // useless in single mesh (since volume is 1 m^3 in single mesh and the default value is 1 too in Cantera)
    // not sure for normal case
    CanteraReactor_.setInitialVolume(V);
    // keep T const before and after sim.advance.
    // this will give you a little improvement
    CanteraReactor_.setEnergy(0);

    CanteraReactorSurf_.setKinetics(surfkin.get());
    CanteraReactorSurf_.setArea(A);//Surface area m^2
    //CanteraReactorSurf_->setReactor(&CanteraReactor_);//void setReactor(ReactorBase* reactor);//use addSurface instead//void addSurface(ReactorSurface* surf);

    CanteraReactor_.addSurface(&CanteraReactorSurf_);
    // 'insert' the gas into the reactor and environment.  Note
    // that it is ok to insert the same gas object into multiple
    // reactors or reservoirs. All this means is that this object
    // will be used to evaluate thermodynamic or kinetic
    // quantities needed.
    CanteraReactor_.insert(gas_Al);

    // create a reactor network
    Cantera::ReactorNet sim;
    sim.addReactor(CanteraReactor_);
    sim.setTolerances(relTol_, absTol_);

    // write down a file
    OFstream OutputAl("Al_surf_gas.csv");

    OutputAl << "t" << "," << "O2" << "," << "AL" << "," << "ALO" << "," << "ALO2" << "," << "AL2O" << "," << "AL2O2" << "," << "AL2O3" << "," << "AL2O3(L)" << "," << "O(S)" << nl;
    //OutputAl << "t" << "," << gas->speciesName(1) << "," << gas->speciesName(5) << "," << gas->speciesName(6) << "," << gas->speciesName(7) << "," << gas->speciesName(8) << "," << gas->speciesName(9) << "," << gas->speciesName(10) << "," << gas->speciesName(11) << "," << surf->speciesName(3) << nl;
    

    // Calculate reaction heat
    scalarField con0(gasSpecie_, 0.0);
    scalarField con1(gasSpecie_, 0.0);

    scalarField surfCon0(surfSpecie_, 0.0);
    scalarField surfCon1(surfSpecie_, 0.0);

    gas -> getConcentrations(con0.begin());
    surf-> getConcentrations(surfCon0.begin());

    scalar hreact0 = sum(con0) * V * gas -> enthalpy_mole() + sum(surfCon0) * A * surf -> enthalpy_mole();


    for(int timeIndex = 0; timeIndex < 401; timeIndex ++)
    {
        scalar time = timeIndex * 0.00001;
        // while(sim.time()<time)
        // {
        sim.advance(time);
        // }
        gas->getMoleFractions(c1.begin());
        surf->getCoverages(surfC1.begin());
        OutputAl << time*1000 << "," << c1[1] << "," << c1[5] << "," << c1[6] << "," << c1[7] << "," << c1[8] << "," << c1[9] << "," << c1[10] << "," << c1[11] << "," << surfC1[3] << nl;
    }

    gas -> getConcentrations(con1.begin());
    surf-> getConcentrations(surfCon1.begin());
    scalar hreact1 = sum(con1) * V * gas -> enthalpy_mole() + sum(surfCon1) * A * surf -> enthalpy_mole();
    Info << "Total Hreact is " << hreact0 - hreact1 << " J." << endl;
    return (hreact0 - hreact1) ;
}



// the main program just calls function simple_demo2 within a 'try' block, and
// catches exceptions that might be thrown
int main(int argc, char* argv[])
{
    const std::string CanteraMechanismFile_ = "Al_surf_gas.yaml";//argv[1];
    const std::string PhaseName1_ = "gasAl"; //argv[2];
    const std::string PhaseName2_ = "bulkAl";
    const std::string PhaseName3_ = "bulkAl2O3";
    const std::string PhaseName4_ = "surfaceAl";

    Info << "Demo programe to call Cantera surface reaction functions" << endl;
    double Hreact = 0.0 ;
    
    try
    {
        Hreact = simple_demo1(CanteraMechanismFile_, PhaseName1_, PhaseName2_, PhaseName3_, PhaseName4_);
    }
    catch (std::exception& err)
    {
        std::cout << err.what() << std::endl;
    }
    
    Info << "End" << endl;
}

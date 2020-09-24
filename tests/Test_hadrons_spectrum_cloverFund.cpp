/*
 * Test_hadrons_spectrum.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;
    
    // run setup ///////////////////////////////////////////////////////////////
    Application              application;
    std::vector<std::string> flavour = {"l"};
    std::vector<double>      mass    = {-0.55};
    double                   csw     = 1.0;
    //double                   csw     = 0.0;
    
    std::cout << "v mass" << mass << std::endl;
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 40;
    globalPar.trajCounter.end   = 149;
    globalPar.trajCounter.step  = 1;

    globalPar.runId             = "b11.0_am-0.45_am-0.45";
    application.setPar(globalPar);
    
    //gauge field
    MIO::LoadNersc::Par gaugePar;
    gaugePar.file = "ckpoint_b11.0_am-0.45_am-0.45_lat";
    application.createModule<MIO::LoadNersc>("gauge", gaugePar);
 
    //source
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);
    // sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    //std::string boundary = "1 1 1 -1";
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";
    
    

    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // actions
        MAction::WilsonClover::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.mass  = mass[i];
        actionPar.boundary = boundary;
        actionPar.csw_r = csw;
        actionPar.csw_t = csw;
        
        actionPar.clover_anisotropy.isAnisotropic= false;
        actionPar.clover_anisotropy.t_direction  = 3    ;   // Explicit for D=4
        actionPar.clover_anisotropy.xi_0         = 1.0  ;
        actionPar.clover_anisotropy.nu           = 1.0  ;

        
          application.createModule<MAction::WilsonClover>("WilsonClover_" + flavour[i], actionPar);
        
        // solvers
        MSolver::RBPrecCG::Par solverPar;
        solverPar.action       = "WilsonClover_" + flavour[i];
        solverPar.residual     = 1.0e-8;
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_" + flavour[i],
                                                    solverPar);
        
        // propagators
        MFermion::GaugeProp::Par quarkPar;
        quarkPar.solver = "CG_" + flavour[i];
        quarkPar.source = "pt";
        application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i], quarkPar);

    }
    for (unsigned int i = 0; i < flavour.size(); ++i)
    for (unsigned int j = i; j < flavour.size(); ++j)
    {
        //std::cout << "FLAVOR SIZE" << flavour.size() << std::endl;
        MContraction::Meson::Par mesPar;
        
        mesPar.output  = "mesons/pt_" + flavour[i] + flavour[j];
       
        mesPar.q1      = "Qpt_" + flavour[i];
        mesPar.q2      = "Qpt_" + flavour[j];
        
       // mesPar.gammas  = "all";
      //  mesPar.gammas = "(Gamma5 Gamma5)(GammaX GammaX)(GammaY GammaY)(GammaZ GammaZ)";
        mesPar.gammas = "(Gamma5 Gamma5)(GammaX GammaX)";
        mesPar.sink    = "sink";
        application.createModule<MContraction::Meson>("meson_pt_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
     
     
    }

    
    // execution
    application.saveParameterFile("spectrum.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}

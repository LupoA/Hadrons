/*******************************************************************************
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: tests/hadrons/Test_hadrons_spectrum.cc
 
 Copyright (C) 2015
 
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 
 See the full license in the file "LICENSE" in the top level distribution
 directory.
 *******************************************************************************/

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

namespace HadronsInputs
{
  class TrajRange : Serializable
  {
  public:
    // -----------------------------------
    GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
				    unsigned int, start,
				    unsigned int, end,
				    unsigned int, step);
  };
  // -----------------------------------
  class GaugeConfigurations : Serializable
  {
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeConfigurations,
				    std::string, gname,
				    std::string, rname);
  };
  // -----------------------------------
  class ValenceFermions : Serializable
  {
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ValenceFermions,
				    double, valencemass,
				    double, coeffcsw);
  };
  // -----------------------------------
}


struct HadronsPar
{
  HadronsInputs::TrajRange           trajRange;
  HadronsInputs::GaugeConfigurations      gaugeConfigurations;
  HadronsInputs::ValenceFermions valenceFermions;
};

int main(int argc, char *argv[])
{
  // parse command line
  std::string parFilename;

  if (argc < 2)
    {
      std::cerr << "usage: " << argv[0] << " <parameter file>";
      std::cerr << std::endl;
      return EXIT_FAILURE;
    }
  parFilename = argv[1];

  // parse parameter file
  HadronsPar par;
  XmlReader reader(parFilename);

  read(reader,        "trajRange",        par.trajRange);
  read(reader,        "gaugeConfigurations",        par.gaugeConfigurations);
  read(reader,        "valenceFermions",        par.valenceFermions);
    
  unsigned int trajStart = par.trajRange.start;
  unsigned int trajEnd   = par.trajRange.end;
  unsigned int trajStep  = par.trajRange.step;
  std::string nameOfGauge = par.gaugeConfigurations.gname;
  std::string nameOfRun = par.gaugeConfigurations.rname;

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
    std::vector<double>      mass    = {par.valenceFermions.valencemass};
    double                   csw     = par.valenceFermions.coeffcsw;

    //Global parameters
    Application::GlobalPar globalPar;
    //globalPar.trajCounter.start = 10;
    //globalPar.trajCounter.end   = 11;
    //globalPar.trajCounter.step  = 1;
    globalPar.trajCounter.start = trajStart;
    globalPar.trajCounter.end   = trajEnd;
    globalPar.trajCounter.step  = trajStep;
    globalPar.runId             = nameOfRun;
    application.setPar(globalPar);

    //Gauge field
    //application.createModule<MGauge::Unit>("gauge");
    MIO::LoadNersc::Par gaugePar;
    //gaugePar.file = "b10_am-0.58_am-0.58";
    gaugePar.file = nameOfGauge;
    //application.createModule<MGauge::Load>("gauge", gaugePar);
    application.createModule<MIO::LoadNersc>("gauge", gaugePar);

    //Update Representation
    MGauge::FundtoTwoIndexAsym::Par GroupPar;
    GroupPar.gaugeconf = "gauge";
    application.createModule<MGauge::FundtoTwoIndexAsym>("gauge_Hirep",GroupPar);

    //Source
    MSource::Point2AS::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point2AS>("pt", ptPar);
   
    //Sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);

/*
    MSink::Point2AS::Par sinkPar;
    sinkPar.mom = "0 0 0" ;
    application.createModule<MSink::Point2AS>("sink", sinkPar);
*/

    //Set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";

    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // actions
        MAction::WilsonClover2AS::Par actionPar;
        actionPar.gauge = "gauge_Hirep";
        actionPar.mass  = mass[i];
        actionPar.boundary = boundary;
        actionPar.csw_r = csw;
        actionPar.csw_t = csw;
        
    actionPar.clover_anisotropy.isAnisotropic= false;
        actionPar.clover_anisotropy.t_direction  = Nd-1 ;
        actionPar.clover_anisotropy.xi_0         = 1.0  ;
        actionPar.clover_anisotropy.nu           = 1.0  ;
    actionPar.boundary = boundary;
        application.createModule<MAction::WilsonClover2AS>("WilsonClover2AS_" + flavour[i], actionPar);
        
        // solvers
        MSolver::RBPrecCG2AS::Par solverPar;
        solverPar.action   = "WilsonClover2AS_" + flavour[i];
        solverPar.residual = 1.0e-8;
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG2AS>("CG_" + flavour[i],
                                                    solverPar);
           
        // propagators
        MFermion::GaugeProp2AS::Par quarkPar;
        quarkPar.solver = "CG_" + flavour[i];
        quarkPar.source = "pt";
        application.createModule<MFermion::GaugeProp2AS>("Qpt_" + flavour[i], quarkPar);
    }
    for (unsigned int i = 0; i < flavour.size(); ++i)
    for (unsigned int j = i; j < flavour.size(); ++j)
    {
        MContraction::Meson2AS::Par mesPar;
        
        mesPar.output  = "Hirep_mesons/pt_" + flavour[i] + flavour[j];
        mesPar.q1      = "Qpt_" + flavour[i];
        mesPar.q2      = "Qpt_" + flavour[j];
        mesPar.gammas  = "(Gamma5 Gamma5)";
        mesPar.sink    = "sink";
        application.createModule<MContraction::Meson2AS>("meson_pt_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
     
    }
   
    // execution
    application.saveParameterFile("Wilson2AS_spectrum.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}

/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./HMC/Mobius2p1fIDSDRprd_40IDpartner.cc

Copyright (C) 2015-2016

Author: Christopher Kelly <ckelly@bnl.gov>
Author: Peter Boyle <pabobyle@ph.ed.ac.uk>


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
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>

using namespace Grid;

//Production binary for the 40ID periodic BC partner ensemble generated for NPR

struct EOFAparameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(EOFAparameters,
				  OneFlavourRationalParams, rat_params,
				  double, action_tolerance,
				  double, md_tolerance,
				  double, relup_freq);

  EOFAparameters() { 
    action_tolerance = 1e-10;
    md_tolerance = 1e-8;
    relup_freq=0.1;

    rat_params.lo = 1.0;
    rat_params.hi = 25.0;
    rat_params.MaxIter  = 50000;
    rat_params.tolerance= 1.0e-9;
    rat_params.degree   = 14;
    rat_params.precision= 50;
  }
};

struct TwoFlavorActionParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(TwoFlavorActionParameters,
				  double, action_tolerance,
				  double, md_tolerance,
				  double, relup_freq);

  TwoFlavorActionParameters() { 
    action_tolerance = 1e-10;
    md_tolerance = 1e-8;
    relup_freq=0.1;
  }
};

struct EvolParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(EvolParameters,
                                  Integer, StartTrajectory,
                                  Integer, Trajectories,
				  Integer, SaveInterval,
				  Integer, Steps,
				  RealD, TrajectoryLength,
                                  bool, MetropolisTest,
				  std::string, StartingType,			  
				  std::vector<TwoFlavorActionParameters>, ratio_l,
				  std::vector<EOFAparameters>, eofa_s,
				  TwoFlavorActionParameters, ratio_DSDR);

  EvolParameters() {
    //For initial thermalization; afterwards user should switch Metropolis on and use StartingType=CheckpointStart
    MetropolisTest    = false;
    StartTrajectory   = 0;
    Trajectories      = 50;
    SaveInterval = 5;
    StartingType      = "ColdStart";
    Steps = 5;
    TrajectoryLength = 1.0;
  }
};

bool fileExists(const std::string &fn){
  std::ifstream f(fn);
  return f.good();
}




struct LanczosParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParameters,
				  double, alpha,
				  double, beta,
				  double, mu,
				  int, ord,
				  int, n_stop,
				  int, n_want,
				  int, n_use,
				  double, tolerance);

  LanczosParameters() {
    alpha = 35;
    beta = 5;
    mu = 0;
    ord = 100;
    n_stop = 10;
    n_want = 10;
    n_use = 15;
    tolerance = 1e-6;
  }
};

template<typename FermionImplPolicy>
void checkEOFA(ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> &EOFA,
	       GridCartesian* FGrid, GridParallelRNG &rng, const LatticeGaugeFieldD &latt){
  std::cout << GridLogMessage << "Starting EOFA action/bounds check" << std::endl;
  typename FermionImplPolicy::FermionField eta(FGrid);
  RealD scale = std::sqrt(0.5);
  gaussian(rng,eta); eta = eta * scale;

  //Use the inbuilt check
  EOFA.refresh(latt, eta);
  EOFA.S(latt);
  std::cout << GridLogMessage << "Finished EOFA upper action/bounds check" << std::endl;
}


template<typename FermionImplPolicy>
class EOFAlinop: public LinearOperatorBase<typename FermionImplPolicy::FermionField>{
  ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> &EOFA;
  LatticeGaugeFieldD &U;
public:
  EOFAlinop(ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> &EOFA, LatticeGaugeFieldD &U): EOFA(EOFA), U(U){}

  typedef typename FermionImplPolicy::FermionField Field;
  void OpDiag (const Field &in, Field &out){ assert(0); }
  void OpDir  (const Field &in, Field &out,int dir,int disp){ assert(0); }
  void OpDirAll  (const Field &in, std::vector<Field> &out){ assert(0); } 

  void Op     (const Field &in, Field &out){ assert(0); }
  void AdjOp  (const Field &in, Field &out){ assert(0); }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ assert(0); }
  void HermOp(const Field &in, Field &out){ EOFA.Meofa(U, in, out); }
};

template<typename FermionImplPolicy>
void upperBoundEOFA(ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> &EOFA,
		    GridCartesian* FGrid, GridParallelRNG &rng, LatticeGaugeFieldD &latt){
  std::cout << GridLogMessage << "Starting EOFA upper bound compute" << std::endl;
  EOFAlinop<FermionImplPolicy> linop(EOFA, latt);
  typename FermionImplPolicy::FermionField eta(FGrid);
  gaussian(rng,eta);
  PowerMethod<typename FermionImplPolicy::FermionField> power_method(1e-8,1000);
  auto lambda_max = power_method(linop,eta);
  std::cout << GridLogMessage << "Upper bound of EOFA operator " << lambda_max << std::endl;
}

//Applications of M^{-1} cost the same as M for EOFA!
template<typename FermionImplPolicy>
class EOFAinvLinop: public LinearOperatorBase<typename FermionImplPolicy::FermionField>{
  ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> &EOFA;
  LatticeGaugeFieldD &U;
public:
  EOFAinvLinop(ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> &EOFA, LatticeGaugeFieldD &U): EOFA(EOFA), U(U){}

  typedef typename FermionImplPolicy::FermionField Field;
  void OpDiag (const Field &in, Field &out){ assert(0); }
  void OpDir  (const Field &in, Field &out,int dir,int disp){ assert(0); }
  void OpDirAll  (const Field &in, std::vector<Field> &out){ assert(0); } 

  void Op     (const Field &in, Field &out){ assert(0); }
  void AdjOp  (const Field &in, Field &out){ assert(0); }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ assert(0); }
  void HermOp(const Field &in, Field &out){ EOFA.MeofaInv(U, in, out); }
};

template<typename FermionImplPolicy>
void lowerBoundEOFA(ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> &EOFA,
		    GridCartesian* FGrid, GridParallelRNG &rng, LatticeGaugeFieldD &latt){
  std::cout << GridLogMessage << "Starting EOFA lower bound compute using power method on M^{-1}. Inverse of highest eigenvalue is the lowest eigenvalue of M" << std::endl;
  EOFAinvLinop<FermionImplPolicy> linop(EOFA, latt);
  typename FermionImplPolicy::FermionField eta(FGrid);
  gaussian(rng,eta);
  PowerMethod<typename FermionImplPolicy::FermionField> power_method(1e-8,1000);
  auto lambda_max = power_method(linop,eta);
  std::cout << GridLogMessage << "Lower bound of EOFA operator " << 1./lambda_max << std::endl;
}


NAMESPACE_BEGIN(Grid);

  template<class FermionOperatorD, class FermionOperatorF, class SchurOperatorD, class  SchurOperatorF> 
  class MixedPrecisionReliableUpdateConjugateGradientOperatorFunction : public OperatorFunction<typename FermionOperatorD::FermionField> {
  public:
    typedef typename FermionOperatorD::FermionField FieldD;
    typedef typename FermionOperatorF::FermionField FieldF;

    using OperatorFunction<FieldD>::operator();

    RealD Tolerance;
    Integer MaxIterations;

    RealD Delta; //reliable update parameter

    GridBase* SinglePrecGrid4; //Grid for single-precision fields
    GridBase* SinglePrecGrid5; //Grid for single-precision fields

    FermionOperatorF &FermOpF;
    FermionOperatorD &FermOpD;;
    SchurOperatorF &LinOpF;
    SchurOperatorD &LinOpD;
    
    MixedPrecisionReliableUpdateConjugateGradientOperatorFunction(RealD tol, 
								  RealD delta,
								  Integer maxit, 
								  GridBase* _sp_grid4, 
								  GridBase* _sp_grid5, 
								  FermionOperatorF &_FermOpF,
								  FermionOperatorD &_FermOpD,
								  SchurOperatorF   &_LinOpF,
								  SchurOperatorD   &_LinOpD): 
      LinOpF(_LinOpF),
      LinOpD(_LinOpD),
      FermOpF(_FermOpF),
      FermOpD(_FermOpD),
      Tolerance(tol), 
      Delta(delta),
      MaxIterations(maxit), 
      SinglePrecGrid4(_sp_grid4),
      SinglePrecGrid5(_sp_grid5)
    { 
    };

    void operator()(LinearOperatorBase<FieldD> &LinOpU, const FieldD &src, FieldD &psi) {

      std::cout << GridLogMessage << " Mixed precision reliable CG update wrapper operator() "<<std::endl;

      SchurOperatorD * SchurOpU = static_cast<SchurOperatorD *>(&LinOpU);
      assert(&(SchurOpU->_Mat)==&(LinOpD._Mat));

      precisionChange(FermOpF.Umu, FermOpD.Umu);

      pickCheckerboard(Even,FermOpF.UmuEven,FermOpF.Umu);
      pickCheckerboard(Odd ,FermOpF.UmuOdd ,FermOpF.Umu);

      ////////////////////////////////////////////////////////////////////////////////////
      // Make a mixed precision conjugate gradient
      ////////////////////////////////////////////////////////////////////////////////////

      ConjugateGradientReliableUpdate<FieldD,FieldF> MPCG(Tolerance,MaxIterations,Delta,SinglePrecGrid5,LinOpF,LinOpD);
      std::cout << GridLogMessage << "Calling mixed precision reliable update Conjugate Gradient" <<std::endl;
      MPCG(src,psi);
    }
  };



NAMESPACE_END(Grid);





int main(int argc, char **argv) {
  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  std::string param_file = "params.xml";
  bool file_load_check = false;

  std::string serial_seeds = "1 2 3 4 5";
  std::string parallel_seeds = "6 7 8 9 10";

  int i=1;
  while(i < argc){
    std::string sarg(argv[i]);
    if(sarg == "--param_file"){
      assert(i!=argc-1);
      param_file = argv[i+1];
      i+=2;
    }else if(sarg == "--read_check"){ //check the fields load correctly and pass checksum/plaquette repro
      file_load_check = true;
      i++;
    }else if(sarg == "--set_seeds"){ //set the rng seeds. Expects two vector args, e.g.  --set_seeds 1.2.3.4 5.6.7.8
      assert(i < argc-2);
      std::vector<int> tmp;
      GridCmdOptionIntVector(argv[i+1],tmp);
      {
	std::stringstream ss;
	for(int j=0;j<tmp.size()-1;j++) ss << tmp[j] << " ";
	ss << tmp.back();
	serial_seeds = ss.str();
      }
      GridCmdOptionIntVector(argv[i+2],tmp);
      {
	std::stringstream ss;
	for(int j=0;j<tmp.size()-1;j++) ss << tmp[j] << " ";
	ss << tmp.back();
	parallel_seeds = ss.str();
      }
      i+=3;
      std::cout << GridLogMessage << "Set serial seeds to " << serial_seeds << std::endl;
      std::cout << GridLogMessage << "Set parallel seeds to " << parallel_seeds << std::endl;
      
    }else{
      i++;
    }
  }

  
  //Read the user parameters
  EvolParameters user_params;
  
  if(fileExists(param_file)){
    std::cout << GridLogMessage << " Reading " << param_file << std::endl;
    Grid::XmlReader rd(param_file);
    read(rd, "Params", user_params);
  }else if(!GlobalSharedMemory::WorldRank){
    std::cout << GridLogMessage << " File " << param_file << " does not exist" << std::endl;
    std::cout << GridLogMessage << " Writing xml template to " << param_file << ".templ" << std::endl;
    {
      Grid::XmlWriter wr(param_file + ".templ");
      write(wr, "Params", user_params);
    }
    std::cout << GridLogMessage << " Done" << std::endl;
    Grid_finalize();
    return 0;
  }

  typedef MobiusEOFAFermionD EOFAfermD;
  typedef MobiusFermionD fermionD;
  typedef typename fermionD::Impl_t FermionImplPolicyD;
  typedef typename fermionD::FermionField FermionFieldD;

  typedef MobiusEOFAFermionF EOFAfermF;
  typedef MobiusFermionF fermionF;
  typedef typename fermionF::Impl_t FermionImplPolicyF;
  typedef typename fermionF::FermionField FermionFieldF;

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IntegratorParameters MD;

  #define USE_FORCE_GRADIENT
  
#ifndef USE_FORCE_GRADIENT
  typedef ConjugateHMCRunnerD<MinimumNorm2> HMCWrapper; //NB: This is the "Omelyan integrator"
  MD.name    = std::string("MinimumNorm2");
#else
  typedef ConjugateHMCRunnerD<ForceGradient> HMCWrapper;
  MD.name    = std::string("ForceGradient");
#endif
  MD.MDsteps = user_params.Steps;
  MD.trajL   = user_params.TrajectoryLength;

  typedef HMCWrapper::ImplPolicy GaugeImplPolicy;
  
  HMCparameters HMCparams;
  HMCparams.StartTrajectory  = user_params.StartTrajectory;
  HMCparams.Trajectories     = user_params.Trajectories;
  HMCparams.NoMetropolisUntil= 0;
  HMCparams.StartingType     = user_params.StartingType;
  HMCparams.MetropolisTest = user_params.MetropolisTest;
  HMCparams.PerformRandomShift = false;
  HMCparams.MD = MD;
  HMCWrapper TheHMC(HMCparams);

  // Grid from the command line arguments --grid and --mpi
  TheHMC.Resources.AddFourDimGrid("gauge"); // use default simd lanes decomposition

  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix    = "ckpoint_rng";
  CPparams.saveInterval  = user_params.SaveInterval;
  CPparams.format        = "IEEE64BIG";
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  //Note that checkpointing saves the RNG state so that this initialization is required only for the very first configuration
  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = serial_seeds;
  RNGpar.parallel_seeds = parallel_seeds;
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  typedef PlaquetteMod<GaugeImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  //////////////////////////////////////////////
  //aiming for ainv=1.723 GeV
  //                                  me         bob
  //Estimated  a(ml+mres) [40ID] = 0.001305    0.00131
  //           a(mh+mres) [40ID] = 0.035910    0.03529
  //Estimate Ls=12, b+c=2  mres~0.0011

  //1/24/2022 initial mres measurement gives mres=0.001,  adjusted light quark mass to 0.0003 from 0.0001
  
  const int Ls      = 12;
  Real beta         = 1.848;
  Real light_mass   = 0.0003;
  Real strange_mass = 0.0342;
  Real pv_mass      = 1.0;
  RealD M5  = 1.8;
  RealD mobius_scale = 2.; //b+c

  RealD mob_bmc = 1.0;
  RealD mob_b = (mobius_scale + mob_bmc)/2.;
  RealD mob_c = (mobius_scale - mob_bmc)/2.;

  std::cout << GridLogMessage
	    << "Ensemble parameters:" << std::endl
	    << "Ls=" << Ls << std::endl
	    << "beta=" << beta << std::endl
	    << "light_mass=" << light_mass << std::endl
	    << "strange_mass=" << strange_mass << std::endl
	    << "mobius_scale=" << mobius_scale << std::endl;
  
  //Setup the Grids
  auto UGridD   = TheHMC.Resources.GetCartesian();
  auto UrbGridD = TheHMC.Resources.GetRBCartesian();
  auto FGridD     = SpaceTimeGrid::makeFiveDimGrid(Ls,UGridD);
  auto FrbGridD   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGridD);

  GridCartesian* UGridF = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexF::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGridF = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridF);
  auto FGridF     = SpaceTimeGrid::makeFiveDimGrid(Ls,UGridF);
  auto FrbGridF   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGridF);

  IwasakiGaugeActionD GaugeAction(beta);

  // temporarily need a gauge field
  LatticeGaugeFieldD Ud(UGridD);
  LatticeGaugeFieldF Uf(UGridF);
 
  //Setup the BCs
  std::vector<Complex> boundary = {1,1,1,-1};
  fermionD::ImplParams Params(boundary);
  
  //Run optional gauge field checksum checker and exit
  if(file_load_check){
    TheHMC.initializeGaugeFieldAndRNGs(Ud);
    std::cout << GridLogMessage << " Done" << std::endl;
    Grid_finalize();
    return 0;
  }


  ////////////////////////////////////
  // Collect actions
  ////////////////////////////////////
  ActionLevel<HMCWrapper::Field> Level1(1); //light quark + strange quark
  ActionLevel<HMCWrapper::Field> Level2(3); //DSDR
  ActionLevel<HMCWrapper::Field> Level3(1); //gauge  2->1 on v8

  /////////////////////////////////////////////////////////////
  // Light two-flavor DWF action
  // have to be careful with the parameters, cf. Test_dwf_gpforce_eofa.cc
  /////////////////////////////////////////////////////////////
  typedef TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicyD> TwoFlavorActionD;
  typedef TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicyF> TwoFlavorActionF;

  typedef SchurDiagMooeeOperator<fermionD,FermionFieldD> SchurOpD;
  typedef SchurDiagMooeeOperator<fermionF,FermionFieldF> SchurOpF;
  typedef MixedPrecisionReliableUpdateConjugateGradientOperatorFunction<fermionD, fermionF, SchurOpD, SchurOpF> relupCG;

  //We will use ml=ms to speed up the ensemble generation. The sea quark masses don't matter for NPR
  std::vector<RealD> lquark_light_masses = { strange_mass,  0.064,  0.112,    0.16,    0.256,  0.442,     0.628    };
  std::vector<RealD> lquark_pv_masses =    { 0.064,         0.112,  0.16,     0.256,   0.442,  0.628,     1.0      };
  int n_light_hsb = 7;

  assert(user_params.ratio_l.size() == n_light_hsb);
  
  for(int i=0;i<n_light_hsb;i++){
    RealD iml = lquark_light_masses[i];
    RealD ipv = lquark_pv_masses[i];

    fermionD* ferm_l_d = new fermionD(Ud, *FGridD, *FrbGridD, *UGridD, *UrbGridD, iml, M5, mob_b, mob_c, Params);
    fermionD* ferm_h_d = new fermionD(Ud, *FGridD, *FrbGridD, *UGridD, *UrbGridD, ipv, M5, mob_b, mob_c, Params);
    fermionF* ferm_l_f = new fermionF(Uf, *FGridF, *FrbGridF, *UGridF, *UrbGridF, iml, M5, mob_b, mob_c, Params);
    fermionF* ferm_h_f = new fermionF(Uf, *FGridF, *FrbGridF, *UGridF, *UrbGridF, ipv, M5, mob_b, mob_c, Params);
  
    SchurOpD* schurop_l_d = new SchurOpD(*ferm_l_d);
    SchurOpD* schurop_h_d = new SchurOpD(*ferm_h_d);

    SchurOpF* schurop_l_f = new SchurOpF(*ferm_l_f);
    SchurOpF* schurop_h_f = new SchurOpF(*ferm_h_f);

    relupCG*  heatbath_solver = new relupCG(user_params.ratio_l[i].action_tolerance, user_params.ratio_l[i].relup_freq, 50000,UrbGridF,FrbGridF,*ferm_h_f,*ferm_h_d,*schurop_h_f,*schurop_h_d);  //inverts PV matrix
    relupCG*  action_solver = new relupCG(user_params.ratio_l[i].action_tolerance, user_params.ratio_l[i].relup_freq, 50000,UrbGridF,FrbGridF,*ferm_l_f,*ferm_l_d,*schurop_l_f,*schurop_l_d);  //inverts light matrix
    relupCG*  deriv_solver = new relupCG(user_params.ratio_l[i].md_tolerance, user_params.ratio_l[i].relup_freq, 50000,UrbGridF,FrbGridF,*ferm_l_f,*ferm_l_d,*schurop_l_f,*schurop_l_d);  //inverts light matrix

    TwoFlavorActionD* pf_action = new TwoFlavorActionD(*ferm_h_d, *ferm_l_d, *deriv_solver, *action_solver, *heatbath_solver);

    Level1.push_back(pf_action);    
  }

  
  /////////////////////////////////////////////////////////////
  // Strange EOFA action
  // have to be careful with the parameters, cf. Test_dwf_gpforce_eofa.cc
  /////////////////////////////////////////////////////////////
  typedef SchurDiagMooeeOperator<EOFAfermD,FermionFieldD> EOFAschuropD;
  typedef SchurDiagMooeeOperator<EOFAfermF,FermionFieldF> EOFAschuropF;
  typedef ExactOneFlavourRatioMixedPrecHeatbathPseudoFermionAction<FermionImplPolicyD, FermionImplPolicyF> EOFAmixPrecPFaction;
  typedef MixedPrecisionReliableUpdateConjugateGradientOperatorFunction<EOFAfermD, EOFAfermF, EOFAschuropD, EOFAschuropF> EOFA_relupCG;

  std::vector<RealD> squark_light_masses = { strange_mass, 0.276,        0.517 };
  std::vector<RealD> squark_pv_masses =    { 0.276,        0.517,        1.0 };
  int n_strange_hsb = 3;
  
  assert(user_params.eofa_s.size() == n_strange_hsb);
  EOFAmixPrecPFaction* EOFA_pfactions[n_strange_hsb];
  for(int i=0;i<n_strange_hsb;i++){
    RealD iml = squark_light_masses[i];
    RealD ipv = squark_pv_masses[i];

    EOFAfermD* LopD = new EOFAfermD(Ud, *FGridD, *FrbGridD, *UGridD, *UrbGridD, iml, iml, ipv, 0.0, -1, M5, mob_b, mob_c, Params);
    EOFAfermF* LopF = new EOFAfermF(Uf, *FGridF, *FrbGridF, *UGridF, *UrbGridF, iml, iml, ipv, 0.0, -1, M5, mob_b, mob_c, Params);
    EOFAfermD* RopD = new EOFAfermD(Ud, *FGridD, *FrbGridD, *UGridD, *UrbGridD, ipv, iml, ipv, -1.0, 1, M5, mob_b, mob_c, Params);
    EOFAfermF* RopF = new EOFAfermF(Uf, *FGridF, *FrbGridF, *UGridF, *UrbGridF, ipv, iml, ipv, -1.0, 1, M5, mob_b, mob_c, Params);

    EOFAschuropD* linopL_D = new EOFAschuropD(*LopD);
    EOFAschuropD* linopR_D = new EOFAschuropD(*RopD);
    
    EOFAschuropF* linopL_F = new EOFAschuropF(*LopF);
    EOFAschuropF* linopR_F = new EOFAschuropF(*RopF);

    EOFA_relupCG* ActionMCG_L = new EOFA_relupCG(user_params.eofa_s[i].action_tolerance, user_params.eofa_s[i].relup_freq, 50000, UGridF, FrbGridF, *LopF, *LopD, *linopL_F, *linopL_D);
    EOFA_relupCG* ActionMCG_R = new EOFA_relupCG(user_params.eofa_s[i].action_tolerance, user_params.eofa_s[i].relup_freq, 50000, UGridF, FrbGridF, *RopF, *RopD, *linopR_F, *linopR_D);

    EOFA_relupCG* DerivMCG_L = new EOFA_relupCG(user_params.eofa_s[i].md_tolerance, user_params.eofa_s[i].relup_freq, 50000, UGridF, FrbGridF, *LopF, *LopD, *linopL_F, *linopL_D);
    EOFA_relupCG* DerivMCG_R = new EOFA_relupCG(user_params.eofa_s[i].md_tolerance, user_params.eofa_s[i].relup_freq, 50000, UGridF, FrbGridF, *RopF, *RopD, *linopR_F, *linopR_D);

    EOFAmixPrecPFaction* EOFA = new EOFAmixPrecPFaction(*LopF, *RopF,
							*LopD, *RopD, 
							*ActionMCG_L, *ActionMCG_R, 
							*ActionMCG_L, *ActionMCG_R, 
							*DerivMCG_L, *DerivMCG_R, 
							user_params.eofa_s[i].rat_params, true);
    EOFA_pfactions[i] = EOFA;
    Level1.push_back(EOFA);
  }

  ///////////////////////////////////
  // DSDR action
  ///////////////////////////////////
  RealD dsdr_mass=-1.8;   
  //Use same DSDR twists as https://arxiv.org/pdf/1208.4412.pdf
  RealD dsdr_epsilon_f = 0.02; //numerator (in determinant)
  RealD dsdr_epsilon_b = 0.5; 
  WilsonTMFermionD Numerator_DSDR_D(Ud, *UGridD, *UrbGridD, dsdr_mass, dsdr_epsilon_f, Params);
  WilsonTMFermionF Numerator_DSDR_F(Uf, *UGridF, *UrbGridF, dsdr_mass, dsdr_epsilon_f, Params);

  WilsonTMFermionD Denominator_DSDR_D(Ud, *UGridD, *UrbGridD, dsdr_mass, dsdr_epsilon_b, Params);
  WilsonTMFermionF Denominator_DSDR_F(Uf, *UGridF, *UrbGridF, dsdr_mass, dsdr_epsilon_b, Params);
 
  typedef SchurDiagMooeeOperator<WilsonTMFermionD,FermionFieldD> SchurOpD_DSDR;
  typedef SchurDiagMooeeOperator<WilsonTMFermionF,FermionFieldF> SchurOpF_DSDR;
  typedef MixedPrecisionReliableUpdateConjugateGradientOperatorFunction<WilsonTMFermionD, WilsonTMFermionF, SchurOpD_DSDR, SchurOpF_DSDR> relupCG_DSDR;

  SchurOpD_DSDR schurop_DSDR_num_d(Numerator_DSDR_D);
  SchurOpD_DSDR schurop_DSDR_den_d(Denominator_DSDR_D);

  SchurOpF_DSDR schurop_DSDR_num_f(Numerator_DSDR_F);
  SchurOpF_DSDR schurop_DSDR_den_f(Denominator_DSDR_F);

  relupCG_DSDR  heatbath_solver_DSDR(user_params.ratio_DSDR.action_tolerance, user_params.ratio_DSDR.relup_freq, 50000,UrbGridF,UrbGridF,
				     Denominator_DSDR_F,Denominator_DSDR_D,schurop_DSDR_den_f,schurop_DSDR_den_d);  //inverts PV matrix (denominator)
  relupCG_DSDR  action_solver_DSDR(user_params.ratio_DSDR.action_tolerance, user_params.ratio_DSDR.relup_freq, 50000,UrbGridF,UrbGridF,
				   Numerator_DSDR_F,Numerator_DSDR_D,schurop_DSDR_num_f,schurop_DSDR_num_d);  //inverts light matrix (numerator)
  relupCG_DSDR  deriv_solver_DSDR(user_params.ratio_DSDR.md_tolerance, user_params.ratio_DSDR.relup_freq, 50000,UrbGridF,UrbGridF,
				  Numerator_DSDR_F,Numerator_DSDR_D,schurop_DSDR_num_f,schurop_DSDR_num_d);  //inverts light matrix
  
  TwoFlavorActionD pf_action_DSDR(Denominator_DSDR_D, Numerator_DSDR_D, deriv_solver_DSDR, action_solver_DSDR, heatbath_solver_DSDR);

  Level2.push_back(&pf_action_DSDR);
  
  /////////////////////////////////////////////////////////////
  // Gauge action
  /////////////////////////////////////////////////////////////
  Level3.push_back(&GaugeAction);

  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  TheHMC.TheAction.push_back(Level3);
  std::cout << GridLogMessage << " Action complete "<< std::endl;


  //Action tuning
  bool 
    check_eofa=false, 
    upper_bound_eofa=false, lower_bound_eofa(false);

  int eofa_which_hsb;

  for(int i=1;i<argc;i++){
    std::string sarg(argv[i]);

    if(sarg == "--check_eofa"){
      assert(i < argc-1);
      check_eofa = true;
      eofa_which_hsb = std::stoi(argv[i+1]); //-1 indicates all hasenbusch
      assert(eofa_which_hsb == -1 || (eofa_which_hsb >= 0 && eofa_which_hsb < n_light_hsb) );
    }
    else if(sarg == "--upper_bound_eofa"){
      assert(i < argc-1);
      upper_bound_eofa = true;
      eofa_which_hsb = std::stoi(argv[i+1]);
      assert(eofa_which_hsb >= 0 && eofa_which_hsb < n_strange_hsb);
    }
    else if(sarg == "--lower_bound_eofa"){
      assert(i < argc-1);
      lower_bound_eofa = true;      
      eofa_which_hsb = std::stoi(argv[i+1]);
      assert(eofa_which_hsb >= 0 && eofa_which_hsb < n_strange_hsb);
    }
  }
  if(check_eofa || upper_bound_eofa || lower_bound_eofa) {
    std::cout << GridLogMessage << "Running checks" << std::endl;
    TheHMC.initializeGaugeFieldAndRNGs(Ud);

    if(check_eofa){
      if(eofa_which_hsb >= 0){
	std::cout << GridLogMessage << "Starting checking EOFA Hasenbusch " << eofa_which_hsb << std::endl;
	checkEOFA(*EOFA_pfactions[eofa_which_hsb], FGridD, TheHMC.Resources.GetParallelRNG(), Ud);
	std::cout << GridLogMessage << "Finished checking EOFA Hasenbusch " << eofa_which_hsb << std::endl;
      }else{
	for(int i=0;i<n_strange_hsb;i++){
	  std::cout << GridLogMessage << "Starting checking EOFA Hasenbusch " << i << std::endl;
	  checkEOFA(*EOFA_pfactions[i], FGridD, TheHMC.Resources.GetParallelRNG(), Ud);
	  std::cout << GridLogMessage << "Finished checking EOFA Hasenbusch " << i << std::endl;
	}
      }
    }	  
    if(upper_bound_eofa) upperBoundEOFA(*EOFA_pfactions[eofa_which_hsb], FGridD, TheHMC.Resources.GetParallelRNG(), Ud);
    if(lower_bound_eofa) lowerBoundEOFA(*EOFA_pfactions[eofa_which_hsb], FGridD, TheHMC.Resources.GetParallelRNG(), Ud);

    std::cout << GridLogMessage << " Done" << std::endl;
    Grid_finalize();
    return 0;
  }


  //Run the HMC
  std::cout << GridLogMessage << " Running the HMC "<< std::endl;
  TheHMC.Run();

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
} // main

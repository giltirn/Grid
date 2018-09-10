    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_cg_prec.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include<bitset>
#include <Grid/Grid.h>

#include <Grid/qcd/action/fermion/CayleyFermion5DmethodImpl.h>
#include <Grid/qcd/action/fermion/CayleyFermion5DcacheMethodImpl.h>
#include <Grid/qcd/action/fermion/WilsonFermion5DmethodImpl.h>
#include <Grid/qcd/action/fermion/WilsonKernelsMethodImpl.h>
#include <Grid/qcd/action/fermion/WilsonKernelsHandMethodImpl.h>
#include <Grid/qcd/action/fermion/WilsonKernelsAsmMethodImpl.h>


#include "Test_compression_operator_fixedpointcomms.h"
#include "Test_compression_operator_bfp16.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template<typename T>
T parse(const std::string &name, std::istream &in){
  std::string p;
  in >> p;
  assert(p==name);
  char eq;
  in >> eq;
  assert(eq == '=');
  T out;
  in >> out;
  return out;
}


//Test preconditioned CG
void testPCG(GridCartesian* UGrid, GridRedBlackCartesian* UrbGrid, GridCartesian* FGrid, GridRedBlackCartesian* FrbGrid,
	     RealD mass, RealD M5, int Ls, LatticeGaugeFieldD &Umu, LatticeFermionD &src_o){

  LatticeFermionD result_o(FrbGrid);
  LatticeFermionD result_o_2(FrbGrid);
  result_o.checkerboard = Odd;
  result_o = zero;
  result_o_2.checkerboard = Odd;
  result_o_2 = zero;

  DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  SchurDiagMooeeOperator<DomainWallFermionD,LatticeFermionD> HermOpEO(Ddwf);
  //DoNothingLinearOperator<LatticeFermionD> Prec;
  //FixedIterConjugateGradientPreconditioner<LatticeFermionD> Prec(HermOpEO, 20);
  SloppyConjugateGradientPreconditioner<LatticeFermionD> Prec(HermOpEO, 1e-2, 1000);

  std::cout << "Preconditioned CG" << std::endl;
  InexactPreconditionedConjugateGradient<LatticeFermionD> pCG(Prec,1.0e-8,10000);
  pCG(HermOpEO,src_o,result_o);

  std::cout << "Starting regular CG" << std::endl;
  ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
  CG(HermOpEO,src_o,result_o_2);

  LatticeFermionD diff_o(FrbGrid);
  RealD diff = axpy_norm(diff_o, -1.0, result_o, result_o_2);

  std::cout << "pCG HermOp applications " << pCG.IterationsToComplete << "(outer) + " << Prec.InnerIterations << "(inner) = " << pCG.IterationsToComplete + Prec.InnerIterations << std::endl;
  std::cout << "CG HermOp applications " << CG.IterationsToComplete << std::endl;
  std::cout << "Diff between results: " << diff << std::endl;
}



void mixedPrecSolve(Integer &inner, Integer &outer, Integer &patchup,
		    LatticeFermionD &result_o, LatticeFermionD &src_o,
		    RealD outer_tol, RealD inner_tol, RealD outer_loop_norm_mult,
		    LinearOperatorBase<LatticeFermionF> &_Linop_f, LinearOperatorBase<LatticeFermionD> &_Linop_d,
		    GridRedBlackCartesian * FrbGrid_f){
  MixedPrecisionConjugateGradient<LatticeFermionD,LatticeFermionF> mCG(outer_tol, 10000, 50, FrbGrid_f, _Linop_f, _Linop_d);
  mCG.InnerTolerance = inner_tol;
  mCG.OuterLoopNormMult = outer_loop_norm_mult;
  mCG(src_o,result_o);
  inner = mCG.TotalInnerIterations; outer = mCG.TotalOuterIterations; patchup = mCG.TotalFinalStepIterations;
}



void relupSolve(Integer &inner, Integer &outer, Integer &patchup,
		LatticeFermionD &result_o, LatticeFermionD &src_o,
		RealD outer_tol, RealD relup_delta,
		LinearOperatorBase<LatticeFermionF> &_Linop_f, LinearOperatorBase<LatticeFermionD> &_Linop_d,
		GridRedBlackCartesian * FrbGrid_f){     
  ConjugateGradientReliableUpdate<LatticeFermionD,LatticeFermionF> relup(outer_tol, 2000, relup_delta, FrbGrid_f, _Linop_f, _Linop_d);
  relup.ErrorOnNoConverge = false;
  relup(src_o,result_o);
  inner = relup.IterationsToComplete; outer = relup.ReliableUpdatesPerformed; patchup = relup.IterationsToCleanup;
}

void precCGsolve(Integer &inner, Integer &outer, 
		 LatticeFermionD &result_o, LatticeFermionD &src_o,
		 RealD outer_tol, Real inner_tol,
		 LinearOperatorBase<LatticeFermionF> &_Linop_f, LinearOperatorBase<LatticeFermionD> &_Linop_d,
		 GridRedBlackCartesian * FrbGrid_f){           
  SloppyConjugateGradientLowerPrecPreconditioner<LatticeFermionD,LatticeFermionF> prec(_Linop_f, FrbGrid_f, inner_tol, 1000);
  InexactPreconditionedConjugateGradient<LatticeFermionD> CG(prec, outer_tol, 100);
  CG(_Linop_d,src_o,result_o);
  inner = prec.InnerIterations; outer = CG.IterationsToComplete;
}


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int Ls=8;
  RealD mass=0.1;
  RealD outer_tol = 1e-8;
  RealD inner_tol_full = 1e-5;
  RealD inner_tol_half = 1e-5;
  RealD inner_tol_16c = 1e-5;
  RealD inner_tol_8c = 1e-5;
  RealD inner_tol_bfp16 = 1e-5;
  RealD outer_loop_norm_mult = 100.;

  RealD relup_delta_full = 0.1;
  RealD relup_delta_half = 0.1;
  RealD relup_delta_16c = 0.1;
  RealD relup_delta_8c = 0.1;
  RealD relup_delta_bfp16 = 0.1;

  std::string config_file = "";

  enum Algorithm { MixedCG, ReliableUpdate, PrecCG };

  Algorithm algorithm = MixedCG;

  for(int i=1;i<argc;i++){
    if(std::string(argv[i]) == "--params"){
      std::ifstream f(argv[i+1]);
      f.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
      Ls = parse<int>("Ls", f);
#define PARSEIT(NM) NM = parse<RealD>(#NM, f)
      PARSEIT(mass);
      PARSEIT(outer_tol);
      PARSEIT(inner_tol_full);
      PARSEIT(inner_tol_half);
      PARSEIT(inner_tol_16c);
      PARSEIT(inner_tol_8c);
      PARSEIT(inner_tol_bfp16);
      PARSEIT(outer_loop_norm_mult);
      PARSEIT(relup_delta_full);
      PARSEIT(relup_delta_half);
      PARSEIT(relup_delta_16c);
      PARSEIT(relup_delta_8c);
      PARSEIT(relup_delta_bfp16);
#undef PARSEIT

      //f >> outer_tol >> inner_tol_full >> inner_tol_half >> inner_tol_16c >> inner_tol_8c;      
    }else if(std::string(argv[i]) == "--config"){
      config_file = argv[i+1];
    }else if(std::string(argv[i]) == "--algorithm"){
      std::string alg = argv[i+1];
      if(alg == "MixedCG") algorithm = MixedCG;
      else if(alg == "ReliableUpdate") algorithm = ReliableUpdate;
      else if(alg == "PrecCG") algorithm = PrecCG;
      else{
	std::cout << "Error: Unknown algorithm " << alg << std::endl;
	exit(-1);
      }      
    }
  }
  
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexD::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  GridCartesian         * UGrid_f   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_f);
  GridCartesian         * FGrid_f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_f);
  GridRedBlackCartesian * FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_f);
  
  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeFermionD    src(FGrid); random(RNG5,src);
  LatticeFermionD result(FGrid); result=zero;
  LatticeGaugeFieldD Umu(UGrid);
  LatticeGaugeFieldF Umu_f(UGrid_f); 

  if(config_file.size() > 0){
    FieldMetaData header;
    NerscIO::readConfiguration(Umu,header,config_file);
  }else{
    SU3::HotConfiguration(RNG4,Umu);
  }

  precisionChange(Umu_f,Umu);
  
  RealD M5=1.8;

  LatticeFermionD    src_o(FrbGrid);
  pickCheckerboard(Odd,src_o,src);

  //testPCG(UGrid, UrbGrid, FGrid, FrbGrid, mass, M5, Ls, Umu, src_o);
  
  //Compare mixed prec restarted single/single internal with same but with single/compressed
  LatticeFermionD result_o_full(FrbGrid);
  LatticeFermionD result_o_half(FrbGrid);
  LatticeFermionD result_o_16(FrbGrid);
  LatticeFermionD result_o_8(FrbGrid);
  LatticeFermionD result_o_bfp16(FrbGrid);
  result_o_full.checkerboard = Odd;
  result_o_full = zero;
  result_o_bfp16 = result_o_16 = result_o_8 = result_o_half = result_o_full;

  //Std
  DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  SchurDiagMooeeOperator<DomainWallFermionD,LatticeFermionD> HermOpEO(Ddwf);

  DomainWallFermionF Ddwf_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);
  SchurDiagMooeeOperator<DomainWallFermionF,LatticeFermionF> HermOpEO_f(Ddwf_f);

  //1/2 prec
  DomainWallFermionFH Ddwfhalf_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);
  SchurDiagMooeeOperator<DomainWallFermionFH,LatticeFermionF> HermOpEOhalf_f(Ddwfhalf_f);

  //16
  DomainWallFermionFixedPointComms16F DdwfC16_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);
  SchurDiagMooeeOperator<DomainWallFermionFixedPointComms16F,LatticeFermionF> HermOpEOC16_f(DdwfC16_f);

  //8
  DomainWallFermionFixedPointComms8F DdwfC8_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);
  SchurDiagMooeeOperator<DomainWallFermionFixedPointComms8F,LatticeFermionF> HermOpEOC8_f(DdwfC8_f);

  //bfp16
  DomainWallFermionBfp16commsF Ddwfbfp16_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);
  SchurDiagMooeeOperator<DomainWallFermionBfp16commsF,LatticeFermionF> HermOpEObfp16_f(Ddwfbfp16_f);



  Integer inner_16, outer_16, patchup_16;
  Integer inner_8, outer_8, patchup_8;
  Integer inner_half, outer_half, patchup_half;
  Integer inner_full, outer_full, patchup_full;
  Integer inner_bfp16, outer_bfp16, patchup_bfp16;

  if(algorithm == MixedCG){
    std::cout << "Starting mixed CG with single/compressed-16 inner\n";
    mixedPrecSolve(inner_16, outer_16, patchup_16, result_o_16, src_o,
		   outer_tol, inner_tol_16c, outer_loop_norm_mult, HermOpEOC16_f, HermOpEO, FrbGrid_f);
      
    std::cout << "Starting mixed CG with single/compressed-8 inner\n";      
    mixedPrecSolve(inner_8, outer_8, patchup_8, result_o_8, src_o,
		   outer_tol, inner_tol_8c, outer_loop_norm_mult, HermOpEOC8_f, HermOpEO, FrbGrid_f);
      
    std::cout << "Starting mixed CG with single/half inner\n";
    mixedPrecSolve(inner_half, outer_half, patchup_half, result_o_half, src_o,
		   outer_tol, inner_tol_half, outer_loop_norm_mult, HermOpEOhalf_f, HermOpEO, FrbGrid_f);
         
    std::cout << "Starting mixed CG with single/single inner\n";    
    mixedPrecSolve(inner_full, outer_full, patchup_full, result_o_full, src_o,
		   outer_tol, inner_tol_full, outer_loop_norm_mult, HermOpEO_f, HermOpEO, FrbGrid_f);

    std::cout << "Starting mixed CG with single/bfp16 inner\n";    
    mixedPrecSolve(inner_bfp16, outer_bfp16, patchup_bfp16, result_o_bfp16, src_o,
		   outer_tol, inner_tol_bfp16, outer_loop_norm_mult, HermOpEObfp16_f, HermOpEO, FrbGrid_f);


  }else if(algorithm == ReliableUpdate){
    std::cout << "Starting relup CG with single/compressed-16 inner\n";    
    relupSolve(inner_16, outer_16, patchup_16, result_o_16, src_o,
	       outer_tol, relup_delta_16c, HermOpEOC16_f, HermOpEO, FrbGrid_f);
      
    std::cout << "Starting relup CG with single/compressed-8 inner\n";
    relupSolve(inner_8, outer_8, patchup_8, result_o_8, src_o,
	       outer_tol, relup_delta_8c, HermOpEOC8_f, HermOpEO, FrbGrid_f);
    
    std::cout << "Starting relup CG with single/half inner\n";      
    relupSolve(inner_half, outer_half, patchup_half, result_o_half, src_o,
	       outer_tol, relup_delta_half, HermOpEOhalf_f, HermOpEO, FrbGrid_f);
     
    std::cout << "Starting relup CG with single/single inner\n";
    relupSolve(inner_full, outer_full, patchup_full, result_o_full, src_o,
	       outer_tol, relup_delta_full, HermOpEO_f, HermOpEO, FrbGrid_f);

    std::cout << "Starting relup CG with single/bfp16 inner\n";
    relupSolve(inner_bfp16, outer_bfp16, patchup_bfp16, result_o_bfp16, src_o,
	       outer_tol, relup_delta_bfp16, HermOpEObfp16_f, HermOpEO, FrbGrid_f);


  }else if(algorithm == PrecCG){
    std::cout << "Starting sloppy pCG with single/compressed-16 inner\n";    
    precCGsolve(inner_16, outer_16, result_o_16, src_o,
		outer_tol, inner_tol_16c, HermOpEOC16_f, HermOpEO, FrbGrid_f);
    
    std::cout << "Starting sloppy pCG with single/compressed-8 inner\n";    
    precCGsolve(inner_8, outer_8, result_o_8, src_o,
		outer_tol, inner_tol_8c, HermOpEOC8_f, HermOpEO, FrbGrid_f);

    std::cout << "Starting sloppy pCG with single/half inner\n";
    precCGsolve(inner_half, outer_half, result_o_half, src_o,
		outer_tol, inner_tol_half, HermOpEOhalf_f, HermOpEO, FrbGrid_f);
     
    std::cout << "Starting sloppy pCG with single/single inner\n";    
    precCGsolve(inner_full, outer_full, result_o_full, src_o,
		outer_tol, inner_tol_full, HermOpEO_f, HermOpEO, FrbGrid_f);

    std::cout << "Starting sloppy pCG with single/bfp16 inner\n";    
    precCGsolve(inner_bfp16, outer_bfp16, result_o_bfp16, src_o,
		outer_tol, inner_tol_bfp16, HermOpEObfp16_f, HermOpEO, FrbGrid_f);

  }else{
    assert(0);
  }


  std::cout << "Ls " << Ls << std::endl;
  std::cout << "Mass " << mass << std::endl;
  std::cout << "Outer tolerance " << outer_tol << std::endl;

  if(algorithm == MixedCG || algorithm == PrecCG){
    std::cout << "Inner tol full " << inner_tol_full << std::endl;
    std::cout << "Inner tol 1/2 prec " << inner_tol_half << std::endl;
    std::cout << "Inner tol compressed-16 " << inner_tol_16c << std::endl;
    std::cout << "Inner tol compressed-8 " << inner_tol_8c << std::endl;
    std::cout << "Inner tol bfp16 " << inner_tol_bfp16 << std::endl;
      
    if(algorithm == MixedCG)
      std::cout << "Outer loop norm mult " << outer_loop_norm_mult << std::endl;

  }else if(algorithm == ReliableUpdate){
    std::cout << "Relup delta full " << relup_delta_full << std::endl;
    std::cout << "Relup delta 1/2 prec " << relup_delta_half << std::endl;
    std::cout << "Relup delta compressed-16 " << relup_delta_16c << std::endl;
    std::cout << "Relup delta compressed-8 " << relup_delta_8c << std::endl;
    std::cout << "Relup delta bfp16 " << relup_delta_bfp16 << std::endl;
  }

  LatticeFermionD diff_o(FrbGrid);
  RealD diff = axpy_norm(diff_o, -1.0, result_o_16, result_o_full);
  std::cout << "Diff between results (s/c16): " << diff << std::endl;

  diff = axpy_norm(diff_o, -1.0, result_o_8, result_o_full);
  std::cout << "Diff between results (s/c8): " << diff << std::endl;

  diff = axpy_norm(diff_o, -1.0, result_o_half, result_o_full);
  std::cout << "Diff between results (s/h): " << diff << std::endl;

  diff = axpy_norm(diff_o, -1.0, result_o_bfp16, result_o_full);
  std::cout << "Diff between results (s/bfp16): " << diff << std::endl;

  if(algorithm == MixedCG || algorithm == ReliableUpdate){
    std::cout << "Iterations (s/c16) inner: " << inner_16 << " outer: " << outer_16 << " patchup: " << patchup_16 << std::endl;
    std::cout << "Iterations (s/c8) inner: " << inner_8 << " outer: " << outer_8 << " patchup: " << patchup_8 << std::endl;
    std::cout << "Iterations (s/h) inner: " << inner_half << " outer: " << outer_half << " patchup: " << patchup_half << std::endl;
    std::cout << "Iterations (s/s) inner: " << inner_full << " outer: " << outer_full << " patchup: " << patchup_full << std::endl;
    std::cout << "Iterations (s/bfp16) inner: " << inner_bfp16 << " outer: " << outer_bfp16 << " patchup: " << patchup_bfp16 << std::endl; 
 }else if(algorithm == PrecCG){
    std::cout << "Iterations (s/c16) inner: " << inner_16 << " outer: " << outer_16 << std::endl;
    std::cout << "Iterations (s/c8) inner: " << inner_8 << " outer: " << outer_8 << std::endl;
    std::cout << "Iterations (s/h) inner: " << inner_half << " outer: " << outer_half << std::endl;
    std::cout << "Iterations (s/s) inner: " << inner_full << " outer: " << outer_full << std::endl;
    std::cout << "Iterations (s/bfp16) inner: " << inner_bfp16 << " outer: " << outer_bfp16 << std::endl;
  }

  Grid_finalize();
}

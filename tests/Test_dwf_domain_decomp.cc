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
#include <Grid/Grid.h>

namespace Grid{

template <class Field>
class ConjugateGradientDomainDecomp : public OperatorFunction<Field> {
 public:
  bool ErrorOnNoConverge;  // throw an assert when the CG fails to converge.
                           // Defaults true.
  RealD Tolerance;
  Integer MaxIterations;
  Integer IterationsToComplete; //Number of iterations the CG took to finish. Filled in upon completion
  Integer CommsFrequency; //Iterations between comms. Use 0 for comms disabled, 1 for every iteration

  ConjugateGradientDomainDecomp(RealD tol, Integer maxit, bool err_on_no_conv = true)
      : Tolerance(tol),
        MaxIterations(maxit),
        ErrorOnNoConverge(err_on_no_conv), CommsFrequency(1){};

  void operator()(LinearOperatorBase<Field> &Linop, const Field &src, Field &psi) {

    psi.checkerboard = src.checkerboard;
    conformable(psi, src);

    RealD cp, c, a, d, b, ssq, qq, b_pred;

    Field p(src);
    Field mmp(src);
    Field r(src);

    // Initial residual computation & set up
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    //Even for comms disabled mode we should compute the initial residual with comms enabled so we actually know where we stand
    Linop.HermOpAndNorm(psi, mmp, d, b);
    

    r = src - mmp;
    p = r;

    a = norm2(p);
    cp = a;
    ssq = norm2(src);

    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient: guess " << guess << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient:   src " << ssq << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient:    mp " << d << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient:   mmp " << b << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient:  cp,r " << cp << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient:     p " << a << std::endl;

    RealD rsq = Tolerance * Tolerance * ssq;

    // Check if guess is really REALLY good :)
    if (cp <= rsq) {
      return;
    }

    if(CommsFrequency == 0) StencilControls::CommsEnabled() = false;

    std::cout << GridLogIterative << std::setprecision(8)
              << "ConjugateGradient: k=0 residual " << cp << " target " << rsq << std::endl;

    GridStopWatch LinalgTimer;
    GridStopWatch MatrixTimer;
    GridStopWatch SolverTimer;

    SolverTimer.Start();
    int k;
    for (k = 1; k <= MaxIterations; k++) {
      c = cp;

      if(CommsFrequency !=0) StencilControls::CommsEnabled() = (k % CommsFrequency == 0); //enable comms only on even multiples of the comms frequency

      MatrixTimer.Start();
      Linop.HermOpAndNorm(p, mmp, d, qq);
      MatrixTimer.Stop();

      LinalgTimer.Start();
      //  RealD    qqck = norm2(mmp);
      //  ComplexD dck  = innerProduct(p,mmp);

      a = c / d;
      b_pred = a * (a * qq - d) / c;

      cp = axpy_norm(r, -a, mmp, r);
      b = cp / c;

      // Fuse these loops ; should be really easy
      psi = a * p + psi;
      p = p * b + r;

      LinalgTimer.Stop();

      std::cout << GridLogIterative << "ConjugateGradient: Iteration " << k
                << " (comms ? " << StencilControls::CommsEnabled() << ") residual " << cp << " target " << rsq << std::endl;
      std::cout << GridLogDebug << "a = "<< a << " b_pred = "<< b_pred << "  b = "<< b << std::endl;
      std::cout << GridLogDebug << "qq = "<< qq << " d = "<< d << "  c = "<< c << std::endl;

      // Stopping condition
      if (cp <= rsq) {
        SolverTimer.Stop();

	StencilControls::CommsEnabled() = true;
        Linop.HermOpAndNorm(psi, mmp, d, qq);
        p = mmp - src;

        RealD srcnorm = sqrt(norm2(src));
        RealD resnorm = sqrt(norm2(p));
        RealD true_residual = resnorm / srcnorm;

        std::cout << GridLogMessage << "ConjugateGradient Converged on iteration " << k << std::endl;
        std::cout << GridLogMessage << "\tComputed residual " << sqrt(cp / ssq)<<std::endl;
	std::cout << GridLogMessage << "\tTrue residual " << true_residual<<std::endl;
	std::cout << GridLogMessage << "\tTarget " << Tolerance << std::endl;

        std::cout << GridLogMessage << "Time breakdown "<<std::endl;
	std::cout << GridLogMessage << "\tElapsed    " << SolverTimer.Elapsed() <<std::endl;
	std::cout << GridLogMessage << "\tMatrix     " << MatrixTimer.Elapsed() <<std::endl;
	std::cout << GridLogMessage << "\tLinalg     " << LinalgTimer.Elapsed() <<std::endl;

        if (ErrorOnNoConverge) assert(true_residual / Tolerance < 10000.0);

	IterationsToComplete = k;	

        return;
      }
    }
    std::cout << GridLogMessage << "ConjugateGradient did NOT converge"
              << std::endl;

    if (ErrorOnNoConverge) assert(0);
    IterationsToComplete = k;

    StencilControls::CommsEnabled() = true;
  }
};

}




using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template<class d>
struct scal {
  d internal;
};

  Gamma::Algebra Gmu [] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT
  };






int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int freq = 1;
  for(int i=1;i<argc;i++){
    if(std::string(argv[i]) == "-freq"){
      std::stringstream ss; ss << argv[i+1]; ss >> freq;
    }
  }

  const int Ls=8;

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
  
  SU3::HotConfiguration(RNG4,Umu);

  precisionChange(Umu_f,Umu);
  
  RealD mass=0.1;
  RealD M5=1.8;
  DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DomainWallFermionF Ddwf_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);

  LatticeFermionD    src_o(FrbGrid);
  LatticeFermionD result_o(FrbGrid);
  LatticeFermionD result_o_2(FrbGrid);
  pickCheckerboard(Odd,src_o,src);
  result_o.checkerboard = Odd;
  result_o = zero;
  result_o_2.checkerboard = Odd;
  result_o_2 = zero;

  SchurDiagMooeeOperator<DomainWallFermionD,LatticeFermionD> HermOpEO(Ddwf);
  SchurDiagMooeeOperator<DomainWallFermionF,LatticeFermionF> HermOpEO_f(Ddwf_f);

  std::cout << "Starting domain decomposed CG" << std::endl;
  ConjugateGradientDomainDecomp<LatticeFermionD> CGdd(1.0e-8,10000);
  CGdd.CommsFrequency = freq;
  CGdd.ErrorOnNoConverge = 0;
  CGdd(HermOpEO,src_o,result_o);

  //MixedPrecisionConjugateGradient<LatticeFermionD,LatticeFermionF> mCG(1.0e-8, 10000, 50, FrbGrid_f, HermOpEO_f, HermOpEO);
  //mCG(src_o,result_o);

  std::cout << "Starting regular CG" << std::endl;
  ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
  CG(HermOpEO,src_o,result_o_2);

  LatticeFermionD diff_o(FrbGrid);
  RealD diff = axpy_norm(diff_o, -1.0, result_o, result_o_2);

  std::cout << "Diff between domain decomp and regular CG: " << diff << std::endl;

  
  Grid_finalize();
}

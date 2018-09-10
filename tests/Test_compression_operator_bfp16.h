#ifndef TEST_COMPRESSION_OPERATOR_BFP16_H_
#define TEST_COMPRESSION_OPERATOR_BFP16_H_

#include <Grid/Grid.h>

#include "Test_compression_compressor_bfp16.h"
#include "Test_compression_wilson_stencil.h"

namespace Grid{
namespace QCD{

template<typename HS,typename S> using WilsonCompressorBfp16 = WilsonCompressorBfp16template<HS,S,WilsonProjector>;

template<class S>
class WilsonBfp16commsImpl: public WilsonImpl<S,FundamentalRepresentation,CoeffReal>{
public:
  typedef WilsonImpl<S,FundamentalRepresentation,CoeffReal> WilsonBase;

#define INHERIT_BASE(TYPE) typedef typename WilsonBase::TYPE TYPE

  INHERIT_BASE(Gimpl);
  INHERIT_GIMPL_TYPES(Gimpl);

  INHERIT_BASE(Coeff_t);
  
  INHERIT_BASE(SiteSpinor);
  INHERIT_BASE(SitePropagator);
  INHERIT_BASE(SiteHalfSpinor);
  INHERIT_BASE(SiteHalfCommSpinor);    
  INHERIT_BASE(SiteDoubledGaugeField);

  INHERIT_BASE(FermionField);
  INHERIT_BASE(PropagatorField);
  INHERIT_BASE(DoubledGaugeField);

  typedef WilsonCompressorBfp16<SiteHalfSpinor, SiteSpinor> Compressor;  

  INHERIT_BASE(ImplParams);
  typedef WilsonStencilBasic<SiteSpinor, SiteHalfSpinor> StencilImpl;

  WilsonBfp16commsImpl(const ImplParams &p = ImplParams()) : WilsonBase(p){}
  
  inline void multLink(SiteHalfSpinor &phi,
		       const SiteDoubledGaugeField &U,
		       const SiteHalfSpinor &chi,
		       int mu,
		       StencilEntry *SE,
		       StencilImpl &St) {
    mult(&phi(), &U(mu), &chi());
  }

#undef INHERIT_BASE
};

typedef WilsonBfp16commsImpl<vComplexF> WilsonBfp16commsImplF;
typedef WilsonBfp16commsImpl<vComplexD> WilsonBfp16commsImplD;

}}

//WilsonBfp16commsImplF
//WilsonBfp16commsImplD

// #define TO_INSTANTIATE \
//   DOIT(WilsonBfp16commsImpl<vComplexF>) \
//   DOIT(WilsonBfp16commsImpl<vComplexD>)

// #include "Test_compression_instantiateImpl.h"

//#undef TO_INSTANTIATE

namespace Grid{
namespace QCD{

  typedef DomainWallFermion<WilsonBfp16commsImplD> DomainWallFermionBfp16commsD;
  typedef DomainWallFermion<WilsonBfp16commsImplF> DomainWallFermionBfp16commsF;

}}


#endif

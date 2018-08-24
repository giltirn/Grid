#ifndef TEST_COMPRESSION_OPERATOR_FIXEDPOINTCOMMS_H_
#define TEST_COMPRESSION_OPERATOR_FIXEDPOINTCOMMS_H_

#include <Grid/Grid.h>

#include "Test_compression_compressor_fixedpoint.h"
#include "Test_compression_wilson_stencil.h"

namespace Grid{
namespace QCD{

template<typename HS,typename S, int packSize> using WilsonCompressorFixedPoint = WilsonCompressorFixedPointTemplate<HS,S,WilsonProjector,packSize>;

template<class S, int packSize = 16>
class WilsonFixedPointCommsImpl: public WilsonImpl<S,FundamentalRepresentation,CoeffReal>{
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

  typedef WilsonCompressorFixedPoint<SiteHalfSpinor, SiteSpinor, packSize> Compressor;  

  INHERIT_BASE(ImplParams);
  typedef WilsonStencilBasic<SiteSpinor, SiteHalfSpinor> StencilImpl;

  WilsonFixedPointCommsImpl(const ImplParams &p = ImplParams()) : WilsonBase(p){}
  
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

typedef WilsonFixedPointCommsImpl<vComplexF,8> WilsonFixedPointComms8ImplF;
typedef WilsonFixedPointCommsImpl<vComplexD,8> WilsonFixedPointComms8ImplD;
typedef WilsonFixedPointCommsImpl<vComplexF,16> WilsonFixedPointComms16ImplF;
typedef WilsonFixedPointCommsImpl<vComplexD,16> WilsonFixedPointComms16ImplD;

}}

#define TO_INSTANTIATE \
  DOIT(WilsonFixedPointComms8ImplF)\
  DOIT(WilsonFixedPointComms8ImplD)\
  DOIT(WilsonFixedPointComms16ImplF)\
  DOIT(WilsonFixedPointComms16ImplD)

#include "Test_compression_instantiateImpl.h"

#undef TO_INSTANTIATE

namespace Grid{
namespace QCD{

typedef DomainWallFermion<WilsonFixedPointComms8ImplD> DomainWallFermionFixedPointComms8D;
typedef DomainWallFermion<WilsonFixedPointComms8ImplF> DomainWallFermionFixedPointComms8F;
typedef DomainWallFermion<WilsonFixedPointComms16ImplD> DomainWallFermionFixedPointComms16D;
typedef DomainWallFermion<WilsonFixedPointComms16ImplF> DomainWallFermionFixedPointComms16F;

}}


#endif

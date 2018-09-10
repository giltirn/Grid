#ifndef TEST_COMPRESSION_COMPRESSOR_BFP16_H
#define TEST_COMPRESSION_COMPRESSOR_BFP16_H

#include <Grid/Grid.h>

#include "Test_compression_utils.h"

namespace Grid{
namespace QCD{

  
//bfp16 is a 16-bit floating point rep with 1 sign bit, 8 exponent bits and 7 mantissa bits

inline int16_t bfp16pack(float val){
  assert(isLittleEndian());
  int16_t out;
  memcpy(&out, ((char const*)&val) + 2, 2);
  return out;
}

inline float bfp16unpack(int16_t val){
  int32_t v32 = int32_t(val) << 16;
  return *( (float const*)&v32 );
}
 

//Compressor that compresses to a single magnitude and Nhs*Dimension fixed point integers of size packSize bits
template<class _Hspinor,class _Spinor, class projector>
class WilsonCompressorBfp16template
{
 public:
  
  int mu,dag;  

  void Point(int p) { mu=p; };

  WilsonCompressorBfp16template(int _dag=0){
    dag = _dag;
  }

  typedef _Spinor         SiteSpinor;
  typedef _Hspinor     SiteHalfSpinor;
  typedef _Hspinor SiteHalfCommSpinor;
  typedef typename SiteHalfSpinor::vector_type     vComplexIn;
  constexpr static int Nw=sizeof(SiteHalfSpinor)/sizeof(vComplexIn); //number of complex numbers in SiteHalfSpinor

  typedef typename SiteHalfSpinor::scalar_object ScalarSiteHalfSpinor;

  constexpr static int Nsimd = vComplexIn::Nsimd();
  constexpr static int Dimension = getHalfSpinorColors<SiteHalfSpinor>::value;
 
  typedef typename ScalarSiteHalfSpinor::scalar_type stype; //std::complex
  typedef typename stype::value_type srtype; //float/double

  typedef uint16_t packedType;

  //Pack and unpack *scalar* SiteHalfSpinor objects
  void packSpinor(void* tov, const ScalarSiteHalfSpinor &from){
    packedType* to = (packedType*)tov;
    packedType p;
    srtype q;
    for(int s=0;s<Nhs;s++)
      for(int c=0;c<Dimension;c++)
	for(int reim=0;reim<2;reim++){
	  q = cmplx_reim( from()(s)(c), reim );
	  *(to++) = bfp16pack(q);
	}
  }

  //Vectorized version, store contiguously
  void packSpinor(void* tov, const SiteHalfSpinor &from){
    packedType* to = (packedType*)tov;
    std::vector<ScalarSiteHalfSpinor> extracted(Nsimd);
    extract(from,extracted);

    static const int incr = Nhs*Dimension*2;

    for(int i=0;i<Nsimd;i++){
      packSpinor((void*)to, extracted[i]);
      to += incr;
    }
  }


  void unpackSpinor(ScalarSiteHalfSpinor &to, void* fromv){
    packedType* from = (packedType*)fromv;

    for(int s=0;s<Nhs;s++)
      for(int c=0;c<Dimension;c++)
	for(int reim=0;reim<2;reim++)
	  cmplx_reim( to()(s)(c), reim ) = bfp16unpack( *(from++) );	
  }

  void unpackSpinor(SiteHalfSpinor &to, void* fromv){
    packedType* from = (packedType*)fromv;
    std::vector<ScalarSiteHalfSpinor> unpacked(Nsimd);

    static const int incr = Nhs*Dimension*2;

    for(int i=0;i<Nsimd;i++){
      unpackSpinor(unpacked[i],(void*)from);
      from += incr;
    }

    merge(to,unpacked);
  }

  inline int CommDatumSize(void) {
    return Nsimd*(  Nhs*Dimension*2*sizeof(packedType) );
  }

  /*****************************************************/
  /* Compress includes precision change if mpi data is not same */
  /*****************************************************/
  void Compress(SiteHalfSpinor *buf,Integer o,const SiteSpinor &in) {
    SiteHalfSpinor hsp;
    projector::Proj(hsp,in,mu,dag);
    
    uint8_t* to = (uint8_t*)buf + o*CommDatumSize();
    packSpinor(to, hsp);
  }

  /*****************************************************/
  /* Exchange includes precision change if mpi data is not same */
  /*****************************************************/
  void Exchange(SiteHalfSpinor *mp,
                       SiteHalfSpinor *vp0,
                       SiteHalfSpinor *vp1,
		       Integer type,Integer o){
    uint8_t* vpp0 = (uint8_t*)vp0 + o*CommDatumSize();
    uint8_t* vpp1 = (uint8_t*)vp1 + o*CommDatumSize();

    SiteHalfSpinor vt0, vt1;
    unpackSpinor(vt0, vpp0);
    unpackSpinor(vt1, vpp1);

    exchange(mp[2*o],mp[2*o+1],vt0,vt1,type);
  }

  /*****************************************************/
  /* Have a decompression step if mpi data is not same */
  /*****************************************************/
  void Decompress(SiteHalfSpinor *out,
			 SiteHalfSpinor *in, Integer o) {    
    uint8_t* hin = (uint8_t*)in + o*CommDatumSize();
    unpackSpinor(out[o],hin);
  }

  /*****************************************************/
  /* Compress Exchange                                 */
  /*****************************************************/
  void CompressExchange(SiteHalfSpinor *out0,
			       SiteHalfSpinor *out1,
			       const SiteSpinor *in,
			       Integer j,Integer k, Integer m,Integer type){
    SiteHalfSpinor temp1, temp2,temp3,temp4;
    projector::Proj(temp1,in[k],mu,dag);
    projector::Proj(temp2,in[m],mu,dag);
    exchange(temp3,temp4,temp1,temp2,type);

    uint8_t* hout0 = (uint8_t*)out0 + j*CommDatumSize();
    uint8_t* hout1 = (uint8_t*)out1 + j*CommDatumSize();
    packSpinor(hout0, temp3);
    packSpinor(hout1, temp4);
  }

  /*****************************************************/
  /* Pass the info to the stencil */
  /*****************************************************/
  inline bool DecompressionStep(void) { return true; }

};


}}

#endif

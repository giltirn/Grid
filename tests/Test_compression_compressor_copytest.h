#ifndef TEST_COMPRESSION_COMPRESSOR_COPY_TEST_H
#define TEST_COMPRESSION_COMPRESSOR_COPY_TEST_H

#include <Grid/Grid.h>

namespace Grid{
namespace QCD{

//Basic copy of WilsonCompressor for demonstration
template<class _Hspinor,class _Spinor, class projector>
class WilsonCompressorCopyTemplate
{
 public:
  
  int mu,dag;  

  void Point(int p) { mu=p; };

  WilsonCompressorCopyTemplate(int _dag=0){
    dag = _dag;
  }

  typedef _Spinor         SiteSpinor;
  typedef _Hspinor     SiteHalfSpinor;
  typedef _Hspinor SiteHalfCommSpinor;
  typedef typename SiteHalfSpinor::vector_type     vComplexIn;
  constexpr static int Nw=sizeof(SiteHalfSpinor)/sizeof(vComplexIn); //number of complex numbers in SiteHalfSpinor

  inline int CommDatumSize(void) {
    return sizeof(SiteHalfCommSpinor);
  }

  /*****************************************************/
  /* Compress includes precision change if mpi data is not same */
  /*****************************************************/
  inline void Compress(SiteHalfSpinor *buf,Integer o,const SiteSpinor &in) {
    projector::Proj(buf[o],in,mu,dag);
  }

  /*****************************************************/
  /* Exchange includes precision change if mpi data is not same */
  /*****************************************************/
  inline void Exchange(SiteHalfSpinor *mp,
                       SiteHalfSpinor *vp0,
                       SiteHalfSpinor *vp1,
		       Integer type,Integer o){
    exchange(mp[2*o],mp[2*o+1],vp0[o],vp1[o],type);
  }

  /*****************************************************/
  /* Have a decompression step if mpi data is not same */
  /*****************************************************/
  inline void Decompress(SiteHalfSpinor *out,
			 SiteHalfSpinor *in, Integer o) {    
    assert(0);
  }

  /*****************************************************/
  /* Compress Exchange                                 */
  /*****************************************************/
  inline void CompressExchange(SiteHalfSpinor *out0,
			       SiteHalfSpinor *out1,
			       const SiteSpinor *in,
			       Integer j,Integer k, Integer m,Integer type){
    SiteHalfSpinor temp1, temp2,temp3,temp4;
    projector::Proj(temp1,in[k],mu,dag);
    projector::Proj(temp2,in[m],mu,dag);
    exchange(out0[j],out1[j],temp1,temp2,type);
  }

  /*****************************************************/
  /* Pass the info to the stencil */
  /*****************************************************/
  inline bool DecompressionStep(void) { return false; }

};

}}

#endif

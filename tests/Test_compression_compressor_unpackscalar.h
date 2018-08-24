#ifndef TEST_COMPRESSION_COMPRESSOR_UNPACK_SCALAR_H
#define TEST_COMPRESSION_COMPRESSOR_UNPACK_SCALAR_H

#include <Grid/Grid.h>

namespace Grid{
namespace QCD{


//Compressor that unpacks vectorized data to scalar
template<class _Hspinor,class _Spinor, class projector>
class WilsonCompressorUnpackScalarTemplate
{
 public:
  
  int mu,dag;  

  void Point(int p) { mu=p; };

  WilsonCompressorUnpackScalarTemplate(int _dag=0){
    dag = _dag;
  }

  typedef _Spinor         SiteSpinor;
  typedef _Hspinor     SiteHalfSpinor;
  typedef _Hspinor SiteHalfCommSpinor;
  typedef typename SiteHalfSpinor::vector_type     vComplexIn;
  constexpr static int Nw=sizeof(SiteHalfSpinor)/sizeof(vComplexIn); //number of complex numbers in SiteHalfSpinor

  typedef typename SiteHalfSpinor::scalar_object ScalarSiteHalfSpinor;

  constexpr static int Nsimd = vComplexIn::Nsimd();

  inline int CommDatumSize(void) {
    return Nsimd*sizeof(ScalarSiteHalfSpinor);
  }

  /*****************************************************/
  /* Compress includes precision change if mpi data is not same */
  /*****************************************************/
  inline void Compress(SiteHalfSpinor *buf,Integer o,const SiteSpinor &in) {
    SiteHalfSpinor hsp;
    projector::Proj(hsp,in,mu,dag);
    
    ScalarSiteHalfSpinor* to = (ScalarSiteHalfSpinor*)buf + o*Nsimd;
    
    std::vector<ScalarSiteHalfSpinor*> extract_args(Nsimd);
    for(int i=0;i<Nsimd;i++) extract_args[i] = to+i;
    extract1(hsp,extract_args,0);
  }

  /*****************************************************/
  /* Exchange includes precision change if mpi data is not same */
  /*****************************************************/
  inline void Exchange(SiteHalfSpinor *mp,
                       SiteHalfSpinor *vp0,
                       SiteHalfSpinor *vp1,
		       Integer type,Integer o){
    ScalarSiteHalfSpinor* vpp0 = (ScalarSiteHalfSpinor*)vp0 + o*Nsimd;
    ScalarSiteHalfSpinor* vpp1 = (ScalarSiteHalfSpinor*)vp1 + o*Nsimd;
    
    std::vector<ScalarSiteHalfSpinor*> merge_args0(Nsimd), merge_args1(Nsimd);
    for(int i=0;i<Nsimd;i++){
      merge_args0[i] = vpp0+i;
      merge_args1[i] = vpp1+i;
    }

    SiteHalfSpinor vt0,vt1;
    merge1(vt0,merge_args0,0);
    merge1(vt1,merge_args1,0);

    exchange(mp[2*o],mp[2*o+1],vt0,vt1,type);
  }

  /*****************************************************/
  /* Have a decompression step if mpi data is not same */
  /*****************************************************/
  inline void Decompress(SiteHalfSpinor *out,
			 SiteHalfSpinor *in, Integer o) {    
    ScalarSiteHalfSpinor* hin = (ScalarSiteHalfSpinor*)in + o*Nsimd;
    std::vector<ScalarSiteHalfSpinor*> merge_args(Nsimd);
    for(int i=0;i<Nsimd;i++) merge_args[i] = hin+i;
    merge1(out[o],merge_args,0);
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
    exchange(temp3,temp4,temp1,temp2,type);

    ScalarSiteHalfSpinor* hout0 = (ScalarSiteHalfSpinor*)out0 + j*Nsimd;
    ScalarSiteHalfSpinor* hout1 = (ScalarSiteHalfSpinor*)out1 + j*Nsimd;

    std::vector<ScalarSiteHalfSpinor*> extract_args0(Nsimd), extract_args1(Nsimd);
    for(int i=0;i<Nsimd;i++){
      extract_args0[i] = hout0+i;
      extract_args1[i] = hout1+i;
    }
    extract1(temp3,extract_args0,0);
    extract1(temp4,extract_args1,0);
  }

  /*****************************************************/
  /* Pass the info to the stencil */
  /*****************************************************/
  inline bool DecompressionStep(void) { return true; }

};


}}

#endif

#ifndef TEST_COMPRESSION_COMPRESSOR_FIXEDPOINT_H_
#define TEST_COMPRESSION_COMPRESSOR_FIXEDPOINT_H_

#include <Grid/Grid.h>

namespace Grid{
namespace QCD{


//Access elements of std::complex
template<typename T>
inline T & cmplx_reim(std::complex<T> &c, const int reim){
  return reinterpret_cast<T(&)[2]>(c)[reim];
}

template<typename T>
inline const T & cmplx_reim(const std::complex<T> &c, const int reim){
  return reinterpret_cast<const T(&)[2]>(c)[reim];
}


//Pack and unpack float/double to fixed point representation of SZ bits
template<int SZ>
struct signedIntMap{};

template<>
struct signedIntMap<8>{ typedef int8_t type; };
template<>
struct signedIntMap<16>{ typedef int16_t type; };


template<typename T, int SZ>
inline typename signedIntMap<SZ>::type packN(T val){
  return typename signedIntMap<SZ>::type( (1<<(SZ-2) ) * val );
}
template<typename T, int SZ>
inline T unpackN(typename signedIntMap<SZ>::type val){
  return T(val)/(1<<(SZ-2));
}

template<typename T>
struct getHalfSpinorColors{
  //template <typename vtype> using iImplHalfSpinor        = iScalar<iVector<iVector<vtype, Dimension>, Nhs> >;
  enum { value = sizeof(typename T::element::element)/sizeof(typename T::element::element::element) };
};

//Compressor that compresses to a single magnitude and Nhs*Dimension fixed point integers of size packSize bits
template<class _Hspinor,class _Spinor, class projector, int packSize = 16>
class WilsonCompressorFixedPointTemplate
{
 public:
  
  int mu,dag;  

  void Point(int p) { mu=p; };

  WilsonCompressorFixedPointTemplate(int _dag=0){
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

  //Pack and unpack *scalar* SiteHalfSpinor objects
  void packSpinor(void* tov, const ScalarSiteHalfSpinor &from){
    uint8_t* to = (uint8_t*)tov;
    typedef typename signedIntMap<packSize>::type packedType;

    srtype max = 0;
    for(int s=0;s<Nhs;s++)
      for(int c=0;c<Dimension;c++)
	for(int reim=0;reim<2;reim++)
	  if(fabs(cmplx_reim( from()(s)(c), reim )) > max )
	    max =  fabs(cmplx_reim( from()(s)(c), reim )) ;
  
    *( (srtype*)to ) = max; //copy the normalization to the buffer
    to += sizeof(srtype);
  
    packedType *top = (packedType*)to;
    packedType p;
    srtype q;
    for(int s=0;s<Nhs;s++)
      for(int c=0;c<Dimension;c++)
	for(int reim=0;reim<2;reim++){
	  q = cmplx_reim( from()(s)(c), reim );
	  if(max != 0.) q /= max;
	  *(top++) = packN<srtype,packSize>(q);
	}
  }

  void packSpinor(void* tov, const SiteHalfSpinor &from){
    uint8_t* to = (uint8_t*)tov;
    std::vector<ScalarSiteHalfSpinor> extracted(Nsimd);
    extract(from,extracted);

    static const int incr = sizeof(srtype) + Nhs*Dimension*2*sizeof(typename signedIntMap<packSize>::type);

    for(int i=0;i<Nsimd;i++){
      packSpinor((void*)to, extracted[i]);
      to += incr;
    }
  }


  void unpackSpinor(ScalarSiteHalfSpinor &to, void* fromv){
    uint8_t* from = (uint8_t*)fromv;
    typedef typename signedIntMap<packSize>::type packedType;

    srtype norm = *( (srtype*)from ); 
    from += sizeof(srtype);

    packedType *fromp = (packedType*)from;
    srtype q;
    for(int s=0;s<Nhs;s++)
      for(int c=0;c<Dimension;c++)
	for(int reim=0;reim<2;reim++){
	  q = unpackN<srtype,packSize>(*(fromp++) );
	  if(norm != 0.) q *= norm;
	  cmplx_reim( to()(s)(c), reim ) = q;
	}
  }

  void unpackSpinor(SiteHalfSpinor &to, void* fromv){
    uint8_t* from = (uint8_t*)fromv;
    std::vector<ScalarSiteHalfSpinor> unpacked(Nsimd);

    static const int incr = sizeof(srtype) + Nhs*Dimension*2*sizeof(typename signedIntMap<packSize>::type);

    for(int i=0;i<Nsimd;i++){
      unpackSpinor(unpacked[i],(void*)from);
      from += incr;
    }

    merge(to,unpacked);
  }

  inline int CommDatumSize(void) {
    return Nsimd*(  sizeof(srtype) + Nhs*Dimension*2*sizeof(typename signedIntMap<packSize>::type) );
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

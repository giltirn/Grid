#ifndef TEST_COMPRESSION_UTILS_H
#define TEST_COMPRESSION_UTILS_H

#include <Grid/Grid.h>

namespace Grid{
namespace QCD{


inline int isLittleEndian(){
  unsigned int x = 1;
  char *c = (char*) &x;
  return (int)*c;
}

template<typename T>
struct getHalfSpinorColors{
  //template <typename vtype> using iImplHalfSpinor        = iScalar<iVector<iVector<vtype, Dimension>, Nhs> >;
  enum { value = sizeof(typename T::element::element)/sizeof(typename T::element::element::element) };
};

//Access elements of std::complex
template<typename T>
inline T & cmplx_reim(std::complex<T> &c, const int reim){
  return reinterpret_cast<T(&)[2]>(c)[reim];
}

template<typename T>
inline const T & cmplx_reim(const std::complex<T> &c, const int reim){
  return reinterpret_cast<const T(&)[2]>(c)[reim];
}

}
}

#endif

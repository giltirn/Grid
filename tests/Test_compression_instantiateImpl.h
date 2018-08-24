#include <Grid/qcd/action/fermion/CayleyFermion5DmethodImpl.h>

namespace Grid{
  namespace QCD{
#define DOIT(A) template class CayleyFermion5D<A>;
    TO_INSTANTIATE;
#undef DOIT
  }
}

#include <Grid/qcd/action/fermion/CayleyFermion5DcacheMethodImpl.h>

namespace Grid{
  namespace QCD{   
#define DOIT(A) INSTANTIATE_DPERP(A)
    TO_INSTANTIATE;
#undef DOIT
  }
}

#include <Grid/qcd/action/fermion/WilsonFermion5DmethodImpl.h>

namespace Grid{
  namespace QCD{
#define DOIT(A) template class WilsonFermion5D<A>;
    TO_INSTANTIATE;
#undef DOIT
  }
}

#include <Grid/qcd/action/fermion/WilsonKernelsMethodImpl.h>

namespace Grid{
  namespace QCD{
#define DOIT(A) template class WilsonKernels<A>;
    TO_INSTANTIATE;
#undef DOIT
  }
}

#include <Grid/qcd/action/fermion/WilsonKernelsHandMethodImpl.h>

namespace Grid{
  namespace QCD{
#define DOIT(A) INSTANTIATE_THEM(A);
    TO_INSTANTIATE;
#undef DOIT
  }
}

#include <Grid/qcd/action/fermion/WilsonKernelsAsmMethodImpl.h>

namespace Grid{
  namespace QCD{
#define DOIT(A) INSTANTIATE_ASM(A);
    TO_INSTANTIATE;
#undef DOIT
  }
}

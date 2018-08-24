#ifndef TEST_COMPRESSION_COMPRESSOR_WILSON_STENCIL_H
#define TEST_COMPRESSION_COMPRESSOR_WILSON_STENCIL_H

#include <Grid/Grid.h>

namespace Grid{
namespace QCD{

  //The Grid standard Wilson stencil has the compressor baked-in in optimized comms routines in order to avoid a switch over direction index. This makes it harder to generalize
  //as we need to generate 8 instantiations of the compressor. This is a simple Wilson stencil that does not perform this optimization

template<class vobj,class cobj>
class WilsonStencilBasic : public CartesianStencil<vobj,cobj> {
public:
  double timer0;
  double timer1;
  double timer2;
  double timer3;
  double timer4;
  double timer5;
  double timer6;
  uint64_t callsi;
  void ZeroCountersi(void)
  {
    timer0=0;
    timer1=0;
    timer2=0;
    timer3=0;
    timer4=0;
    timer5=0;
    timer6=0;
    callsi=0;
  }
  void Reporti(int calls)
  {
    if ( timer0 ) std::cout << GridLogMessage << " timer0 (HaloGatherOpt) " <<timer0/calls <<std::endl;
    if ( timer1 ) std::cout << GridLogMessage << " timer1 (Communicate)   " <<timer1/calls <<std::endl;
    if ( timer2 ) std::cout << GridLogMessage << " timer2 (CommsMerge )   " <<timer2/calls <<std::endl;
    if ( timer3 ) std::cout << GridLogMessage << " timer3 (commsMergeShm) " <<timer3/calls <<std::endl;
    if ( timer4 ) std::cout << GridLogMessage << " timer4 " <<timer4 <<std::endl;
  }


  std::vector<int> same_node;
  std::vector<int> surface_list;

  WilsonStencilBasic(GridBase *grid,
		int npoints,
		int checkerboard,
		const std::vector<int> &directions,
		const std::vector<int> &distances)  
    : CartesianStencil<vobj,cobj> (grid,npoints,checkerboard,directions,distances) ,
    same_node(npoints)
  { 
    ZeroCountersi();
    surface_list.resize(0);
  };

  void BuildSurfaceList(int Ls,int vol4){

    // find same node for SHM
    // Here we know the distance is 1 for WilsonStencil
    for(int point=0;point<this->_npoints;point++){
      same_node[point] = this->SameNode(point);
    }
    
    for(int site = 0 ;site< vol4;site++){
      int local = 1;
      for(int point=0;point<this->_npoints;point++){
	if( (!this->GetNodeLocal(site*Ls,point)) && (!same_node[point]) ){ 
	  local = 0;
	}
      }
      if(local == 0) { 
	surface_list.push_back(site);
      }
    }
  }


  template < class compressor>
  void HaloExchangeOpt(const Lattice<vobj> &source,compressor &compress) 
  {
    std::vector<std::vector<CommsRequest_t> > reqs;
    this->HaloExchangeOptGather(source,compress);
    double t1=usecond();
    this->Communicate();
    double t2=usecond(); timer1 += t2-t1;
    this->CommsMerge(compress);
    double t3=usecond(); timer2 += t3-t2;
    this->CommsMergeSHM(compress);
    double t4=usecond(); timer3 += t4-t3;
  }
  
  template <class compressor>
  void HaloExchangeOptGather(const Lattice<vobj> &source,compressor &compress){
    this->Prepare();
    double t0=usecond();
    this->HaloGatherOpt(source,compress);
    double t1=usecond();
    timer0 += t1-t0;
    callsi++;
  }

  template <class compressor>
  void HaloGatherOpt(const Lattice<vobj> &source,compressor &compress)
  {
    this->halogtime-=usecond();
    this->HaloGather(source,compress);
    this->halogtime+=usecond();
  }
};



}}

#endif

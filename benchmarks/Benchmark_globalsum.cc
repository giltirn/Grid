    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_comms.cc

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

using namespace std;
using namespace Grid;

void stats(double &mean, double &std, double &median, double &min, double &max, double &iqr_lo, double &iqr_hi, const std::vector<double> &vals){
  std::vector<double> vs(vals);
  std::sort(vs.begin(),vs.end());

  size_t sz = vs.size();
  median = vs[sz/2];
  min = vs[0];
  max = vs[sz-1];

  iqr_lo = vs[ int(floor( 0.1587*sz )) ];
  iqr_hi = vs[ int(ceil( 0.8413*sz )) ];
  
  mean = 0;
  std = 0;
  for(double v: vals){
    mean += v;
    std += v*v;
  }
  mean = mean / sz;
  std = sqrt( std/sz - mean*mean );
}

void normalMPIreduction(double* ptr, size_t N, bool data_on_device, bool allow_acc_aware_mpi, CartesianCommunicator &comm){
  bool enable_acc_aware_mpi = false;
#ifdef ACCELERATOR_AWARE_MPI
  enable_acc_aware_mpi = allow_acc_aware_mpi;
#endif
  double* buf = ptr;
  size_t bytes = N*sizeof(double);

  if(data_on_device && !enable_acc_aware_mpi){
    buf = (double*)acceleratorAllocHost(bytes);
    acceleratorCopyFromDevice(ptr, buf, bytes);
  }

  MPI_Allreduce(MPI_IN_PLACE, buf, N, MPI_DOUBLE, MPI_SUM, comm.communicator);

  if(data_on_device && !enable_acc_aware_mpi){
    acceleratorCopyToDevice(buf, ptr, bytes);
    acceleratorFreeHost(buf);
  }
}

template<typename Implementation>
void benchmarkAllReduce(const std::vector<int> &sizes_MB, CartesianCommunicator &comm, bool device_ptr, bool allow_acc_aware_mpi, int nrpt){
  Implementation impl;
  impl.setup(comm);
  
  { //test the full reduction
    int ranks = comm.ProcessorCount();

    size_t N = 4*1024*1024; //32MB
    size_t bytes = N*sizeof(double);
    
    std::vector<double> real_expect(N,0.);
    hostVector<double> data_host(N);
    for(size_t i=0;i<N;i++){
      data_host[i] = sin(0.1 + i*comm.ThisRank());
      for(int r=0;r<ranks;r++){
	real_expect[i] += sin(0.1 + i*r);
      }
    }    
    hostVector<double> got_host(N), expect_host(N);

    if(device_ptr){
      deviceVector<double> got_device(N), expect_device(N);
      acceleratorCopyToDevice(data_host.data(), got_device.data(),bytes);
      acceleratorCopyToDevice(data_host.data(), expect_device.data(),bytes);

      impl.reduce(got_device.data(), N, true, allow_acc_aware_mpi);
      normalMPIreduction(expect_device.data(), N, true, allow_acc_aware_mpi, comm);
      
      acceleratorCopyFromDevice(got_device.data(), got_host.data(), bytes);
      acceleratorCopyFromDevice(expect_device.data(), expect_host.data(), bytes);
    }else{
      got_host = data_host;
      expect_host = data_host;

      impl.reduce(got_host.data(), N, false, allow_acc_aware_mpi);
      normalMPIreduction(expect_host.data(), N, false, allow_acc_aware_mpi, comm);
    }

    bool fail =false;
    for(int i=0;i<N;i++){
      if ( fabs(expect_host[i]-got_host[i])/expect_host[i] > 1e-7 ){
	std::cout << "Check " << i << " got " << got_host[i] << " expect " << expect_host[i] << " diff " << got_host[i]-expect_host[i] << " real expect " << real_expect[i] << std::endl;
	fail =true;
      }
    }
    assert(!fail);
  }

  for(int size_MB : sizes_MB){
    size_t bytes = size_t(size_MB) * 1024*1024;
    size_t ndouble = bytes / sizeof(double);

    assert( (size_t) ((int)ndouble) == ndouble ); //make sure it can be downcast to int

    double* buf_base;
    if(device_ptr)
      buf_base = (double*)acceleratorAllocDevice(bytes);
    else
      buf_base = (double*)malloc(bytes);

    std::vector<double> time_base(nrpt);
    {      
      for(int n=0;n<nrpt;n++){
	double dt = -usecond();			       
	normalMPIreduction(buf_base, ndouble, device_ptr, allow_acc_aware_mpi, comm);
	dt += usecond();
	time_base[n] = dt / 1000000;
      }
    }

    if(device_ptr)
      acceleratorFreeDevice(buf_base);
    else free(buf_base);

    double mu_base, std_base, med_base, min_base, max_base, iqr_lo_base, iqr_hi_base;
    stats(mu_base, std_base, med_base, min_base, max_base, iqr_lo_base, iqr_hi_base, time_base);
       
    double MB_per_s_base = size_MB / mu_base;
    double MB_per_s_base_err = size_MB/mu_base/mu_base * std_base;
    
    double* buf_new;
    if(device_ptr)
      buf_new = (double*)acceleratorAllocDevice(bytes);
    else
      buf_new = (double*)malloc(bytes);

    std::vector<double> time_new(nrpt);
    {
      for(int n=0;n<nrpt;n++){
	double dt = -usecond();
	impl.reduce(buf_new,ndouble,device_ptr,allow_acc_aware_mpi);
	dt += usecond();
	time_new[n] = dt / 100000;	
      }
    }

    double mu_new, std_new, med_new, min_new, max_new, iqr_lo_new, iqr_hi_new;
    stats(mu_new, std_new, med_new, min_new, max_new, iqr_lo_new, iqr_hi_new, time_new);

    double MB_per_s_new = size_MB / mu_new;
    double MB_per_s_new_err = size_MB/mu_new/mu_new * std_new;
    
    if(device_ptr)
      acceleratorFreeDevice(buf_new);
    else free(buf_new);

    std::cout << "Size " << size_MB << " MB" << std::endl;
    std::cout << "Base time mu=" << mu_base << " std=" << std_base << " med=" << med_base << " iqr_lo=" << iqr_lo_base << " iqr_hi=" << iqr_hi_base << " min=" << min_base << " max=" << max_base << "  rate=" << MB_per_s_base << "+-" << MB_per_s_base_err << std::endl;
    std::cout << "New  time mu=" << mu_new << " std=" << std_new << " med=" << med_new << " iqr_lo=" << iqr_lo_new << " iqr_hi=" << iqr_hi_new << " min=" << min_new << " max=" << max_new << "  rate=" << MB_per_s_new << "+-" << MB_per_s_new_err << std::endl << std::endl;   
  }
}


struct ImplRing{
  CartesianCommunicator* comm;

  void setup(CartesianCommunicator &_comm){
    comm = &_comm; 
  }
  void reduce(double* p, size_t N, bool on_device, bool allow_acc_aware_mpi){
    comm->GlobalSumVectorRing(p, N, on_device, comm->communicator, allow_acc_aware_mpi);
  }
};

void benchmarkAllReduceRing(const std::vector<int> &sizes_MB, CartesianCommunicator &comm, bool device_ptr, bool allow_acc_aware_mpi, const int nrpt){
  benchmarkAllReduce<ImplRing>(sizes_MB, comm, device_ptr, allow_acc_aware_mpi, nrpt);
}
 
struct ImplSharedRing{
  CartesianCommunicator* comm;

  void setup(CartesianCommunicator &_comm){
    comm = &_comm;    
  }
  void reduce(double* p, size_t N, bool on_device, bool allow_acc_aware_mpi){
    comm->GlobalSumVectorRingShared(p, N, on_device, comm->communicator, allow_acc_aware_mpi);
  }
};

void benchmarkAllReduceSharedRing(const std::vector<int> &sizes_MB, CartesianCommunicator &comm, bool device_ptr, bool allow_acc_aware_mpi, const int nrpt){
  benchmarkAllReduce<ImplSharedRing>(sizes_MB, comm, device_ptr, allow_acc_aware_mpi, nrpt);
};


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  std::vector<int> sizes_MB({1,10,50,100,150,200}); //reduction sizes in MB
  int nrpt = 50; //how many measurements to take

  std::string arg;
  if( GridCmdOptionExists(argv,argv+argc,"--sizesMB") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--mpi");
    GridCmdOptionIntVector(arg,sizes_MB);
  }    
  if( GridCmdOptionExists(argv,argv+argc,"--nrpt") ){
    std::vector<int> tmp(0);
    arg= GridCmdOptionPayload(argv,argv+argc,"--nrpt");
    GridCmdOptionIntVector(arg,tmp);
    assert(tmp.size()==1);
    nrpt = tmp[0];
  }

#ifdef ACCELERATOR_AWARE_MPI
  std::cout<<GridLogMessage << "Accelerator-aware MPI is available" << std::endl;
#endif

  CartesianCommunicator comm(mpi_layout);

  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking ring reduction on host"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  
  benchmarkAllReduceRing(sizes_MB, comm, false, true, nrpt);

  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking ring reduction on device with accelerator-aware MPI enabled"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;

  benchmarkAllReduceRing(sizes_MB, comm, true, true, nrpt);

  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking ring reduction on device with accelerator-aware MPI disabled"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;

  benchmarkAllReduceRing(sizes_MB, comm, true, false, nrpt);


  if(Enable_shared_mem_buffer){
    std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
    std::cout<<GridLogMessage << "= Benchmarking shared ring reduction on host"<<std::endl;
    std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  
    benchmarkAllReduceSharedRing(sizes_MB, comm, false, true, nrpt);

    std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
    std::cout<<GridLogMessage << "= Benchmarking shared ring reduction on device with accelerator-aware MPI enabled"<<std::endl;
    std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;

    benchmarkAllReduceSharedRing(sizes_MB, comm, true, true, nrpt);

    std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
    std::cout<<GridLogMessage << "= Benchmarking shared ring reduction on device with accelerator-aware MPI disabled"<<std::endl;
    std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;

    benchmarkAllReduceSharedRing(sizes_MB, comm, true, false, nrpt);
  }else{
    std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
    std::cout<<GridLogMessage << "= Skipped benchmarking shared ring reduction because the shared memory buffer is disabled"<<std::endl;
    std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  }


  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= All done; Bye Bye"<<std::endl;
  std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;

  Grid_finalize();
}

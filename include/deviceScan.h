#define DSORT LIBP_DIR"/libs/core/"

#ifndef DEVICE_SCAN
#define DEVICE_SCAN 1

// was 512
#define SCAN_BLOCK_SIZE 1024


//template <class T>
class deviceScan_t{
private:
  occa::kernel blockShflScanKernel;
  occa::kernel finalizeScanKernel;
  occa::kernel findStartsKernel;
  occa::kernel trashCompactorKernel;
public:

#if 1
  // constructor
  deviceScan_t(occa::device &device, const char *entryType, const char *entryMap, occa::properties props);

  dlong blockCount(dlong entries);

  void scan(const dlong    entries,// might need to upgrade
	    occa::memory &o_list,
	    occa::memory &o_tmp,
	    dlong *h_tmp,
	    occa::memory &o_scan);

  void mallocTemps(occa::device &device, dlong entries, occa::memory &o_tmp, dlong **h_tmp);

  dlong trashCompactor(occa::device &device,
		       const dlong entries,
		       const int entrySize,
		       occa::memory &o_list,
		       occa::memory &o_compactedList);

  
#else

  
  dlong blockCount(const dlong entries){
    return ((entries+SCAN_BLOCK_SIZE-1)/SCAN_BLOCK_SIZE);  
  }
  
  void scan(const int    entries,
	    occa::memory &o_list,
	    occa::memory &o_tmp,
	    T *h_tmp,
	    occa::memory &o_scan) {

    dlong Nblocks = blockCount(entries);
    
    // 1. launch DEVICE block scan kernel
    blockShflScanKernel(entries, o_list, o_scan, o_tmp);
    
    // 2. copy offsets back to HOST
    o_tmp.copyTo(h_tmp);

    // 3. scan offsets on HOST
    for(int n=1;n<Nblocks;++n){
      h_tmp[n] += h_tmp[n-1];
    }
  
    // 4. copy scanned offsets back to DEVCE
    o_tmp.copyFrom(h_tmp);

    // 5. finalize scan
    finalizeScanKernel(entries, o_tmp, o_scan);
  
  }

  void  mallocTemps(occa::device &device, dlong entries, occa::memory &o_tmp, T **h_tmp){
  
    size_t sz =  sizeof(T)*blockCount(entries);
    o_tmp = device.malloc(sz);
    *h_tmp = (T*) malloc(sz);
  }


  deviceScan_t(occa::device &device, const char *entryType, const char *entryMap, occa::properties props){

    // Compile the kernel at run-time
    occa::settings()["kernel/verbose"] = true;

    occa::properties kernelInfo = props;
  
    kernelInfo["includes"] += entryType; // "entry.h";
    kernelInfo["includes"] += entryMap;  // "compareEntry.h";
    kernelInfo["defines/SCAN_BLOCK_SIZE"] = (int)SCAN_BLOCK_SIZE;
  
    blockShflScanKernel = device.buildKernel(LIBCORE_DIR "/okl/blockShflScan.okl",
					     "blockShflScan", kernelInfo);

    finalizeScanKernel = device.buildKernel(LIBCORE_DIR "/okl/blockShflScan.okl",
					    "finalizeScan", kernelInfo);
  

  }
#endif
};
  
#endif


/* The MIT License (MIT)
 *
 * Copyright (c) 2014-2018 David Medina and Tim Warburton
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 */

#include "mesh.hpp"

// was 512
#define SCAN_BLOCK_SIZE 1024


//template <class T>
dlong deviceScan_t::blockCount(const dlong entries){
  return ((entries+SCAN_BLOCK_SIZE-1)/SCAN_BLOCK_SIZE);  
}

void deviceScan_t::scan(const dlong   entries,
			occa::memory &o_list,
			occa::memory &o_tmp,
			dlong *h_tmp,
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

dlong deviceScan_t::trashCompactor(platform_t &platform,
				   const dlong entries,
				   const int entrySize,
				   occa::memory &o_list,
				   occa::memory &o_compactedList){


  // set up some temprs
  dlong  *h_tmp;
  occa::memory o_tmp;
  mallocTemps(platform, entries, o_tmp, &h_tmp);
  
  // scan
  occa::memory o_scan = platform.device.malloc(entries*sizeof(dlong));
  scan(entries, o_list, o_tmp, h_tmp, o_scan);

  // find number of unique values
  dlong Nstarts = 0;
  (o_scan+sizeof(dlong)*(entries-1)).copyTo(&Nstarts);
  ++Nstarts;

  // find starts
  occa::memory o_starts = platform.device.malloc(Nstarts*sizeof(dlong));
  findStartsKernel(entries, o_scan, o_starts);

  // compactify duplicate entries  
  o_compactedList = platform.device.malloc(Nstarts*entrySize);
  trashCompactorKernel(entries, Nstarts, o_starts, o_list, o_compactedList);

  // tidy up
  free(h_tmp);
		      
  return Nstarts;
		      
}
			     

void  deviceScan_t::mallocTemps(platform_t &platform, dlong entries, occa::memory &o_tmp, dlong **h_tmp){
  
  size_t sz =  sizeof(dlong)*blockCount(entries);
  o_tmp = platform.device.malloc(sz);
  *h_tmp = (dlong*) malloc(sz);
}


deviceScan_t::deviceScan_t(platform_t &platform, const char *entryType, const char *entryMap, occa::properties props){

  // Compile the kernel at run-time
  occa::settings()["kernel/verbose"] = true;

  occa::properties kernelInfo = props;
  
  kernelInfo["includes"] += entryType; // "entry.h";
  kernelInfo["includes"] += entryMap;  // "compareEntry.h";
  kernelInfo["defines/SCAN_BLOCK_SIZE"] = (int)SCAN_BLOCK_SIZE;
  
  blockShflScanKernel = platform.buildKernel(LIBCORE_DIR "/okl/blockShflScan.okl",
				      "blockShflScan", kernelInfo);

  finalizeScanKernel = platform.buildKernel(LIBCORE_DIR "/okl/blockShflScan.okl",
				      "finalizeScan", kernelInfo);
  
  findStartsKernel = platform.buildKernel(LIBCORE_DIR "/okl/blockShflScan.okl",
					"findStarts", kernelInfo);

  
  trashCompactorKernel = platform.buildKernel(LIBCORE_DIR "/okl/blockShflScan.okl",
					    "trashCompactor", kernelInfo);


  
}


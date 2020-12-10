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


class platform_t;

#define DSORT LIBP_DIR"/libs/core/"

#ifndef DEVICE_SCAN
#define DEVICE_SCAN 1

// was 512
//#define SCAN_BLOCK_SIZE 256


//template <class T>
class deviceScan_t{
private:
  occa::kernel blockShflScanKernel;
  occa::kernel finalizeScanKernel;
  occa::kernel findStartsKernel;
  occa::kernel segmentedReductionKernel;
public:

  // constructor
  deviceScan_t(platform_t &platform, const char *entryType, const char *entryMap, occa::properties props);

  dlong blockCount(dlong entries);

  void scan(const dlong    entries,// might need to upgrade
	    occa::memory &o_list,
	    occa::memory &o_tmp,
	    dlong *h_tmp,
	    occa::memory &o_scan);

  void mallocTemps(platform_t &platform, dlong entries, occa::memory &o_tmp, dlong **h_tmp);

  dlong segmentedReduction(platform_t &platform,
			   const dlong entries,
			   const int entrySize,
			   const int includeLast,
			   occa::memory &o_list,
			   occa::memory &o_compactedList);

  
};
  
#endif


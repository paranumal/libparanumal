#define DSORT LIBP_DIR"/libs/core/"

class deviceSort_t{
 private:
  occa::kernel bitonicSortSharedKernel;
  occa::kernel bitonicSwapGlobalKernel;
  occa::kernel bitonicMergeGlobalKernel;
  occa::kernel bitonicMergeSharedKernel;
 public:
  // constructor
  deviceSort_t(occa::device &device, const char *entryType, const char *entryHeader, occa::properties props);

  // sort
  void sort(const int entries, occa::memory o_list);
};

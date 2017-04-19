

dcsr *newDCSR(almond_t *almond, csr *B);

void dcsrHaloExchangeStart(dcsr *A, size_t Nbytes, void *sendBuffer, void *recvBuffer);

void dcsrHaloExchangeFinish(dcsr *A);



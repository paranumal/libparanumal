
int compareEntry(const entry_t *a, const entry_t *b){

  if(a->row < b->row) return -1;
  if(a->row > b->row) return +1;

  if(a->col < b->col) return -1;
  if(a->col > b->col) return +1;

  return 0;
}


int entryMap(const entry_t a, const entry_t b){

  if(a.row != b.row || a.col != b.col) return 1;
  
  return 0;
}

template <class T>
int setupAide::getArgs(string key, T& t){
  matrix<T> m;

  getArgs(key,m);

  if(m.size()){
    t = m[1];

    return 1;
  }

  printf("Failed to find [%s].\n", key.c_str());
  return 0;
}

template <class T>
int setupAide::getArgs(string key, matrix<T>& m){
  stringstream args;
  vector<T> argv;
  int argc;
  T input;

  args.str( getArgs(key) );

  while(args >> input)
    argv.push_back(input);

  argc = argv.size();

  if(!argc){
    printf("Failed to find [%s].\n", key.c_str());
    return 0;
  }

  m.resize(argc,1);

  for(int i=1; i<=argc; i++)
    m[i] = argv[i-1];

  return 1;
}

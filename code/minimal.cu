# include "uammd.cuh" /* UAMMD */

using namespace uammd; /* UAMMD */
using std::make_shared;
using std::endl;
using std::cout;

int main(int argc, char *argv[]){

  auto sys = make_shared<System>(argc,argv); /* UAMMD */

  cout<<endl<<"--> Hello, UAMMD! <--"<<endl<<endl;

  sys->finish(); /* UAMMD */

  return 0;
}

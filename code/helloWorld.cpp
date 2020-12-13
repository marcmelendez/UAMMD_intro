# include <iostream>

using std::cout;
using std::endl;

int add(int a, int b) {
  return a + b;
}

void goodbye() {
  cout<<"Goodbye!"<<endl;
  return;
}
int main(int argc, char * argv[])
{
  cout<<"Hello, World!"<<endl;
  int printEverynSteps = 13;

  for(int i = 1; i <= 100; ++i) {
    if(printEverynSteps > 0
       && i % printEverynSteps == 0) {
        cout<<"Iteration: "<<i<<endl;
    }
  }
  cout<<"3 + 5 = "<<add(3,5)<<endl;
  goodbye();
	return 0;
}

#include <iomanip> 
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <cstdio>
#include <new>


using namespace std;

int main( int argc, char *argv[]) 
{
  
 int rsize, psize;
 
 sscanf ( argv[1], "%i", &rsize);
 sscanf ( argv[2], "%i", &psize);
 
  double s[200][200];
  
  int j1, j2; 
  
 ifstream coeff("./tmp/mat.tmp"); 
 for(j1=1;j1<=50;j1++){for(j2=1;j2<=50;j2++){coeff>>s[j1][j2];}}
 coeff.close();
 
int k1, k2, i1, i2; 
double K1, K2, I1, I2; 
double f1, f2, c;  

for(k1=0;k1<rsize;k1++){
for(k2=0;k2<psize;k2++){

c=0;
K1=double(k1);
K2=double(k2);

for(i1=0;i1<50;i1++){
for(i2=0;i2<50;i2++){
I1=double(i1);
I2=double(i2);
f1=cos( (3.14/3000)*(K1+1/2)*(I1+1/2));
f2=cos( (3.14/750 )*(K2+1/2)*(I2+1/2));
 c=c+s[i1+1][i2+1]*f1*f2; 
 
} 
                    }
cout<<c<<" ";
                    }
cout<<" "<<endl;
                    }


}

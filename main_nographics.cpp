/* how to compile under Windows 64bit:

make -f Makefile_nographics.win
(this optimizes for this CPU and 64bit)

this creates no_graphics.exe which can be executed from the command line with parameters as usual

*/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
using namespace std;
#define SQR(x)  (x)*(x)
#include "functions.h"
#include "rod_shaped.h"

extern vector <Rod_shaped*> bact ;
int sample ;

void err(char *reason)
{
  cout <<reason<<endl ; 
  exit(0) ;
}

void err(char *reason, int a)
{
  cout <<reason<<": "<<a<<endl ; 
  exit(0) ;
}

void err(char *reason, double a)
{
  cout <<reason<<": "<<a<<endl ; 
  exit(0) ;
}

extern long long int xdr48 ;
extern long long tinit ; // real time at start

int main(int argc, char *argv[])
{
  NUM=new char[256] ; strcpy(NUM,"test_dir") ;
  int nsam=1,stopcond=0 ;
  start_file=NULL ;
  float max_time ;
  int max_size=0, stc ;
  if (argc<3) { err(" Error:: at least two arguments needed, option and value, for example:\n"
                    "n [directory name]\n"
                    "r [RNG seed]\n"
                    "N [no. of samples]\n"
                    "S [stop when: 1=only one type in the population]\n"
                    "\t[4 = after reaching size (next argument)]\n"
					"\t[8 = after given time (next arg)]. Conditions can be combined e.g. S 1 S 8 1000]\n"                    
                    "i [start filename]\n"
                    "E [elastic modulus]\n"
                    "w or h [width or height]\n"
                    "z [zeta friction coefficient]\n"
                    "Z [zeta_motile wrt to growing]\n"
                    "d [death_rate]\n"
                    "v [velocity_motile]\n"
                    "f [ini_fraction_motile]\n"
                    "Program terminated. \n"); } 
  else { 
    for (int i=1;i<argc;i++) {
      switch (argv[i][0]) {
        case 'E': Rod_shaped::E_bact=atof(argv[++i]) ; cout <<"E_bact="<<Rod_shaped::E_bact<<" " ; break ;
        case 'i': start_file=argv[++i] ; cout <<"startfile="<<start_file<<" " ;break ;
        case 'h': height=atoi(argv[++i]) ; cout <<"height="<<height<<" " ; break ;
        case 'n': NUM=argv[++i] ; cout <<"NUM="<<NUM<<" " ; break ;
        case 'N': nsam=atoi(argv[++i]) ; cout <<"nsam="<<nsam<<" " ; break ;
        case 'r': RAND=atoi(argv[++i]) ; cout <<"RAND="<<RAND<<" " ; break ;
        case 'S': stc=atoi(argv[++i]) ; cout <<"stopcond="<<stc<<" " ; 
          if (stc==4) max_size=atoi(argv[++i]) ;
          if (stc==8) max_time=atof(argv[++i]) ;
          stopcond+=stc ; cout<<"total stopcond="<<stopcond<<" " ;
          break ;
        case 'w': width=atoi(argv[++i]) ; cout <<"width="<<width<<" " ; break ;
        case 'z': Rod_shaped::zeta=atof(argv[++i]) ; cout <<"zeta="<<Rod_shaped::zeta<<" " ; break ;
        case 'Z': zeta_motile_wrt_growing=atof(argv[++i]) ; cout <<"zeta_motile_wrt_growing="<<zeta_motile_wrt_growing<<" " ; break ;
        case 'd': death_rate=atof(argv[++i]) ; cout <<"death_rate="<<death_rate<<" " ; break ;
        case 'v': velocity_motile=atof(argv[++i]) ; cout <<"velocity_motile="<<velocity_motile<<" " ; break ;
        case 'f': ini_fraction_motile=atof(argv[++i]) ; cout <<"ini_fraction_motile="<<ini_fraction_motile<<" " ; break ;
        default: err("unrecognized option.") ; break ;
      }
    }
    cout <<endl ;
  }
  _srand48(RAND) ;
  init();
  
  for (sample=0;sample<nsam;sample++) {
    if (sample>0) { // reset if necessary for more than 1 sample
      reset(false) ; 
    } 
    for (long long int i=0;;i++) {
      run() ;
      if (i%1000==0) {
        cout <<1.*(clock()-tinit)/CLOCKS_PER_SEC<<"\ts"<<sample<<": N="<<bact.size()<<" T="<<tt<<endl ; 
      }
      if ((stopcond&1) && is_only_one_type()==true) {
      	cout <<"only_one_type reached\n" ;
        break ;
      }
      if ((stopcond&4) && bact.size()>=max_size) {
      	cout <<"max_size reached\n";
        break ;
      }
      if ((stopcond&8) && tt>=max_time) {
      	cout <<"max_time reached\n";
        break ;
      }
    }
  }
	return 0 ;
}

void Rod_shaped::draw()
{    
}



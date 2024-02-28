#define _USE_MATH_DEFINES //AM
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vector>
using namespace std;
#include "functions.h"
#include "vecs.h"
#include "rod_shaped.h"
#include <time.h>

#define __MAIN
#include "params.h"


// -----------variables used in the main loop--------------------
// cc[][] is the concentration of food
// tt is time

// BW - diffusion is deleted in this version of the code

char *start_file ;
double dt ;
double d0=2*r0 ; // width of the cell
char *NUM ;
int RAND=1 ;
int gr_ok=1 ; // if 1 then growth is simulated
double tt ;

int ttt;  // number of simulation steps so far
int ttt_last_save, ttt_last_dt_update ; 
int ttt_last_save_2;//AM
int how_many_steps_save ;
int how_many_steps_save_2;//AM 

//AM vector <Bacterium*> bact ; // vector of cells
vector <Rod_shaped*> bact ; //AM relying on what BW did in one of the visual verions of the code

// these we don't save:

long long tinit ; // real time at start
int noncriterr ; // if set, it means that a non-critical error has occured and the last saved state has to be loaded 

int noclosest ;
double avclosest ; // average number of closest neightbours
int nbox, boxwidth=int(10*d0) ;
vector <vector <vector <int> > > boxes ; // table of boxes (vectors) holding the cells

vector <int> number_past;
vector <vec2> locations_past;
int frame=0;
int division_register=0;
double initial_d_s=0;
int n_to_kill=0;

ofstream secs, params, locs ;

#ifdef __linux
#include <unistd.h>
typedef unsigned int DWORD ;
int memory_taken() // return memory available in MB
{
  long long int rss = 0L;
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
		return (size_t)0L;		/* Can't open? */
	if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
	{
		fclose( fp );
		return (size_t)0L;		/* Can't read? */
	}
	fclose( fp );
	long long int mem=((size_t)rss * (size_t)sysconf( _SC_PAGESIZE)) ;
	return (int) (mem/(1<<20));
}
#include <sys/sysinfo.h>
unsigned int freemem() // returns available memory in MB
{
  struct sysinfo s ;
  sysinfo(&s) ;
  return ((s.freeram)>>20) ;
}
#else
#include <windows.h>
#include <psapi.h>
int memory_taken() // return memory available in MB
{
	//PROCESS_MEMORY_COUNTERS info;
	//GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	//return (int) (info.WorkingSetSize/(1<<20));
  return 0;
}
#endif

// random number generator - like drand48() under Linux, but faster
// works only with compilers which have long long int!
long long unsigned int xdr48=0x000100010001LL, mndr48=0x0005deece66dLL, doddr48=0xbLL ;
double _drand48(void)  
{
  xdr48=mndr48*xdr48+doddr48 ; xdr48&=0xffffffffffffLL ;
  return (xdr48/281474976710656.0) ;
}

void _srand48(unsigned int x0)
{
	xdr48=(x0<<16)+0x330e ; 
}

double gauss48(){

  static int i=0;
  double phi,r;
  static double x,y;

  if(i==0){
    phi = 2.0*M_PI*_drand48();
    r = sqrt(-2.0*log(_drand48()));
    x = cos(phi)*r;
    y = sin(phi)*r;
    i=1;
    return x;
  }
  else{
    i=0;
    return y;
  }
} 

double sign(double x) { if (x<0) return -1 ; else return 1 ; }

inline int xtobox(double x) { return int(nbox/2+x/boxwidth) ; }
inline int ytobox(double y) { return int(nbox/2+y/boxwidth) ; }

void   populate_boxes()
{
  int i,j;
  for (i=0;i<nbox;i++) for (j=0;j<nbox;j++) boxes[i][j].clear() ;
  for (i=0;i<bact.size();i++) {
    bact[i]->bi=boxes[xtobox(bact[i]->x())][ytobox(bact[i]->y())].size() ;
    boxes[xtobox(bact[i]->x())][ytobox(bact[i]->y())].push_back(i) ;
    bact[i]->bx=xtobox(bact[i]->x()) ;
    bact[i]->by=ytobox(bact[i]->y()) ;
    if (bact[i]->bx<0 || bact[i]->by<0 || bact[i]->bx>=nbox || bact[i]->by>=nbox) err("bxby outside nobx") ;
  }
}

void update_box(Rod_shaped *b, int i)
{
  while (b->x()>width/2) b->add_x(-width) ; //p.b.c
  while (b->x()<-width/2) b->add_x(width) ;
  while (b->y()>height/2) b->add_y(-height) ;//p.b.c.
  while (b->y()<-height/2) b->add_y(height) ;
  int bx=xtobox(b->x()) ;
  int by=ytobox(b->y()) ;
  if (bx<0 || by<0 || bx>=nbox || by>=nbox) err("bxby outside nobx") ;
  if (b->bi==-1) { // freshly made bacterium
    b->bx=bx ; b->by=by ;
    b->bi=boxes[bx][by].size() ;
    boxes[bx][by].push_back(i) ;
  } else if (b->bx!=bx || b->by!=by) {
    if (boxes[b->bx][b->by].size()==0) err("zero",b->by) ;
    int j=boxes[b->bx][b->by].back() ;  // index of bacterium at the end of the box
    if (j<0 || j>=bact.size()) err("j=",j) ;
    boxes[b->bx][b->by].pop_back() ;
    if (j!=i) { boxes[b->bx][b->by][b->bi]=j ; bact[j]->bi=b->bi ; }
    
    b->bi=boxes[bx][by].size() ; 
    b->bx=bx ; b->by=by ;
    boxes[bx][by].push_back(i) ;
  } 
}

void just_mutated(int i) 
{
  float yp=bact[i]->y(),ym=yp ;
  for (int j=0;j<bact.size();j++) {
    if (i!=j && fabs(bact[i]->x()-bact[j]->x())<2*d0 && bact[j]->y()>yp) yp=bact[j]->y() ;
    if (i!=j && fabs(bact[i]->x()-bact[j]->x())<2*d0 && bact[j]->y()<ym) ym=bact[j]->y() ;
  }
}

int Rod_shaped::tot_number_of_bacteria = 0 ; // initialize tot. bacterium counter

void reset(bool reset_rand)
{
  dt=dtmax ; //set up time step
  int i;
  char txt[256] ;
  for (i=0;i<bact.size();i++) delete bact[i] ; //delete all bacteria
  bact.clear() ;

  float xc=-width/2+d0,yc=-height/2+d0 ;
  do {
    Rod_shaped *b=new Rod_shaped(0) ; //i%2) ; 
    b->r.x=xc ; b->r.y=yc ; 
    xc+=1.1*d0 ; if (xc>width/2) { xc-=width ; yc+=(2*Rod_shaped::motile_length+2.3)*d0 ; }
    b->phi=M_PI/2+(_drand48()-0.5)*0.2 ; 
    b->recalc_cspt() ;
    b->type=0 ; //if (yc>-height/3) b->type=1 ;
    if (_drand48()<ini_fraction_motile)
    {
   	 b->type=1 ;
         b->d_final=b->d=Rod_shaped::motile_length ;
    }
    bact.push_back(b) ;
  } while (yc<height/2-d0) ;
  populate_boxes() ;
  bact[0]->init() ; // initialize some tables inside "Rod_shaped" and find closest n.n.

  sprintf(txt,"%s/par(t)_%d_%d.dat",NUM,RAND,sample) ; params.open(txt) ; params.close() ;
  sprintf(txt, "%s/loc_and_v(t)_%d_%d.dat", NUM, RAND, sample); locs.open(txt) ; locs << width << "\t" << height << endl; locs.close() ;

  long long unsigned int temp ;   
  if (!reset_rand) temp=xdr48 ; 
  if (start_file!=NULL) { printf("loading file %s",start_file) ; fflush(stdout) ; load_all(start_file) ; printf(" done\n") ; fflush(stdout) ; }
  if (frame < 3) 
  {
      for (int i = 0; i < bact.size(); i++)
      {
          if (bact[i]->type == 1)
          {
              number_past.push_back(bact[i]->number);
              locations_past.push_back(bact[i]->r);
          }
      }
  }

  if (!reset_rand) xdr48=temp ; // make sure that after resetting and reloading init. conf. we do not restart RNG
  ttt_last_save=-(1<<30) ; ttt_last_dt_update=0 ; ttt=0 ; tt=0 ;
}

// initialize everything: allocate memory, set initial concentrations of food
// throw some cells at random
void init()
{
  if (int(width/boxwidth)*boxwidth!=width) err("width must be the multiple of 10*d0") ;  //chceck if 
  dt=dtmax ; //set up initial time step
  nbox=2+MAX(int(width/boxwidth),int(height/boxwidth)) ; // number of boxes along one side of the simulation board
  nbox=int(nbox/2)*2 ; // to make it even
  boxes.resize(nbox) ; //make number of boxes along one dimension of the board equal to nbox
  for (int i=0;i<nbox;i++) boxes[i].resize(nbox) ; //make number of boxes along the other dimension of the board equal to nbox

  char txt[256];
  sprintf(txt, "mkdir %s", NUM);
  system(txt); //make a folder named using NUM variable

  sprintf(txt, "del /Q %s", NUM);
  system(txt); //delete the contents of that folder; makes sense if the folder was created during a previous run of the program

  reset(false) ;
  tinit=clock() ; noncriterr=0 ;
}

#define FW(x) fwrite(&(x),sizeof((x)),1,all) ;
#define FR(x) fread(&(x),sizeof((x)),1,all) ;

void save_all(char *name)
{
  FILE *all;
  int i;
  all=fopen(name,"wb") ;
  i=bact.size() ;
  FW(i) ; //FW(tabx) ; FW(taby) ; FW(heightone) ; 
	FW(width) ; FW(height) ;
  FW(tt) ; FW(xdr48) ;

  for (i=0;i<bact.size();i++) bact[i]->save_data(all) ;
  fclose(all) ;
}

void load_all(char *name)
{
  int i,nn;
  int tabxn, tabyn, heightonen ;
  long long int newxdr48 ;
  FILE *all;
  all=fopen(name,"rb") ;
  FR(nn) ; //FR(tabxn) ; FR(tabyn) ;FR(heightonen) ; 
	FR(width) ; FR(height) ;
  FR(tt) ; FR(newxdr48) ;
    for (i=0;i<bact.size();i++) delete bact[i] ;
  bact.clear() ;
  for (i=0;i<nn;i++) { 
    bact.push_back(new Rod_shaped) ; bact[i]->load_data(all) ;
  }
     
  fclose(all) ; 
  populate_boxes() ;
  xdr48=newxdr48 ;
  _maxdisp=0 ; 
}

void save_data(char *name, int cc_all=0)
{
  FILE *pos=fopen(name,"w") ;
  if (pos==NULL) err(name) ;
  fprintf(pos,"%d\t%d %d\n",bact.size(),width,height) ;

  for (int i=0;i<bact.size();i++) {
    char txt[256] ;
    bact[i]->to_string(txt) ;
    fprintf(pos,txt) ;
  }
  fclose(pos) ;
}

double globalOrder() { //output the global order parameter

     double Q11=0,Q12=0;

     for (int i=0;i<bact.size();i++) {
         Q11+=SQR(bact[i]->cos_phi())-0.5;
         Q12+=bact[i]->cos_phi() * bact[i]->sin_phi();
     }

     Q11 /= bact.size(); Q12 /= bact.size();
     return 2*sqrt(SQR(Q11) + SQR(Q12));
}

int count_types(int t)
{
  int n=0 ;
  for (int i=0;i<bact.size();i++) if (bact[i]->type==t) n++ ;
  return n ;
}

bool is_only_one_type()
{
  for (int i=1;i<bact.size();i++) if (bact[i]->type!=bact[0]->type) return false ;
  return true ; 
}

void save_how_many_bact_each_type(char *name)
{
  int maxt=0 ;
  for (int i=0;i<bact.size();i++) if (bact[i]->type+1>maxt) maxt=bact[i]->type+1 ;
  int *types=new int[maxt], *act_types=new int[maxt] ;
  for (int i=0;i<maxt;i++) types[i]=act_types[i]=0 ;
  for (int i=0;i<bact.size();i++) {
    types[bact[i]->type]++ ;
    if (bact[i]->active==1) act_types[bact[i]->type]++ ;
  }

  FILE *f ;
  if (ttt==0) f=fopen(name,"w") ; else f=fopen(name,"a") ;
  if (f==NULL) err("cannot open file for save_how_many_bact_each_type") ;
  for (int i=0;i<maxt;i++) if (types[i]>0) fprintf(f,"%d %d %d %d\n",ttt,i,types[i],act_types[i]) ; // time_step  type  #bacteria
  fclose(f) ;
  delete [] types ; delete [] act_types ;
}

int active_cells()
{
  int nn0=0 ;
  for (int i=0;i<bact.size();i++) if (bact[i]->active==1) nn0++ ;
  return nn0 ;
}  

int growing_cells()
{
  float max_gr_rate=0 ;
  int nn=0 ;
  for (int i=0;i<bact.size();i++) if (bact[i]->elong_rate()>max_gr_rate) max_gr_rate=bact[i]->elong_rate() ;
  for (int i=0;i<bact.size();i++) if (bact[i]->elong_rate()>0.5*max_gr_rate) nn++ ;
  return nn ;
}

double tot_area()
{
	double a=0;
  for (int i=0;i<bact.size();i++) a+=bact[i]->area() ;
  return a ;
}

void end()
{
// if any files open, close them
}

int count_type(int tp, int i)
{
	int count=0;
  for (int j=0;j<bact[i]->closest.size();j++)
  	if (bact[bact[i]->closest[j]]->type==tp) count++ ;
  return count ;
}

//------------cluster identification--------------------------

vector <int> vis;
vector <vector <int> > clusters ;
int maxclsize ;
double orderinmaxcl,surfmaxcl, vxmaxcl,vymaxcl ;

void count_clusters() //new code from Bartek
{
  vector <int> heap;
  int i, j, k,v,n;
	int no_cluster=0 ;

	vis.clear() ; vis.resize(bact.size(),-1) ;
  clusters.clear() ;

  for (k=0;k<bact.size();k++) {
    if (vis[k]==-1 && bact[k]->type==1) {
      heap.push_back(k) ; n=1 ;
      vis[k]=no_cluster ; 
			vector <int> nv ;
			nv.push_back(k) ;
			clusters.push_back(nv) ;
      do {
        n-- ;
				v=heap[n] ; 
        heap.pop_back();
        for (j=0;j<bact[v]->closest.size();j++) {
        	int v2=bact[v]->closest[j] ;
          if (bact[v2]->type==1 && vis[v2]==-1 && Rod_shaped::find_distance_ij(v, v2)<2*r0) { heap.push_back(v2) ; vis[v2]=no_cluster ; clusters[no_cluster].push_back(v2) ; n++;}
				}
      } while (n>0) ;
      heap.clear();
      no_cluster++ ;
    }
  }
  //printf("%d\n",int(clusters.size())) ;
}
//		max_cluster_size(&maxclsize,&orderinmaxcl,&velinmaxcl,&surfmaxcl) ;

void max_cluster_size(int *maxsize,double *order,double *vx, double *vy,double *surf) 
{
	count_clusters() ;
	//return clusters.size() ; // gives the number of clusters
	int which=-1 ;
	(*maxsize)=0 ;
	for (int i=0;i<clusters.size();i++) if (clusters[i].size()>(*maxsize)) { which=i ; (*maxsize)=clusters[i].size() ; }
	// marking cells in the max cluster

  double Q11=0,Q12=0;
  (*vx)=(*vy)=0 ;
  for (int j=0;j<clusters[which].size();j++) {
  	int i=clusters[which][j] ;
    Q11+=SQR(bact[i]->cos_phi())-0.5;
    Q12+=bact[i]->cos_phi() * bact[i]->sin_phi();
    (*vx)+=bact[i]->vx() ; (*vy)+=bact[i]->vy() ;
  }

  Q11 /= clusters[which].size(); Q12 /= clusters[which].size();
  (*order)=2*sqrt(SQR(Q11) + SQR(Q12));

	(*vx)/= clusters[which].size(); (*vy)/= clusters[which].size();

}


// main loop - does a single step of time evolution
void run()
{
  if (_stop || bact.size()==0) return ;

  int i, j;

  how_many_steps_save=int(0.5/dt) ; 
  how_many_steps_save_2=int(1./dt); //AM  
  _maxdisp=0 ; noncriterr=0 ;
  noclosest=0 ;
  
  for (i=0;i<bact.size();i++) {
    bact[i]->recalc_mass_and_I() ;
    bact[i]->clear_forces() ;
    bact[i]->n_closest_of_closest=0;//20220531 
  }

  for (i=0;i<bact.size();i++) bact[i]->find_acceleration(i) ;  
  for (i=0;i<bact.size();i++) if (!noncriterr && bact[i]->active>=0) bact[i]->one_Euler_step(i) ; 
  
  if (!noncriterr) { // if everything is ok then....
    int nn0=bact.size() ;
    for (i=0;i<nn0;i++) {
			if (bact[i]->type==0 && gr_ok) bact[i]->one_growth_step(i) ; 
    }
    for (i=0;i<bact.size();i++) update_box(bact[i],i) ;

    avclosest=0 ;
    for (i=0;i<bact.size();i++) {
      if (bact[i]->nn_needs_update && bact[i]->active>=0) bact[i]->find_closest(i) ;
      avclosest+=bact[i]->closest.size() ;
    }
    avclosest/=bact.size() ;
  } else if ((save_cells&1) || backup) { // but if not, load the last configuration and decrease dt
    char txt[256] ; 
    if (backup) sprintf(txt,"%s/temp_%d.dat",NUM,RAND) ; 
    if ((save_cells&1)) sprintf(txt,"%s/temp_%d_%d.dat",NUM,RAND,sample) ; 
    printf("_maxdisp=%lf, dt decreased at t=%lf to %lf, reloading configuration %s...",_maxdisp,tt,dt/2,txt) ; fflush(stdout) ;
    load_all(txt) ; printf(" done\n") ; fflush(stdout) ;
    dt/=2 ; noncriterr=0 ;
    if (dt<1e-6) err("dt too small!") ;
    // BW data files should be modified accordingly!
    return ;    
  } else {//err("dt too large. Cannot reload, config. non-existent") ;
    // reinitialize if nothing else can be done, with smaller dt
    printf("cannot decrease dt because temp.dat does not exists. Restarting...\n") ; fflush(stdout) ;
    reset(true) ;
    dt/=2 ; noncriterr=0 ;
    return ;
  }
  
  if (ttt%10==0) {
    int some_change ;
    some_change=0 ;
  // change state of cells

		if (tot_area()>width*height*confluence) {
			//err("n",n_to_kill) ;
			for (;n_to_kill>0;n_to_kill--) {
				do { i=int(bact.size()*_drand48()) ; } while (bact[i]->type!=0 || bact[i]->active==-1) ;
				some_change=1 ; bact[i]->active=-1 ;
			}
		} else n_to_kill=0 ;

    ////sort(ints.begin(), ints.end());
    for (i=0;i<bact.size();i++) {
      if (bact[i]->active==-1) { // delete bacterium 
        // first save its position etc. into the file.......
        // .......and then delete it:

        j=boxes[bact[i]->bx][bact[i]->by].back() ; 
        if (j<0 || j>bact.size()) err("er. j=",j) ;
        boxes[bact[i]->bx][bact[i]->by].pop_back() ;
        if (j!=i) { boxes[bact[i]->bx][bact[i]->by][bact[i]->bi]=j ; bact[j]->bi=bact[i]->bi ; }
  
        delete bact[i] ;
        if (i<bact.size()-1) {
          bact[i]=bact.back() ; 
          boxes[bact[i]->bx][bact[i]->by][bact[i]->bi]=i ;           
        }
        bact.pop_back() ;
        i-- ; some_change=1 ;
      }
      if (i>=bact.size()-1) break ;
    }
// BW this is necessary only if bacteria are removed, but can slow down the program by 50% (?)
    if (some_change) for (i=0;i<bact.size();i++) if (bact[i]->active>=0) bact[i]->find_closest(i) ;
    }
  
  char txt[256];
if (ttt>ttt_last_save_2+how_many_steps_save_2){
ttt_last_save_2=ttt;
//if (frame>1000){
//sprintf(txt,"%s/tempx_%d_%d.dat",NUM,ttt,sample) ; save_all(txt) ;
//}

frame++;
}
  if (ttt>ttt_last_save+how_many_steps_save) {
    ttt_last_save=ttt ; 
        double m_ncoc=0;
        double m_d=0;
        double m_d_final=0;
    for (i=0; i<bact.size();i++)
    {
      if (bact[i]->type==0)
      {m_ncoc=m_ncoc+bact[i]->n_closest_of_closest; m_d=m_d+bact[i]->d;m_d_final=m_d_final+bact[i]->d_final;}
    }
    m_ncoc=m_ncoc/count_types(0);
    m_d=m_d/count_types(0);
    m_d_final=m_d_final/count_types(0);
    initial_d_s=initial_d_s/division_register;

    if ((save_cells&1)) { sprintf(txt,"%s/temp_%d_%d.dat",NUM,RAND,sample) ; save_all(txt) ; load_all(txt) ; }
    else if (backup) { sprintf(txt,"%s/temp_%d.dat",NUM,RAND) ; save_all(txt) ; load_all(txt) ; }
//    load_all(txt) ;
    if ((save_cells&2)) { sprintf(txt,"%s/data_%d_%d.dat",NUM,RAND,sample) ; save_data(txt) ; }

    if ((save_cells&16)) { sprintf(txt,"%s/types_%d_%d.dat",NUM,RAND,sample) ; save_how_many_bact_each_type(txt) ; }

    sprintf(txt,"%s/par(t)_%d_%d.dat",NUM,RAND,sample) ; params.open(txt,ios::app) ;

		max_cluster_size(&maxclsize,&orderinmaxcl,&vxmaxcl,&vymaxcl,&surfmaxcl) ;
    // 1.time   2.iterations  3.size  4.no.of_type_0  5. max_cluster_size
    params<<tt<<"\t"<<ttt<<"\t"<<bact.size()<<"\t"<<count_types(0)<<"\t"<<maxclsize ;
		// 5.order_in_max_cluster  6. global_order 7. surface_max_cluster  8. vx_in_max_cluster  9.vy_in_max_cluster
    params<<"\t"<<orderinmaxcl<<"\t"<<globalOrder()<<"\t"<<surfmaxcl<<"\t"<<vxmaxcl<<"\t"<<vymaxcl ;

    //  time[s]  memory_taken[MB]
    params<<"\t"<<(clock()-tinit)/CLOCKS_PER_SEC<<"\t"<<memory_taken()<<"\t"<<m_ncoc<<"\t"<<division_register<<"\t"<<m_d<<"\t"<<initial_d_s<<"\t"<<m_d_final<<"\t"<<tot_area()/6400.<<endl ;
    params.close() ;
    if (params.fail()) params.clear() ;
    //sprintf(txt, "%s/cluster_sizes(t)_%d_%d.dat", NUM, RAND, sample); params.open(txt, ios::app);
    //params << tt << "\t" << ttt;
    //int sizee;
    //int i_size;
    //i_size = int(clusters.size());
    //for (int i = 0; i < i_size; i++) {
    //    //sizee = clusters[i].size();
    //            params << "\t" << clusters[i].size();
    //                }
    //params << endl;
    //params.close();
    //if (params.fail()) params.clear();
    division_register=0;
    initial_d_s=0;
    //sprintf(txt, "%s/loc_and_v(t)_%d_%d.dat", NUM, RAND, sample); locs.open(txt, ios::app);
    //locs << tt << "\t" << ttt<<"\t"<<bact.size();
    //for (int i = 0; i < bact.size(); i++) {
    //  locs <<"\t"<< bact[i]->number <<"\t"<< bact[i]->type << "\t" << bact[i]->r.x << "\t" << bact[i]->r.y <<"\t" << bact[i]->v.x << "\t" << bact[i]->v.y << "\t" << bact[i]->n.x << "\t" << bact[i]->n.y << "\t" << bact[i]->d ;
    //}
    //locs<<endl<<flush ;
    //locs.close() ;
    //if (locs.fail()) locs.clear() ;

 }



  ttt++ ;
  tt+=dt ;

}

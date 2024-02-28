#define _USE_MATH_DEFINES //AM
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;
#include "functions.h"
#include "vecs.h"
#include "rod_shaped.h"
#include "params.h"

extern vector <Rod_shaped*> bact ;
extern vector <vector <vector <int> > > boxes ;
extern int division_register;
extern double initial_d_s;
extern int n_to_kill;
//inline int xtotab(double x) { return int(tabx/2+x/dx0) ; }
//inline int ytotab(double y) { return int(taby/2+y/dx0) ; }

int Rod_shaped::powmax15=1000, Rod_shaped::powmax05=1000 ;
int Rod_shaped::rods=0 ;
double *Rod_shaped::pow15=new double[Rod_shaped::powmax15+1], *Rod_shaped::pow05=new double[Rod_shaped::powmax05+1] ;
double max_max_length=MAX(Rod_shaped::max_length+0.5, Rod_shaped::motile_length);

typedef vec2 vec ;

void Rod_shaped::init()
{
  int i;
  for (i=0;i<=powmax15;i++) pow15[i]=(4./3)*E_bact*sqrt(r0*d0)*d0*pow(1.*i/powmax15,1.5) ; ; // force give by pow15 is in Pa um^2 = pN 
  for (i=0;i<=powmax05;i++) pow05[i]=pow(1.*i/powmax05,0.5) ;  
  for (int i=0;i<bact.size();i++) bact[i]->find_closest(i) ; 
}

/*bool Rod_shaped::inside_bact(vec p) // returns true if point p lies inside the bacterium
{
  double dd,t ;
  find_sqr_distance_to_rod(p,dd,t) ;
  if (dd<r0*r0) return true ;
  return false ;
}*/

void per_diffxy(vec &r1, vec &r2, vec &dr) // calculates the different between r1 and r2 assuming pbc in x and y
{
  dr=r1-r2 ; 
  if (dr.x<-width/2) dr.x+=width ; 
  else if (dr.x>width/2) dr.x-=width ; 
  if (dr.y<-height/2) dr.y+=height; 
  else if (dr.y>height/2) dr.y-=height ; 
}

const int kx[9]={0,1,1,0,-1,-1,-1,0,1},ky[9]={0,0,1,1,1,0,-1,-1,-1} ;  
inline int xtobox(double x) { return int(nbox/2+x/boxwidth) ; }
inline int ytobox(double y) { return int(nbox/2+y/boxwidth) ; }

void pbc_location(vec& r)
{
    while (r.x > width / 2) r.x -= width;
    while (r.x < -width / 2) r.x += width;
    while (r.y > height / 2) r.y -= height;
    while (r.y < -height / 2) r.y += height;
}

void Rod_shaped::find_closest(int thisi) 
{
  noclosest++ ;
  vec dr ;
  rold1=r+n*d ; rold2=r-n*d ; dold=d ;
  closest.clear() ;
  pbc_location(this->r);
  int bxn, byn, *box ;  
  if (bi==-1 || bx==-1 || by==-1) { bxn=xtobox(r.x) ; byn=ytobox(r.y) ; } else { bxn=bx ; byn=by ; }
  for (int k=0;k<=8;k++) {
    int xx=int(width/boxwidth) ; //err("x",xx) ;
    int yy=int(height/boxwidth) ;
    box=&(boxes[int(nbox/2-xx/2+(bxn+kx[k]-nbox/2+3*xx/2)%xx)][int(nbox/2-yy/2+(byn+ky[k]-nbox/2+3*yy/2)%yy)][0]) ;
    int nmax=boxes[int(nbox/2-xx/2+(bxn+kx[k]-nbox/2+3*xx/2)%xx)][int(nbox/2-yy/2+(byn+ky[k]-nbox/2+3*yy/2)%yy)].size() ;
    for (int n=0;n<nmax;n++) {
      int i=box[n] ;
      per_diffxy(this->r, bact[i]->r, dr) ;
      double r2=squared(dr) ;
      if (i!=thisi && r2<(SQR(2*max_max_length*d0+d0))) {
        closest.push_back(i) ; //err("x",i) ;
        int j ;
        for (j=0;j<bact[i]->closest.size();j++) if (bact[i]->closest[j]==thisi) break ;
        if (j==bact[i]->closest.size()) bact[i]->closest.push_back(thisi) ;
      }
    }  
  }
  nn_needs_update=0 ;
}

inline void Rod_shaped::find_sqr_distance_to_rod(vec2 &rr, double &dist2, double &t)
{
  double drx=rr.x-r.x, dry=rr.y-r.y ;
  t=(drx*n.x+dry*n.y) ;   
  if (t>d) t=d ;
  if (t<-d) t=-d ;  
  dist2=(SQR(drx-t*n.x)+SQR(dry-t*n.y)) ;
}

void find_forces_SS(vec &r1, vec &r2,  vec &f1)
{
  const double d20=d0*d0 ;
  vec n=r2 ; n-=r1 ;
  double n2=squared(n) ;
  if (n2>d20) { f1.zero() ; return ; }
  double ln=sqrt(n2) ;
  n/=ln ;
  double fn=Rod_shaped::pow15[int(Rod_shaped::powmax15*(d0-ln)/d0)] ; //pow(d0-n,1.5) ;
  f1=n*(-fn) ;
}

void Rod_shaped::update_acceleration_i()
{
  vec f1,f2, r1=n*d,r2=n*(-d) ;
  f+=f1 ; f+=f2 ; 
  torque+=cross(r1,f1) ;
  torque+=cross(r2,f2) ;
}

void find_sqr_distance_lines(Rod_shaped *bi,Rod_shaped *bj,double &dij,double &ti, double &tj) 
{
  vec rij= bi->r, ni=bi->n, nj=bj->n, r1,r2,f1 ;
  rij-=bj->r ;
  double ninj=scalar(ni,nj), w=ninj*ninj-1 ;
  if (w==0) { ti=tj=0 ; }
  else {
    double snirij=scalar(ni,rij), snjrij=scalar(nj,rij) ;
    ti=(snirij-ninj*snjrij)/w ;
    tj=(-snjrij+ninj*snirij)/w ;
    if (ti>bi->d) ti=bi->d ;
    if (ti<-bi->d) ti=-bi->d ;
    if (tj>bj->d) tj=bj->d ;
    if (tj<-bj->d) tj=-bj->d ;
  }
  ni*=ti ; nj*=tj ;
  vec d=rij ; d+=ni ; d-=nj ;
  dij=scalar(d,d) ;    
}

void Rod_shaped::add_abs_forces(vec f) 
{
  fmag+=scalar(f,n) ;  // axias stress only
}

// increses accelerations of interacting rod-shaped cells i,j
//void Rod_shaped::update_accelerations_ij_rods(int i,int j) 
double Rod_shaped::update_accelerations_ij_rods(Rod_shaped *bi,Rod_shaped *bj) 
{
  static vec r1,r2,f1,ri1,ri2,rj1,rj2,bin,bjn,dr1,dr2,v1,v2 ;
  bin=bi->n ; bin*=bi->d ;
  bjn=bj->n ; bjn*=bj->d ;
  ri1=bi->r ; ri1+= bin ;
  ri2=bi->r ; ri2-= bin ;
  rj1=bj->r ; rj1+= bjn ;
  rj2=bj->r ; rj2-= bjn;

  double di1, di2, dj1, dj2, ti1, ti2, tj1, tj2, dij, ti, tj, dmin ;
  bj->find_sqr_distance_to_rod(ri1,di1,ti1) ; //if (di1<d02) goto g1 ;
  bj->find_sqr_distance_to_rod(ri2,di2,ti2) ; //if (di2<d02) goto g1 ; 
  bi->find_sqr_distance_to_rod(rj1,dj1,tj1) ; //if (dj1<d02) goto g1 ; 
  bi->find_sqr_distance_to_rod(rj2,dj2,tj2) ; //if (dj2<d02) goto g1 ; 
  find_sqr_distance_lines(bi,bj,dij,ti,tj) ;

  if (dij<1e-6) return 1e-6 ; // BW this should not happen!
  int sgn=0 ;
  Rod_shaped *b1=NULL, *b2=NULL ;
  if (di1<=di2 && di1<=dj1 && di1<=dj2 && di1<=dij) { dmin=di1 ; b1=bi ; b2=bj ; r1=ri1 ; r2=bj->n ; r2*=ti1 ; r2+=bj->r ; sgn=1 ; goto g2 ; }
  if (di2<=di1 && di2<=dj1 && di2<=dj2 && di2<=dij) { dmin=di2 ; b1=bi ; b2=bj ; r1=ri2 ; r2=bj->n ; r2*=ti2 ; r2+=bj->r ; sgn=-1 ; goto g2 ; }
  if (dj1<=di1 && dj1<=di2 && dj1<=dj2 && dj1<=dij) { dmin=dj1 ; b1=bj ; b2=bi ; r1=rj1 ; r2=bi->n ; r2*=tj1 ; r2+=bi->r ; sgn=1 ; goto g2 ; }
  if (dj2<=di1 && dj2<=di2 && dj2<=dj1 && dj2<=dij) { dmin=dj2 ; b1=bj ; b2=bi ; r1=rj2 ; r2=bi->n ; r2*=tj2 ; r2+=bi->r ; sgn=-1 ; goto g2 ; }
  if (dij<=di1 && dij<=di2 && dij<=dj1 && dij<=dj2) { dmin=dij ; b1=bi ; b2=bj ; r1=bi->r ; r1+=bi->n*ti ; r2=bj->r ; r2+=bj->n*tj ; goto g2 ; } 

  if (b1==NULL || b2==NULL || r1==r2) {
    char txt[100] ;
    txt[0]='!' ; txt[1]=0 ;
    for (int i=0;i<bact.size();i++) if (bj==bact[i]) {
        if (isnan(bj->r.x) || isnan(bj->r.y)) err("isnan_bj") ;
        if (isnan(bi->r.x) || isnan(bi->r.y)) err("isnan_bi") ;
  	    sprintf(txt, "xyi=(%lf %lf) xyj=(%lf %lf) %lf %lf %lf %lf %lf xxx %lf %lf %lf %lf %lf %lf",
          bj->r.x,bj->r.y,bi->r.x,bi->r.y,di1,di2,dj1,dj2,dij,ti1,ti2,tj1,tj2,ti,tj);
      }
	  err(txt) ;
  }
g2:
  dr1=r1 ; dr1-=b1->r ; dr2=r2 ; dr2-= b2->r ;
  find_forces_SS(r1, r2, f1) ;

  if (f1.iszero()) return dmin;

  b1->f+=f1 ; b1->add_abs_forces(f1) ;
  b2->f-=f1 ; b2->add_abs_forces(f1) ;

  add_cross(b1->torque,dr1,f1) ;
  subs_cross(b2->torque,dr2,f1) ;

  return dmin ;
}


double Rod_shaped::find_distance_ij(int i, int j) 
{
  Rod_shaped *bi=bact[i], *bj=bact[j] ;
  static vec r1,r2,f1,ri1,ri2,rj1,rj2,bin,bjn ;
  bin=bi->n ; bin*=bi->d ;
  bjn=bj->n ; bjn*=bj->d ;
  ri1=bi->r ; ri1+= bin ;
  ri2=bi->r ; ri2-= bin ;
  rj1=bj->r ; rj1+= bjn ;
  rj2=bj->r ; rj2-= bjn;

  if (bj->r.x-bi->r.x>width/2) { rj1.x-=width ; rj2.x-=width ; }
  if (bj->r.x-bi->r.x<-width/2) { rj1.x+=width ; rj2.x+=width ; } // periodic b.c.
  if (bj->r.y-bi->r.y>height/2) { rj1.y-=height ; rj2.y-=height ; }
  if (bj->r.y-bi->r.y<-height/2) { rj1.y+=height ; rj2.y+=height ; } // periodic b.c.

  double di1, di2, dj1, dj2, ti1, ti2, tj1, tj2, dij, ti, tj, dmin ;
  bj->find_sqr_distance_to_rod(ri1,di1,ti1) ; 
  bj->find_sqr_distance_to_rod(ri2,di2,ti2) ; 
  bi->find_sqr_distance_to_rod(rj1,dj1,tj1) ; 
  bi->find_sqr_distance_to_rod(rj2,dj2,tj2) ; 
  find_sqr_distance_lines(bi,bj,dij,ti,tj) ;

  if (di1<=di2 && di1<=dj1 && di1<=dj2 && di1<=dij) dmin=di1 ;  
  if (di2<=di1 && di2<=dj1 && di2<=dj2 && di2<=dij) dmin=di2 ;  
  if (dj1<=di1 && dj1<=di2 && dj1<=dj2 && dj1<=dij) dmin=dj1 ;  
  if (dj2<=di1 && dj2<=di2 && dj2<=dj1 && dj2<=dij) dmin=dj2 ; 
  if (dij<=di1 && dij<=di2 && dij<=dj1 && dij<=dj2) dmin=dij ; 
  return dmin ;
}


void Rod_shaped::recalc_mass_and_I() 
{
  double q=1.*r0*r0*M_PI ;
  //if (type==1) q*=motile_mass_wrt_growing_mass ; // motile particles have different mass
  fmagav+=dt*(fmag-fmagav) ;
  fmag=0 ;
  mass=q*(2*d0*d+(4/3.)*r0*d0) ; // mass in pg (1e-15 kg)
  double a=1.5*d*r0*r0+0.533*r0*r0*r0+1.333*d*d*r0+0.666*d*d*d, b=r0*r0*(d+(8./15)*r0) ;
  mom_inertia=a*q ;
}

void Rod_shaped::clear_forces() 
{
  f.zero() ; 
  torque=0 ;
}

void Rod_shaped::find_acceleration(int i)
{
  int j,k,m;
  if (active<0) return ; // do not calculate forces for inactive cells (which do not grow and are deep in the colony)
  this->update_acceleration_i() ;
  for (m=0;m<closest.size();m++) {
    j=closest[m] ; 
    if (j>i) { 

      double xjold=bact[j]->r.x ; //err("x",i) ;
      double yjold=bact[j]->r.y ; //err("x",i) ;
      if (bact[j]->r.x-bact[i]->r.x>width/2) bact[j]->r.x=bact[j]->r.x-width ; 
      if (bact[j]->r.x-bact[i]->r.x<-width/2) bact[j]->r.x=width+bact[j]->r.x ;  // periodic b.c.
      if (bact[j]->r.y-bact[i]->r.y>height/2) bact[j]->r.y=bact[j]->r.y-height ; 
      if (bact[j]->r.y-bact[i]->r.y<-height/2) bact[j]->r.y=height+bact[j]->r.y ;  // periodic b.c.

      double dmin=Rod_shaped::update_accelerations_ij_rods(bact[i],bact[j]) ;
      if (dmin<SQR(d0*2.0)) { //20220531  BW 24/06/2022 changed from 1.2 to 1.0
        n_closest_of_closest++;
        bact[j]->n_closest_of_closest++;
      } 
      if (dmin>SQR(2*max_max_length*d0+d0)) { // remove neighbour if calculated minimal distance too large
        closest[m]=closest[closest.size()-1] ; closest.pop_back() ; m-- ;
      }
      bact[j]->r.x=xjold ; 
      bact[j]->r.y=yjold ; 
    }
  }
}

void Rod_shaped::one_Euler_step(int thisi) {
  if (active<1) { f.zero() ; v.zero() ; torque=0 ; } // no acceleration

	zetabact=zeta;
	if (type==1) zetabact=zeta*zeta_motile_wrt_growing ;// change in zeta for motile cells
  v=f/(mass*zetabact) ; // um/h

	if (type==1) { v+=n*velocity_motile ;} // motility

  vec dr=v*dt ;
  vec rrr=f ; 
  r+=dr ;  
  
  if (manhattan(dr)>_maxdisp) _maxdisp=manhattan(dr) ;
  if (isnan(r.x) || isnan(r.y)) { noncriterr=1 ; return ; } //err("isnan_euler") ;

  double dphi = dt*torque/(zetabact*mom_inertia) ;
  if (d*dphi>_maxdisp) _maxdisp=d*dphi ;
  phi+=dphi ; 
  phi=fmod(phi,2*M_PI) ;
  n.x=csp=cos(phi) ; n.y=snp=sin(phi) ;

  if (_maxdisp>d0) { noncriterr=1 ; return ; } //err("jump!") ;

  pbc_location(this->r);

  dr=r+n*d - rold1 ; if (squared(dr)>SQR(0.3*d0)) nn_needs_update=1 ; //{ this->find_closest(thisi) ; return ; }
  dr=r-n*d - rold2 ; if (squared(dr)>SQR(0.3*d0)) nn_needs_update=1 ; //this->find_closest(thisi) ; 

}


void Rod_shaped::one_growth_step(int thisi) 
{
   //double *c=conc(r.x,r.y) ;
   //*uptakefun((*c)) ; if ((*c)<cutoff_c) gr=0 ; 
   //uptake=eat_rate*fitness ; if (uptake<0) uptake=0 ; 
   //(*c)-= uptake*uptakefun((*c))*dt/(dx0*dx0) ; if ((*c)<0) (*c)=0 ; // decrease concetration of food

  if (active>=1) { // only active cells can grow and replicate

    vd=0.5*growth_rate ;
    d+=dt*vd ; // increase the length of the rod
        
    if (tt>0 && d>=d_final) { // birth of a new cell
      // position of a new cell is calculated, the old cell is shifted
      // so that after replication the two cells touch each other with circular
      // caps as if the mother cell divided in half
      Rod_shaped *nb=new Rod_shaped(type) ;
      bact.push_back(nb) ;
      division_register+=1;
      n_to_kill+=1;
      nb->r=r-n*(r0+d)*0.5 ;
      r+=n*(r0+d)*0.5 ;
      nb->v=v ; 
      d=(d-r0)/2.0 ; nb->d=d ; if (d<0) err("d<0 in one_growth_step()") ;
      initial_d_s=initial_d_s+d;
      nb->phi=phi ; 
      nb->active=1 ; //nb->prev_number=number ;
      nb->vd=vd ;

      double qq=0.002 ;      
      phi+=(_drand48()-0.5)*qq ; nb->phi+=(_drand48()-0.5)*qq ; // some randomness introduced
      //phi=_drand48()*M_PI*2 ; nb->phi=_drand48()*M_PI*2 ; // randomized direction

      this->recalc_cspt() ;
      nb->recalc_cspt() ;
      
      this->nn_needs_update=1 ;
      nb->nn_needs_update=1 ;
      
    }
  } 
}

void Rod_shaped::recalc_cspt() 
{
  csp=cos(phi) ; snp=sin(phi) ; n=vec(csp,snp) ; 
}



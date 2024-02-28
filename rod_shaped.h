#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;

#include "functions.h"
#include "vecs.h"

#ifndef rod_shaped
#define rod_shaped

#include "params.h"

typedef vec2 vec ;

extern double d0, dt ;
extern double tt, r0;

class Rod_shaped {
  public:
    static int tot_number_of_bacteria;
      // these need to be saved    
    int type, number, active, nn_needs_update; // type = "species" of cell (colour), number = unique tag for each newly produced bacterium
      // active = 1 if the cell moves, =0 if it does not move but is seen by other cells, 
      // =-1 if it does not interact at all but can still consume nutrient
      // nn_needs_update if set to >0 means that find_closest() must be executed to update the list of nearest neighbours
    vector <int> closest; // closest neighbours
    int bx, by; // coordinates of the box
    int bi; // number of the bacterium in the box

      // these are updated upon loading the file

    double fmag, fmagav; // total absolute sum of forces
    // these need to be saved    
    vec r,v,rold1,rold2,f,n ; // r = position, v = velocity, rold1,2=old position when closest[] determined
    double d, dold, vd, d_final, phi, csp,snp, mass, zetabact ;
    //double motile_length;
    double torque, mom_inertia;
// d is 1/2 length of the rod (without caps)
// vd is the speed at which d grows
// d_final is 1/2 length of cell beyond which it splits

// csp, snp, cst, snt are cosines and sines of the current value of phi, theta, 
// used for speeding up some calculations
    int n_closest_of_closest; //20220531 number of the closest of the closest neighbours

    // these are updated upon loading the file
    vec fext1,fext2 ;

    static double *pow15, *pow05 ;
    static int powmax15, powmax05 ;
    static double E_bact, growth_rate, max_length, zeta, motile_length ;
    static int rods ;

#define FW(x) fwrite(&(x),sizeof((x)),1,filename) ;
#define FR(x) fread(&(x),sizeof((x)),1,filename) ;
    
    void save_data(FILE *filename) { 
      FW(type) ; FW(active) ; FW(nn_needs_update) ; 
      FW(r) ; FW(v) ; FW(rold1) ; FW(rold2) ; FW(f) ; FW(n) ; 
      FW(d) ; FW(dold) ; FW(vd) ; FW(d_final) ; FW(phi) ; FW(csp) ; FW(snp) ; FW(mass) ; FW(zetabact) ;
      FW(torque); FW(mom_inertia) ; FW(number); //AM number
      int j=closest.size() ; FW(j) ; for (int i=0;i<j;i++) FW(closest[i]) ;
    }  
    void load_data(FILE *filename) { 
      FR(type) ; FR(active) ; FR(nn_needs_update) ;
      FR(r) ; FR(v) ; FR(rold1) ; FR(rold2) ; FR(f) ; FR(n) ; 
      FR(d) ; FR(dold) ; FR(vd) ; FR(d_final) ; FR(phi) ; FR(csp) ; FR(snp) ; FR(mass) ; FR(zetabact) ;
      FR(torque); FR(mom_inertia); FR(number); //AM number
      int j ; FR(j) ; closest.resize(j) ; for (int i=0;i<closest.size();i++) FR(closest[i]) ; 
    }  

    Rod_shaped() {
      rods++ ; number=(tot_number_of_bacteria++) ; 
      closest.clear() ;
    }

    Rod_shaped(int typ) {
      rods++ ; active=1 ; type=typ ; bx=by=bi=-1 ; number=(tot_number_of_bacteria++) ; 
      d_final=max_length*d0*(1+gauss48()*0.15) ; if (d_final<0.5*d0) d_final=0.5*d0 ; // this ensures the minimal doubling length is enough to create separate circles
      d=dold=d_final*(1-_drand48()*0.5) ; 
      if (type==1)
      {
      d=dold=d_final=0.05;
      }
      vd=0 ; 
      zetabact=zeta*(0.75+0.5*_drand48()) ;
      r.zero() ; phi=0 ; v.zero() ; csp=cos(phi) ; snp=sin(phi) ; 
      
      fext1.zero() ; fext2.zero() ;
      fmag=fmagav=0 ;

      n=vec2(csp,snp) ;
      rold1=r+n*d ; rold2=r-n*d ; 
      closest.clear() ;
    }
    
    ~Rod_shaped() { rods-- ; }

    float elong_rate() { return vd ; }
    double cos_phi() { return csp ; }
    double sin_phi() { return snp ; }
    double x() { return r.x ; }
    double y() { return r.y ; }
    double area() {return d*r0*4+M_PI*r0*r0;}
    void lowest_y(vec *rr) { vec r1=r+n*d, r2=r-n*d ; if (r1.y<r2.y) (*rr)=r1 ; else (*rr)=r2 ; } ;
    void highest_y(vec *rr) { vec r1=r+n*d, r2=r-n*d ; if (r1.y>r2.y) (*rr)=r1 ; else (*rr)=r2 ; } ;
    void fill_lowest_y(float *h, int w) { 
      int j ; 
      for (float s=-d-0.5*d0;s<d+0.5*d0;s++) {
        vec rs=r+n*s ;
        j=int(rs.x+w/2+w)%w ;
        if (r.y<h[j]) h[j]=r.y ; 
      }
    } 
    void fill_highest_y(float *h, int w) {
      int j ; 
      for (float s=-d-0.5*d0;s<d+0.5*d0;s++) {
        vec rs=r+n*s ;
        j=int(rs.x+w/2+w)%w ;
        if (r.y>h[j]) h[j]=r.y ; 
      }      
    }

    double z() { return 0 ; }
    void to_string(char *str) { // this function is used to make a human-readable
    // string containing data about position, velocity, orientation etc.
      sprintf(str,"%d %d 1 %lf %lf 0  %lf %lf 0  %lf 0 %lf  %le\n",
        type,active, r.x,r.y, v.x,v.y, phi,d,  fmag*1e-12) ; 
        // fmag is saved in [N]
    }
    double vx() { return v.x ; }
    double vy() { return v.y ; }
    vec velocity(void) { return v ; }
    void add_y(double dy) { r.y+=dy ; }
    void add_x(double dx) { r.x+=dx ; }

    void init() ;
    void clear_forces() ;
    void find_closest(int) ;
    static double update_accelerations_ij_rods(Rod_shaped *bi,Rod_shaped *bj) ;   
	  static double find_distance_ij(int i, int j) ; 
    inline void find_sqr_distance_to_rod(vec &rr, double &dist2, double &t) ;
    void recalc_mass_and_I() ;
    void update_acceleration_i() ;
    void find_acceleration(int i) ;
    void one_Euler_step(int) ;
    void one_growth_step(int) ;
    void recalc_cspt() ;
    void draw() ;
    //bool inside_bact(vec) ;
    
    void add_abs_forces(vec f) ;
    
};

#endif	



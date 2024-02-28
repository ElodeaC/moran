//---------------------parameters-----------------------------

#define WHITE_LINES_INSIDE // draw thin white lines inside rod-shaped cells
#define AUTOSCALE // scale so that whole colony is always visible. Ignored for TUBE

#define VERBOSE // display additional parameters on the screen in the graphics mode

#ifdef __MAIN

#define MORAN // if defined, keeps the number of growing cells constant
double confluence=0.846 ;

// ------------ parameters of the simulation ----------------------

bool backup=true ; // store old configuration in case the simulation crashes and dt needs to be decreased
int save_cells=1+2 ; // if non-zero, save all and save positions etc. of all cells. 
                    //1==save_all | 2==save_data | 4==save_profile | 8==save_envelope
                    //16==save numbers of each type of cells
                    //32==save positions of all cells every time frame

int width=80, height=80 ; // size [um] of the simulation box for BOX, otherwise not used
// width must be a multiple of boxwidth=10*d0

double dtmax=1./(1<<12) ; // max. integration step, [h]

double death_rate=40;
double zeta_motile_wrt_growing=1 ; 
double velocity_motile=5 ;
double ini_fraction_motile=0.85 ;

double r0=0.5 ; // radius of the rod [um] (bacteria are rod-shaped with 
              // semi-circular caps on both ends, each having also radius r0 )            
double Rod_shaped::E_bact=1e+5 ; //3.75e+5 ; // elastic modulus of the bacterium, Pa
double Rod_shaped::growth_rate = 2. ; // elongation rate in [um/h]  
double Rod_shaped::max_length=1; //0.5 ; // dimensionless. 1/2 x maximal length[um] of the rod (without caps) before division - 2==E. Coli, 0.5=circles
double Rod_shaped::zeta=50 ; // dynamic friction, Pa*h. 500 for fast simulations in which the density however increases towards the centre. 
double Rod_shaped::motile_length=0.5 ; //1.0 ;

                                // 10 for slow simulation with little overlap between cells

// -------------------------------------------------------------

// this is used only for WINDOWS and graphics
int winx=1000,winy=1000 ; // bitmap size. Don't change
int _stop=0 ; // if ==1, the simulation is stopped but still displayed
double _maxdisp ; // maximal displacement in one step, used to monitor convergence

#endif


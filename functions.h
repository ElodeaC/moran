#define SQR(x)  (x)*(x)
#define MAX(x,y) (x>y?x:y)

void load_all(char *) ;
void main_proc() ;
void err(char *reason) ;
void err(char *reason, int x) ;
void err(char *reason, double x) ;
double _drand48(void) ;
void _srand48(unsigned int x0) ;
double gauss48() ;
double sign(double x) ;
void run() ;
void init() ;
void reset(bool) ;
void end() ;
double globalOrder() ;
int count_types(int t);
bool is_only_one_type() ;
double uptakefun(double c) ;
void just_mutated(int) ;
void end_capture() ;
int max_cluster_size() ;

extern char *NUM,*start_file ;
extern double avclosest, zeta_motile_wrt_growing, velocity_motile,death_rate,ini_fraction_motile ;
extern double orderinmaxcl, surfmaxcl ;
extern double vxmaxcl,vymaxcl ;
//extern vector <Bacterium*> bact ;
extern float fit_adv, _flatten, ext_torque, ext_torque_angle ;
extern int ttt, noncriterr, noclosest, RAND, gr_ok ;

extern int sample, video, plot_diff, winx,winy,_stop,shift_cm_y, width, height,maxclsize ;
extern float scale ;
extern double shiftx, shifty, r0, dt, dtmax, diff_rate, c0, _maxdisp, dx0;
//extern double** cc, ** sinkcc, ** slime, * ccone;
extern int ****tab, maxn ;
extern int nbox, boxwidth ;
extern bool backup ;

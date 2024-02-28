#ifndef vectors_already_defined
#define vectors_already_defined

class vec2  // class of 2d vectors
{
public:
	double x, y;

	vec2( double InX, double InY) : x( InX ), y( InY )
		{
		}
	vec2( ) : x(0), y(0)
		{
		}
	inline void operator= ( const vec2& V2 ) 
		{
  		x = V2.x; y = V2.y; 
		}
	inline bool operator== (const vec2& V2) const 
		{
		return (x == V2.x && y == V2.y);
		}

	inline vec2 operator+ (const vec2& V2) const 
		{
		return vec2( x + V2.x,  y + V2.y);
		}
	inline vec2 operator- (const vec2& V2) const
		{
		return vec2( x - V2.x,  y - V2.y);
		}
	inline vec2 operator- ( ) const
		{
		return vec2(-x, -y);
		}

	inline vec2 operator/ (double S ) const
		{
		double fInv = 1.0f / S;
		return vec2 (x * fInv , y * fInv);
		}
/*	inline vec3 operator* (const vec3& V2) const
		{
		return vec3 (x * V2.x,  y * V2.y,  z * V2.z);
		}*/
	inline vec2 operator* (double S) const
		{
		return vec2 (x * S,  y * S);
		}

	inline void operator+= ( const vec2& V2 )
		{
		x += V2.x;
		y += V2.y;
		}
	inline void operator-= ( const vec2& V2 )
		{
		x -= V2.x;
		y -= V2.y;
		}
	inline vec2 operator*= (double S) 
		{
		x*=S ; y*=S ; return *this ;
		}
	inline vec2 operator/= (double S) 
		{
      double S2=1./S ;
		  x*=S2 ; y*=S2 ; return *this ;
		}
  void zero() { x=y=0 ; }
  int iszero() { return (x==0 && y==0) ; }
};

inline double norm(vec2 &a)
{
  return sqrt(a.x*a.x+a.y*a.y) ;
}

inline double manhattan(vec2 &a)
{
  return fabs(a.x)+fabs(a.y) ; 
}

inline double squared(vec2 &a)
{
  return (a.x*a.x+a.y*a.y) ;
}

inline double scalar(vec2 &a, vec2 &b) {
  return a.x*b.x+a.y*b.y ; 
}

inline void normalize(vec2 &a)
{
  double l=sqrt(a.x*a.x+a.y*a.y) ;
  a.x/=l ; a.y/=l ;  
}

inline double cross(vec2 &a,vec2 &b)
{
  return (a.x*b.y-a.y*b.x) ;
}

inline void add_cross(double &c, vec2 &a,vec2 &b)
{
  c+=a.x*b.y-a.y*b.x ;
}

inline void subs_cross(double &c, vec2 &a,vec2 &b)
{
  c-=a.x*b.y-a.y*b.x ;
}


#endif

/*============================================================================

 Function      ncdf

 Usage         ncdf(q)

 Arguments     q: constant REAL representing a quantile of the normal
                  density.
                   
 Description   This is a CUDA device function that computes the normal CDF at
               q.

 Return value  REAL, representing the value of the normal CDF.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================*/

// normal cdf
__device__ REAL ncdfgpu(const REAL q)
{

 // quadrature nodes
  const REAL x1 = 0.9739065285;
  const REAL x2 = -0.9739065285;
  const REAL x3 = 0.8650633666;
  const REAL x4 = -0.8650633666;
  const REAL x5 = 0.6794095682;
  const REAL x6 = -0.6794095682;
  const REAL x7 = 0.4333953941;
  const REAL x8 = -0.4333953941;
  const REAL x9 = 0.1488743389;
  const REAL x10 = -0.1488743389;

  // quadrature weights
  const REAL w1 = 0.06667134430;
  const REAL w2 = 0.06667134430;
  const REAL w3 = 0.1494513491;
  const REAL w4 = 0.1494513491;
  const REAL w5 = 0.2190863625;
  const REAL w6 = 0.2190863625;
  const REAL w7 = 0.2692667193;
  const REAL w8 = 0.2692667193;
  const REAL w9 = 0.2955242247;
  const REAL w10 = 0.2955242247;

  // symmetry
  REAL c;
  if(q < 0){
    c = -q;
  } else {
    c = q;
  }

  REAL integral;
  // bound
  if(q >= 5.2){
    integral = 1.0;
  } else if(q <= -5.2) {
    integral = 0.0;
  } else {

    // function values
    const REAL pi = 3.14159265;
    const REAL mult = 1/sqrt(2*pi);
    const REAL f1 = mult*exp(-0.5*pow(0.5*(x1+1)*c, 2));
    const REAL f2 = mult*exp(-0.5*pow(0.5*(x2+1)*c, 2));
    const REAL f3 = mult*exp(-0.5*pow(0.5*(x3+1)*c, 2));
    const REAL f4 = mult*exp(-0.5*pow(0.5*(x4+1)*c, 2));
    const REAL f5 = mult*exp(-0.5*pow(0.5*(x5+1)*c, 2));
    const REAL f6 = mult*exp(-0.5*pow(0.5*(x6+1)*c, 2));
    const REAL f7 = mult*exp(-0.5*pow(0.5*(x7+1)*c, 2));
    const REAL f8 = mult*exp(-0.5*pow(0.5*(x8+1)*c, 2));
    const REAL f9 = mult*exp(-0.5*pow(0.5*(x9+1)*c, 2));
    const REAL f10 = mult*exp(-0.5*pow(0.5*(x10+1)*c, 2));

    if(q < 0){
      integral = 0.5 - 0.5*c*(w1*f1+w2*f2+w3*f3+w4*f4+w5*f5+w6*f6+w7*f7+w8*f8+w9*f9+w10*f10);
    } else {
      integral = 0.5 + 0.5*c*(w1*f1+w2*f2+w3*f3+w4*f4+w5*f5+w6*f6+w7*f7+w8*f8+w9*f9+w10*f10);
    }

    // enforce bounds
    if(integral > 1) integral = 1;
    if(integral < 0) integral = 0;
  }

  return integral;
}

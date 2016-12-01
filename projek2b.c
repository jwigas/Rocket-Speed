//Jediyanu wigas tu'u , rocket integration solution, succeed, with checked result
// notes : if you using gauss quadrature legendre 2pts function, the integral result is 0.000050 
//         but if you using gauss quadrature from numerical recipes book, the integral result is 0.000001
//	   I have checked the exact result using my "own hand" (manual), and the exact integral result is +/ 0.000001
//	   and the final result is 1.999814 equal with use of gauss quadrature legendre 2pts

#include <stdio.h>
#include <math.h>

#define k 0.0025
#define g 9.81  /*gravitation const*/

static const double x = 5.77350269189625764507e-01;

static const double A =  1.0;

float f(float x) { 				//function f(x)
  return log(100/(100*(1-(k*x))));
}

float qgaus( float (*f)(float) ,float a, float b) // this is gauss quadrature legendre 2pts function //
{						  
   float c = 0.5 * (b - a);
   float d = 0.5 * (b + a);
   float dum = c * x;

   return c *  A * ((*f)(d - dum) + (*f)(d + dum));
}

/*   // this is gauss quadrature numerical integration function from numerical recipes in C book //
  float qgaus(float (*f)(float),float a,float b){
  int j;
  float xr,xm,dx,s;
  static float x[]={0.0,0.1488743389,0.4333953941, 0.6794095682, 0.8650633666, 0.9739065285};
  static float w[]={0.0, 0.2955242247, 0.2692667193, 0.2190863625, 0.1494513491, 0.0666713443};
  
  xm = 0.5*(b+a);
  xr = 0.5*(b-a);
  s = 0;
    for (j=1;j<=5;j++){
      dx = xr*x[j];
      printf("dx %i : %f\n",j,dx);
      s += w[j]*((*f)(xm+dx)+(*f)(xm-dx));
      printf("iteration result %i : %f\n",j,s);
    }
    return s;
}
*/

int main(){
  float vs,v0, integral, t,x_t;
  
  printf("k constan %f and gravitation constant %f \n",k,g);
  printf("insert vs : ");scanf("%f",&vs);
  printf("insert v0 : ");scanf("%f",&v0);
  printf("insert  t : ");scanf("%f",&t); // in this time variable, i'm use 0.02 s according to instruction
  
  integral = qgaus(f,0,t);
  printf("the integral result is : %f\n",integral);
  
  x_t = (v0*t) - (0.05)*(g*t*t) + (vs*integral);

  printf("the v0*t result is : %f \n",(v0*t));
  printf("the (1/2)*(g*t*t) result is : %f \n",((0.05)*(g*t*t)));
  printf("the vs*integral result is : %f \n",(vs*integral));
  printf("the result is : %f \n",x_t);

  return 0;
}


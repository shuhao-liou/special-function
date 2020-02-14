#ifndef _ytlu_class_ORTHOPOLY_
#define _ytlu_class_ORTHOPOLY_
#include <stdio.h>
#include <math.h>
#define XCRIT (1.0e-12)
#define MAX_ITER  200
#define ORTHOPOLY laguerre

class ORTHOPOLY
{   public:
        int     order;          //order of Laguerre
        int     max_iter;   // maximum iterations for Aitken method.
        int     n_root;     // number of root found
        double *root;       // store the roots
        double *weit;
        double  xcrit;      // tolerance for root convergence.
        int    *ncount;    // no of interation in finding each root.
        int     n1_newton;  // pre-newton loops before Aitken.
        int     n2_newton;  // pro-newton looops after Aitken.

   
                  ORTHOPOLY()      {;}
                  ORTHOPOLY(int nn) {this->order = nn;}
   inline void    n(const int nn)  {this->order = nn;}
   inline int     n() const       { return this->order;}           
          double  f(const double) const;
          double  df(const double) const;
          void    setparameters();
          void    findroot(double);
          void    printroot();
          double  onestep(double);
};

double ORTHOPOLY::f(const double x) const
{
    double f0, f1, f2;
    int i;
    
    f0 = 1.0;        if(this->order == 0) return f0;
    f1 = 1.0 - x;    if(this->order == 1) return f1;
    for (i=2; i <= this->order; i++)
        {
            f2 = 2.0*f1 - f0 - ((1.0+x)*f1 -f0) / (double)(i);
            f0 = f1;
            f1 = f2;
        }    
    return f2;        
}    
double ORTHOPOLY::df(const double x) const
{
    double f0, f1, f2;
    int i;
    
    if (fabs(x)== 0.0) return (double)(-this->order); // Ln(0) = -n
    if(this->order == 0) return 0.0;                  // L0(x) = 0
    if(this->order == 1) return (-1.0);               // L1(x) = -1
    f0 = 1.0;
    f1 = 1.0 - x;
    for (i=2; i<= this->order; i++)
        {
            f2 = 2.0*f1 - f0 - ((1.0+x)*f1 -f0) / (double)(i);
            f0 = f1;
            f1 = f2;
        }    
    return ((double)(this->order) * (f1-f0) / x);  // L'n= n[Ln-L(n-1)]/x      
}

void ORTHOPOLY::findroot(double xinit)
{ double dx, x0, x1, x2, sum, gamma;
  int i, n_iter;

  this->setparameters();
  x0 = xinit;
  do {  n_iter = 0;
        x0 = x0 + 0.0001;
        for (i=0; i<5; i++) x0 = this->onestep(x0);
        do {
             x1 = this->onestep(x0);
             x2 = this->onestep(x1);
             gamma  = (fabs(x1-x0) < 1.0e-13)? gamma = 0.0 : 
                      (x2 - x1) / (x1 - x0);
             dx = gamma * (x2 - x1) / (1.0 - gamma);
             x0 = x2 + dx;
             n_iter++;
             if (n_iter > this->max_iter) return;
            } while (fabs(dx) > this->xcrit);
         this->root[this->n_root] = x0;
// compute the weit function for each root x0
// w(k) = 1.0/(x * Ln'(xk)^2)
// ytlu @ mail.ncku.edu.tw    2004/10/27
// 
         x1 = this->df(x0);
         this->weit[this->n_root] = 1.0/x0/x1/x1;
         this->ncount[n_root] = n_iter;
         n_root++;   
    } while (this->n_root < this->order);     
   return;   
}
double ORTHOPOLY::onestep(double x)
{ double dx, sum, fx, dfx;
  int i;

   fx = this->f(x);
   dfx = this->df(x);
   sum = 0.0;
   for (i=0; i < this->n_root; i++) sum += (1.0/(x - this->root[i]));
   dx = fx / (dfx - fx * sum);
   x -= dx;
   return x;   
}
void ORTHOPOLY::setparameters()
{
  xcrit = XCRIT;
  max_iter = MAX_ITER;
  n_root = 0;
  this->root = new double [this->order];
  this->weit = new double [this->order];
  this->ncount = new int [this->order];
  this->n1_newton = 5;
  this->n2_newton = 5;
  return;

}
void ORTHOPOLY::printroot()
{ int i;
   printf(" n = %d\n",this->n_root);
   printf("   k             x(k)                         w(k)\n");
   for (i=0; i<n_root; i++)
      { printf(" %3d  %25.16le    %25.16le\n",(i+1),this->root[i],
           this->weit[i]);
      }    
  return;    
}    
#endif

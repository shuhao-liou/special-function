#ifndef _ytlu_class_laguerre_
#define _ytlu_class_laguerre_
class laguerre
{ int nlag;

   public:
       laguerre(){;}
       laguerre(int n) {this->nlag = n;}
   inline void n(const int n){this->nlag = n;}
   inline int n() const { return this->nlag;}           
       double f(const double) const;
       double df(const double) const;
};

double laguerre::f(const double x) const
{
    double f0, f1, f2;
    int i;
    
    f0 = 1.0;        if(this->nlag == 0) return f0;
    f1 = 1.0 - x;    if(this->nlag == 1) return f1;
    for (i=2; i <= this->nlag; i++)
        {
            f2 = 2.0*f1 - f0 - ((1.0+x)*f1 -f0) / (double)(i);
            f0 = f1;
            f1 = f2;
        }    
    return f2;        
}    
double laguerre::df(const double x) const
{
    double f0, f1, f2;
    int i;
    
    if (fabs(x) < 1.0e-15) return (double)(-this->nlag);
    if(this->nlag == 0) return 0.0;
    if(this->nlag == 1) return (-1.0);
    f0 = 1.0;
    f1 = 1.0 - x;
    for (i=2; i<= this->nlag; i++)
        {
            f2 = 2.0*f1 - f0 - ((1.0+x)*f1 -f0) / (double)(i);
            f0 = f1;
            f1 = f2;
        }    
    return ((double)(this->nlag) * (f1-f0) / x);        
}
#endif

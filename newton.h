#ifndef _ytlu_struct_twodim_
#define _ytlu_struct_twodim_
  struct twodim
    { double x, y;

      twodim(){};
      twodim(double r1, double r2) //constructor
            { x= r1;   y = r2;     }  
    };
#endif

#ifndef _ytlu_interp_newton_
#define _ytlu_interp_newton_
class newton
{     int nintp;
      double *xintp, *yintp, *acoef;
  
  public:
      
      // constructor & destructor
      
      newton( int, double*, double*);
      ~newton()
      { delete [] xintp;
        delete [] yintp;
        delete [] acoef;
      }    
      
      //member function
      
      twodim f(const double) const;
 //     double df(const double) const;
      
}; //end of class newton

newton::newton(int n, double *x, double *y)
{ int i, j;
  double xtp, wn, sum;
  
    this->nintp = n;
    this->xintp = new double [n];
    this->yintp = new double [n];
    this->acoef = new double [n];
    for (i=0; i<n; i++) this->xintp[i] = x[i];
    for (i=0; i<n; i++) this->yintp[i] = y[i];
// compute the divided-difference coefficients.
    this->acoef[0] = y[0];
    for (i=1; i<n; i++)
      { wn = 1.0;
        sum = 0.0;
        for (j=i; j>0; j--)
            { xtp = x[i] - x[j-1];
              wn = wn * xtp;
              sum = sum * xtp + this->acoef[j-1];
            }
        this->acoef[i] = (y[i] - sum) / wn;
      }
}
/* 
double newton::f(const double xt) const
{ double sum;
  int i;
    sum = 0.0;
    for (i=nintp; i>0; i--)
       { sum = sum * (xt - this->xintp[i-1]) + this->acoef[i-1];
       }
    return(sum);    
}
*/    
twodim newton::f(const double xt) const
{ double sum, dsum;
  int i;
    sum = 0.0;
    dsum = 0.0;
    for (i=nintp; i>0; i--)
       { 
         dsum = dsum * (xt - xintp[i-1]) + sum;
         sum = sum * (xt - xintp[i-1]) + acoef[i-1];
       }
    return twodim(sum, dsum);    
}//end of class newton.      
#endif            

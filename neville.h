#ifndef _ytlu_struct_twodim_
#define _ytlu_struct_twodim_
  struct twodim
    { double x, y;

      twodim(){};
      twodim(double r1, double r2) //constructor
            { x= r1;   y = r2;     }  
    };
#endif

#ifndef _ytlu_interp_neville_
#define _ytlu_interp_neville_
class neville
{     int nintp;
      double *xintp, *yintp;
  
  public:
      
      // constructor
      
      neville( int, double*, double*);
      //member function
      
      twodim f(const double) const;
//      double df(const double) const;
}; //end of class interp

neville::neville(int n, double *x, double *y)
{ int i;
    this->nintp = n;
    this->xintp = new double [n];
    this->yintp = new double [n];
    for (i=0; i<n; i++) this->xintp[i] = x[i];
    for (i=0; i<n; i++) this->yintp[i] = y[i];
} 
/*
double neville::f(const double xi) const
{ double x1, x2, f1, f2, *ytt;
  int i, j;
  
  ytt = new double [nintp];
  for (i=0; i<nintp; i++) ytt[i] = yintp[i];
  for (i=1; i<nintp; i++)
     { for (j=0; j<(nintp-i); j++)
          { x1 = xintp[j];
            x2 = xintp[j+i];
            f1 = ytt[j];
            f2 = ytt[j+1];
            ytt[j] = (f1* (xi-x2) - f2*(xi-x1)) / (x1-x2);
           } 
     }
  f2 = ytt[0];
  delete [] ytt;
  return(f2);
}
*/
twodim neville::f(const double xi) const
{ double x1, x2, f1, f2, *ytt;
  double df1, df2, *dff;
  int i, j;

  
  ytt = new double [nintp];
  dff = new double [nintp];
  for (i=0; i<nintp; i++) ytt[i] = yintp[i];
  for (i=0; i<nintp; i++) dff[i] = 0.0;
  for (i=1; i<nintp; i++)
     { for (j=0; j<(nintp-i); j++)
          { x1 = xintp[j];
            x2 = xintp[j+i];
            f1 = ytt[j];
            f2 = ytt[j+1];
            df1 = dff[j];
            df2 = dff[j+1];
            ytt[j] = (f1* (xi-x2) - f2*(xi-x1)) / (x1-x2);
            dff[j] = (f1 + df1*(xi-x2) - f2 - df2*(xi-x1))/(x1-x2);
           } 
     }
  f2 = ytt[0];
  df2 = dff[0];
  delete [] ytt;
  delete [] dff;
  return  twodim(f2, df2);
}
//end of class neville.      
#endif

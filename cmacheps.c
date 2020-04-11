/*

Purpose:  approximate [machine+GCC compiler] epsilon

Build:  gcc -o cmacheps64lebeau cmacheps.c
Build:  gcc -o cmacheps32hogan cmacheps.c

Run:  ./macheps64lebeau

*/
#include <stdio.h>

int
main() {
double  dbl_one = 1.0;
float  flt_one = 1.0;
double dbl, dbl_eps;
float flt, flt_eps;

  dbl_eps = 1.0;
  dbl = 1.5;
  while (dbl != dbl_one) {
    dbl_eps *= 0.5;
    dbl = 1.0 + (dbl_eps*0.5);
  }
  fprintf(stdout,"%.3e ~ Double epsilon, variable comparison\n", dbl_eps);

  flt_eps = 1.0;
  flt = 1.5;
  while (flt != flt_one) {
    flt_eps *= 0.5;
    flt = 1.0 + (flt_eps*0.5);
  }
  fprintf(stdout,"%.3e ~ Float epsilon, variable comparison\n", flt_eps);



  dbl_eps = 1.0;
  while (((double)(1.0+(0.5*dbl_eps))) != dbl_one) {
    dbl_eps *= 0.5;
  }
  fprintf(stdout,"%.3e ~ Double epsilon, cast double expression comparison\n", dbl_eps);

  flt_eps = 1.0;
  while (((float)(1.0+(0.5*flt_eps))) != flt_one) {
    flt_eps *= 0.5;
  }
  fprintf(stdout,"%.3e ~ Float epsilon, cast float expression comparison\n", flt_eps);



  dbl_eps = 1.0;
  while ((1.0+(0.5*dbl_eps)) != dbl_one) {
    dbl_eps *= 0.5;
  }
  fprintf(stdout,"%.3e ~ Double epsilon, expression comparison\n", dbl_eps);

  flt_eps = 1.0;
  while ((1.0+(0.5*flt_eps)) != flt_one) {
    flt_eps *= 0.5;
  }
  fprintf(stdout,"%.3e ~ Float epsilon, expression comparison\n", flt_eps);

  return 0;
}


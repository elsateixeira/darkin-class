#include "common.h"
#include "evolver_ndf15.h"

static int nonfinite_derivs(
                            double x,
                            double * y,
                            double * dy,
                            void * parameters_and_workspace,
                            ErrorMsg error_message
                            ) {
  int index;
  int neq = *((int *)parameters_and_workspace);

  (void)x;
  (void)y;
  (void)error_message;

  for (index=0; index<neq; index++) {
    dy[index] = NAN;
  }

  return _SUCCESS_;
}

int main(void) {
  struct jacobian jac;
  struct numjac_workspace nj_ws;
  ErrorMsg error_message;
  double y[3] = {0.,1.,2.};
  double fval[3] = {0.,1.,1.};
  int neq = 2;
  int nfe = 0;
  int status;

  error_message[0] = '\0';

  if (initialize_jacobian(&jac,neq,error_message) == _FAILURE_) {
    fprintf(stderr,"Could not initialize Jacobian: %s\n",error_message);
    return _FAILURE_;
  }
  if (initialize_numjac_workspace(&nj_ws,neq,error_message) == _FAILURE_) {
    fprintf(stderr,"Could not initialize numjac workspace: %s\n",error_message);
    uninitialize_jacobian(&jac);
    return _FAILURE_;
  }

  status = numjac(nonfinite_derivs,
                  0.,
                  y,
                  fval,
                  &jac,
                  &nj_ws,
                  1.e-15,
                  neq,
                  &nfe,
                  &neq,
                  error_message);

  uninitialize_numjac_workspace(&nj_ws);
  uninitialize_jacobian(&jac);

  if (status != _FAILURE_) {
    fprintf(stderr,"numjac accepted a non-finite derivative column.\n");
    return _FAILURE_;
  }
  if (strstr(error_message,"non-finite derivative") == NULL) {
    fprintf(stderr,"numjac returned the wrong diagnostic: %s\n",error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;
}

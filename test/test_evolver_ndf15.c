#include "common.h"
#include "evolver_ndf15.h"

struct root_test_parameters {
  double root;
  double upper_bound;
  short quadratic;
};

static int bounded_root_function(
                                 double * x,
                                 int x_size,
                                 void * parameters,
                                 double * output,
                                 ErrorMsg error_message
                                 ) {
  struct root_test_parameters * test = parameters;

  class_test(x_size != 1,
             error_message,
             "expected a one-dimensional root test, got %d",
             x_size);
  class_test((isfinite(x[0]) == 0) ||
             (x[0] <= 0.) ||
             (x[0] > test->upper_bound),
             error_message,
             "root-test point %.16e is outside (0,%.16e]",
             x[0],
             test->upper_bound);

  if (test->quadratic == _TRUE_)
    output[0] = x[0]*x[0]-test->root*test->root;
  else
    output[0] = x[0]-test->root;

  return _SUCCESS_;
}

static int check_root_solver(
                             double initial_x,
                             double root,
                             double upper_bound,
                             short quadratic,
                             double tolx,
                             double expected_tolerance,
                             const char * label
                             ) {
  struct root_test_parameters test = {root,upper_bound,quadratic};
  ErrorMsg error_message;
  double x = initial_x;
  double dxdF = 1.;
  double output;
  int fevals = 0;

  error_message[0] = '\0';
  if (fzero_Newton(bounded_root_function,
                   &x,
                   &dxdF,
                   1,
                   tolx,
                   1.e-10,
                   &test,
                   &fevals,
                   error_message) == _FAILURE_) {
    fprintf(stderr,"%s failed: %s\n",label,error_message);
    return _FAILURE_;
  }
  if (bounded_root_function(&x,1,&test,&output,error_message) == _FAILURE_ ||
      fabs(x-root) > expected_tolerance ||
      fabs(output) > 1.e-10) {
    fprintf(stderr,
            "%s returned x=%.16e, residual=%.16e after %d evaluations\n",
            label,x,output,fevals);
    return _FAILURE_;
  }

  return _SUCCESS_;
}

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

  /* The nominal positive Jacobian probe crosses the upper domain boundary;
     fzero_Newton must shrink it before taking the linear step to the root. */
  if (check_root_solver(0.1,0.05,0.1002,_FALSE_,1.e-12,1.e-10,
                        "Jacobian-probe retry") == _FAILURE_)
    return _FAILURE_;

  /* The full Newton step from x=0.1 for x^2=1 is outside the allowed domain;
     the backtracking update must recover the positive root. */
  if (check_root_solver(0.1,1.,2.,_TRUE_,1.e-12,1.e-8,
                        "backtracked Newton step") == _FAILURE_)
    return _FAILURE_;

  /* The central seed is invalid, while the deterministic positive restart is
     computable and must converge. */
  if (check_root_solver(0.,0.1,1.,_FALSE_,1.e-12,1.e-10,
                        "seeded restart") == _FAILURE_)
    return _FAILURE_;

  /* With a deliberately loose step tolerance, the first Newton update is
     small enough to satisfy tolx but still has a large residual. */
  if (check_root_solver(2.,1.,10.,_TRUE_,1.,1.e-8,
                        "small-step residual guard") == _FAILURE_)
    return _FAILURE_;

  return _SUCCESS_;
}

template<class Type>
struct diffusionkernel_sp {
  vector<Type> times;              
  matrix<Type> locs;
  Type beta;
  vector<Type> w;
  matrix<Type> xyloc;
  matrix<Type> Qbase;
  diffusionkernel_sp(vector<Type> times_, matrix<Type> locs_, Type beta_,
		     vector<Type> w_, matrix<Type> xyloc_, matrix<Type> Qbase_)   // Constructor of integrand
    : times(times_), locs(locs_), beta(beta_), w(w_), xyloc(xyloc_), Qbase(Qbase_) {}       // Initializer list
  vector<Type> ratesep_sp(Type t){
    using namespace density;
    using atomic::tiny_ad::isfinite;
    matrix<Type> ans(times.size(), w.size());
    ans.setZero();
    for (int i = 0; i < times.size(); ++i)
      if (t > times[i]){
        matrix<Type> Q = Qbase * (t - times[i]);
        MVNORM_t<Type> bivnorm_sp(Q);
        for (int j = 0; j < w.size(); ++j){
          vector<Type> loci = xyloc.row(j) - locs.row(i);
          Type temp = exp(-beta * (t - times[i]) - bivnorm_sp(loci));
          if (isfinite(temp)) ans(i, j) += temp;
        }
      }
    return ans * w;
  }
  Type operator()(Type t){  // Evaluate integrand
    return sum(this->ratesep_sp(t));
  }
};

// checking whether the point is in any triangle specified by mesh$graph$tv from R side.
// mesh$graph$tv stores the index information about triangulation.
// if a given point is in any of the triangulation, it can be written as a positive 
// linear combination of these three points and the weights sum up to 1.
// This is a simple way to check, and it is not the same as the one in INLA.
template<class Type>
int pointinSPbare_sp(vector<Type> loci, matrix<Type> xyloc, matrix<int> tv){
  matrix<Type> A(tv.cols(), tv.cols());
  vector<Type> b(loci.size() + 1);
  int ind = -1;
  b << loci, Type(1.);
  for (int i = 0; i < tv.rows(); ++i){
    A.col(0) << xyloc(tv(i, 0) - 1, 0), xyloc(tv(i, 0) - 1, 1), Type(1.);
    A.col(1) << xyloc(tv(i, 1) - 1, 0), xyloc(tv(i, 1) - 1, 1), Type(1.);
    A.col(2) << xyloc(tv(i, 2) - 1, 0), xyloc(tv(i, 2) - 1, 1), Type(1.);
    vector<Type> temp = atomic::matinv(A) * b;
    if ((temp.array() >= 0).all()){
      ind = i;
      break;
    }
  }
  
  return ind; // output the index of tv, if the given point is not in any of the triangulation, output -1;
}

template<class Type>
bool pointinSP_sp(vector<Type> loci, matrix<Type> xyloc, matrix<int> tv){
  int ind = pointinSPbare_sp(loci, xyloc, tv);
  if (ind >= 0)
    return true;
  else
    return false;
}

// generate bivariate normal sample with mean loci and covariance matrix sigma.
// This function uses cholesky decomposition to decompose the covariance matrix.
template<class Type>
vector<Type> rbivnorm_sp(vector<Type> loci, matrix<Type> sigma){
  vector<Type> ans(2);
  density::MVNORM_t<Type>(sigma).simulate(ans);
  ans += loci;
  return ans;
}

// obtain the predicator vector for a single observation
// checked with INLA. The results are not ordered according to x;
template<class Type>
vector<Type> predproj(vector<Type> loci, matrix<Type> xyloc, matrix<int> tv, int & ind){
  using namespace Eigen;
  // run this after running pointinSPbare_sp.
  vector<Type> x1 = xyloc.row(tv(ind, 0) - 1), x2 = xyloc.row(tv(ind, 1) - 1), x3 = xyloc.row(tv(ind, 2) - 1);
  vector<Type> ans(3); // temp = (k, b);
  // first index
  ans[0] = Hyperplane<Type, 2>::Through(x2, x3).absDistance(loci) * sqrt((x2[0] - x3[0]) * (x2[0] - x3[0]) + (x2[1] - x3[1]) * (x2[1] - x3[1])) * Type(.5);
  // second index
  ans[1] = Hyperplane<Type, 2>::Through(x1, x3).absDistance(loci) * sqrt((x1[0] - x3[0]) * (x1[0] - x3[0]) + (x1[1] - x3[1]) * (x1[1] - x3[1])) * Type(.5);
  // third index
  ans[2] = Hyperplane<Type, 2>::Through(x1, x2).absDistance(loci) * sqrt((x2[0] - x1[0]) * (x2[0] - x1[0]) + (x2[1] - x1[1]) * (x2[1] - x1[1])) * Type(.5);
  return ans / sum(ans);
}

// generate sample from non-homogenous poisson plane
template<class Type>
vector<Type> rpoisplane_sp(matrix<Type> xyloc, matrix<int> tv, vector<Type> x, vector<Type> w){
  vector<Type> ans(2);
  Type xmin = xyloc.col(0).minCoeff(), xmax = xyloc.col(0).maxCoeff();
  Type ymin = xyloc.col(1).minCoeff(), ymax = xyloc.col(1).maxCoeff();
  int ind;
  bool mflag = false;
  Type envlp = exp(x.maxCoeff()), u, temp;
  while (!mflag){
    ans[0] = runif(xmin, xmax); ans[1] = runif(ymin, ymax); // generate uniformly over [min(x coord of mesh nodes), max(x coord of mesh nodes)] * [min(y coord of mesh nodes), max(y coord of mesh nodes)]
    ind = pointinSPbare_sp(ans, xyloc, tv); // check whether the generated point is in the area of interest.
    // rejection sampling:
    // The one point generated above is uniform, but in spde case the rate varies (random effect x)
    // use runif(0, exp(max(x))) to perform rejection sampling.
    if (ind >= 0){
      u = runif(Type(0.), envlp);
      vector<Type> lincomb = predproj(ans, xyloc, tv, ind); // get the row in predicator matrix
      temp = exp(x[tv(ind, 0) - 1] * lincomb[0] + x[tv(ind, 1) - 1] * lincomb[1] + x[tv(ind, 2) - 1] * lincomb[2]); // calculate the rate given this generated point.
      if (u < temp){ // rejection step.
        mflag = true;
      }
    }
  }
  return ans;
}

#ifndef spde_hawkes_hpp
#define spde_hawkes_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj


template<class Type>
Type spde_hawkes(objective_function<Type>* obj) {
  using namespace R_inla; // Where Q_spde is defined.
  using namespace density; // this where the structure for GMRF and AR_t is defined
  using namespace Eigen;  // probably for sparseness class
  // vector of time
  DATA_VECTOR(times);
  DATA_MATRIX(locs);
  DATA_VECTOR(w);
  DATA_MATRIX(xyloc);
  DATA_SCALAR(tmax);
  DATA_INTEGER(simple); // If 1, use time-independent Gaussian fields
  // parameters of the hawkes process
  DATA_MATRIX(designmat); // first column is 1s, the rest are spatial covariates 
  PARAMETER_VECTOR(coefs); // first entry is log_mu, the rest are coefficients
  PARAMETER(logit_abratio);
  PARAMETER(log_beta);
  PARAMETER(log_xsigma);
  PARAMETER(log_ysigma);
  PARAMETER(atanh_rho);
  Type log_mu = coefs[0];
  Type mu = exp(coefs[0]);
  Type beta = exp(log_beta);
  Type alpha = exp(logit_abratio + log_beta) / (Type(1.) + exp(logit_abratio)); // enforcing 0 <= alpha <= beta;
  Type xsigma = exp(log_xsigma);
  Type ysigma = exp(log_ysigma);
  Type rho = tanh(atanh_rho);

  Type nll = 0.0;

  DATA_SPARSE_MATRIX(lmat); // the predicator matrix.
  DATA_STRUCT(spde,spde_t); // this as structure of spde object is defined in r-inla hence using that namespace
  PARAMETER(log_kappa); 
  PARAMETER(log_tau);
  PARAMETER_VECTOR(x); // random field.
  Type kappa = exp(log_kappa);
  Type tau = exp(log_tau);


  SparseMatrix<Type> Q = Q_spde(spde,kappa) * pow(tau, 2); // create the precision matrix from the spde model for the GMRF
  SIMULATE {
    // simulate the random field
    rnorm_fill(x);
    x = GMRF_t<Type>(Q).sqrt_cov_scale(x);
  }
  nll += GMRF(Q)(x); // field the random effect is a GMRF with precision Q

  // For the definition of the explanation, see Reinhart (2018).

  // term-1
  // An APPROXIMATION of the integral of the field over the area.
  // int_{t, s} mu(s).
  vector<Type> L = exp(x + designmat * coefs);
  nll += (L * w).sum() * tmax;

  // term 2
  // sum_{i = 1}^n log(lambda(s_i, t_i)), where lambda(s, t) = mu(s_i) + \sum_{i: t_i < t} g(s - s_i, t - t_i)
  matrix<Type> Qbase(2, 2), Q2(2, 2);
  Qbase << exp(Type(2.) * log_xsigma), rho * exp(log_ysigma + log_xsigma), rho * exp(log_ysigma + log_xsigma), exp(Type(2.) * log_ysigma);
  vector<Type> loci(2);
  vector<Type> A(times.size());
  A.setZero();
  if (simple == 0){
    for (int j = 0; j < times.size(); ++j)
      for (int i = 0; i < j; ++i)
	if (times[j] - times[i] > 0){
	  Q2 = Qbase * (times[j] - times[i]);
	  loci = locs.row(j) - locs.row(i);
	  A[j] += exp(-beta * (times[j] - times[i]) - MVNORM(Q2)(loci));
	}
  } else {
    MVNORM_t<Type> bivnorm_sp(Qbase);
    for (int j = 1; j < times.size(); ++j)
      for (int i = 0; i < j; ++i){
        loci = locs.row(j) - locs.row(i);
        A[j] += exp(-beta * (times[j] - times[i]) - bivnorm_sp(loci));
      }
  }
  vector<Type> C = log(lmat * L + alpha * A);
  nll -= sum(C);

  // term 3
  if (simple == 0){
    diffusionkernel_sp<Type> diffker(times, locs, beta, w, xyloc, Qbase);
    nll += alpha * romberg::integrate(diffker, Type(0.), tmax);
  } else {
    MVNORM_t<Type> bivnorm_sp2(Qbase);
    vector<Type> marks(times.size());
    matrix<Type> ans(times.size(), w.size());
    for (int j = 0; j < w.size(); ++j)
      for (int k = 0; k < times.size(); ++k){
        vector<Type> loci = xyloc.row(j) - locs.row(k);
        ans(k, j) = exp(-bivnorm_sp2(loci));
      }
    // For the purposes of integrating lambda, can treat like a marked model. 
    // The mark is the volume of the Gaussian within the domain (0<=V<=1)
    marks =  ans * w; 
    vector<Type> B = vector<Type>::Zero(times.size());
    
    for(int i = 1; i < times.size(); ++i){
      B[i] = exp(-beta * (times[i] - times[i - 1])) * (marks[i - 1] + B[i - 1]);
    }
    
    nll += (alpha/beta) * Type(sum(marks) - marks.template tail<1>()[0] - B.template tail<1>()[0]);
  }

  SIMULATE {
    // This simulation process follows Algorithm 4 in Section 3.3 of Reinhart (2018).
    DATA_IMATRIX(tv);
    /* 
       Need the triangulation to check whether points generated are in the area of interest.

       For the points generated by non-homogeneous Poisson intensity, the area is chosen and then,
       the point is generated by bivariate normal distribution with mean at corresponding mesh point and 
       very small covariance matrix, see notes on why a random sample is needed.

       Notes: generated points can not be on the mesh point. This will cause issue with evaluating
       lambda_X(t) in Section 3.3. Basically we need to avoid the density evaluation at the mean.
    */

    locs.setZero();
    times.setConstant(INFINITY);
    // // try to implement as exactly done in the paper.
    int i = 0;
    Type temp, urate, D = (w * exp(log_mu + x)).array().sum();
    vector<Type> locibase(2);
    vector<Type> lambdaXsasep(w.size());
    // Step 1.
    Type Ub = runif(Type(0.), Type(1.)), gammac = D, ua = -log(Ub) / gammac, sa = ua;
    // Step 2.
    times[0] = ua;
    // Step 7 for the first time. 
    loci = rpoisplane_sp(xyloc, tv, x, w);
    locs.row(0) = loci;
    // Step 10 for the first time.
    // std::cout<<times[i]<<" "<<locs(i, 0)<<" "<<locs(i, 1)<<" "<<gammac<<std::endl;
    i++;
    // Repeat from Step 3
    while (i < times.size()){
      // Step 3.
      Ub = runif(Type(0.), Type(1.));
      ua = -log(Ub) / gammac;
      // Step 4.
      sa += ua;
      Ub = runif(Type(0.), Type(1.));
      // Step 5.
      diffusionkernel_sp<Type> diffker2(times, locs, beta, w, xyloc, Qbase);
      lambdaXsasep = diffker2.ratesep_sp(sa) * alpha;
      if (Ub > (sum(lambdaXsasep) + D) / gammac){
        gammac = sum(lambdaXsasep) + D;
        continue;
      }
      // Step 6.
      times[i] = sa;
      urate = runif(Type(0.), sum(lambdaXsasep) + D);
      if (urate > D){
        // Step 8
        temp = D;
        for (int j = 0; j < lambdaXsasep.size(); ++j){
          temp += lambdaXsasep[j];
          if (temp > urate) {
            locibase = locs.row(j);
            Q2 = Qbase * (times[i] - times[j]);
            loci = rbivnorm_sp(locibase, Q2);
            if (pointinSP_sp(loci, xyloc, tv)){
              // Step 10
              locs.row(i) = loci;
              // std::cout<<times[i]<<" "<<locs(i, 0)<<" "<<locs(i, 1)<<" "<<gammac<<std::endl;
              i++;
            }
            break;
          }
        }
        // else Step 9.
      }else{
        // Step 7
        loci = rpoisplane_sp(xyloc, tv, x, w);
        locs.row(i) = loci;
        // Step 10.
        // std::cout<<times[i]<<" "<<locs(i, 0)<<" "<<locs(i, 1)<<" "<<gammac<<std::endl;
        i++;
      }
      gammac = D + sum(lambdaXsasep);
    }

    tmax = max(times);
    REPORT(times);
    REPORT(locs);
    REPORT(tmax);
  }

  ADREPORT(mu);
  ADREPORT(coefs);
  ADREPORT(alpha);
  ADREPORT(beta);
  ADREPORT(xsigma);
  ADREPORT(ysigma);
  ADREPORT(rho);
  ADREPORT(kappa);
  ADREPORT(tau);
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif

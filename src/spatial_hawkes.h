#ifndef spatial_hawkes_hpp
#define spatial_hawkes_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
struct diffusionkernel {
  vector<Type> times;              
  matrix<Type> locs;
  Type beta;
  vector<Type> w;
  matrix<Type> xyloc;
  matrix<Type> Qbase;
  diffusionkernel(vector<Type> times_, matrix<Type> locs_, Type beta_,
    vector<Type> w_, matrix<Type> xyloc_, matrix<Type> Qbase_)   // Constructor of integrand
  : times(times_), locs(locs_), beta(beta_), w(w_), xyloc(xyloc_), Qbase(Qbase_) {}       // Initializer list
  vector<Type> ratesep(Type t){
    using namespace density;
    matrix<Type> ans(times.size(), w.size());
    ans.setZero();
    //MVNORM_t<Type> bivnorm(Qbase);
    for (int i = 0; i < times.size(); ++i)
      if (t > times[i]){
        matrix<Type> Q = Qbase * (t - times[i]);
        MVNORM_t<Type> bivnorm(Q);
        for (int j = 0; j < w.size(); ++j){
          vector<Type> loci = xyloc.row(j) - locs.row(i);
          ans(i, j) = exp(-beta * (t - times[i]) - bivnorm(loci));
        }
      }
    return ans * w;
  }
  Type operator()(Type t){  // Evaluate integrand
    return sum(this->ratesep(t));
  }
};

// checking whether the point is in any triangle specified by mesh$graph$tv from R side.
// mesh$graph$tv stores the index information about triangulation.
// and use triangle interior to check whether a given point is in any of the triangulation.
// There may be better way to do this.
template<class Type>
vector<Type> pointinSPbare(vector<Type> loci, matrix<Type> xyloc, matrix<int> tv, int & ind){
  matrix<Type> A(tv.cols(), tv.cols());
  vector<Type> b(loci.size() + 1), ans(tv.cols());
  ind = -1;
  b << loci, Type(1.);
  for (int i = 0; i < tv.rows(); ++i){
    A.col(0) << xyloc(tv(i, 0) - 1, 0), xyloc(tv(i, 0) - 1, 1), Type(1.);
    A.col(1) << xyloc(tv(i, 1) - 1, 0), xyloc(tv(i, 1) - 1, 1), Type(1.);
    A.col(2) << xyloc(tv(i, 2) - 1, 0), xyloc(tv(i, 2) - 1, 1), Type(1.);
    vector<Type> temp = atomic::matinv(A) * b;
    if ((temp.array() >= 0).all()){
      ans = temp;
      ind = i;
      break;
    }
  }
  
  return ans;
}

// wrap it up. Check whether loci is in any triangulation specified by tv.
template<class Type>
bool pointinSP(vector<Type> loci, matrix<Type> xyloc, matrix<int> tv){
  int ind;
  pointinSPbare(loci, xyloc, tv, ind);
  if (ind >= 0)
    return true;
  else
    return false;
}

// generate bivariate normal sample with mean loci and covariance matrix sigma.
// This function uses cholesky decomposition to decompose the covariance matrix.
template<class Type>
vector<Type> rbivnorm(vector<Type> loci, matrix<Type> sigma){
  vector<Type> ans(2);
  density::MVNORM_t<Type>(sigma).simulate(ans);
  ans += loci;
  return ans;
}

// generate sample from homogenous poisson plane (in this case uniform distribution)
template<class Type>
vector<Type> rpoisplane(matrix<Type> xyloc, matrix<int> tv){
  vector<Type> ans(2);
  // get the range of the total field.
  Type xmin = xyloc.col(0).minCoeff(), xmax = xyloc.col(0).maxCoeff();
  Type ymin = xyloc.col(1).minCoeff(), ymax = xyloc.col(1).maxCoeff();
  bool mflag = false;
  while (!mflag){
    ans << runif(xmin, xmax), runif(ymin, ymax);
    // make rejections to the observation if it is outside of the inner triangulation.
    // This may be in sufficient if the shape of the field is weird.
    if (pointinSP(ans, xyloc, tv)){
      mflag = true;
    }
  }
  return ans;
}




template<class Type>
Type spatial_hawkes(objective_function<Type>* obj) {
  using namespace R_inla; // Where Q_spde is defined.
  using namespace density; // this where the structure for GMRF and AR_t is defined
  using namespace Eigen;  // probably for sparseness class
  // vector of time
  DATA_VECTOR(times);
  DATA_MATRIX(locs);
  DATA_VECTOR(w);
  DATA_SPARSE_MATRIX(lmat); // the predicator matrix
  DATA_MATRIX(xyloc);
  DATA_SCALAR(tmax);
  DATA_INTEGER(simple); // If 1, use time-independent Gaussian fields
  // parameters of the Hawkes process
  DATA_MATRIX(designmat); // first column is 1s, the rest are spatial covariates 
  PARAMETER_VECTOR(coefs); // first entry is log_mu, the rest are coefficients
  PARAMETER(logit_abratio);
  PARAMETER(log_beta);
  PARAMETER(log_xsigma);
  PARAMETER(log_ysigma);
  PARAMETER(atanh_rho);
  Type mu = exp(coefs[0]);
  Type beta = exp(log_beta);
  Type alpha = exp(logit_abratio + log_beta) / (Type(1.) + exp(logit_abratio)); // enforcing 0 <= alpha <= beta;
  Type xsigma = exp(log_xsigma);
  Type ysigma = exp(log_ysigma);
  Type rho = tanh(atanh_rho);

  Type nll = 0.0;

  // For the definition of the explanation, see equation 8 of Reinhart (2018).
  //https://projecteuclid.org/journals/statistical-science/volume-33/issue-3/A-Review-of-Self-Exciting-Spatio-Temporal-Point-Processes-and/10.1214/17-STS629.full

  // term-1
  // An APPROXIMATION of the integral of the field over the area.
  // int_{t, s} mu(s).
  vector<Type> L = exp(designmat * coefs);
  nll += (L * w).sum() * tmax;

  // term 2
  // sum_{i = 1}^n log(lambda(s_i, t_i)), where lambda(s, t) = mu(s_i) + \sum_{i: t_i < t} g(s - s_i, t - t_i)
  matrix<Type> Qbase(2, 2), Q2(2, 2);
  Qbase << xsigma * xsigma, rho * xsigma * ysigma, rho * xsigma * ysigma, ysigma * ysigma;
  vector<Type> loci(2);
  vector<Type> A(times.size());
  A.setZero();
  if (simple == 0) {
    for (int j = 1; j < times.size(); ++j)
      for (int i = 0; i < j; ++i){
        Q2 = Qbase * (times[j] - times[i]);
        loci = locs.row(j) - locs.row(i);
        A[j] += exp(-beta * (times[j] - times[i]) - MVNORM(Q2)(loci)); // MVNORM returns -log of density
      }
  } else {
    MVNORM_t<Type> bivnorm(Qbase);
    for (int j = 1; j < times.size(); ++j)
      for (int i = 0; i < j; ++i){
        loci = locs.row(j) - locs.row(i);
        A[j] += exp(-beta * (times[j] - times[i]) - bivnorm(loci));
        }
  }
  vector<Type> C = log(lmat * L + alpha * A);
  nll -= sum(C);

  // term 3
  if (simple == 0){
    diffusionkernel<Type> diffker(times, locs, beta, w, xyloc, Qbase);
    nll += alpha * romberg::integrate(diffker, Type(0.), tmax);
  } else {
    MVNORM_t<Type> bivnorm2(Qbase);
    vector<Type> marks(times.size());
    matrix<Type> ans(times.size(), w.size());
    for (int j = 0; j < w.size(); ++j)
      for (int k = 0; k < times.size(); ++k){
        vector<Type> loci = xyloc.row(j) - locs.row(k);
        ans(k, j) = exp(-bivnorm2(loci));
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

  SIMULATE { // Only for constant background without covariates
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
    Type temp, urate, D = mu * sum(w);
    vector<Type> locibase(2);
    vector<Type> lambdaXsasep(w.size());
    // Step 1.
    Type Ub = runif(Type(0.), Type(1.)), gammac = D, ua = -log(Ub) / gammac, sa = ua;
    // Step 2.
    times[0] = ua;
    // Step 7 for the first time. 
    loci = rpoisplane(xyloc, tv);
    locs.row(0) = loci;
    // Step 10 for the first time.
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
      diffusionkernel<Type> diffker2(times, locs, beta, w, xyloc, Qbase);
      lambdaXsasep = diffker2.ratesep(sa) * alpha;
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
            loci = rbivnorm(locibase, Q2);
            if (pointinSP(loci, xyloc, tv)){
              // Step 10
              locs.row(i) = loci;
              i++;
            }
            break;
          }
        }
        // else Step 9.
      }else{
        // Step 7
        loci = rpoisplane(xyloc, tv);
        locs.row(i) = loci;
        // Step 10.
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
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif

#include <fstream>
#include <TMB.hpp>
#include <vector>
#include <iostream>
#include <numeric>
#include <math.h>
// #include "init.h"
#include "hawkes.h"
#include "lgcp.h"
#include "marked_lgcp.h"
#include "custom_hawkes.h"
#include "neg_alpha_custom_hawkes.h"
#include "neg_alpha_hawkes.h"
#include "spatial_hawkes.h"
#include "spde_hawkes.h"
#include "multi_hawkes.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(model_type);
  if (model_type == "hawkes"){
    return hawkes(this);
  } else
  if (model_type == "lgcp") {
      return lgcp(this);
  } else
  if (model_type == "marked_lgcp") {
    return marked_lgcp(this);
  } else
  if (model_type == "custom_hawkes") {
    return custom_hawkes(this);
  } else
  if (model_type == "neg_alpha_custom_hawkes") {
    return neg_alpha_custom_hawkes(this);
  } else
  if (model_type == "neg_alpha_hawkes") {
    return neg_alpha_hawkes(this);
  } else
  if (model_type == "spde_hawkes") {
    return spde_hawkes(this);
  } else
  if (model_type == "spatial_hawkes") {
    return spatial_hawkes(this);
  } else 
  if (model_type == "multi_hawkes") {
    return multi_hawkes(this);
  } else {
 	Rf_error("Unknown model.");
   }
   return 0;
}

  

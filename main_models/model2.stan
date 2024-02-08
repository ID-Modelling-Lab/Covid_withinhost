functions {
  real[] OneSys(real t, real[] y, real[] params, real[] x_r, int[] x_i) {
  real dydt[4]; 

  dydt[1] = 0.14 * (8 * 10^7 - y[1]) - params[1] * y[1] * y[3];
  dydt[2] = params[1] * y[1] * y[3] - 0.14 * y[2] - params[2] * y[4] * y[2];
  dydt[3] = 1.12 * 10^4 * y[2] - params[3] * y[3];
  dydt[4] = params[4] * y[2] * y[4];
  
  return dydt;
  }
}


data {
  int<lower = 0> L; //no. of rows in the dataset
  int<lower = 0> VGtsv; // no. of vacgr + time since vac split
  int<lower = 0> AGtsv; // no. of age groups + time since vac split
  int l_tsolve; //length of vector of times to solve
  real logvirus[L]; //swab data
  real t0; //initial value for t
  real ts[L]; //all days (arranged by individual)
  int ts_ind[L]; //index of days (arranged by individual)
  real t_solve[l_tsolve]; //vector of times to solve
  int pos_VGtsvofAGtsv[AGtsv]; // vector of which vgtsv group the agtsv belongs to
  int pos_AGofAGtsv[AGtsv]; // vector of which age group (1-7) the agtsv belongs to
  int pos_agtsv[AGtsv]; // vector of indices of first swab for each AGtsv group
  int count_agtsv[AGtsv]; // vector of no. of swabs in each AGtsv group
  real rel_tol;
  real abs_tol;
  int max_num_steps;
  int which_fixed; //which parameter to fix as reference group
  real fixed_value; //value to fix reference group at
}


transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  real<lower = 10^-12, upper = 10^-9> beta_par;
  real<lower = 0.01, upper = 10> c_par;
  real<lower = 0, upper = 10> gamma_par[VGtsv];
  real<lower = 0> rho;
  real<lower = 0> eta;
  real<lower = 0, upper = 1> omega_par[VGtsv];
  real<lower = 0> phi;
  real<lower = 0> lambda;
  real<lower = 0> theta_mult[5]; // total 6 values for theta, theta[1] = 1.0
}

transformed parameters {
  real predVal[L];
  real params[4];
  real y0[4];
  real<lower = 0> theta[6];
  // real<lower = 0> gamma_mod[AGtsv];
  real<lower = 0> omega_mod[AGtsv];
  
  {
    y0[1] = 8*10^7;
    y0[2] = 0;
    y0[3] = 1;
    y0[4] = 10;
  
    params[1] = beta_par;
    params[3] = c_par;
    
    // checks
    print("beta: ", params[1]);
    print("c: ", params[3]);
    //
    
    // theta multiplier for model variants 2-4
    theta[which_fixed] = fixed_value;
    for (i in 2:6) {
      theta[i] = theta[which_fixed] * theta_mult[i-1];
    }
    //

    // iterate over total age groups
    for (ag in 1:AGtsv) {
      real yout[l_tsolve, 4];
      
      // model 2, gamma vg omega ag
      omega_mod[ag] = theta[pos_AGofAGtsv[ag]] * omega_par[pos_VGtsvofAGtsv[ag]];
      
      params[2] = gamma_par[pos_VGtsvofAGtsv[ag]];
      params[4] = omega_mod[ag];

      
      // // model 3, gamma ag omega vg
      // gamma_mod[ag] = theta[pos_AGofAGtsv[ag]] * gamma_par[pos_VGtsvofAGtsv[ag]];
      // 
      // params[2] = gamma_mod[ag];
      // params[4] = omega_par[pos_VGtsvofAGtsv[ag]];
      // 
      // // model 4, gamma ag omega ag
      // gamma_mod[ag] = theta[pos_AGofAGtsv[ag]] * gamma_par[pos_VGtsvofAGtsv[ag]];
      // omega_mod[ag] = theta[pos_AGofAGtsv[ag]] * omega_par[pos_VGtsvofAGtsv[ag]];
      // 
      // params[2] = gamma_mod[ag];
      // params[4] = omega_mod[ag];
      

      yout = integrate_ode_bdf(OneSys, y0, t0, t_solve, params, x_r, x_i, rel_tol, abs_tol, max_num_steps);
      predVal[pos_agtsv[ag] : (pos_agtsv[ag] + count_agtsv[ag] - 1)] = yout[ts_ind[pos_agtsv[ag] : (pos_agtsv[ag] + count_agtsv[ag] - 1)], 3];
      if (ag == 5) {
      print("yout[1] ", log10(yout[1,3]));
      print("gamma: ", params[2]);
      print("theta: ", theta[pos_AGofAGtsv[ag]]);
      print("omega_par: ", omega_par[pos_VGtsvofAGtsv[ag]]);
      print("omega_mod: ", omega_mod[ag]);
      
      }

    }
    
  }
  
}


model {
  omega_par ~ beta(phi, lambda);
  log10(gamma_par) ~ beta(rho, eta);
  
  for (i in 1:L) {
    logvirus[i] ~ normal(log(predVal[i]), 1);
  }
}

generated quantities {
  vector[L] log_lik;
  real R0_imm[VGtsv];
  real sumloglike;
  real R0;
  sumloglike = 0;
  
  for (i in 1:L) {
    log_lik[i] = normal_lpdf(logvirus[i] | log(predVal[i]), 1);
    sumloglike += log_lik[i];
  }
  
  for (i in 1:VGtsv) {
    R0_imm[i] = (beta_par * (1.12 * 10^4) * (8 * 10^7)) / (c_par *(0.14 + gamma_par[i] * 10));
  }
  
  R0 = (beta_par * (1.12 * 10^4) * (8 * 10^7)) / (c_par *0.14);
}

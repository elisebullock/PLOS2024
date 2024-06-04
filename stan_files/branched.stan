functions{
#include functions.stan

  real intplus(int j, vector y, int numofki67int, int numofbrduint, int start) 
  {
      // Calculating a jumping index (calculates indexes down a ki67 strand)
      array[numofki67int+2] int indices;
      for (i in 1:(numofki67int+2)) {
          indices[i]=start+j+(numofbrduint+2)*(i-1);
        }
      real summed = sum(y[indices]);
      
      return summed;   
  }
  
  vector KIAB_labeling(real t,  vector y,
    // unpack parameters
     real alpha_A,
     real alpha_B,
     real delta_A,
     real delta_B,
     real beta,
     real fs,
     real Source,
     real mu,
     real eff,
     int numofki67int,
     int numofbrduint
     ) 
   { // BRDU labeling of A and B

        real numki67comp = numofki67int+1;
        int switch1 = (numofbrduint+2)*(numofki67int+2);
        vector[switch1*2] dy;
           
        // compute derivatives
        // Directionality of the array 
        // A++ Aint+ A-+ A+int Aintint A-int A+- Aint- A--
          
        //A++
        dy[1] = Source*fs - (delta_A+numki67comp*beta+alpha_A)*y[1] + 2*alpha_A*eff*sum(y[1:switch1]);
        
        //Aint+
        for (j in 0:(numofbrduint-1)) {
            dy[2+j] = - (delta_A+numki67comp*beta+alpha_A)*y[2+j]+2*alpha_A*(1-eff)*intplus(j,y,numofki67int,numofbrduint,1);
        }
        
        //A-+
        dy[2+numofbrduint] = -(delta_A+alpha_A+numki67comp*beta)*y[2+numofbrduint] +2*alpha_A*(1-eff)*(intplus(numofbrduint,y,numofki67int,numofbrduint,1)+intplus(numofbrduint+1,y,numofki67int,numofbrduint,1));
            
        //A(+ int -)int 
        for (i in 0:(numofki67int-1)) {
            for (j in 0:(numofbrduint+1)) {
                dy[(2+numofbrduint)*(1+i)+j+1] = -(delta_A+alpha_A+numki67comp*beta)*y[(2+numofbrduint)*(1+i)+j+1]+numki67comp*beta*y[(2+numofbrduint)*(i)+j+1];
            }
        }
            
        //A+-
        dy[(1+numofki67int)*(2+numofbrduint)+1] = -(delta_A+alpha_A)*y[(1+numofki67int)*(2+numofbrduint)+1]+numki67comp*beta*y[(numofki67int)*(2+numofbrduint)+1];
        
        //Aint-
        for (j in 0:(numofbrduint-1)) {
            dy[(1+numofki67int)*(2+numofbrduint)+2+j] = -(delta_A+alpha_A)*y[(1+numofki67int)*(2+numofbrduint)+2+j]+numki67comp*beta*y[(numofki67int)*(2+numofbrduint)+2+j];
        }
        
        //A--
        dy[switch1]  = -(delta_A+alpha_A)*y[switch1]+numki67comp*beta*y[(2+numofbrduint)*(1+numofki67int)];
         
           
        //B++
        dy[switch1+1] =  Source*(1-fs) -(delta_B+numki67comp*beta+alpha_B)*y[switch1+1] + 2*alpha_B*eff*sum(y[(switch1+1):]); 
        
        //Bint+
        for (j in 0:(numofbrduint-1)) {
            dy[switch1+2+j] = - (delta_B+numki67comp*beta+alpha_B)*y[switch1+2+j] + 2*alpha_B*(1-eff)*intplus(j,y,numofki67int,numofbrduint,switch1+1);
         }
         
        //B-+
        dy[switch1+2+numofbrduint] = -(delta_B+numki67comp*beta+alpha_B)*y[switch1+2+numofbrduint]+2*alpha_B*(1-eff)*(intplus(numofbrduint,y,numofki67int,numofbrduint,switch1+1)+intplus(numofbrduint+1,y,numofki67int,numofbrduint,switch1+1));
        
        //B(+ int -)int
        for (i in 0:(numofki67int-1)) {
            for (j in 0:(numofbrduint+1)) {
                dy[switch1+(2+numofbrduint)*(1+i)+j+1] = -(delta_B+alpha_B+numki67comp*beta)*y[switch1+(2+numofbrduint)*(1+i)+j+1]+numki67comp*beta*y[switch1+(2+numofbrduint)*(i)+j+1];
            }
        }
        
        //B+-
        dy[switch1+(1+numofki67int)*(2+numofbrduint)+1] = -(delta_B+alpha_B)*y[switch1+(1+numofki67int)*(2+numofbrduint)+1] + numki67comp*beta*y[switch1+(numofki67int)*(2+numofbrduint)+1];
        
        //Bint-
        for (j in 0:(numofbrduint-1)) {
            dy[switch1+(1+numofki67int)*(2+numofbrduint)+2+j] = -(delta_B+alpha_B)*y[switch1+(1+numofki67int)*(2+numofbrduint)+2+j]+numki67comp*beta*y[switch1+(numofki67int)*(2+numofbrduint)+2+j];
        }
        
        //B--
        dy[2*switch1] = -(delta_B+alpha_B)*y[2*switch1]+numki67comp*beta*y[(2+numofbrduint)*(1+numofki67int)+switch1];
           
    return dy;
  }
  
  vector KIAB_unlabeling(real t,  vector y,
    // unpack parameters
     real alpha_A,
     real alpha_B,
     real delta_A,
     real delta_B,
     real beta,
     real fs,
     real Source,
     real mu,
     real eff,
     int numofki67int,
     int numofbrduint,
     real switchtime
     ) 
  { // Unlabeling of BRDU
  
          int switch1 = (numofbrduint+2)*(numofki67int+2);
          real exp_decay_value = exp(-mu*(t-switchtime));
          real numki67comp = numofki67int+1;
          vector[switch1*2] dy;
          
          // compute derivatives
          // Directionality of the array 
          // A++ Aint+ A-+ A+int Aintint A-int A+- Aint- A--
                  
            //A++
            dy[1] = Source*fs*exp_decay_value - (delta_A+numki67comp*beta+alpha_A)*y[1];
            
            //Aint+
            for (j in 0:(numofbrduint-1)) {
                dy[2+j] = - (delta_A+numki67comp*beta+alpha_A)*y[2+j] + 2*alpha_A*intplus(j,y,numofki67int,numofbrduint,1);
            }
           
            //A-+
            dy[2+numofbrduint] = Source*fs*(1-exp_decay_value)-(delta_A+alpha_A+numki67comp*beta)*y[2+numofbrduint] +2*alpha_A*(intplus(numofbrduint,y,numofki67int,numofbrduint,1)+intplus(numofbrduint+1,y,numofki67int,numofbrduint,1));
            
            //A(+ int -)int 
            for (i in 0:(numofki67int-1)) {
                for (j in 0:(numofbrduint+1)) {
                    dy[(2+numofbrduint)*(1+i)+j+1] = -(delta_A+alpha_A+numki67comp*beta)*y[(2+numofbrduint)*(1+i)+j+1]+numki67comp*beta*y[(2+numofbrduint)*(i)+j+1];
                }
            }
              
            //A+-
            dy[(1+numofki67int)*(2+numofbrduint)+1] = -(delta_A+alpha_A)*y[(1+numofki67int)*(2+numofbrduint)+1]+numki67comp*beta*y[(numofki67int)*(2+numofbrduint)+1];
                  
            //Aint-
            for (j in 0:(numofbrduint-1)) {
                dy[(1+numofki67int)*(2+numofbrduint)+2+j] = -(delta_A+alpha_A)*y[(1+numofki67int)*(2+numofbrduint)+2+j]+numki67comp*beta*y[(numofki67int)*(2+numofbrduint)+2+j];
            }
                  
            //A--
            dy[switch1]  = -(delta_A+alpha_A)*y[switch1]+numki67comp*beta*y[(2+numofbrduint)*(1+numofki67int)];

            //B++
            dy[switch1+1] = Source*(1-fs)*exp_decay_value -(delta_B+numki67comp*beta+alpha_B)*y[switch1+1];
            
            //Bint+
            for (j in 0:(numofbrduint-1)) {
                dy[switch1+2+j] = - (delta_B+numki67comp*beta+alpha_B)*y[switch1+2+j]  + 2*alpha_B*intplus(j,y,numofki67int,numofbrduint,switch1+1);
            }
            
            //B-+
            dy[switch1+2+numofbrduint] = Source*(1-fs)*(1-exp_decay_value) -(delta_B+numki67comp*beta+alpha_B)*y[switch1+2+numofbrduint]+2*alpha_B*(intplus(numofbrduint,y,numofki67int,numofbrduint,switch1+1)+intplus(numofbrduint+1,y,numofki67int,numofbrduint,switch1+1)); 
              
             //B(+ int -)int
             for (i in 0:(numofki67int-1)) {
                 for (j in 0:(numofbrduint+1)) {
                     dy[switch1+(2+numofbrduint)*(1+i)+j+1] = -(delta_B+alpha_B+numki67comp*beta)*y[switch1+(2+numofbrduint)*(1+i)+j+1]+numki67comp*beta*y[switch1+(2+numofbrduint)*(i)+j+1];
                 }
             }
             
            //B+-
            dy[switch1+(1+numofki67int)*(2+numofbrduint)+1] = -(delta_B+alpha_B)*y[switch1+(1+numofki67int)*(2+numofbrduint)+1] + numki67comp*beta*y[switch1+(numofki67int)*(2+numofbrduint)+1];
             
            //Bint-
            for (j in 0:(numofbrduint-1)) {
                dy[switch1+(1+numofki67int)*(2+numofbrduint)+2+j] = -(delta_B+alpha_B)*y[switch1+(1+numofki67int)*(2+numofbrduint)+2+j]+numki67comp*beta*y[switch1+(numofki67int)*(2+numofbrduint)+2+j]; 
            }
             
            //B--
            dy[2*switch1] = -(delta_B+alpha_B)*y[2*switch1] + numki67comp*beta*y[(2+numofbrduint)*(1+numofki67int)+switch1];
            
    return dy;
  }
 
  vector equilibrium_ode(real t,  vector y,
      // unpack parameters
       real alpha_A,
       real alpha_B,
       real delta_A,
       real delta_B,
       real beta,
       real fs,
       real Source,
       real mu,
       real eff,
       int numofki67int,
       int numofbrduint
       ) 
    { // Calculating the steady state values
        
        real numki67comp = numofki67int+1;
        int switch1 = (numofbrduint+2)*(numofki67int+2);
        vector[switch1*2] dy;
        
        // compute derivatives
        // Directionality of the array 
        // A++ Aint+ A-+ A+int Aintint A-int A+- Aint- A--
                        
        //A++
        dy[1] = 0;
        
        //Aint+
        for (j in 0:(numofbrduint-1)) {
            dy[2+j] = 0;
        }
            
        //A-+
        dy[2+numofbrduint] = Source*fs - (delta_A+alpha_A+numki67comp*beta)*y[2+numofbrduint] + 2*alpha_A*sum(y[:(switch1)]);
        
        
        //A(+ int -)int 
        for (i in 0:(numofki67int-1)) { 
            for (j in 0:(numofbrduint+1)) { 
                if(j==(numofbrduint+1)) {
                    dy[(2+numofbrduint)*(1+i)+j+1] = -(delta_A+alpha_A+numki67comp*beta)*y[(2+numofbrduint)*(1+i)+j+1]+numki67comp*beta*y[(2+numofbrduint)*(i)+j+1];
                    }
                else{
                    dy[(2+numofbrduint)*(1+i)+j+1] = 0;
                }
            }
        }
            
        //A+-
        dy[(1+numofki67int)*(2+numofbrduint)+1] = 0;
           
        //Aint-
        for (j in 0:(numofbrduint-1)) {
            dy[(1+numofki67int)*(2+numofbrduint)+2+j] = 0;
        }
        
        //A--
        dy[switch1]  = -(delta_A+alpha_A)*y[switch1]+numki67comp*beta*y[(2+numofbrduint)*(1+numofki67int)];
               
        //B++
        dy[switch1+1] = 0;
        
        //Bint+
        for (j in 0:(numofbrduint-1)) {
            dy[switch1+2+j] = 0;
        }
         
        //B-+
        dy[switch1+2+numofbrduint] = Source*(1-fs) -(delta_B+numki67comp*beta+alpha_B)*y[switch1+2+numofbrduint]+2*alpha_B*sum(y[(switch1+1):]);
        
        //B(+ int -)int
        for (i in 0:(numofki67int-1)) {
            for (j in 0:(numofbrduint+1)) {
                if(j==(numofbrduint+1)) {
                    dy[switch1+(2+numofbrduint)*(1+i)+j+1] = -(delta_B+alpha_B+numki67comp*beta)*y[switch1+(2+numofbrduint)*(1+i)+j+1]+numki67comp*beta*y[switch1+(2+numofbrduint)*(i)+j+1];
                }
                else {
                    dy[switch1+(2+numofbrduint)*(1+i)+j+1] = 0;
                }
            }
        }
        
        //B+-
        dy[switch1+(1+numofki67int)*(2+numofbrduint)+1] = 0;
        
        //Bint-
        for (j in 0:(numofbrduint-1)) {
            dy[switch1+(1+numofki67int)*(2+numofbrduint)+2+j] = 0;
        }
                
        //B--
        dy[2*switch1] = -(delta_B+alpha_B)*y[2*switch1] + numki67comp*beta*y[(2+numofbrduint)*(1+numofki67int)+switch1];
        
    return  dy;
  }
    
    vector lamda(real t,  vector y,
      // unpack parameters
       real alpha_A,
       real alpha_B,
       real delta_A,
       real delta_B
       ) 
    { // Calculating the steady state values
        
        vector[2] dy;
        
        // compute derivatives
        // Directionality of the array 
        // A++ Aint+ A-+ A+int Aintint A-int A+- Aint- A--
                        
        //Fast
        dy[1] = -(delta_A-alpha_A)*y[1];
   
        //A-+
        dy[2] = -(delta_B-alpha_B)*y[2];
        
        
    return  dy;
  }   
      
  array[] vector solveode(vector init_cond, array[] real time_index, array[] real time_index_equilibrium,
        real alpha_A,
        real alpha_B,
        real delta_A,
        real delta_B,
        real beta,
        real fs,
        real Source,
        real mu,
        real eff,
        int numofki67int,
        int numofbrduint,
        int switchtime,
        int switch1,
        int numObs,
        int index
        )
    {
          
        // code for either non-swtiching or for real data
        array[numObs] vector[2*switch1] k_hat; //= ode_rk45(KIAB_labeling, init_conditions[1], 0.0, time_index, alpha_A, alpha_B, delta_A, delta_B, beta, fs, Source, mu,eff, numofki67int,numofbrduint);
                
        k_hat[1:index] = ode_rk45(KIAB_labeling, init_cond, 0.0, time_index[1:index], alpha_A, alpha_B, delta_A, delta_B, beta, fs, Source, mu,eff, numofki67int,numofbrduint);
        k_hat[(index+1):numObs] = ode_rk45(KIAB_unlabeling, k_hat[index], switchtime, time_index[(index+1):numObs], alpha_A, alpha_B, delta_A, delta_B, beta, fs, Source, mu,eff, numofki67int,numofbrduint, switchtime);
        
    return k_hat;
  }
    
}

data{
    int<lower  = 1>    numObs;        // number of observations for donor fractions may be different that cell counts
    int switch1;                      // Switching array point between A and B cells
    array[numObs] real time_index;    // time array corresponding to data points
    array[numObs] real kcounts;       // Counts of all cells
    array[numObs] int brdhi_kihi_c;
    array[numObs] int brdhi_kilo_c;
    array[numObs] int brdlo_kilo_c;
    array[numObs] int brdlo_kihi_c;
    array[100] real fulltime;
    
    int switchtime;                   // Time where BRDU is stopped, unlabeling starts
    vector[2*switch1] init_cond_obs;  // Initial conditions used for equilibrium setting
    vector[2] init_lamdas;  // Initial conditions used for equilibrium setting
    int numofki67int;                 // Num of ki67 intermediates
    int numofbrduint;                 // Num of BRDU intermediates
    array[1] real time_index_equilibrium; // time array for steady state equili
    int counts_t0;
    real scalescale;
    real stdcounts;
    int indexswitch;
    
    real alpha_A_data;
    real alpha_B_data;
    real delta_A_data;
    real delta_B_data;
    real gamma_data;
    real fs_data;
    real beta_data;
    real Source_data;
    real eff_data;
    real mu_data;
     
   real Source_hi;
   real Source_lo;
 
 //real Source;                      // Naive pool
 //real fs;
 //real eff;                         // Efficiency of BRDU uptake upon division
 //real mu;                          // Degradation of BRDU during unlabeling phase
 //real alpha_B;
 //real delta_A;
 //real delta_B;
 //real gamma;
 //real beta;

 
  }

transformed data{
  array[numObs,4] int data_fractional;

  // data transformations
  for (i in 1:numObs){
      data_fractional[i,1] = brdhi_kilo_c[i];
      data_fractional[i,2] = brdhi_kihi_c[i];
      data_fractional[i,3] = brdlo_kilo_c[i];
      data_fractional[i,4] = brdlo_kihi_c[i];
  }
}



parameters{
    // priors for parameters
    //real<lower = gamma_data/10,     upper = gamma_data*100>    gamma;
    real<lower = delta_A_data/2,    upper = delta_A_data*2>   delta_A; //death of fast
    real<lower = alpha_A_data/2,    upper = delta_A>          alpha_A; //division of fast
    real<lower = alpha_B_data/2,    upper = alpha_A>          alpha_B; //division of slow
    real<lower = alpha_B,           upper = delta_A>          delta_B; //death of slow
    real<lower = fs_data/2,         upper = 0.95>             fs;
    real<lower = beta_data/2,       upper = beta_data*2>      beta;
    real<lower = Source_lo,         upper = Source_hi>    Source;
    real<lower = eff_data/2,        upper = 1>                eff;
    real<lower = mu_data/2,         upper = mu_data*2>        mu;
 
    
    // priors for data likelihoods
    real<lower = 0>           scale;
    real<lower = 0>           phi_inv;
  }



  transformed parameters{
    vector<lower = 0>[numObs]   c_total;        // Total cell counts
    vector<lower = 0,upper=1>[numObs]  f_brdhi_kihi;   // fraction of BRDuHi ki67Hi
    vector<lower = 0,upper=1>[numObs]  f_brdhi_kilo;   // fraction of BRDuHi ki67Lo 
    vector<lower = 0,upper=1>[numObs]  f_brdhi;        // fraction of BRDUHi
    vector<lower = 0,upper=1>[numObs]  f_kihi;         // fraction of ki67Hi 
  
      
    vector<lower = 0>[numObs]  c_kihi;        // fraction of ki67Hi 
    vector<lower = 0>[numObs]  c_kilo;        // fraction of ki67Hi 
    vector<lower = 0>[numObs]  c_brdlo;       // fraction of ki67Hi 
    vector<lower = 0>[numObs]  c_brdhi_kihi;  // fraction of ki67Hi
    vector<lower = 0>[numObs]  c_brdhi_kilo;  // fraction of ki67Hi
    vector<lower = 0>[numObs]  c_brdlo_kilo;  // fraction of ki67Hi  
        
    array[numObs] vector[4] norm_c_total;    // Simplex to store the normalized quadrant counts

    array[1] vector[2*switch1] init_conditions  = ode_bdf(equilibrium_ode, init_cond_obs, 0.0, time_index_equilibrium, alpha_A, alpha_B, delta_A, delta_B, beta, fs, Source, mu,eff, numofki67int,numofbrduint);

    array[numObs] vector[2*switch1] k_hat = solveode(init_conditions[1], time_index, time_index_equilibrium,
      alpha_A,
         alpha_B,
         delta_A,
         delta_B,
         beta,
         fs,
         Source,
         mu,
         eff,
         numofki67int,
         numofbrduint,
         switchtime,
         switch1,
         numObs,
         indexswitch); 

    for (i in 1:numObs){
        // calculating
        c_brdlo_kilo[i] =  k_hat[i,switch1]+k_hat[i,2*switch1];
        
        c_kihi[i] = sum(k_hat[i,1:((2+numofbrduint)*(1+numofki67int))])+sum(k_hat[i,(1+switch1):(switch1+(2+numofbrduint)*(1+numofki67int))]);
        c_brdlo[i] = intplus(numofbrduint+1, k_hat[i,:], numofki67int,  numofbrduint, 1) + intplus(numofbrduint+1, k_hat[i,:], numofki67int,  numofbrduint, switch1+1);
        
        // total counts
        c_total[i] = sum(k_hat[i, 1:2*switch1]);
        

        c_kilo[i] = c_total[i] - c_kihi[i];
        c_brdhi_kihi[i] = c_total[i] - c_brdlo[i] - sum(k_hat[i,((2+numofbrduint)*(1+numofki67int)+1):(switch1-1)])-sum(k_hat[i,((2+numofbrduint)*(1+numofki67int)+1+switch1):(2*switch1-1)]);

        // fractions of BRDUhi cells in ki67Hi
        f_brdhi_kihi[i] = (c_brdhi_kihi[i])/c_kihi[i];
      
        // fractions of BRDUhi cells in KI67lo
        c_brdhi_kilo[i] = (sum(k_hat[i,(1+(2+numofbrduint)*(1+numofki67int)):(switch1-1)])+sum(k_hat[i,(switch1+1+(2+numofbrduint)*(1+numofki67int)):(2*switch1-1)]));
        f_brdhi_kilo[i] = c_brdhi_kilo[i]/c_kilo[i];
        
        // fractions of BRDUHi
        f_brdhi[i] = (c_total[i]-c_brdlo[i])/c_total[i];
        
        // fractions of ki67Hi
        f_kihi[i] = c_kihi[i]/c_total[i];
                
        // Simplex formation of quadrants normalized to total
        norm_c_total[i,1] = c_brdhi_kilo[i]/c_total[i];
        norm_c_total[i,2] = c_brdhi_kihi[i]/c_total[i];
        norm_c_total[i,3] = c_brdlo_kilo[i]/c_total[i];
        norm_c_total[i,4] = (c_total[i]-c_brdhi_kilo[i]-c_brdhi_kihi[i]-c_brdlo_kilo[i])/c_total[i];

    }
    
    
}



model{
    //priors for liklihoods
    scale ~ lognormal(log(stdcounts), scalescale);
    phi_inv ~ exponential(5);

    //priors based on 4cm data
    alpha_A ~ lognormal(log(alpha_A_data),  0.5);
    alpha_B ~ lognormal(log(alpha_B_data),  1);
    delta_A ~ lognormal(log(delta_A_data),  0.5);
    delta_B ~ lognormal(log(delta_B_data),  1);
    beta    ~ lognormal(log(beta_data),     0.1);
    //gamma   ~ lognormal(log(gamma_data),    1);
    fs      ~ lognormal(log(fs_data),       0.5);
    mu      ~ lognormal(log(mu_data),       0.1);
    eff     ~ lognormal(log(eff_data),      0.1);
    Source  ~ lognormal(log(Source_data),   0.5);

    // prior for total counts at t=0
    counts_t0 ~ lognormal(log(sum(init_conditions[1,:])), stdcounts*10/numObs);

    for (i in 1:numObs) {
        data_fractional[i,:] ~ dirichlet_multinomial((norm_c_total[i,:]/phi_inv));  
    }

}



generated quantities{
    //posteriors
    real ppd_alpha_A = lognormal_rng(log(alpha_A_data), 0.5);
    real ppd_alpha_B = lognormal_rng(log(alpha_B_data), 1);
    real ppd_delta_A = lognormal_rng(log((delta_A_data)), 0.1);
    real ppd_delta_B = lognormal_rng(log(delta_B_data), 1);
    real ppd_beta = lognormal_rng(log(beta_data), 0.5);
    real ppd_gamma = lognormal_rng(log(gamma_data+0.001), 1);
    real ppd_fs = lognormal_rng(log(fs_data), 0.5);

    real ppd_Source = lognormal_rng(log(Source_data), 0.1);
    real ppd_mu = lognormal_rng(log(mu_data), 0.5);
    real ppd_eff = lognormal_rng(log(eff_data), 0.1);

    vector[numObs] log_lik;
    vector[numObs] data_counts_hat;
    array[numObs,4] int  norm_f_total_hat;
    
    for (i in 1:numObs) {
    	log_lik[i] = dirichlet_multinomial_lpmf(data_fractional[i,:] | norm_c_total[i,:]/phi_inv);
                  
        norm_f_total_hat[i,:] = dirichlet_multinomial_rng(norm_c_total[i,:]/phi_inv,counts_t0); 
        
        data_counts_hat[i] = lognormal_rng(log(c_total[i]),scale);

        }
    
       array[100] vector[2*switch1] k_hat_calc = solveode(init_conditions[1], fulltime, time_index_equilibrium, 
      alpha_A,
        alpha_B,
         delta_A,
         delta_B,
         beta,
         fs,
         Source,
         mu,
         eff,
         numofki67int,
         numofbrduint,
         switchtime,
         switch1,
         100,
         50); 
         
            vector<lower = 0>[100]   c_total_calc;        // Total cell counts
            vector<lower = 0,upper=1>[100]  f_brdhi_kihi_calc;   // fraction of BRDuHi ki67Hi
            vector<lower = 0,upper=1>[100]  f_brdhi_kilo_calc;   // fraction of BRDuHi ki67Lo 
            vector<lower = 0,upper=1>[100]  f_brdhi_calc;        // fraction of BRDUHi
            vector<lower = 0,upper=1>[100]  f_kihi_calc;         // fraction of ki67Hi 
          
            vector<lower = 0>[100]  c_kihi_calc;        // fraction of ki67Hi 
            vector<lower = 0>[100]  c_kilo_calc;        // fraction of ki67Hi 
            vector<lower = 0>[100]  c_brdlo_calc;       // fraction of ki67Hi 
            vector<lower = 0>[100]  c_brdhi_kihi_calc;  // fraction of ki67Hi
            vector<lower = 0>[100]  c_brdhi_kilo_calc;  // fraction of ki67Hi
            vector<lower = 0>[100]  c_brdlo_kilo_calc;  // fraction of ki67Hi  

            array[100] vector[4] norm_c_total_calc;    // Simplex to store the normalized quadrant counts
            vector[100] clonal;
            vector[2] init_lamdas_=init_lamdas;
            init_lamdas_[1] = init_lamdas_[1]*fs;
            init_lamdas_[2] = init_lamdas_[1]*(1-fs);
            
            array[100] vector[2] lamdas  = ode_rk45(lamda,init_lamdas_, 0.0, linspaced_array(100, 1, 1000), alpha_A, alpha_B, delta_A, delta_B);

            for (i in 1:100){
                // calculating
                c_brdlo_kilo_calc[i] =  k_hat_calc[i,switch1]+k_hat_calc[i,2*switch1];
                
                c_kihi_calc[i] = sum(k_hat_calc[i,1:((2+numofbrduint)*(1+numofki67int))])+sum(k_hat_calc[i,(1+switch1):(switch1+(2+numofbrduint)*(1+numofki67int))]);
                c_brdlo_calc[i] = intplus(numofbrduint+1, k_hat_calc[i,:], numofki67int,  numofbrduint, 1) + intplus(numofbrduint+1, k_hat_calc[i,:], numofki67int,  numofbrduint, switch1+1);
                
                // total counts
                c_total_calc[i] = sum(k_hat_calc[i, 1:2*switch1]);
                

                c_kilo_calc[i] = c_total_calc[i] - c_kihi_calc[i];
                c_brdhi_kihi_calc[i] = c_total_calc[i] - c_brdlo_calc[i] - sum(k_hat_calc[i,((2+numofbrduint)*(1+numofki67int)+1):(switch1-1)])-sum(k_hat_calc[i,((2+numofbrduint)*(1+numofki67int)+1+switch1):(2*switch1-1)]);

                // fractions of BRDUhi cells in ki67Hi
                f_brdhi_kihi_calc[i] = (c_brdhi_kihi_calc[i])/c_kihi_calc[i];
              
                // fractions of BRDUhi cells in KI67lo
                c_brdhi_kilo_calc[i] = (sum(k_hat_calc[i,(1+(2+numofbrduint)*(1+numofki67int)):(switch1-1)])+sum(k_hat_calc[i,(switch1+1+(2+numofbrduint)*(1+numofki67int)):(2*switch1-1)]));
                f_brdhi_kilo_calc[i] = c_brdhi_kilo_calc[i]/c_kilo_calc[i];
                
                // fractions of BRDUHi
                f_brdhi_calc[i] = (c_total_calc[i]-c_brdlo_calc[i])/c_total_calc[i];
                
                // fractions of ki67Hi
                f_kihi_calc[i] = c_kihi_calc[i]/c_total_calc[i];
                
                norm_c_total_calc[i,1] = c_brdhi_kilo_calc[i]/c_total_calc[i];
                norm_c_total_calc[i,2] = c_brdhi_kihi_calc[i]/c_total_calc[i];
                norm_c_total_calc[i,3] = c_brdlo_kilo_calc[i]/c_total_calc[i];
                norm_c_total_calc[i,4] = (c_total_calc[i]-c_brdhi_kilo_calc[i]-c_brdhi_kihi_calc[i]-c_brdlo_kilo_calc[i])/c_total_calc[i];  
    
                clonal[i] = (lamdas[i,1]+ lamdas[i,2])/1000;            
            }

}
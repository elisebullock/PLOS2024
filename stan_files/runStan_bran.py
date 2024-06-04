#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 14:55:02 2022

@author: elise
"""

import matplotlib.pyplot as plt
import scipy.stats as sts
import numpy as np
import pandas as pd
import cmdstanpy ## import stan interface for Python
from scipy.special import logit,expit

import sys
sys.path.append("..")
    
np.random.seed(2101)

import paras
    
sm = cmdstanpy.CmdStanModel(stan_file=paras.stanloc+"branched.stan")

data_dict = {
    "numObs" : paras.numObs,
    "time_index" : paras.time_index,
    "kcounts" : (paras.counts_log),
    "brdhi_kihi_c" : paras.brdhi_kihi_c,
    "brdhi_kilo_c" : paras.brdhi_kilo_c,
    "brdlo_kilo_c" : paras.brdlo_kilo_c,
    "brdlo_kihi_c" : paras.brdlo_kihi_c,
    "switchtime" : paras.switchtime,
    "init_cond_obs" : paras.y0_est,
    "init_lamdas":[1000,0],
    "numofki67int" : paras.numofki67int,
    "numofbrduint" : paras.numofbrduint,
    "switch1" : (paras.numofbrduint+2)*(paras.numofki67int+2),
    "time_index_equilibrium" : paras.time_index_equilibrium,
    "counts_t0": np.int_(np.round_(np.mean(paras.dfslice.Qall))),
    "scalescale":paras.scalescale,
    "stdcounts":(np.std(np.log(paras.counts_log))),
    "fulltime": np.concatenate((np.logspace(-1,np.log10(paras.switchtime),num=int(100/2)), (np.logspace(-1,np.log10(paras.tmax-paras.switchtime),num=int(100/2)))+paras.switchtime)),
    "indexswitch":paras.indexswitch,
    
    "alpha_A_data": paras.alpha_A,
    "alpha_B_data": paras.alpha_B,
    "delta_A_data" : paras.delta_A,
    "delta_B_data" : paras.delta_B,
    "gamma_data" : paras.gammaB,
    "fs_data" : paras.fsB,
    "beta_data" : paras.beta,
    "Source_data" : paras.SourceB,
    "eff_data" : paras.eff,
    "mu_data" : paras.mu,

    "Source_hi":paras.SourceB_hi,
    "Source_lo":paras.SourceB_lo,
        
    }

init_d = {

    "Source" : paras.SourceB,
    "mu" : paras.mu,
    "alpha_A": paras.alpha_A,
    "alpha_B": paras.alpha_B,
    "delta_A" : paras.delta_A,
    "delta_B" : paras.delta_B,
    "fs" : paras.fsB,
    "beta" : paras.beta,
    "eff" : paras.eff,
    #"data_total": total
    }


cmdstanpy.write_stan_json(paras.filelocation+"branched.json", data_dict)

sam = sm.sample(
    data=paras.filelocation+"branched.json", 
    chains=5, 
    parallel_chains=5,
    refresh=5,
    inits=init_d,
    output_dir=paras.filelocation,
    show_progress=True,
    show_console=False,
    seed=74654,
    iter_warmup=100,
    iter_sampling=100,
    max_treedepth=20,
    adapt_delta=.95
   
)

import pickle
with open(paras.filelocation+"branched.pkl", 'wb') as f:
    pickle.dump(sam, f)

#--------------------------------------------------------------------------------#
#- This script replicate the experiments of Graham (2016ECMA)                   -#
#--------------------------------------------------------------------------------#

# Xin Zheng
# March 2017

from __future__ import division

import sys
sys.path.append('/Users/xin/')

import netrics as netrics
import time
import multiprocessing as mp
import numpy as np
import scipy as sp
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------------#
#- Define main simulation                                                       -#
#--------------------------------------------------------------------------------#
def runSimulationDesign(m):

    # Step 0: set the random state for reproducibility

    randst = np.random.RandomState(m*361)

    K          = MonteCarloDesign[m][0]
    Ashape1    = MonteCarloDesign[m][1]
    Ashape2    = MonteCarloDesign[m][2]
    SuppLength = MonteCarloDesign[m][3]
    AMean0     = MonteCarloDesign[m][4]
    AMean1     = MonteCarloDesign[m][5]
    beta       = MonteCarloDesign[m][6]

    SimulationResults_TL  = np.zeros((B,6))
    SimulationRssults_JFE = np.zeros((B,5))
    SimulationResults_BC  = np.zeros((B,5))
    NetworkProperties     = np.zeros((B,5))

    for b in xrange(0,B):

        # Step 1: Simulate network

        W_nnk = []
        W     = []
        x_bar = np.zeros((N,1))

        for k in range(0,K):

            X = 2 * (randst.binomial(1, 1/2, (N,1)) - 1/2)

            x_bar += X/K
            W_k   =  X * X.T - np.eye(N)
            W_nnk.append(W_k)

            W_k   =  W_k[ij_LowTri].reshape((-1,1))
            W.append(W_k)


        W = np.column_stack(W)

        A_i          = (AMean0 + AMean1)/2 + (AMean1 - (AMean0 + AMean1)/2) * X_bar \
                 + SuppLength*(randst.beta(AShape1,AShape2, (N,1)) - AShape1/(AShape1+AShape2))
        A            = A_i + A_i.T - 2*np.diag(np.ravel(A_i))
        A            = A[ij_LowTri].reshape((-1,1))
        p            = np.exp(np.dot(W,np.ones((K,1))*beta) + A) / (1 + np.exp(np.dot(W,np.ones((K,1))*beta) + A))
        D            = np.zeros((N,N), dtype='int8')
        D[ij_LowTri] = np.ravel(randst.uniform(0, 1, (n,1)) <= p)
        D            = D + D.T

        del X, X_bar, W_k, A, p

        G = nx.Graph(D)
        deg_seq = nx.degree(G).values()
        NetworkProperties[b,0] = nx.density(G)
        NetworkProperties[b,1] = nx.transitivity(G)
        NetworkProperties[b,2] = np.mean(deg_seq)
        NetworkProperties[b,3] = np.std(deg_seq)
        NetworkProperties[b,4] = len(max(nx.connected_components(G), key=len))/N

        del G, deg_seq

        # Step 3: Compute tetrad logit estimate of beta

        try:

            [beta_TL, vcov_beta_TL, terad_frac_TL, success_TL] = netrics.tetrad_logit(D, W_nk, dtcon=fixed_dtcon, \
W_names=None)

            SimulationResults_TL[b,0] = tetrad_frac_TL
            SimulationResults_TL[b,1] = True
            SimulationResults_TL[b,2] = beta_TL[0]
            SimulationResults_TL[b,3] = np.sqrt(vcov_beta_TL[0,0])
            SimulationResults_TL[b,4] = (beta_TL[0] - 1.96*np.sqrt(vcov_beta_TL[0,0]) <= beta <= beta_TL[0] \
                                         +1.96*np.sqrt(vcov_beta_TL[0,0]))
            SimulationResults_TL[b,5] = (beta_TL[0] - 1.645*np.sqrt(vcov_beta_TL[0,0]) <= beta <= beta_TL[0] \
                                         +1.645*np.sqrt(vcov_beta_TL[0,0]))
        except Exception, e;

            SimulationResults_TL[b,1] = False

        # Step 4: Compute JFE/BC logit esitmate of beta

        try:

            [beta_JFE, beta_JFE_BC, vcov_beta_JFE, A_JFE, success_JFE] = netrics.dyad_jfe_logit(D, W_nnk, T, \
                                                                                                silent=True, W_names=None,\
                                                                                                beta_sv=beta*np.ones((K,)))

            SimulationResults_JFE[b,0] = True
            SimulationResults_JFE[b,1] = beta_JFE[0]
            SimulationResults_BC[b,0]  = beta_JFE_BC[0]
            SimulationResults_JFE[b,2] = np.sqrt(vcov_beta_JFE[0,0])
            SimulationResults_BC[b,1]  = np.sqrt(vcov_beta_JFE[0,0])
            SimulationResults_JFE[b,3] = (beta_JFE[0] - 1.96*np.sqrt(vcov_beta_JFE[0,0]) <= beta <= beta_JFE[0] \
                                          +1.96*np.sqrt(vcov_beta_JFE[0,0]))
            SimulationResults_BC[b,2] = 
                                                

#!/usr/bin/env python
# coding: utf-8

# In[2590]:


import numpy as np
import jax
import jax.numpy as jnp
import jax.random as jxr
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
from functools import partial
from tqdm import tqdm

import sys
sys.path.append('../')
import utils
import inference
from models import CLDS, WeightSpaceGaussianProcess, ParamsCLDS

import pandas as pd
import pynapple as nap
import requests, math, os
import tqdm
import sklearn
import seaborn as sns

from scipy.io import savemat
from scipy.io import loadmat


# In[2591]:
whichDate   = '260303'
howManyBhv  = 1
for howManyC in range(6):
    for whichSeed in range(10):
        mat_data       = loadmat(f'C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code/TMLROut/forCLDS_saveFMT_260303_seed{whichSeed+1:02d}_TMLR_b{howManyBhv}_c{howManyC+1}.mat')
        filenameToSave = f'{whichDate}_seed{whichSeed+1:02d}_TMLR_b{howManyBhv}_c{howManyC+1}.mat'
        # mat_data      = loadmat('saveFMT_forCLDS.mat')
        print(f'\n{howManyC}\n')
        print(f'\n{whichSeed}\n')

        # In[2592]:


        dFF           = mat_data['dFF3d'] # np.array(mat_data['dFF'])
        behavior      = mat_data['bhv3d']# np.array(mat_data['behaviorData'])
        gtstates      = mat_data['gts3d']# np.array(mat_data['groundTruthStates'])
        trueFC        = mat_data['fcMat']


        # In[2593]:


        type(gtstates)


        # In[2594]:


        trueFC.shape


        # In[2595]:


        gtstates


        # In[2596]:


        print(f"Is arr NaN? {np.isnan(gtstates)}")


        # In[2597]:


        partition = int(0.8 * dFF.shape[2])

        dFF1  = np.swapaxes(dFF,0,2)
        dFF2  = np.swapaxes(dFF1,1,2)

        gt1   = np.swapaxes(gtstates,0,2)
        gt2   = np.swapaxes(gt1,1,2)

        bhv1  = np.swapaxes(behavior,0,2)
        bhv2  = np.swapaxes(bhv1,1,2)

        X     = jnp.array(dFF2) #jnp.array(gt2)
        Y     = jnp.array(dFF2)
        U     = jnp.array(bhv2)

        # U = U[:, 1, :]  # Use only one behavior readout

        # X_train, Y_train = X[:,:,:partition], Y[:,:,:partition]
        # X_test, Y_test = X[:,:,partition:], Y[:,:,partition:]

        # U_train = U[:,:partition]
        # U_test  = U[:,partition:]

        U = U[:, :, 1]  # Use only one behavior readout

        X_train, Y_train = X[:partition], Y[:partition]
        X_test, Y_test = X[partition:], Y[partition:]

        U_train, U_test = U[:partition], U[partition:]


        # In[2598]:


        dFF1  = np.swapaxes(dFF,0,2)
        dFF2  = np.swapaxes(dFF1,1,2)


        # In[2599]:


        dFF2.shape


        # In[2600]:


        U[0,:]


        # In[2601]:


        U.shape


        # In[2602]:


        plt.plot(U[0,:])
        # plt.show()


        # In[2603]:


        U.shape


        # In[2604]:


        X.shape


        # In[2605]:


        Y.shape


        # In[2606]:


        X_train.shape


        # In[2607]:


        U_train.shape


        # In[2608]:


        plt.plot(X_train[0,:,:])
        # plt.show()


        # In[2609]:


        plt.plot(Y_train[0,:,:])
        # plt.show()


        # In[2610]:


        print(f"Is arr NaN? {np.isnan(U)}")


        # In[2611]:


        trial_time_length = 3000


        # # Fit model

        # In[2612]:


        # Define model

        latent_dim = 8
        n_neurons = Y.shape[2]

        # sigma = 1 too high, results in NaNs
        _sigma, _kappa, c_period = 0.75, 0.15, 2*jnp.pi #1.0, 0.15, 2*jnp.pi
        t_period = trial_time_length + 6 * _kappa
        T2_basis_funcs = utils.Tm_basis(5, M_conditions=1, sigma=_sigma, kappa=_kappa, period=jnp.array([t_period, c_period]))
        T1_basis_funcs = utils.T1_basis(5, _sigma, _kappa, c_period)

        A_prior = WeightSpaceGaussianProcess(T1_basis_funcs, D1=latent_dim, D2=latent_dim)
        b_prior = WeightSpaceGaussianProcess(T1_basis_funcs, D1=latent_dim, D2=1)
        m0_prior = WeightSpaceGaussianProcess(T1_basis_funcs, D1=latent_dim, D2=1)
        C_prior = WeightSpaceGaussianProcess(T1_basis_funcs, D1=n_neurons, D2=latent_dim)
        model = CLDS(
            wgps={
                'A': A_prior, 
                'b': b_prior,
                'C': C_prior,
                'm0': m0_prior,
                }, 
            state_dim=latent_dim, 
            emission_dim=n_neurons,
            )



        # In[2613]:


        t_period


        # In[2614]:


        n_neurons


        # In[2615]:


        fig, axs = plt.subplots(figsize=[24,6], ncols=1,nrows=6) #2)
        for i in range(3):
            print(i)
            axs[i*2 + 0].plot(A_prior.sample(jxr.PRNGKey(3), U[i])[:,0,0], c='tab:blue');
            # axs[i*3 + 1].plot(U[i,:], A_prior.sample(jxr.PRNGKey(3), U[i])[:,0,1], c='tab:orange')
            axs[i*2 + 1].plot(U[i], c='tab:green')
            # axs[1].plot(U[i][:,1], A_prior.sample(jxr.PRNGKey(3), U[i])[:,0,1], '.', c='tab:orange')


        # In[2616]:


        # # Initialize A and b weights as ring attractor

        # def ring_weights(kappa, sigma):
        #     def weight_space_coefficients(m):
        #         return jnp.sqrt(utils.squared_exponential_spectral_measure(m, sigma, kappa))

        #     # A weights

        #     ring_wA_weights = jnp.zeros((len(A_prior.basis_funcs), latent_dim, latent_dim))
        #     # A[0,0] = 1/2 - cos(2 t)/2
        #     ring_wA_weights = ring_wA_weights.at[11,0,0].set(1/weight_space_coefficients(0) * 1/2)
        #     ring_wA_weights = ring_wA_weights.at[15,0,0].set(-1/weight_space_coefficients(2) * 1/2)

        #     # A[1,1] = 1/2 + cos(2 t)/2
        #     ring_wA_weights = ring_wA_weights.at[11,1,1].set(1/weight_space_coefficients(0) * 1/2)
        #     ring_wA_weights = ring_wA_weights.at[15,1,1].set(1/weight_space_coefficients(2) * 1/2)

        #     # A[0,1] = -sin(2 t)/2
        #     ring_wA_weights = ring_wA_weights.at[14,1,0].set(-1/weight_space_coefficients(2) * 1/2)
        #     ring_wA_weights = ring_wA_weights.at[14,0,1].set(-1/weight_space_coefficients(2) * 1/2)

        #     # b weights

        #     ring_wb_weights = jnp.zeros((len(b_prior.basis_funcs), 2, 1))

        #     ring_wb_weights = jnp.zeros((len(b_prior.basis_funcs), b_prior.D1, b_prior.D2)) # len_basis D1 D2
        #     ring_wb_weights = ring_wb_weights.at[12,0,0].set(1/weight_space_coefficients(1))
        #     ring_wb_weights = ring_wb_weights.at[13,1,0].set(1/weight_space_coefficients(1))

        #     return ring_wA_weights, ring_wb_weights


        # In[2617]:


        num_timesteps = Y.shape[1]

        # ring_wA_weights, ring_wb_weights = ring_weights(_kappa, _sigma)

        # Initialize
        seed = 2 # use 0 for plotting
        A_key, b_key, C_key, m0_key = jxr.split(jxr.PRNGKey(seed), 4)
        initial_params = ParamsCLDS(
            dynamics_gp_weights =  A_prior.sample_weights(A_key),# ring_wA_weights, 
            Q = jnp.eye(latent_dim),
            R = jnp.eye(n_neurons),
            m0 = jnp.zeros(latent_dim),
            m0_gp_weights = m0_prior.sample_weights(m0_key),
            S0 = jnp.eye(latent_dim),
            emissions_gp_weights = C_prior.sample_weights(C_key),
            bias_gp_weights = b_prior.sample_weights(b_key),
            Cs = jnp.tile(jxr.normal(C_key, (n_neurons, latent_dim)), (num_timesteps, 1, 1)), #C_prior.sample(C_key, U[0]), #jnp.tile(jxr.normal(C_key, (n_neurons, latent_dim)), (num_timesteps, 1, 1)), #C_prior.sample(C_key, U[0]),
            bs = None, #jnp.zeros((num_timesteps-1, latent_dim)),# None, #b_prior.sample(b_key, U[0]).squeeze(), #jnp.zeros((num_timesteps-1, latent_dim)), #b_prior.sample(b_key, conditions).squeeze(),
        )


        # In[2618]:


        num_timesteps


        # In[2619]:


        # Fit model
        params, log_probs = inference.fit_em(model, initial_params, emissions=Y_train, conditions=U_train, num_iters=200)
        test_ll = model.marginal_log_lik(params, emissions=Y_test, conditions=U_test)
        print(f"Test log likelihood: {test_ll}")

        # Show results
        fig, ax = plt.subplots(figsize=[4,3])
        ax.plot(log_probs)
        ax.set_ylabel('log prob')
        ax.set_xlabel('Epochs');


        # In[2620]:


        fig, axs = plt.subplots(ncols=latent_dim, nrows=latent_dim, figsize=(10, 8), constrained_layout=True);

        for i in range(8):
            for j in range(8):
                ax = axs[i,j]
                # true_As, true_bs, true_Cs = dynamics(theta[0], omega, switchtimeTL[0], switchtimeBR[0])
                ax.plot(trueFC[i,j,:,0], 'k--', label='True')
                ax.plot(A_prior(initial_params.dynamics_gp_weights, U[100][:-1])[:,i,j], label='init')
                ax.plot(A_prior(params.dynamics_gp_weights, U[49][:-1])[:,i,j], label='EM')
                # print(A_prior(params.dynamics_gp_weights, theta[0][:-1]))
                # print(A_prior(params.dynamics_gp_weights, theta[0][:-1]).shape)
                # # ax.set_title(f'A[{i},{j}]')
        axs[3,0].set_xlabel('Time')
        axs[3,3].legend()

        for ax in axs.flatten():
            sns.despine(ax=ax, trim=True, offset=True)
        for ax in axs[0]:
            ax.set_xticklabels([])

        fig.suptitle('Dynamics matrix $A$');


        # dictSimulation = {'true_As':true_As, 'true_bs':true_bs, 'true_Cs':true_Cs}
        # savemat('dictSimulation.mat',dictSimulation)


        # In[2621]:


        U.shape


        # In[2622]:


        U[200]


        # In[2623]:


        U[200].shape


        # In[2624]:


        (U.shape)[0]


        # In[2625]:


        cldsAMats   = np.zeros((8,8,(U.shape)[1],(U.shape)[0]))
        # permutOrder = [1, 2, 0]
        for i in range((U.shape)[0]):
            thisA = A_prior(params.dynamics_gp_weights, U[i][:])
            cldsAMats[:,:,:,i] = np.transpose(thisA,(1, 2, 0)) #[thisA[j] for j in permutOrder]

        dictCLDS   = {'trueFC': trueFC[:,:,:,:],
                    'clds': cldsAMats
                    }
        savemat(filenameToSave,dictCLDS)
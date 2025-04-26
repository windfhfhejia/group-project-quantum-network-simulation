#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 23:33:41 2025

@author: chenyj
"""



"""
Improved QKD Network Implementation with Wavelength Division Multiplexing
"""

import SingleLink as sl
import numpy as np
from typing import Tuple, List, Dict
from scipy.optimize import differential_evolution
import logging
import sys
import random
import csv
import os
import matplotlib.pyplot as plt

class QKDNetwork:
    def __init__(self, num_users: int = 6):
        self.num_users = num_users
        
        
        self.wavelength_matrix1 = np.array([
            [5, -5, -3, 5, -5, -3],
            [4, 3, -4, 4, 3, -4],
            [2, 1, 1, -2, -1, -1]])
        
        self.wavelength_matrix2 = np.array(
           [[6,-6,-4,6,-6,-4],
            [5,4,-5, 5,4,-5],
            [3,2,1,-3,-2,-1]])
        
        
        self.wavelength_matrix3=np.array(
            [[1,-1,6,-6,1,-1],
             [2,4,-2,-5,2,4],
             [3,5,-4,-3,3,5],
             [7,7,8,8,-7,-7]])

        self.wavelength_matrix4=np.array(
            [[1,-1,-2,-3,-4,-5],
             [2,6,-6,-7,-8,-9],
             [3,7,10,-10,-11,-12],
             [4,8,11,13,-13,-14],
             [5,9,12,14,15,-15]])
            
        
        #self.current_wavelength_allocation = self.wavelength_matrix1
        
        self.performance_metrics = {
            'average_key_rate': 0,
            'min_key_rate': 0,
            'channel_efficiency': 0,
        }
    
    #check if the wavelength_allocation is valid
    def validate_wavelength_allocation(self,wavelength_allocation: np.ndarray)-> np.ndarray: 
        sender=[]
        receiver=[]
        num_cols = wavelength_allocation.shape[1]  
        sender, receiver = np.triu_indices(num_cols, k=1)  # Get upper triangle indices

        for i in range(len(sender)):
            chan1 = wavelength_allocation[:, sender[i]].tolist()  
            chan2 = wavelength_allocation[:, receiver[i]].tolist()
            #check if wavelenth for all users are unique
            if len(chan1)!=len(np.unique(chan1)) or len(chan2)!=len(np.unique(chan2)):
                sys.exit("Error: Invalid wavelength allocation detected.")
        #check if all users have channels to use
        if not any(-item in chan2 for item in chan1):
            sys.exit("Error: Invalid wavelength allocation detected.")
        else:
            return wavelength_allocation
    
    def generate_fiber_length(self, wa:int)->np.ndarray:
        if wa==0:
            #fiber_length=[random.randint(1, 30) for _ in range(6)]
            fiber_length=[25,19,4,25,24,23]
            self.fb_variable=fiber_length
        if wa>0:
            fiber_length=self.fb_variable
        #fiber_length=[26, 30, 2, 27, 22, 8]
        #print(fiber_length)
        return fiber_length
        
        
        
    def generate_random_coincidence_window_matrix(self)-> np.ndarray: 
        coincidence_matrix=np.zeros((self.num_users, self.num_users), dtype=int)
        for row in range(self.num_users):
            for col in range(self.num_users):
                if row==col:
                    coincidence_matrix[row][col]=0
                else:
                    coincidence_matrix[row][col]=random.randint(40,800)
                    coincidence_matrix[col][row]=coincidence_matrix[row][col]
        return coincidence_matrix
                
                

    '''
    (1)assume no counts from other users and detector efficiency is 70%, total pair 
    generation rate=10^6
    (2)input: users: order of this user in the sequence E.g. Alice->0, Chloe->2
    (3)output: ideal_single_link: counts at the link correspond to the user assume counts
    from other user=0
    pgr_per_channel: pair generation rate per cahnnel corresponds to wavelength allocation
    '''

    def ideal_single_counts(self, wavelength_allocation: np.ndarray, user: int,\
                          fiber_length: np.ndarray, pair_generation_rate:int, bs:int) -> Tuple[int, int]:
 
        ideal_single_link = np.zeros(self.num_users, dtype=int)
        pgr_per_channel = (pair_generation_rate/bs) / np.max(wavelength_allocation)
        
        singlecounts=sl.Calc_QKD_Parameters_For_Single_Link\
        (PairGenerationRate_per_second=pgr_per_channel, FibreLength_in_km = \
        [fiber_length[user],0],SingleCountsFromOtherUsers_per_second=[0,0,0,0])
        
        ideal_single_link = singlecounts['Singles rate for Alice'].HA / 0.7
        return ideal_single_link, pgr_per_channel
    
    
    
    '''
    (1)This function calculates the number of useless channel (noise) and 
    effective channelfor the specific node 
    (2) input: row and col: sender and receiver 
    (3)output: channels_at_nodes: availble channels at both nodes for the effective link
               effective_channel: which channels can be used for effective communication
               num_eff_chan: number of effective channel
               num_noise_chan: number of noise channel at the sender node
    '''

    def num_channel_from_other_users(self, wavelength_allocation: np.ndarray,
                                   row: int, col: int) -> Tuple[List, List, int, int]:
        channels_at_nodes = []
        effective_channel = []
        count_dict = {}
        
        
        num_full_eff = 0
        num_half_eff = 0
        num_full_noise = 0
        num_half_noise = 0

            
        for x in range(len(wavelength_allocation)):
            for y in range(len(wavelength_allocation[0])):
                if y == row:
                    channels_at_nodes.append(wavelength_allocation[x][y])

        for x in range(len(wavelength_allocation)):
            for y in range(len(wavelength_allocation[0])):
                 if y == col:
                     channels_at_nodes.append(wavelength_allocation[x][y])
        
        for num in channels_at_nodes:
            if -num in channels_at_nodes and num>0:
                effective_channel.append(num)
        
        #record the number of wavelength occured at the wavelength assignment    
        for x in wavelength_allocation:
            for y in x:
                count_dict[y]=count_dict.get(y, 0) + 1       
        for element in effective_channel:
            if count_dict.get(element, 0)==1:
                num_full_eff += 1
            if count_dict.get(element, 0)==2:
                num_half_eff += 1
        for element in channels_at_nodes[0:len(wavelength_allocation)]:
            if count_dict.get(-element, 0)==1:
                num_full_noise += 1
            if count_dict.get(-element, 0)==2:
                num_half_noise += count_dict.get(element, 0)
        num_full_noise=num_full_noise-num_full_eff
        num_half_noise=num_half_noise-num_half_eff
        
        if num_full_noise+num_half_noise>4:
            num_full_noise=num_full_noise + 1
            num_half_noise=num_half_noise-2
        
        return num_full_eff, num_half_eff, num_full_noise, num_half_noise
    
    
    
    '''
    calculate the total counts at the nodes from other users(noise) given with effective link
    '''
    def pgr_noise_counts(self, wavelength_allocation: np.ndarray, fiber_length: np.ndarray, 
                          pair_generation_rate:int, row: int, col: int) -> Tuple[int, int]:
     
        num_full_eff, num_half_eff, num_full_noise, num_half_noise = self.num_channel_from_other_users(
            wavelength_allocation, row, col
            )
        single_counts_full, pgr_full = self.ideal_single_counts(
            wavelength_allocation, row, fiber_length, pair_generation_rate, 1
            )
        single_counts_half, pgr_half = self.ideal_single_counts(
            wavelength_allocation, row, fiber_length, pair_generation_rate, 2
            )
        link_pgr = pgr_full * num_full_eff + pgr_half * num_half_eff
        node_noise = single_counts_full * num_full_noise + single_counts_half * num_half_noise

        return link_pgr, node_noise


    '''
    calculate the new secure key rate of each link base on the counts from other user
    and the number of effective channel 
    '''
    def multi_secure_key_rate(self, wavelength_allocation: np.ndarray,fiber_length:np.ndarray,
                 num_users: int, coincidence_matrix: np.ndarray, pair_generation_rate: int) -> Tuple[List, List]:
        
        secure_key_rate_matrix = np.zeros((self.num_users, self.num_users), dtype=int)
        QBER_matrix=np.zeros((self.num_users, self.num_users), dtype=float)
        HE_matrix=np.zeros((self.num_users, self.num_users), dtype=float)
        for row in range(self.num_users):
            for col in range(self.num_users):
                if row < col:
                        
                    new_pgr, node_noise1=self.pgr_noise_counts(wavelength_allocation, 
                                                fiber_length, pair_generation_rate, row, col)
                    _, node_noise2=self.pgr_noise_counts(wavelength_allocation, 
                                                fiber_length, pair_generation_rate, col, row)
                    new_single_link = sl.Calc_QKD_Parameters_For_Single_Link(
                        PairGenerationRate_per_second=new_pgr,
                        CoincidenceWindow_in_ps=coincidence_matrix[row][col],
                        FibreLength_in_km=[fiber_length[row], fiber_length[col]],
                        SingleCountsFromOtherUsers_per_second=[
                            node_noise1, node_noise1,
                            node_noise2, node_noise2
                            ]
                    )
                
                    secure_key_rate_matrix[row][col] = new_single_link['Total secure key rate']
                    #make sure all the channel can be used and minimum skr for the channel is assumed to be 10
                    if secure_key_rate_matrix[row][col]<3:
                        return 0, secure_key_rate_matrix, QBER_matrix 
                
                    QBER_matrix[row][col]=new_single_link['QBER between Alice and Bob']
                    HE_matrix[row][col]=new_pgr/(((new_single_link['Singles rate for Alice'].HA*2)*
                                                 (new_single_link['Singles rate for Bob'].HA*2))**0.5)
        sum_of_skr=np.sum(secure_key_rate_matrix)
        return sum_of_skr,secure_key_rate_matrix, QBER_matrix
    
    def list_to_matrix(self, params: np.ndarray) -> np.ndarray:
        matrix = np.zeros((self.num_users, self.num_users))
        idx = 0
        for i in range(self.num_users):
            for j in range(i + 1, self.num_users):
                matrix[i,j] = matrix[j,i] = params[idx]
                idx += 1
        return matrix
    
    def optimize_network_parameters(self, wavelength_allocation: np.ndarray, num_users:int, wa:int) -> Dict:
        

        num_channel = (self.num_users * (self.num_users - 1)) // 2 
        fiber_length=self.generate_fiber_length(wa)
        print('fiber length：',fiber_length)
        #set boundary for coincidence window and pair generation rate
        bounds = [(40, 800)] * num_channel + [(4e5, 1e7)]
        
        
        wrapped_objective = lambda x: -self.multi_secure_key_rate(wavelength_allocation,fiber_length, num_users,
                                                                  self.list_to_matrix(x[:-1]), x[-1] )[0]
        
        
        result = differential_evolution(
            wrapped_objective, 
            bounds,
            maxiter=100,
            popsize=30,
            updating='immediate',  
            workers=1,  
            disp=True  
        )
        
        
        if not result.success:
            logging.warning(f"Optimization did not converge: {result.message}")
        
        
        best_params = result.x
        n_cw = (self.num_users * (self.num_users - 1)) // 2
        best_cw_params = best_params[:n_cw]
        best_pgr = best_params[n_cw]
        
        
        best_cw = np.zeros((self.num_users, self.num_users))
        idx = 0
        for i in range(self.num_users):
            for j in range(i + 1, self.num_users):
                best_cw[i,j] = best_cw[j,i] = int(best_cw_params[idx])
                idx += 1


        
        sum_skr,skr_matrix,qber_matrix = self.multi_secure_key_rate(wavelength_allocation, fiber_length,
                     num_users, best_cw, best_pgr)
        average_qber=np.sum(qber_matrix)/(self.num_users*(self.num_users-1)/2)   
        print(skr_matrix)
        #print(Heralding_matrix)
        
        csv_store = [
            fiber_length,  
            int(sum_skr),  
            round(float(average_qber), 4),  
            [int(param) for param in best_params[:-1]],  
            round(best_params[-1] / 1000) * 1000, 
            int(result.nit),  
            bool(result.success)  
            ]
        with open("output3.csv", mode="a", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(csv_store)
        print("New data added to CSV file")
        
        return {
            'coincidence_window': best_cw,
            'pair_generation_rate': best_pgr,
            'sum of secure key rate': sum_skr,
            'secure key rate matrix': skr_matrix,
            'average qber': average_qber,
            'QBER matrix': qber_matrix,
            'optimization_success': result.success,
            'optimization_message': result.message,
            'number_of_iterations': result.nit,
            'function_evaluations': result.nfev
        }
    
    def plot_results(self, results: Dict):
        
        plt.figure(figsize=(15, 10))
        
        # Plot SKR matrix
        plt.subplot(2, 2, 1)
        plt.imshow(results['secure key rate matrix'])
        plt.colorbar(label='Secure Key Rate')
        plt.title('Secure Key Rate Matrix')
        
        # Plot QBER matrix
        plt.subplot(2, 2, 2)
        plt.imshow(results['QBER matrix'])
        plt.colorbar(label='QBER')
        plt.title('QBER Matrix')
        
        # Plot coincidence window matrix
        plt.subplot(2, 2, 3)
        plt.imshow(results['coincidence_window'])
        plt.colorbar(label='Coincidence Window')
        plt.title('Coincidence Window Matrix')
        
        plt.tight_layout()
        plt.show()
        
    def check_and_create_csv(self,filename='output4.csv'):
        data = ["fiber_length", "total secure key rate","average QBER", "coincidence window", "pair generation rate", "number of iterations","optimization_success"]
        if not os.path.exists(filename):
            with open(filename, mode='w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(data)
                    
    



def main():

    network = QKDNetwork(num_users=6)
    

    #x=sl.Calc_QKD_Parameters_For_Single_Link(FibreLength_in_km=[1.3,3],SingleCountsFromOtherUsers_per_second=[1803501,1803501,1811788,1811788])
    #single_counts, pgr_per_channel=network.ideal_single_counts(network.wavelength_matrix1, 2,network.fiber_length)
    '''
    for m in range(6):
        for n in range(6):
            if m<n:
                a,b,c,d=network.num_channel_from_other_users(network.wavelength_matrix2,m,n)
                print(m,n,a,b,c,d)
    '''
    #node_others=network.node_noise_counts(network.wavelength_matrix2)
    #test,_=network.multi_secure_key_rate(network.wavelength_matrix1, network.fiber_length, 6, 150, 1E6)
    #validated_allocation=network.validate_wavelength_allocation(network.wavelength_matrix2)
    #network.check_and_create_csv('output1.csv')
    
    
    wavelength_allocation_variable=[network.wavelength_matrix1, network.wavelength_matrix2, network.wavelength_matrix4]
    for m in range(3):
        i=m%3
        network.check_and_create_csv('output4.csv')
        wavelength_allocation = network.validate_wavelength_allocation(
            wavelength_allocation_variable[i]
        )
        
        
        logging.info("Starting network optimization...")
        results = network.optimize_network_parameters(wavelength_allocation, network.num_users,i)
        
        
        logging.info("\nOptimization Results:")
        logging.info(f"Best pair generation rate: {results['pair_generation_rate']:.2e}")
        logging.info(f"secure key rate matrix: \n {results['secure key rate matrix']}")
        logging.info(f"total secure key rate: {results['sum of secure key rate']:.2f}")
        logging.info(f"Average QBER: {results['average qber']:.4f}")
        logging.info(f"Number of iterations: {results['number_of_iterations']}")
        
        
        #network.save_results(results)
        
        
        network.plot_results(results)
    
    
    
    #print('total counts per second from other users: \n', node_others)
    #print(a,b,c,d)
    #print(test)
    #print('single counts per seconds for specific users: \n',single_counts)
    #print(x)#test for the single link file
    
   
    '''
    new_allocation = network.dynamic_wavelength_allocation()
    if not np.array_equal(new_allocation, network.current_wavelength_allocation):
        print('New wavelength allocation suggested:\n', new_allocation)
    '''
    
if __name__ == "__main__":
    main()
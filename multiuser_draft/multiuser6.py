#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 23:56:47 2025

@author: chenyj
"""

"""
Improved QKD Network Implementation with Wavelength Division Multiplexing
"""

import SingleLink as sl
import numpy as np
from typing import Tuple, List, Dict
import sys
import random
import csv
import os

class QKDNetwork:
    def __init__(self, num_users: int = 6):
        self.num_users = num_users
        
        
        self.wavelength_matrix1 = np.array([
            [4, -4, -2, 4, -4, -2],
            [3, 2, -3, 3, 2, -3],
            [1, 1, 1, -1, -1, -1]
        ])
        
        self.wavelength_matrix2 = np.array([
            [5, -5, -3, 5, -5, -3],
            [4, 3, -4, 4, 3, -4],
            [2, 1, 1, -2, -1, -1]
        ])
        
        self.wavelength_matrix3=np.array(
            [[5,-5,5,5,5,5],
             [4,4,-4,4,4,4],        
             [3,3,3,-3,3,3],
             [2,2,2,2,-2,2],
             [1,1,1,1,1,-1]])

        self.wavelength_matrix4=np.array(
            [[6,-6,-4,6,-6,-4],
             [5,4,-5, 5,4,-5],
             [3,2,1,-3,-2,-1]])
        
        self.current_wavelength_allocation = self.wavelength_matrix1
        
        self.performance_metrics = {
            'average_key_rate': 0,
            'min_key_rate': 0,
            'channel_efficiency': 0,
        }
    
    #check if the wavelength_allocation is valid
    def validate_wavelength_allocation(self,wavelength_allocation)-> np.ndarray: 
        sender=[]
        receiver=[]
        num_cols = wavelength_allocation.shape[1]  
        sender, receiver = np.triu_indices(num_cols, k=1)  # Get upper triangle indices

        for i in range(len(sender)):
            chan1 = wavelength_allocation[:, sender[i]].tolist()  
            chan2 = wavelength_allocation[:, receiver[i]].tolist()  

        if not any(-item in chan2 for item in chan1):
            sys.exit("Error: Invalid wavelength allocation detected.")
        else:
            return wavelength_allocation
    
    def generate_fiber_length(self)->np.ndarray:
        #fiber_length=[random.randint(1, 30) for _ in range(6)]
        fiber_length=[8, 23, 25, 4, 8, 28]
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
                          fiber_length: np.ndarray, pair_generation_rate:int) -> Tuple[int, float]:
 
        ideal_single_link = np.zeros(self.num_users, dtype=int)
        pgr_per_channel = pair_generation_rate / np.max(wavelength_allocation)
        
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
        
        for x in range(len(wavelength_allocation)):
            for y in range(len(wavelength_allocation[0])):
                if y == row or y == col:
                    channels_at_nodes.append(wavelength_allocation[x][y])
        
        num_eff_chan = 0
        if row > col:
            for num in channels_at_nodes[0::2]:
                if -num in channels_at_nodes:
                    effective_channel.append(num)    
                    num_eff_chan += 1
        else:
            for num in channels_at_nodes[1::2]:
                if -num in channels_at_nodes:
                    effective_channel.append(num)    
                    num_eff_chan += 1
        
        num_noise_chan = 0
        for x in range(len(wavelength_allocation)):
            for y in range(len(wavelength_allocation[0])):
                if wavelength_allocation[x][y] * (-1) in effective_channel:
                    num_noise_chan += 1
        
        num_noise_chan = num_noise_chan - len(effective_channel)
        return channels_at_nodes, effective_channel, num_eff_chan, num_noise_chan
    
    
    
    '''
    calculate the total counts at the nodes from other users(noise) given with effective link
    '''
    def node_noise_counts(self, wavelength_allocation: np.ndarray, fiber_length: np.ndarray, 
                          pair_generation_rate:int) -> np.ndarray:
     
        node_others = np.zeros((self.num_users, self.num_users), dtype=int)
        
        for row in range(self.num_users):
            for col in range(self.num_users):
                single_counts, _ = self.ideal_single_counts(
                    wavelength_allocation, row, fiber_length, pair_generation_rate
                )
                _, _, _, num_noise_chan = self.num_channel_from_other_users(
                    wavelength_allocation, row, col
                )
                node_others[row][col] = single_counts * num_noise_chan
                
        return node_others


    '''
    calculate the new secure key rate of each link base on the counts from other user
    and the number of effective channel 
    '''
    def multi_secure_key_rate(self, wavelength_allocation: np.ndarray,fiber_length:np.ndarray,
                 num_users: int, coincidence_matrix: np.ndarray, pair_generation_rate: int) -> Tuple[List, List]:
        
        secure_key_rate_matrix = np.zeros((self.num_users, self.num_users), dtype=int)
        QBER_matrix=np.zeros((self.num_users, self.num_users), dtype=float)
        node_others = self.node_noise_counts(wavelength_allocation, fiber_length, pair_generation_rate)
        
        for row in range(self.num_users):
            for col in range(self.num_users):
                if row == col:
                    continue
                    
                _, pgr_per_channel = self.ideal_single_counts(
                    wavelength_allocation, row, fiber_length, pair_generation_rate
                )
                _, _, num_eff_chan, _ = self.num_channel_from_other_users(
                    wavelength_allocation, row, col
                )
                
                new_single_link = sl.Calc_QKD_Parameters_For_Single_Link(
                    PairGenerationRate_per_second=pgr_per_channel * num_eff_chan,
                    CoincidenceWindow_in_ps=coincidence_matrix[row][col],
                    FibreLength_in_km=[fiber_length[row], fiber_length[col]],
                    SingleCountsFromOtherUsers_per_second=[
                        node_others[row][col], node_others[row][col],
                        node_others[col][row], node_others[col][row]
                    ]
                )
                
                secure_key_rate_matrix[row][col] = new_single_link['Total secure key rate']
                QBER_matrix[row][col]=new_single_link['QBER between Alice and Bob']
                
        return secure_key_rate_matrix, QBER_matrix
    
    
    
    '''
    (1)find the best combination of coincidence window and pair generation rate to achieve the largest total 
    secure key rate
    (2)cd_break: break the cd for loop when coincidence window increases but total secure 
    rate decrease
    (3)pgr_break: break the 'factor' for loop when any link has 0 secure key rate (When
     QBER exceed security threshold)
    (4)optimized_parameter:1.sum of total secure key rate 2. average QBER 3. coincidence window 
       4. pair generation rate
    '''
    
    def optimize_coincidence_window(self, wavelength_allocation: np.ndarray, 
                                    num_users:int) -> Tuple[np.ndarray, np.ndarray, list]:
        
        optimized_secure_key_rate_matrix = np.zeros((self.num_users, self.num_users), dtype=int)
        optimized_QBER_matrix=np.zeros((self.num_users, self.num_users), dtype=float)
        csv_store=np.zeros((3,5), dtype=list)
        optimized_parameter=[0,0,0,0];
        fiber_length=self.generate_fiber_length()
        total_sample=1000
        for num_sample in range(1,total_sample):
            cw_matrix=self.generate_random_coincidence_window_matrix()
            pgr_range1 = np.arange(4e5, 1e6, 1e5)
            pgr_range2= np.arange(1e6, 1e7, 1e6)
            pgr_range= list(pgr_range1) + list(pgr_range2)
            for pgr in pgr_range:
                sum_of_secure_key_rate=0;
                skr_matrix, qber_matrix=self.multi_secure_key_rate(wavelength_allocation, fiber_length, num_users, cw_matrix,pgr)
                for row in range(len(skr_matrix)):
                    for col in range(len(skr_matrix)):
                        if row!=col and skr_matrix[row][col]==0:
                            sum_of_secure_key_rate=0
                        else:
                            sum_of_secure_key_rate=sum_of_secure_key_rate+skr_matrix[row][col]
                if sum_of_secure_key_rate>optimized_parameter[0]:
                    average_qber=np.sum(qber_matrix)/(num_users*(num_users-1))
                    #average_skr=sum_of_secure_key_rate/(num_users*(num_users-1))
                    optimized_parameter=[sum_of_secure_key_rate,average_qber,cw_matrix,pgr]
                    optimized_secure_key_rate_matrix=skr_matrix
                    optimized_QBER_matrix=qber_matrix
                    optimized_coincidence_window=cw_matrix
                        
                    '''
                    store the top 3 secure key rate combination to csv with given fiber length
                    '''
                    csv_store[2]=csv_store[1]
                    csv_store[1]=csv_store[0]
                    csv_secure_key_rate=optimized_parameter[0]
                    cd_array=np.array(optimized_coincidence_window)
                    csv_coincidece_matrix=cd_array[np.tril_indices(cd_array.shape[0], k=-1)]
                    csv_pair_generation_rate=int(optimized_parameter[3])
                    csv_store[0]=[fiber_length, csv_secure_key_rate, csv_coincidece_matrix,csv_pair_generation_rate,total_sample ]
            
            for i in range(10):
                if num_sample==(i+1)*total_sample/10:
                    print(f'loading: {int((i + 1) * 10)}%')
        
        with open("output1.csv", mode="a", newline="") as file:
            writer = csv.writer(file)
            for i in range(3):
                writer.writerow(csv_store[i])
        print("New data added to CSV file")
        
        return optimized_secure_key_rate_matrix, optimized_QBER_matrix, optimized_parameter, optimized_coincidence_window, total_sample
        
        
    def check_and_create_csv(self,filename='output1.csv'):
        data = ["fiber_length", "total secure key rate", "coincidence window", "pair generation rate", "number of cd sample"]
        if not os.path.exists(filename):
            with open(filename, mode='w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(data)
                    
                    

    def dynamic_wavelength_allocation(self) -> np.ndarray:
       
        current_performance = self.evaluate_network_performance(
            self.multi_secure_key_rate(self.current_wavelength_allocation, 150)
        )
        

        test_allocations = [self.wavelength_matrix1, self.wavelength_matrix2]
        best_allocation = self.current_wavelength_allocation
        best_performance = current_performance['average_key_rate']
        
        for allocation in test_allocations:
            test_skr = self.multi_secure_key_rate(allocation, 150)
            test_performance = np.mean(test_skr[test_skr > 0])
            
            if test_performance > best_performance:
                best_performance = test_performance
                best_allocation = allocation
        
        return best_allocation

def main():

    network = QKDNetwork(num_users=6)
    

    #x=sl.Calc_QKD_Parameters_For_Single_Link(FibreLength_in_km=[1.3,3],SingleCountsFromOtherUsers_per_second=[1803501,1803501,1811788,1811788])
    #single_counts, pgr_per_channel=network.ideal_single_counts(network.wavelength_matrix1, 2,network.fiber_length)
    #a,b,c,d=network.num_channel_from_other_users(network.wavelength_matrix3,4,0)
    #node_others=network.node_noise_counts(network.wavelength_matrix2)
    #test,_=network.multi_secure_key_rate(network.wavelength_matrix1, network.fiber_length, 6, 150, 1E6)
    validated_allocation=network.validate_wavelength_allocation(network.wavelength_matrix2)
    network.check_and_create_csv('output1.csv')
    for i in range(1):
        opt_skr_matrix, opt_qber_matrix, opt_parameter,_,_=network.optimize_coincidence_window(validated_allocation, network.num_users)
        print(i)
        
    '''
    print('optimized coincidence window (ps) matrix:\n',opt_parameter[2])
    print('optimized pair generation rate:\n', "{:.0E}".format(opt_parameter[3]))
    print('optimized secure key rate for each links :\n', opt_skr_matrix)
    print('total secure key rate: \n', opt_parameter[0])
    print('optimized quantum BER for each links: \n', opt_qber_matrix)
    print('average quantum BER: \n', opt_parameter[1])
    
    #print('total counts per second from other users: \n', node_others)
    #print(a,b,c,d)
    #print(test)
    #print('single counts per seconds for specific users: \n',single_counts)
    #print(x)#test for the single link file
    
   
    ‘’‘
    new_allocation = network.dynamic_wavelength_allocation()
    if not np.array_equal(new_allocation, network.current_wavelength_allocation):
        print('New wavelength allocation suggested:\n', new_allocation)
    '''
    
if __name__ == "__main__":
    main()
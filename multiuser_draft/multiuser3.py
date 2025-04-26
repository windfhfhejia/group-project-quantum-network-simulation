#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 00:37:32 2025

@author: chenyj
"""

import SingleLink as sl
import numpy as np

num_users=6
fblength=[0.1,0.2,0.3,0.4,0.5,0.6]#6 users and each user connectes with QNSP

#wavelength allocation
wavelength_matrix1=[[4,-4,-2,4,-4,-2],
                    [3, 2,-3,3,2, -3],
                    [1, 1, 1,-1,-1,-1]]

wavelength_matrix2=[[5,-5,-3,5,-5,-3],
                    [4,3,-4, 4,3,-4],
                    [2,1, 1,-2,-1,-1]]

wavelength_matrix3=[[5,-5,5,5,5,5],
                    [4,4,-4,4,4,4],
                    [3,3,3,-3,3,3],
                    [2,2,2,2,-2,2],
                    [1,1,1,1,1,-1]]

wavelength_matrix4=[[6,-6,-4,6,-6,-4],
                     [5,4,-5, 5,4,-5],
                     [3,2,1,-3,-2,-1]]

#def check_wavelength_allocation(wavelength_allocation, num_users):
   



'''
(1)assume no counts from other users and detector efficiency is 70%, total pair 
generation rate=10^6
(2)input: users: order of this user in the sequence E.g. Alice->0, Chloe->2
(3)output: ideal_single_link: counts at the link correspond to the user assume counts
from other user=0
pgr_per_channel: pair generation rate per cahnnel corresponds to wavelength allocation
'''
def ideal_single_counts(wavelength_allocation, num_users, user, fiber_length,pair_generation_rate):
    ideal_single_link=np.zeros(num_users, dtype=int)
    pgr_per_channel=pair_generation_rate/np.max(wavelength_allocation)
    singlecounts=sl.Calc_QKD_Parameters_For_Single_Link\
    (PairGenerationRate_per_second=pgr_per_channel, FibreLength_in_km = \
    [fiber_length[user],0],SingleCountsFromOtherUsers_per_second=[0,0,0,0])
    ideal_single_link=singlecounts['Singles rate for Alice'].HA/0.7
       
    return ideal_single_link, pgr_per_channel



'''
(1)This function calculates the number of channel (noise) for the specific node and effective channel
(2) input: row and col: correspond node and the effective channel in the table
(3)output: channels_at_nodes: availble channels at both nodes for the effective link
           effective_channel: which channels are used for effective communication
           num_eff_chan: number of effective channel
           num_noise_chan: number of noise channel at the sender node
'''
def num_channel_from_other_users(wavelength_allocation,sender,receiver):
    channels_at_nodes=[]
    effective_channel=[]
    for x in range(len(wavelength_allocation)):
        for y in range(len(wavelength_allocation[0])):
            if y==sender or y==receiver:
                channels_at_nodes.append(wavelength_allocation[x][y])
    num_eff_chan=0
    if sender>receiver:
        for num in channels_at_nodes[0::2]:
            if -num in channels_at_nodes:
                effective_channel.append(num)    
                num_eff_chan+=1
    else:
         for num in channels_at_nodes[1::2]:
             if -num in channels_at_nodes:
                 effective_channel.append(num)    
                 num_eff_chan+=1
    num_noise_chan=0
    for x in range(len(wavelength_allocation)):
        for y in range(len(wavelength_allocation[0])):
            if wavelength_allocation[x][y]*(-1)in effective_channel:
                num_noise_chan+=1
    num_noise_chan=num_noise_chan-len(effective_channel)
    return channels_at_nodes,effective_channel,num_eff_chan,num_noise_chan
      
       
      
'''
calculate the total counts at the nodes from other users(noise) given with effective link
'''
def node_noise_counts(wavelength_allocation, fiber_length, num_users,pair_generation_rate):
    node_others=np.zeros((num_users, num_users), dtype=int)
    for row in range(num_users):
        for col in range(num_users):
            single_counts,_=ideal_single_counts(wavelength_allocation, num_users,row, fiber_length,pair_generation_rate)
            _, _, _, num_noise_chan = num_channel_from_other_users(wavelength_allocation, row, col)
            node_others[row][col]=single_counts*num_noise_chan
    return node_others



'''
calculate the new secure key rate of each link base on the counts from other user
and the number of effective channel 
'''
def multi_secure_key_rate(wavelength_allocation, fiber_length, num_users, coincidence_window, pair_generation_rate):
    secure_key_rate_matrix=np.zeros((num_users, num_users), dtype=int)
    QBER_matrix=np.zeros((num_users, num_users), dtype=float)
    node_others=node_noise_counts(wavelength_allocation, fiber_length, num_users,pair_generation_rate)
    for row in range(len(secure_key_rate_matrix)):
        for col in range(len(secure_key_rate_matrix)):
            _,pgr_per_channel= ideal_single_counts(wavelength_allocation, num_users, row, fiber_length,pair_generation_rate)
            _,_,num_eff_chan,_=num_channel_from_other_users(wavelength_allocation,row,col)
            if row==col:
                secure_key_rate_matrix[row][col]=0
            else:
                new_single_link=sl.Calc_QKD_Parameters_For_Single_Link\
                (PairGenerationRate_per_second=pgr_per_channel*num_eff_chan, \
                 CoincidenceWindow_in_ps = coincidence_window, \
                 FibreLength_in_km = [fiber_length[row],fiber_length[col]], \
                 SingleCountsFromOtherUsers_per_second=[node_others[row][col],\
                node_others[row][col],node_others[col][row],node_others[col][row]])
                
                secure_key_rate_matrix[row][col]=new_single_link['Total secure key rate']
                QBER_matrix[row][col]=new_single_link['QBER between Alice and Bob']
    return secure_key_rate_matrix,QBER_matrix


'''
(1)find the best combination of coincidence window and pair generation rate to achieve the largest total 
secure key rate
(2)cd_break: break the cd for loop when coincidence window increases but total secure 
rate does not increase
(3)pgr_break: break the 'factor' for loop when any link has 0 secure key rate (When
 QBER exceed security threshold)
(4)optimized_parameter:1.sum of total secure key rate 2. average QBER 3. coincidence window 
   4. pair generation rate
'''
def optimize_coincidence_window(wavelength_allocation, fiber_length, num_users):
    optimized_secure_key_rate_matrix=np.zeros((num_users, num_users), dtype=int)
    optimized_QBER_matrix=np.zeros((num_users, num_users), dtype=float)
    optimized_parameter=[0,0,0,0];
    cd_step=5
    for cd in range(40,320,cd_step):
        for factor in range(5,8):
            for num in range(1,10):
                pgr=num*10**factor
                sum_of_secure_key_rate=0;
                skr_matrix, qber_matrix=multi_secure_key_rate(wavelength_allocation, fiber_length, num_users, cd,pgr)
                for row in range(len(skr_matrix)):
                    for col in range(len(skr_matrix)):
                        if row!=col and skr_matrix[row][col]==0:
                            sum_of_secure_key_rate=0
                        else:
                            sum_of_secure_key_rate=sum_of_secure_key_rate+skr_matrix[row][col]
                if sum_of_secure_key_rate>optimized_parameter[0]:
                    average_qber=np.sum(qber_matrix)/(num_users*(num_users-1))
                    optimized_parameter=[sum_of_secure_key_rate,average_qber,cd,pgr]
                    optimized_secure_key_rate_matrix=skr_matrix
                    optimized_QBER_matrix=qber_matrix
                
    return optimized_secure_key_rate_matrix, optimized_QBER_matrix, optimized_parameter
    

if __name__ == "__main__":
    
    #x=sl.Calc_QKD_Parameters_For_Single_Link(FibreLength_in_km=[1.3,3],SingleCountsFromOtherUsers_per_second=[1803501,1803501,1811788,1811788])
    #single_counts, pgr_per_channel=ideal_single_counts(wavelength_matrix1, num_users,2,fblength, 1E6)
    #a,b,c,d=num_channel_from_other_users(wavelength_matrix3,4,0)
    #node_others=node_noise_counts(wavelength_matrix2, fblength, num_users,1E6)
    #test,_=multi_secure_key_rate(wavelength_matrix1, fblength, 6, 150, 1E6)
    opt_skr_matrix, opt_qber_matrix, opt_parameter=optimize_coincidence_window(wavelength_matrix1, fblength, num_users)
    print('optimized coincidence window (ps):\n',opt_parameter[2])
    print('optimized pair generation rate:\n', "{:.0E}".format(opt_parameter[3]))
    print('optimized secure key rate for each links :\n', opt_skr_matrix)
    print('optimized total secure key rate: \n', opt_parameter[0])
    print('optimized quantum BER for each links: \n', opt_qber_matrix)
    print('average quantum BER: \n', opt_parameter[1])
    #print('total counts per second from other users: \n', node_others)
    #print(a,b,c,d)
    #print(test)
    #print('single counts per seconds for specific users: \n',single_counts)
    #print(x)#test for the single link file
    
    
    
    
    
    
    
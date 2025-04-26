#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 00:37:32 2025

@author: chenyj
"""

import SingleLink as sl
import numpy as np

num_users=6
fblength=[2, 4, 10, 11, 12, 30]#6 users and each user connectes with QNSP

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


'''
(1)assume no counts from other users and detector efficiency is 70%, total pair 
generation rate=10^6
(2)input: users: order of this user in the sequence E.g. Alice->0, Chloe->2
(3)output: ideal_single_link: counts at the link correspond to the user assume counts
from other user=0
pgr_per_channel: pair generation rate per cahnnel corresponds to wavelength allocation
'''
def ideal_single_counts(wavelength_allocation, num_users, user, fiber_length):
    ideal_single_link=np.zeros(num_users, dtype=int)
    pgr_per_channel=1E6/np.max(wavelength_allocation)
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
def num_channel_from_other_users(wavelength_allocation,row,col):
    channels_at_nodes=[]
    effective_channel=[]
    for x in range(len(wavelength_allocation)):
        for y in range(len(wavelength_allocation[0])):
            if y==row or y==col:
                channels_at_nodes.append(wavelength_allocation[x][y])
    num_eff_chan=0
    if row>col:
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
def node_noise_counts(wavelength_allocation, fiber_length, num_users):
    node_others=np.zeros((num_users, num_users), dtype=int)
    for row in range(num_users):
        for col in range(num_users):
            single_counts,_=ideal_single_counts(wavelength_allocation, num_users,row, fiber_length)
            _, _, _, num_noise_chan = num_channel_from_other_users(wavelength_allocation, row, col)
            node_others[row][col]=single_counts*num_noise_chan
    return node_others



'''
calculate the new secure key rate of each link base on the counts from other user
and the number of effective channel 
'''
def multi_secure_key_rate(wavelength_allocation, fiber_length, num_users, coincidence_window):
    secure_key_rate_matrix=np.zeros((num_users, num_users), dtype=int)
    node_others=node_noise_counts(wavelength_allocation, fiber_length, num_users)
    for row in range(len(secure_key_rate_matrix)):
        for col in range(len(secure_key_rate_matrix)):
            _,pgr_per_channel= ideal_single_counts(wavelength_allocation, num_users, row, fiber_length)
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
    return secure_key_rate_matrix


'''
find the best coincidence window for this scenerio to achieve the largest total 
secure key rate
'''
def optimize_coincidence_window(wavelength_allocation, fiber_length, num_users):
    optimized_secure_key_rate_matrix=np.zeros((num_users, num_users), dtype=int)
    optimized_cd=[0,0];
    for cd in range(10,800,10):
        sum_of_secure_key_rate=0;
        skr_matrix=multi_secure_key_rate(wavelength_allocation, fiber_length, num_users, cd)
        for row in range(len(skr_matrix)):
            for col in range(len(skr_matrix)):
                sum_of_secure_key_rate=sum_of_secure_key_rate+skr_matrix[row][col]
        if sum_of_secure_key_rate>optimized_cd[0]:
            optimized_cd[0]=sum_of_secure_key_rate
            optimized_cd[1]=cd
            optimized_secure_key_rate_matrix=skr_matrix
    return optimized_secure_key_rate_matrix, optimized_cd[1], optimized_cd[0]
    

if __name__ == "__main__":
    
    #x=sl.Calc_QKD_Parameters_For_Single_Link(FibreLength_in_km=[1.3,3],SingleCountsFromOtherUsers_per_second=[1803501,1803501,1811788,1811788])
    #single_counts, pgr_per_channel=ideal_single_counts(wavelength_matrix1, num_users,2,fblength)
    #a,b,c,d=num_channel_from_other_users(wavelength_matrix3,0,4)
    #node_others=node_noise_counts(wavelength_matrix2, fblength, num_users)
    opt_skr_matrix, opt_cd, opt_tskr=optimize_coincidence_window(wavelength_matrix2, fblength, num_users)
    
    print('optimized coincidence window:\n',opt_cd)
    print('optimized secure key rate for each links :\n', opt_skr_matrix)
    print('optimized total secure key rate: \n', opt_tskr)
    #print('total counts per second from other users: \n', node_others)
    #print(a,b,c,d)
    #print('single counts per seconds for specific users: \n',single_counts)
    #print(x)#test for the single link file
    
    
    
    
    
    
    
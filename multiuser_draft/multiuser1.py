#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 23:15:12 2025

@author: chenyj
"""

import SingleLink as sl
import numpy as np

#fiber length of all the links in km
def fiberlength(num_users):
    fblength=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5]#15 fiber length for 6 users
    
    fbmatrix = np.zeros((num_users, num_users), dtype=float)
    array_position=0
    # Fill the matrix with values based on positions
    for i in range(num_users):
        for j in range(i, num_users):
            if i==j:
                fbmatrix[i, j] = 0
            else:
                fbmatrix[i, j] = fblength[array_position]
                fbmatrix[j, i] = fblength[array_position]
                array_position=array_position+1

    return fbmatrix


#record single counts correspond to the fiber length from the previous function
#assume no counts from other users and detector efficiency is 70%
def all_single_counts(fbmatrix, num_users):
    total_counts=0;
    single_detection_matrix=np.zeros((num_users, num_users), dtype=int)
    node_counts=np.zeros((num_users), dtype=int)
    for row in range(len(fbmatrix)):
        for col in range(len(fbmatrix)):
            if row==col:
                single_detection_matrix[row][col]=0
            else:
                singlecounts=sl.Calc_QKD_Parameters_For_Single_Link(FibreLength_in_km = \
                [fbmatrix[row][col],fbmatrix[row][col]],SingleCountsFromOtherUsers_per_second=[0,0,0,0])
                single_detection_matrix[row][col]=singlecounts['Singles rate for Alice'].HA/0.7
            total_counts=total_counts+single_detection_matrix[row][col]
        node_counts[row]=total_counts
        total_counts=0
    return single_detection_matrix, node_counts


#calculate the counts from other users at each nodes and each links
#counts from other users=total counts at the node-counts from its link
def node_other_counts(single_detection_matrix, node_counts):
    node_others=np.zeros((num_users, num_users), dtype=int)
    for row in range(len(node_others)):
        node_others[row]=(node_counts-single_detection_matrix[row])
        for col in range(len(node_others)):
            if row==col:
                node_others[row][col]=0
    return node_others
                                       
'''
calculate the new secure key rate of each link base on the counts from other user
node_others:counts from other users
fbmatrix: fiber length of each links
'''
def multi_secure_key_rate(node_others, fbmatrix, coincidence_window):
    secure_key_rate_matrix=np.zeros((num_users, num_users), dtype=int)
    for row in range(len(secure_key_rate_matrix)):
        for col in range(len(secure_key_rate_matrix)):
            if row==col:
                secure_key_rate_matrix[row][col]=0
            else:
                new_single_link=sl.Calc_QKD_Parameters_For_Single_Link\
                (CoincidenceWindow_in_ps = coincidence_window, \
                 FibreLength_in_km = [fbmatrix[row][col],fbmatrix[row][col]], \
                 SingleCountsFromOtherUsers_per_second=[node_others[row][col],\
                node_others[row][col],node_others[col][row],node_others[col][row]])
                
                secure_key_rate_matrix[row][col]=new_single_link['Total secure key rate']
    return secure_key_rate_matrix

'''
find the best coincidence window for this scenerio to achieve the largest total 
secure key rate
'''
def optimize_coincidence_window():
    optimized_secure_key_rate_matrix=np.zeros((num_users, num_users), dtype=int)
    optimized_cd=[0,0];
    for cd in range(10,400,1):
        sum_of_secure_key_rate=0;
        skr_matrix=multi_secure_key_rate(node_others, fbmatrix,cd)
        for row in range(len(skr_matrix)):
            for col in range(len(skr_matrix)):
                sum_of_secure_key_rate=sum_of_secure_key_rate+skr_matrix[row][col]
        if sum_of_secure_key_rate>optimized_cd[0]:
            optimized_cd[0]=sum_of_secure_key_rate
            optimized_cd[1]=cd
            optimized_secure_key_rate_matrix=skr_matrix
    return optimized_secure_key_rate_matrix, optimized_cd[1], optimized_cd[0]
    

if __name__ == "__main__":
    
    x=sl.Calc_QKD_Parameters_For_Single_Link(FibreLength_in_km=[1.3,1.3],SingleCountsFromOtherUsers_per_second=[1803501,1803501,1811788,1811788])
    print(x)#test for the single link file
    num_users=6
    fbmatrix = fiberlength(num_users)
    single_detection_matrix, node_counts = all_single_counts(fbmatrix, num_users)
    node_others=node_other_counts(single_detection_matrix, node_counts)
    secure_key_rate=multi_secure_key_rate(node_others, fbmatrix,106)
    optimized_secure_key_rate_matrix, optimized_coincidence_window, optimized_total_secure_rate=optimize_coincidence_window()
    #print(x)#test for the single link file
    print('fiber length between all the nodes（km):\n' ,fbmatrix)
    print('single counts per second at all the links (assume counts=0 for other users): \n',single_detection_matrix)
    print('total counts per second at the nodes:\n', node_counts)
    print('total counts per second from other users: \n', node_others)
    #print('secure key rate at each links(adjust coincidence window manually): \n',secure_key_rate)
    print('optimized coincidence window:\n',optimized_coincidence_window )
    print('optimized secure key rate for each links :\n', optimized_secure_key_rate_matrix)
    #print('optimized total secure key rate: \n', optimized_total_secure_rate)
    

    
    
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 00:03:46 2025

@author: chenyj
"""

import multiuser_final as mu

import numpy as np
import matplotlib.pyplot as plt

netplot=mu.QKDNetwork(num_users=6)


def plot_dynamic_user(netplot):
    dynamic_user_WA_set = [netplot.wa_fullmesh, netplot.wa_partial_mesh1,
    netplot.wa_partial_mesh2]
    dynamic_skr_list = []
    for m in range(3):
        i=m%3
        wavelength_allocation = netplot.validate_wavelength_allocation(
            dynamic_user_WA_set[i]
        )
        results = netplot.optimize_network_parameters(wavelength_allocation, netplot.num_users,i,0, 0)
        dynamic_skr_list.append(netplot.matrix_to_list(results['secure key rate matrix']))
    labels = ['AB', 'AC', 'AD', 'AF', 'AG','BC', 'BD', 'BF', 'BG', 'CD', 'CF', 'CG', 'DF', 'DG', 'FG']
    x = np.arange(len(labels))  
    bar_width = 0.25  

    
    plt.figure(figsize=(12, 6))
    plt.bar(x - bar_width, dynamic_skr_list[0], width=bar_width, label='full mesh', color='blue')
    plt.bar(x, dynamic_skr_list[1], width=bar_width, label='partial mesh1', color='red')
    plt.bar(x + bar_width, dynamic_skr_list[2], width=bar_width, label='partial mesh2', color='green')


    plt.xticks(ticks=x, labels=labels, rotation=45)  # Set x-axis labels and rotate
    plt.xlabel("Names (D)")
    plt.ylabel("Values")
    plt.title(f'dynamic topologies with fiber length=: {results["fiber_length"]}')
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.savefig('dynamic_topologies.png')
    plt.show()



def plot_wavelength_multiplex(netplot):
    dynamic_user_WA_set = [netplot.wa_subnet1, netplot.wa_subnet2,
    netplot.wa_fullmesh]
    dynamic_skr_list = []
    for m in range(3):
        i=m%3
        wavelength_allocation = netplot.validate_wavelength_allocation(
            dynamic_user_WA_set[i]
        )    
        results = netplot.optimize_network_parameters(wavelength_allocation, netplot.num_users,i,1,1)
        print(results['sum of secure key rate'])
        dynamic_skr_list.append(netplot.matrix_to_list(results['secure key rate matrix']))
    labels = ['AB', 'AC', 'AD', 'AF', 'AG','BC', 'BD', 'BF', 'BG', 'CD', 'CF', 'CG', 'DF', 'DG', 'FG']
    x = np.arange(len(labels))  
    bar_width = 0.25  

    
    plt.figure(figsize=(12, 6))
    plt.bar(x - bar_width, dynamic_skr_list[0], width=bar_width, label='5 pairs wavelength', color='blue')
    plt.bar(x, dynamic_skr_list[1], width=bar_width, label='6 pairs wavelength', color='red')
    plt.bar(x + bar_width, dynamic_skr_list[2], width=bar_width, label='15 pairs wavelength', color='green')


    plt.xticks(ticks=x, labels=labels, rotation=45)  # Set x-axis labels and rotate
    plt.xlabel("Names (D)")
    plt.ylabel("SKR (bps)")
    plt.title(f'Relationship between number of wavelength channels and SKR, fiber length=: {results["fiber_length"]}')
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.savefig('num_wavelength_channel.png')
    plt.show()
    
    
    
def plot_additional_channels(netplot):
    dynamic_user_WA_set = [netplot.wa_fullmesh, netplot.wa_add_channel1, netplot.wa_add_channel2]
    dynamic_skr_list = []
    for m in range(3):
        i=m%3
        wavelength_allocation = netplot.validate_wavelength_allocation(
            dynamic_user_WA_set[i]
        )
        results = netplot.optimize_network_parameters(wavelength_allocation, netplot.num_users,i,1,1)
        dynamic_skr_list.append(netplot.matrix_to_list(results['secure key rate matrix']))
    labels = ['AB', 'AC', 'AD', 'AF', 'AG','BC', 'BD', 'BF', 'BG', 'CD', 'CF', 'CG', 'DF', 'DG', 'FG']
    x = np.arange(len(labels))  
    bar_width = 0.25  

    
    plt.figure(figsize=(12, 6))
    plt.bar(x - bar_width, dynamic_skr_list[0], width=bar_width, label='AB with 1 pair', color='blue')
    plt.bar(x, dynamic_skr_list[1], width=bar_width, label='AB with 2 pairs', color='red')
    plt.bar(x + bar_width, dynamic_skr_list[2], width=bar_width, label='AB with 3 pairs', color='green')


    plt.xticks(ticks=x, labels=labels, rotation=45)  # Set x-axis labels and rotate
    plt.xlabel("Names (D)")
    plt.ylabel("SKR (bps)")
    plt.title(f'Effects of additional wavelength channel to SKR, fiber length=: {results["fiber_length"]}')
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.savefig('additional_wavelength_channel.png')
    plt.show()


def plot_optimization_priority(netplot):
    dynamic_user_WA_set = [netplot.wa_fullmesh, netplot.wa_fullmesh]
    dynamic_skr_list = []
    for m in range(2):
        i=m%2
        wavelength_allocation = netplot.validate_wavelength_allocation(
            dynamic_user_WA_set[i]
        )
        results = netplot.optimize_network_parameters(wavelength_allocation, netplot.num_users,i,0,2*i)
        dynamic_skr_list.append(netplot.matrix_to_list(results['secure key rate matrix']))
    labels = ['AB', 'AC', 'AD', 'AF', 'AG','BC', 'BD', 'BF', 'BG', 'CD', 'CF', 'CG', 'DF', 'DG', 'FG']
    x = np.arange(len(labels))  
    bar_width = 0.25  

    
    plt.figure(figsize=(12, 6))
    plt.bar(x - bar_width, dynamic_skr_list[0], width=bar_width, label='Total SKR priority', color='blue')
    plt.bar(x, dynamic_skr_list[1], width=bar_width, label='Minimum SKR priority', color='red')
    


    plt.xticks(ticks=x, labels=labels, rotation=45)  # Set x-axis labels and rotate
    plt.xlabel("Names (D)")
    plt.ylabel("SKR (bps)")
    plt.title(f'comparison between different optimization mode, fiber length=: {results["fiber_length"]}')
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.savefig('optimization_mode.png')
    plt.show()


def plot_matrix(netplot, wavelength_allocation):
    
    results = netplot.optimize_network_parameters(wavelength_allocation, netplot.num_users, 0, 0, 0)
    labels = ['Alice', 'Bob','Chloe', 'Dave', 'Feng', 'Gopi']
    x = np.arange(len(labels)) 
    
    
    plt.figure(figsize=(15, 10))
    fig = plt.figure(figsize=(15, 10))
    fig.suptitle(f'Fiber length for this simulation: {results["fiber_length"]}',fontsize=14)

    # Plot SKR matrix
    plt.subplot(2, 2, 1)
    plt.xticks(ticks=x, labels=labels, rotation=45)
    plt.yticks(ticks=x, labels=labels)
    plt.imshow(results['secure key rate matrix'])
    plt.colorbar(label='Secure Key Rate')
    plt.title(f'Secure Key Rate Matrix for: {results["wavelength_allocation"]}')
    
    # Plot QBER matrix
    plt.subplot(2, 2, 2)
    plt.xticks(ticks=x, labels=labels, rotation=45)
    plt.yticks(ticks=x, labels=labels)
    plt.imshow(results['QBER matrix'])
    plt.colorbar(label='QBER')
    plt.title(f'QBER Matrix for: {results["wavelength_allocation"]}')
    '''
    # Plot coincidence window matrix
    plt.subplot(2, 2, 3)
    plt.xticks(ticks=x, labels=labels, rotation=45)
    plt.yticks(ticks=x, labels=labels)
    plt.imshow(results['coincidence_window'])
    plt.colorbar(label='Coincidence Window')
    plt.title(f'Coincidence Window Matrix for: {results["wavelength_allocation"]}')
    '''
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.savefig('network_parameter_matrix.jpg')
    plt.show()

def clean_line(line):
    # Remove unsupported characters (non-latin-1)
    return line.encode('latin-1', errors='ignore').decode('latin-1')
    
    
def main():
    plot_dynamic_user(netplot)
    #plot_wavelength_multiplex(netplot)
    #plot_additional_channels(netplot)
    #plot_optimization_priority(netplot)
    #plot_matrix(netplot,netplot.wa_fullmesh)
    
    '''
    from fpdf import FPDF

    out = 'multiuser.pdf'     # output pdf file

    # Create PDF instance
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Helvetica", size=10)

    with open('multiuser_final.py', 'r', encoding='utf-8') as file:
        for line in file:
            cleaned = clean_line(line)
            pdf.cell(200, 5, txt=cleaned.strip(), ln=True)

    pdf.output(out)

    print("PDF is Saved Successfully")
    '''

if __name__ == "__main__":
    main()
    
    
    
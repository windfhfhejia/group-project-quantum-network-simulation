# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 15:51:14 2019

@author: xv18766
"""

import os
import csv
import LinkModel5 as lm
import numpy as np
import SplittingTableGenerator as STG
from scipy import optimize


def pair_rate_vs(network, output_file_prefix, nodes, node_pairs,
                 pair_rate_limit='secure_key_rate death',
                 pair_rate_step_size=1E4):

    node_names = {node: node.__str__() for node in nodes}
    node_pair_names = {node_pair: node_pair[0].__str__() + ' and '
                       + node_pair[1].__str__() for node_pair in node_pairs}
    variables = ['Pair rate of sources',
                 *['Singles rate for ' + name 
                   for name in node_names.values()],
                 *['Pair rate between ' + names
                   for names in node_pair_names.values()],
                 *['Shifted key rate between ' + names
                   for names in node_pair_names.values()],
                 *['Accidentals between ' + names
                   for names in node_pair_names.values()],
                 *['QBER between ' + names
                   for names in node_pair_names.values()],
                 *['Secure key rate between ' + names
                   for names in node_pair_names.values()],
                 'Total secure key rate',
                ]

    def get_data(pair_rate):
        for source in network.splitting_table.sources:
            source.pair_rate = pair_rate
        data = {}
        data['Pair rate of sources'] = pair_rate
        for node, name in node_names.items():
            data['Singles rate for ' + name] = (
                    network.total_singles_rate(node))
        for node_pair, names in node_pair_names.items():
            data['Pair rate between ' + names] = (
                    network.total_pair_rate(node_pair))
            data['Shifted key rate between ' + names] = (
                    network.shifted_key_rate(node_pair))
            data['Accidentals between ' + names] = (
                    network.accidentals(node_pair))
            data['QBER between ' + names] = (
                    network.QBER(node_pair))
            data['Secure key rate between ' + names] = (
                    network.secure_key_rate(node_pair))
        data['Total secure key rate'] = (network.total_secure_key_rate())
        return data

    with open(output_file_prefix + '.csv', 'w', newline = '') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=variables)
        writer.writeheader()
        if pair_rate_limit == 'secure_key_rate death':
            pair_rate = 0
            i = 1
            writer.writerow(get_data(pair_rate))
            while (max([network.secure_key_rate(node_pair)
                    for node_pair in node_pairs]) > 0 or i < 10):
                pair_rate += pair_rate_step_size
                data = get_data(pair_rate)
                writer.writerow(data)
                i += 1
        else:
            pair_rate_limit = float(pair_rate_limit)
            for pair_rate in np.linspace(0, pair_rate_limit, 100):
                data = get_data(pair_rate)
                writer.writerow(data)

#    fixed_vars = {}
#    with open(output_file_prefix + '_fixed_variables.csv',
#              'r', newline = '') as csvfile:
#        writer = csv.DictWriter(csvfile)
#        writer.writeheader()
#        writer.writerow(fixed_vars)
    
    return variables


def get_negative_secure_key_rate(pair_rate, network, node_pair):
    for source in network.splitting_table.sources:
        source.pair_rate = pair_rate
    return -network.secure_key_rate(node_pair)


def max_secure_key_rate(network, node_pair):
    f = get_negative_secure_key_rate
    maximum = optimize.minimize_scalar(f, bounds=(0, 1E8),
                                       args=(network, node_pair),
                                       method='bounded')
    if maximum.x > 9*1E7:
        print("I'm worried about how close the value is to the bound.")
    return (maximum.x, -maximum.fun)


def jitter_vs_max_SKR(network, output_file_prefix, node_pairs,
              jitter_limit=1E3):

    node_pair_names = {node_pair: node_pair[0].__str__() + ' and '
                       + node_pair[1].__str__() for node_pair in node_pairs}
    variables = ['Detector jitter',
                 *['Optimal pair rate for ' + names
                   for names in node_pair_names.values()],
                 *['Max secure key rate between ' + names
                   for names in node_pair_names.values()],
                 'Total secure key rate'
                ]

    def get_max_secure_key_rate(jitter):
        data = {}
        data['Detector jitter'] = jitter
        for node_pair, names in node_pair_names.items():
            for node in node_pair:
                for detector in node.detectors:
                    detector.jitter = jitter
        # This line just makes the network recognise the new jitter
        for k in network.coincidence_probabilities:
            network.coincidence_probabilities[k] = (
                    network.calc_cp_between_nodes(k))
        for node_pair, names in node_pair_names.items():
            optimum = max_secure_key_rate(network, node_pair)
            data['Optimal pair rate for ' + names] = optimum[0]
            data['Max secure key rate between ' + names] = optimum[1]
        data['Total secure key rate'] = network.total_secure_key_rate()
        return data

    with open(output_file_prefix + '.csv', 'w', newline = '') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=variables)
        writer.writeheader()
        for jitter in np.linspace(1, jitter_limit, 100):
            data = get_max_secure_key_rate(jitter)
            writer.writerow(data)

    return variables
#
#
#def negative_SKR_given_coinc_window(coincidence_window, network, node_pair):
#    network.set_coincidence_windows(node_pair, coincidence_window)
#    return -network.secure_key_rate(node_pair, nonnegative=False)
#
#def optimal_coinc_window(network, node_pair):
#    f = negative_SKR_given_coinc_window
#    maximum = optimize.minimize_scalar(
#        f, bounds=(0, 1E12), args=(network, node_pair,), method='bounded')
#    if maximum.x > 9*1E11:
#        print("The coincidence window is nearly one second wide.")
#    return (maximum.x, -maximum.fun)
#

def jitter_vs_optimal_coincidence_window(
        network, output_file_prefix, node_pairs,
        jitter_limit='secure_key_rate death',
        ubound="default", jitter_step_size=25):
# jitter_limit = 10**3.34 works well for 4p WDM only
    node_pair_names = {node_pair: node_pair[0].__str__() + ' and '
                       + node_pair[1].__str__() for node_pair in node_pairs}
    indexed_sources = {node_pair: (list(enumerate(list(network.sources[node_pair]))))
                        for node_pair in node_pairs}
    variables = ['Detector jitter',
                 *['Optimal coincidence window for ' + names +
                   ', Source ' + str(i)
                   for node_pair, names in node_pair_names.items()
                   for i, source in indexed_sources[node_pair]],
                 *['Secure key rate between ' + names
                   for names in node_pair_names.values()],
                 'Total secure key rate',
                ]

    def get_optimal_coinc_window(jitter):
        data = {}
        data['Detector jitter'] = jitter
        for node_pair, names in node_pair_names.items():
            for node in node_pair:
                for detector in node.detectors:
                    detector.jitter = jitter
        for node_pair, names in node_pair_names.items():
            for i, source in indexed_sources[node_pair]:
                optimum = network.optimize_coincidence_window(
                            node_pair, source=source, ubound=ubound)
                data['Optimal coincidence window for ' + names +
                     ', Source ' + str(i)] = optimum[0]
            data['Secure key rate between ' + names] = optimum[1]
        data['Total secure key rate'] = network.total_secure_key_rate()
        return data

    with open(output_file_prefix + '.csv', 'w', newline = '') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=variables)
        writer.writeheader()
        if jitter_limit == 'secure_key_rate death':
            jitter = 1
            data = get_optimal_coinc_window(jitter)
            writer.writerow(data)
            while (max([network.secure_key_rate(node_pair)
                    for node_pair in node_pairs]) > 0):
                jitter += jitter_step_size
                data = get_optimal_coinc_window(jitter)
                writer.writerow(data)
        else:
            jitter_limit = float(jitter_limit)
            for jitter in np.linspace(1, jitter_limit, 100):
                data = get_optimal_coinc_window(jitter)
                writer.writerow(data)

    return variables


def pair_rate_vs_optimal_coincidence_window(
        network, output_file_prefix, nodes, node_pairs,
        pair_rate_limit='secure_key_rate death',
        ubound="default", pair_rate_step_size=1E4):

    node_names = {node: node.__str__() for node in nodes}
    node_pair_names = {node_pair: node_pair[0].__str__() + ' and '
                       + node_pair[1].__str__() for node_pair in node_pairs}
    indexed_sources = {node_pair: (list(enumerate(list(network.sources[node_pair]))))
                        for node_pair in node_pairs}
    variables = ['Pair rate of sources',
                 *['Optimal coincidence window for ' + names +
                   ', Source ' + str(i)
                   for node_pair, names in node_pair_names.items()
                   for i, source in indexed_sources[node_pair]],
                 *['Singles rate for ' + name 
                   for name in node_names.values()],
                 *['Pair rate between ' + names
                   for names in node_pair_names.values()],
                 *['Shifted key rate between ' + names
                   for names in node_pair_names.values()],
                 *['Accidentals between ' + names
                   for names in node_pair_names.values()],
                 *['QBER between ' + names
                   for names in node_pair_names.values()],
                 *['Secure key rate between ' + names
                   for names in node_pair_names.values()],
                 'Total secure key rate',
                ]

    def get_data(pair_rate):
        data = {}
        data['Pair rate of sources'] = pair_rate
        for source in network.splitting_table.sources:
            source.pair_rate = pair_rate
        for node, name in node_names.items():
            data['Singles rate for ' + name] = (
                    network.total_singles_rate(node))
        for node_pair, names in node_pair_names.items():
            data['Pair rate between ' + names] = (
                    network.total_pair_rate(node_pair))
            data['Shifted key rate between ' + names] = (
                    network.shifted_key_rate(node_pair))
            data['Accidentals between ' + names] = (
                    network.accidentals(node_pair))
            data['QBER between ' + names] = (
                    network.QBER(node_pair))
            for i, source in indexed_sources[node_pair]:
                optimum = network.optimize_coincidence_window(
                            node_pair, source=source, ubound=ubound)
                data['Optimal coincidence window for ' + names +
                     ', Source ' + str(i)] = optimum[0]
            data['Secure key rate between ' + names] = optimum[1]
        data['Total secure key rate'] = network.total_secure_key_rate()
        return data

    with open(output_file_prefix + '.csv', 'w', newline = '') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=variables)
        writer.writeheader()
        if pair_rate_limit == 'secure_key_rate death':
            pair_rate = 0
            i = 1
            writer.writerow(get_data(pair_rate))
            while (max([network.secure_key_rate(node_pair)
                    for node_pair in node_pairs]) > 0 or i < 10):
                pair_rate += pair_rate_step_size
                data = get_data(pair_rate)
                writer.writerow(data)
                i += 1
        else:
            pair_rate_limit = float(pair_rate_limit)
            for pair_rate in np.linspace(0, pair_rate_limit, 100):
                data = get_data(pair_rate)
                writer.writerow(data)

    return variables


def transmission_vs(network, nodes, node_pairs,
                    node_transmission):
    """

    Note that this function only considers transmission at nodes, NOT
    in the splitting table. Set the fibre transmissions to 1 if you
    want the full range of possible transmissions!

    Only the transmissions of the 'nodes' will be changed, not
    necessarily all the nodes in the 'node_pairs'.
    """

    node_names = {node: node.__str__() for node in nodes}
    node_pair_names = {node_pair: node_pair[0].__str__() + ' and '
                       + node_pair[1].__str__() for node_pair in node_pairs}
    variables = [
                *['Singles rate (/s) for ' + name 
                   for name in node_names.values()],
                *['Shifted key rate (/s) between ' + names
                   for names in node_pair_names.values()],
                *['Accidentals (/s) between ' + names
                   for names in node_pair_names.values()],
                *['QBER between ' + names
                   for names in node_pair_names.values()],
                *['Secure key rate (/s) between ' + names
                   for names in node_pair_names.values()],
                'Total secure key rate (/s)',
                ]

    def get_data(transmission):
        for node in nodes:
            node.transmissions = (
                    transmission, transmission, transmission, transmission)
        data = {}
        for node, name in node_names.items():
            data['Singles rate for ' + name] = (
                    network.detection_rates(node))
        for node_pair, names in node_pair_names.items():
            data['Shifted key rate between ' + names] = (
                    network.shifted_key_rate(node_pair))
            data['Accidentals between ' + names] = (
                    network.accidentals(node_pair))
            data['QBER between ' + names] = (
                    network.QBER(node_pair))
            data['Secure key rate between ' + names] = (
                    network.secure_key_rate(node_pair))
        data['Total secure key rate'] = network.total_secure_key_rate()
        return data
    # ~ if output_file_prefix != "":
        # ~ with open(output_file_prefix + '.csv', 'w', newline = '') as csvfile:
            # ~ writer = csv.DictWriter(csvfile, fieldnames=variables)
            # ~ writer.writeheader()
            # ~ if min_transmission_in_dB == 'secure_key_rate death':
                # ~ transmission = 1
                # ~ writer.writerow(get_data(transmission))
                # ~ while (max([network.secure_key_rate(node_pair)
                        # ~ for node_pair in node_pairs]) > 0):
                    # ~ transmission *= transmission_ratio
                    # ~ data = get_data(transmission)
                    # ~ writer.writerow(data)
            # ~ else:
                # ~ min_transmission_in_dB = float(min_transmission_in_dB)
                # ~ for transmission in np.logspace(
                        # ~ min_transmission_in_dB/10, 0, num=100):
                    # ~ data = get_data(transmission)
                    # ~ writer.writerow(data)
                    
                    
    data=get_data(node_transmission) #Note all iteration over teh various transmision values has been removed.
    
    return data


def transmission_vs_optimal_coincidence_window(
        network, output_file_prefix, nodes, node_pairs,
        min_transmission_in_dB='secure_key_rate death',
        ubound="default", transmission_ratio=0.95):
    """

    Note that this function only considers transmission at nodes, NOT
    in the splitting table. Set the fibre transmissions to 1 if you
    want the full range of possible transmissions!

    Only the transmissions of the 'nodes' will be changed, not
    necessarily all the nodes in the 'node_pairs'.
    """

    node_names = {node: node.__str__() for node in nodes}
    node_pair_names = {node_pair: node_pair[0].__str__() + ' and '
                       + node_pair[1].__str__() for node_pair in node_pairs}
    indexed_sources = {node_pair: (list(enumerate(list(network.sources[node_pair]))))
                        for node_pair in node_pairs}
    variables = ['Transmission at nodes',
                 *['Optimal coincidence window for ' + names +
                   ', Source ' + str(i)
                   for node_pair, names in node_pair_names.items()
                   for i, source in indexed_sources[node_pair]],
                 *['Singles rate for ' + name 
                   for name in node_names.values()],
                 *['Pair rate between ' + names
                   for names in node_pair_names.values()],
                 *['Shifted key rate between ' + names
                   for names in node_pair_names.values()],
                 *['Accidentals between ' + names
                   for names in node_pair_names.values()],
                 *['QBER between ' + names
                   for names in node_pair_names.values()],
                 *['Secure key rate between ' + names
                   for names in node_pair_names.values()],
                 'Total secure key rate',
                ]

    def get_data(transmission):
        for node in nodes:
            node.transmissions = (
                    transmission, transmission, transmission, transmission)
        data = {}
        data['Transmission at nodes'] = transmission
        for node, name in node_names.items():
            data['Singles rate for ' + name] = (
                    network.total_singles_rate(node))
        for node_pair, names in node_pair_names.items():
            data['Pair rate between ' + names] = (
                    network.total_pair_rate(node_pair))
            data['Shifted key rate between ' + names] = (
                    network.shifted_key_rate(node_pair))
            data['Accidentals between ' + names] = (
                    network.accidentals(node_pair))
            data['QBER between ' + names] = (
                    network.QBER(node_pair))
            for i, source in indexed_sources[node_pair]:
                optimum = network.optimize_coincidence_window(
                            node_pair, source=source, ubound=ubound)
                data['Optimal coincidence window for ' + names +
                     ', Source ' + str(i)] = optimum[0]
            data['Secure key rate between ' + names] = optimum[1]
        data['Total secure key rate'] = network.total_secure_key_rate()
        return data

    with open(output_file_prefix + '.csv', 'w', newline = '') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=variables)
        writer.writeheader()
        if min_transmission_in_dB == 'secure_key_rate death':
            transmission = 1
            writer.writerow(get_data(transmission))
            while (max([network.secure_key_rate(node_pair)
                    for node_pair in node_pairs]) > 0):
                transmission *= transmission_ratio
                data = get_data(transmission)
                writer.writerow(data)
        else:
            min_transmission_in_dB = float(min_transmission_in_dB)
            for transmission in np.logspace(
                    min_transmission_in_dB/10, 0, num=100):
                data = get_data(transmission)
                writer.writerow(data)

    return variables



def coincidence_window_vs(network, output_file_prefix, nodes, node_pairs,
                 coincidence_window_limit='secure_key_rate death',
                 coincidence_window_step_size=25,
                 pair_rate=1E5):

    node_names = {node: node.__str__() for node in nodes}
    node_pair_names = {node_pair: node_pair[0].__str__() + ' and '
                       + node_pair[1].__str__() for node_pair in node_pairs}
    variables = ['Pair rate of sources',
                 'Coincidence window',
                 *['Singles rate for ' + name 
                   for name in node_names.values()],
                 *['Pair rate between ' + names
                   for names in node_pair_names.values()],
                 *['Shifted key rate between ' + names
                   for names in node_pair_names.values()],
                 *['Accidentals between ' + names
                   for names in node_pair_names.values()],
                 *['QBER between ' + names
                   for names in node_pair_names.values()],
                 *['Secure key rate between ' + names
                   for names in node_pair_names.values()],
                 'Total secure key rate',
                ]

    for source in network.splitting_table.sources:
        source.pair_rate = pair_rate

    def get_data(coincidence_window):
        for node_pair in node_pairs:
            network.set_coincidence_windows(node_pair, value=coincidence_window)
        data = {}
        data['Pair rate of sources'] = pair_rate
        data['Coincidence window'] = coincidence_window
        for node, name in node_names.items():
            data['Singles rate for ' + name] = (
                    network.total_singles_rate(node))
        for node_pair, names in node_pair_names.items():
            data['Pair rate between ' + names] = (
                    network.total_pair_rate(node_pair))
            data['Shifted key rate between ' + names] = (
                    network.shifted_key_rate(node_pair))
            data['Accidentals between ' + names] = (
                    network.accidentals(node_pair))
            data['QBER between ' + names] = (
                    network.QBER(node_pair))
            data['Secure key rate between ' + names] = (
                    network.secure_key_rate(node_pair))
        data['Total secure key rate'] = (network.total_secure_key_rate())
        return data

    with open(output_file_prefix + '.csv', 'w', newline = '') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=variables)
        writer.writeheader()
        if coincidence_window_limit == 'secure_key_rate death':
            coincidence_window = 0
            i = 1
            writer.writerow(get_data(coincidence_window))
            while (max([network.secure_key_rate(node_pair)
                    for node_pair in node_pairs]) > 0 or i < 10):
                coincidence_window += coincidence_window_step_size
                data = get_data(coincidence_window)
                writer.writerow(data)
                i += 1
        else:
            coincidence_window_limit = float(coincidence_window_limit)
            for coincidence_window in np.linspace(
                    0, coincidence_window_limit, 100):
                data = get_data(coincidence_window)
                writer.writerow(data)

#    fixed_vars = {}
#    with open(output_file_prefix + '_fixed_variables.csv',
#              'r', newline = '') as csvfile:
#        writer = csv.DictWriter(csvfile)
#        writer.writeheader()
#        writer.writerow(fixed_vars)
    
    return variables










def convolve_Jitters(*jitters):
    sum_of_squares = 0
    for j in jitters:
        sum_of_squares = sum_of_squares + j*j
    convolvedJitter = np.sqrt(sum_of_squares)
    return (convolvedJitter)



def distance_vs_optimal_coincidence_window(
        network, output_file_prefix, nodes, node_pairs,
        ubound=1000, lossperkm = 0.3, chromaticdispersion = 10, bandwidth = 0.8,Distance_step_size = 1):
    """

    Note that this function only considers transmission at nodes, NOT
    in the splitting table. Set the fibre transmissions to 1 if you
    want the full range of possible transmissions!

    Only the transmissions of the 'nodes' will be changed, not
    necessarily all the nodes in the 'node_pairs'.
    
    This function accepts new parameters:
     lossperkm which is in dB 
     chromatic dispersion in ps per nm per km
     bandwidth in nm
     
    """

    node_names = {node: node.__str__() for node in nodes}
    node_pair_names = {node_pair: node_pair[0].__str__() + ' and '
                       + node_pair[1].__str__() for node_pair in node_pairs}
    indexed_sources = {node_pair: (list(enumerate(list(network.sources[node_pair]))))
                        for node_pair in node_pairs}
    variables = ['Distance (km)', 'chromatic dispersion (ps)', 'Transmission at nodes',
                 *['Optimal coincidence window for ' + names +
                   ', Source ' + str(i)
                   for node_pair, names in node_pair_names.items()
                   for i, source in indexed_sources[node_pair]],
                 *['Singles rate for ' + name 
                   for name in node_names.values()],
                 *['Pair rate between ' + names
                   for names in node_pair_names.values()],
                 *['Shifted key rate between ' + names
                   for names in node_pair_names.values()],
                 *['Accidentals between ' + names
                   for names in node_pair_names.values()],
                 *['QBER between ' + names
                   for names in node_pair_names.values()],
                 *['Secure key rate between ' + names
                   for names in node_pair_names.values()],
                 'Total secure key rate',
                ]
    #Save jitter values
    base_jitter= []
    for node_pair, names in node_pair_names.items():
            for node in node_pair:
                for detector in node.detectors:
                    base_jitter.append(detector.jitter)
                    
    
    def get_data(distance):
        transmission =  10**(-lossperkm *distance /10)
        dispersion = distance*bandwidth*chromaticdispersion
        
        i=0
        for node_pair, names in node_pair_names.items():
            for node in node_pair:
                for detector in node.detectors:
                    detector.jitter = convolve_Jitters(base_jitter[i],dispersion) #always calculate based on saved values not on previous values.
                    i = i + 1
        for node_pair, names in node_pair_names.items():
            for i, source in indexed_sources[node_pair]:
                optimum = network.optimize_coincidence_window(
                            node_pair, source=source, ubound=ubound)
        
        
        for node in nodes:
            node.transmissions = (
                    transmission, transmission, transmission, transmission)
        data = {}
        data['Distance (km)'] = distance
        data['chromatic dispersion (ps)'] = dispersion
        data['Transmission at nodes'] = transmission
        for node, name in node_names.items():
            data['Singles rate for ' + name] = (
                    network.total_singles_rate(node))
        for node_pair, names in node_pair_names.items():
            data['Pair rate between ' + names] = (
                    network.total_pair_rate(node_pair))
            data['Shifted key rate between ' + names] = (
                    network.shifted_key_rate(node_pair))
            data['Accidentals between ' + names] = (
                    network.accidentals(node_pair))
            data['QBER between ' + names] = (
                    network.QBER(node_pair))
            for i, source in indexed_sources[node_pair]:
                optimum = network.optimize_coincidence_window(
                            node_pair, source=source, ubound=ubound)
                data['Optimal coincidence window for ' + names +
                     ', Source ' + str(i)] = optimum[0]
            data['Secure key rate between ' + names] = optimum[1]
        data['Total secure key rate'] = network.total_secure_key_rate()
        return data

    with open(output_file_prefix + '.csv', 'w', newline = '') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=variables)
        writer.writeheader()
        distance = 0
        writer.writerow(get_data(distance))
        while (max([network.secure_key_rate(node_pair)
                for node_pair in node_pairs]) > 0):
            distance = distance + Distance_step_size
            data = get_data(distance)
            writer.writerow(data)
 
    #restore daved jitter values
    i=0
    for node_pair, names in node_pair_names.items():
            for node in node_pair:
                for detector in node.detectors:
                    detector.jitter = base_jitter[i] 
                    i=i+1
    
    
    
    return variables






































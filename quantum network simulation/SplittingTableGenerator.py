# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 18:14:14 2019

@author: xv18766
"""

import CSVReader
import sys
import os
import numpy as np


class SplittingDict(dict):
    """ """

    @property
    def number_of_indices(self):
        return self._number_of_indices

    @number_of_indices.setter
    def number_of_indices(self, value):
        if int(value) == value:
            self._number_of_indices = value
        else:
            raise ValueError("value should be an integer")
        return self._number_of_indices

    def count_number_of_indices(self):
        indices = set()
        for list_of_tuples in self.values():
            for tuple in list_of_tuples:
                indices.add(abs(tuple[0]))
        return len(indices)


def noWDM_splitting_dict(node_transmission_dict, splitting=1, index=1):
    """ """
    # Check the inputs look vaguely correct
    if not 0 < splitting <= 1:
        raise ValueError("splitting must be >0 and <=1.")
    if not isinstance(index, int):
        raise TypeError("index must be a strictly positive integer.")
    if not index > 0:
        raise ValueError("index must be a strictly positive integer.")
    
    # Make tuples of the required form.
    # Divide the initial level of splitting equally between all nodes.
    def source_tuples(transmission):
        return [(index, splitting/len(node_transmission_dict), transmission),
                (-index, splitting/len(node_transmission_dict), transmission)]

    # Final dictionary, populated with pairs of 3-tuples
    splitting_dict = SplittingDict({
            node: source_tuples(node_transmission_dict[node])
            for node in node_transmission_dict
            })
    splitting_dict.number_of_indices = 1
    return splitting_dict


def WDM_splitting_dict(node_transmission_dict, splitting=1, start_index=1):
    """Return a splitting dictionary and number of indices used.
    ({node: set((source_index, splitting, transmission))}, indices)

    node_transmission_dict must be a dictionary of the form
    {node: transmission}. The transmission will be used as the
    default fibre transmission (3rd entry) of sources to that node.

    splitting is the same for all sources. It can be less than one, in
    which case the splitting table will be invalid by itself, but could
    form part of a larger splitting table.

    start_index is the first index with which the table is populated.
    It must be a strictly positive integer.
    """

    # Check the inputs look vaguely correct
    if not 0 < splitting <= 1:
        raise ValueError("splitting must be >0 and <=1.")
    if not isinstance(start_index, int):
        raise TypeError("start_index must be a strictly positive integer.")
    if not start_index > 0:
        raise ValueError("start_index must be a strictly positive integer.")

    # Make tuples of the required form.
    def source_tuple(index, transmission):
        return (index, splitting, transmission)

    nodes = list(node_transmission_dict.keys())
    # Final dictionary, to be populated with sets of 3-tuples
    splitting_dict = SplittingDict({node: list() for node in nodes})
    # Auxiliary dictionary, keeps track of which nodes have unused indices
    unused_indices = {node: list() for node in nodes}
    index = start_index
    for i, node in enumerate(nodes):
        transmission = node_transmission_dict[node]
        for previous_node in nodes[0:i]:
            # Pop an index from the previous_node to form a connection
            negative = -unused_indices[previous_node].pop()
            splitting_dict[node].append(source_tuple(negative, transmission))
        for next_node in nodes[i+1:]:
            # Add a tuple to the dictionary with the current index
            splitting_dict[node].append(source_tuple(index, transmission))
            # Note that the index has been used with this node
            unused_indices[node].append(index)
            index += 1
        # Reverse the list to prepare for popping
        unused_indices[node].reverse()
    number_of_indices = index - start_index
    splitting_dict.number_of_indices = number_of_indices
    return splitting_dict


def signal_split_splitting_dict(node_transmission_dict,
                             splitting=1, start_index=1):
    """ """
    # Check the inputs look vaguely correct
    if not 0 < splitting <= 1:
        raise ValueError("splitting must be >0 and <=1.")
    if not isinstance(start_index, int):
        raise TypeError("start_index must be a strictly positive integer.")
    if not start_index > 0:
        raise ValueError("start_index must be a strictly positive integer.")

    # Create the set of indices which will be used in the table
    number_of_nodes = len(node_transmission_dict)
    indices = set(range(start_index, start_index+number_of_nodes))

    # 
    def source_tuples(index, transmission):
        tuples = [(i, splitting/(number_of_nodes-1), transmission)
                  for i in indices.difference({index})]
        tuples.append((-index, splitting, transmission))
        return tuples

    index = start_index - 1
    splitting_dict = SplittingDict({})
    for i, item in enumerate(node_transmission_dict.items()):
        index += 1
        splitting_dict[item[0]] = source_tuples(index, item[1])
    splitting_dict.number_of_indices = number_of_nodes
    return splitting_dict

def log2_N_splitting_dict(node_transmission_dict,
                             splitting=1, start_index=1):
    """Return a splitting dict using the log_2-architecture."""
    # Check the inputs look vaguely correct
    if not 0 < splitting <= 1:
        raise ValueError("splitting must be >0 and <=1.")
    if not isinstance(start_index, int):
        raise TypeError("start_index must be a strictly positive integer.")
    if not start_index > 0:
        raise ValueError("start_index must be a strictly positive integer.")
    # Check that the number of nodes is a power of 2
    number_of_nodes = len(node_transmission_dict)
    log2_N = np.log2(number_of_nodes)
    if not log2_N == int(log2_N):
        raise ValueError("The number of nodes is not a power of 2")
    log2_N = int(log2_N)

    # Set the splitting (accounting for the base splitting)
    splitting = splitting/(number_of_nodes/2)

    # Make tuples of the required form.
    def source_tuple(index, transmission):
        return (index, splitting, transmission)

    list_of_nodes = list(node_transmission_dict.keys())
    index = start_index - 1
    splitting_dict = SplittingDict({node: [] for node in list_of_nodes})
    for m in range(log2_N):
        index += 1
        # The maths here is convoluted but works
        size_of_two_subgroups = 2**(log2_N - m)
        for i, node in enumerate(list_of_nodes):
            trans = node_transmission_dict[node]
            if (i/size_of_two_subgroups)%1 < 0.5:
                splitting_dict[node].append(source_tuple(index, trans))
            else:
                splitting_dict[node].append(source_tuple(-index, trans))
    splitting_dict.number_of_indices = log2_N
    return splitting_dict

def dummy_node_transmission_dict(size, transmission=1):
    """Return a node_transmission_dict of the desired size."""
    return {i: transmission for i in range(1, size+1)}

def splitting_table_from_splitting_dict(splitting_dict):
    table = []
    for k, v in splitting_dict.items():
        table.append([k, *v])
    return table


def import_splitting_dict(csv_file, nodes):
    sd = SplittingDict(CSVReader.import_splitting_dict(csv_file, nodes))
    sd.number_of_indices = sd.count_number_of_indices()
    return sd
def import_FAKE_splitting_dict(nodes):
    sd = SplittingDict(CSVReader.import_FAKE_splitting_dict(nodes))
    sd.number_of_indices = sd.count_number_of_indices()
    return sd


def import_splitting_table(csv_file, nodes):
    d = import_splitting_dict(csv_file, nodes)
    t = splitting_table_from_splitting_dict(d)
    return t


def csv_from_splitting_dict(csv_file, splitting_dict):
    r = CSVReader.write_dict_to_csv(csv_file, splitting_dict)
    return r

def dir_and_csv_from_splitting_dict(csv_file, splitting_dict):
    os.mkdir(csv_file)
    os.chdir(csv_file)
    csv_file = csv_file + '.csv'
    CSVReader.write_dict_to_csv(csv_file, splitting_dict)
    os.chdir('..')
    return csv_file


if __name__ == "__main__":
#    nodes = ['Alice', 'Bob', 'Chloe', 'Dave']
#    ntdict = dict.fromkeys(nodes, 0.2) 
#    s = WDM_splitting_dict(ntdict)
#    for k, v in s.items():
#        print(k, v)
#    print("\n")
#
#    ns = noWDM_splitting_dict(ntdict, index=3)
#    for k, v in ns.items():
#        print(k, v)
#    print("\n")
#
#    os = signal_split_splitting_dict(ntdict, start_index=2)
#    for k, v in os.items():
#        print(k, v)
#    print("\n")
#
#    sd = import_splitting_dict('4p WDM only.csv', nodes)
#    for k, v in sd.items():
#        print(k, v)
#    print(sd.number_of_indices)

    numarg=len(sys.argv)
    if numarg <2:
        print ("Usage: \n SplittingTableGenerator -f [Input file name with full path] \n creates a spiltting table from the file and passes it to XXXXXX.py")
    else:
        for i in range(numarg):
            if sys.argv[i] == '-f' or '-F':
                try:
                    inputfile = sys.argv[i+1]
                except:
                    break;
        print (inputfile)    

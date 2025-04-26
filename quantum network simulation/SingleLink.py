# -*- coding: utf-8 -*-
"""
This is a single function that takes a long list of input parameters and returns the key rate under those circumstances. 
The program is massively overkill for this purpouse but it has been adapted from a complete network simulator.
"""
import os
import csv
import SplittingTableGenerator as STG #Self written script for teh network simulation
import LinkModel5 as lm #Self written script for teh network simulation
import DataGenerator5 as DG #Self written script for teh network simulation
import GraphGenerator as GG #Self written script for teh network simulation
import numpy as np
from scipy import optimize
import itertools
import re
import math

def Calc_QKD_Parameters_For_Single_Link( 
DetectorDark_per_second = [1000,1000,1000,1000], #the dark counts of all 4 detectors. the first two for the first user
DetectorEfficiency = [0.7,0.7,0.7,0.7], #efficiency of all 4 detectors. 
NodeAdditionalDetectorNoise_HA_VD_per_second = [0,0,0,0],
NodeBasisTransmission_HAVD = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5],
NodeBasisChoiceProbability = [0.5,0.5],
TimeTaggerJitter_in_ps = [40,40], 
FiberLoss_in_dB_per_Km = [0.2,0.2], 
FiberDispersion_in_ps_per_nm_per_km = [16.7,16.7], 
CompensatedDispersion_in_ps = [0,0], 
FibreLength_in_km = [0.02,0.02], 
OtherTransmission = [0.75,0.75],
OtherLoss_in_dB = [6,6], #OtherTransmission and other loss in dB are equivalent. 
SingleCountsFromOtherUsers_per_second = [30000,30000,30000,30000],  #Additional single counts on each of the two detectors HA and VD for both users
SingleCountsFromClassicalCoexistance_per_second = [0,0,0,0], #When using classical co-existance in the same fibre the additional counts caused in each of the detectors due to this. If not using classical quantum co-existance this is 0 for all 4 values
BeamSplittingRatio = [1,1], 
DetectorParalyzable = [True, True],
DetectorJitter_in_ps=40, #Although there is a seperate jitter for each pair of detectors, we model a single link assuming a single value.
TimeSyncJitter_in_ps = 60,
CoincidenceWindow_in_ps = 120,#change this value
PairGenerationRate_per_second = 1E6, #2E5 to 5E6
SourceHeraldingEfficiency = 0.2, 
centralWavelength_in_nm = 1550.12, 
SignalBandwidth_in_nm = 0.8,
SourceEntangledStateParameterTheta = math.pi/4, #The source will produce a state sin(theta) |HH> + e^(i Phi) cos(Theta) |VV>. Set theta =pi/4 and phi =0 for a perfect Bell state. Alter this for imperfect entanglememt.
SourceEntangledStateParameterPhi = 0
):
    """NOTE:
    In any realistic scenario several conditions must be manually met for the input of this function to be physically meaningful
        1) SingleCountsFromOtherUsers_per_second should accurately reflect the number of other links connected to the same users. As the source pair rate is increased or decreased this should change by a similar factor. When a single user has highly lossy links and normal links then computing these values can be involved.
        2) Compensating for more jitter than actually caused can be a very bad thing.
        3) the BeamSplittingRatio is only implemented here as an additional loss and checks to ensure that the sum of all splitting ratios for each wavelength channel are disabled by default.
        4) This simplified simulation no longer keeps track of the wavelength channels being used. It only uses 1 and -1 for all calculations.
        5) The detector jitter is actually the g2 jitter and is assumed to be the same for all 4 pairs of detectors between any two given nodes.
        6) There are no checks on the coincidence window used. Please make sure that teh values are reasonable. Like the jitter, the coincidence window is a single value. assumed to be the same for all 4 pairs of detectors.
        7) There are no checks to ensure that key is generated from both basis choices.  
    """
    extraNoise = np.zeros(len(NodeAdditionalDetectorNoise_HA_VD_per_second))
    for i in range(len(NodeAdditionalDetectorNoise_HA_VD_per_second)):
        extraNoise[i]=NodeAdditionalDetectorNoise_HA_VD_per_second[i] + SingleCountsFromClassicalCoexistance_per_second[i] + SingleCountsFromOtherUsers_per_second[i]
    
    number_of_nodes = 2
    
    dispersion = np.zeros(number_of_nodes)
    linkTransmission = np.zeros(number_of_nodes)
    p1 = [(0,0,0)] #will be used later on to ammend the splitting table
    p2 = [(0,0,0)]
    for i in range (number_of_nodes):
        dispersion[i] = FiberDispersion_in_ps_per_nm_per_km[i] * SignalBandwidth_in_nm * FibreLength_in_km[i] - CompensatedDispersion_in_ps[i]
        linkTransmission[i] = 10**(-FiberLoss_in_dB_per_Km[i]/10 * FibreLength_in_km[i])*OtherTransmission[i] *10**(-OtherLoss_in_dB[i]/10)
    p1 = [(1,float(BeamSplittingRatio[0]),float(linkTransmission[0]))]
    p2 = [(-1,float(BeamSplittingRatio[1]),float(linkTransmission[1]))]
        
    coherenceTime = ((centralWavelength_in_nm*1E-9)**2 )/(299792458*SignalBandwidth_in_nm*1E-9)
    jitter = convolve_Jitters(DetectorJitter_in_ps, *TimeTaggerJitter_in_ps, TimeSyncJitter_in_ps, coherenceTime, *dispersion)
    
    
    detectors = [(lm.Detector(paralyzable=DetectorParalyzable[i],dark_rate=DetectorDark_per_second[2*i], jitter=jitter,base_efficiency=DetectorEfficiency[2*i]),
                      lm.Detector(paralyzable=DetectorParalyzable[i],dark_rate=DetectorDark_per_second[2*i+1], jitter=jitter,base_efficiency=DetectorEfficiency[2*i+1]))
                        for i in range(number_of_nodes)]
                        
                        
    names = ['Alice', 'Bob', 'Chloe', 'Dave', 'Fay', 'Gopi', 'Heidi', 'Ivan',
                 'Jess', 'Kevin', 'Lisa', 'Mat', 'Nora', 'Omar', 'Puja', 'Quinn',] #The names are a legacy thing from the network simulation. 
    nodes = [lm.Node(detectors=detectors[i], transmissions=(NodeBasisTransmission_HAVD[4*i],NodeBasisTransmission_HAVD[4*i+1],NodeBasisTransmission_HAVD[4*i+2],NodeBasisTransmission_HAVD[4*i+3]),
                 noise_rate=(extraNoise[2*i] , extraNoise[2*i+1]), basis_choice_HVtoAD=NodeBasisChoiceProbability[i], name=name)
                 for i, name in enumerate(names[0:number_of_nodes])]
    # Since this is a legacy from the network simulation lets just import a 2 node network. The csvfile is manually coaded for a 2 user network with only 1 wavelength channel.
    indices_of_interesting_nodes = range(number_of_nodes)
    splitting_dict = STG.import_FAKE_splitting_dict(nodes)
    #Now I replace the values from the CSV file with real ones!
    node1 = {nodes[0]:p1}
    splitting_dict.update(node1)
    node2 = {nodes[1]:p2}
    splitting_dict.update(node2)
    
    number_of_sources = splitting_dict.number_of_indices
    
    
    sources = [lm.Source(theta=SourceEntangledStateParameterTheta, phi=SourceEntangledStateParameterPhi,
                 pair_rate=PairGenerationRate_per_second, heralding_efficiency=SourceHeraldingEfficiency)
                for _ in range(number_of_sources)]
    
    
    
    splitting_table = STG.splitting_table_from_splitting_dict(splitting_dict)
    splitting_table = lm.SplittingTable(splitting_table, sources,VerifySplittingRatios =False)
    network = lm.Network(splitting_table)
    
    
    # another net sim legacy. just choose the first two nodes
    node_pairs = [(nodes[0], nodes[1])]
    all_node_pairs = [(nodes[x], nodes[y]) for x, y in 
                      itertools.combinations(range(number_of_nodes), 2)]
    node_pair_names = {node_pair: node_pair[0].__str__() + ' and '
                               + node_pair[1].__str__() for node_pair in node_pairs}
    nodes_of_interest = list(nodes[i] for i in indices_of_interesting_nodes)
    node_names = {node: node.__str__() for node in nodes_of_interest}
    #### This section saves parameters so that they can be reset after each output

    pair_rates = {source: source.pair_rate for source in sources}
    jitters = {detector: detector.jitter
               for pair in detectors
               for detector in pair}
    transmissions = {node: node.transmissions
                     for node in nodes_of_interest}
    
    network.set_coincidence_windows(node_pairs[0],value = CoincidenceWindow_in_ps)
    
    #Very simply return the key rate without optimising the coinc Window.
    if True:
            data = DG.transmission_vs(
                    network, 
                    nodes_of_interest, all_node_pairs,
                    node_transmission= 1) 
            # Reset values
            for node in nodes_of_interest:
                node.transmissions = transmissions[node]
            
            #print('QBER between Alice and Bob:',data['QBER between Alice and Bob'])
            #print('Shifted key rate between Alice and Bob:', data['Shifted key rate between Alice and Bob'])
            #print("Secure key rate:",data['Total secure key rate'])  
    data.update({'Jitter (ps)':jitter})
    data.update({'Link transmission (dB)': (10*np.log10(linkTransmission[0]),10*np.log10(linkTransmission[1]))})
    data.update({'Dispersion (ps)':(dispersion[0],dispersion[1])})
    return data

def convolve_Jitters(*jitters):
        sum_of_squares = 0
        for j in jitters:
            sum_of_squares = sum_of_squares + j*j
        convolvedJitter = np.sqrt(sum_of_squares)
        return (convolvedJitter)    





if __name__ == "__main__":
    a=Calc_QKD_Parameters_For_Single_Link(FibreLength_in_km = [0.02,0.02], BeamSplittingRatio = [0.25,0.25])
    print(a)
    max_secure_rate=[0,0]
    for cw in range(10,2000,1):
        #a=Calc_QKD_Parameters_For_Single_Link(CoincidenceWindow_in_ps = cw )
        #print(a['Total secure key rate'])
        if a['Total secure key rate']>max_secure_rate[0]:
            max_secure_rate[0]=a['Total secure key rate']
            max_secure_rate[1]=cw
    #print('max secure key rate:',max_secure_rate[0],'correspond conincidence window(ps):',max_secure_rate[1] )
    #print(Calc_QKD_Parameters_For_Single_Link())#max_secure_rate[1]))
    
    
    # ~ print (Calc_QKD_Parameters_For_Single_Link(FibreLength_in_km = [0100.02,0.02], BeamSplittingRatio = [1,1], 
    # ~ SingleCountsFromClassicalCoexistance_per_second = [0,0,0,0],PairGenerationRate_per_second = 6e6,OtherLoss_in_dB=[8,4],
    # ~ CompensatedDispersion_in_ps =[200,0]))
    singlecounts=Calc_QKD_Parameters_For_Single_Link\
    ( FibreLength_in_km = [10,10],PairGenerationRate_per_second = 1E6,
     SingleCountsFromOtherUsers_per_second=[0,0,0,0], 
    OtherLoss_in_dB = [9,9])
    
    ideal_single_link = singlecounts['Singles rate for Alice'].HA / 0.7
    #print(ideal_single_link)
    






























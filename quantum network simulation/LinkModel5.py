# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 16:53:04 2019

@author: xv18766 Alex Qui worte the original script as part of the network sim for project A modifications for long distance by Siddarth
"""

import numpy as np
import scipy.integrate as integrate
import math
import cmath
import matplotlib.pyplot as plt
import abstractlinks as al
import itertools
from scipy import optimize


class Source(al.AbstractSource):
    """General Source object.

    Public methods:
    state_probs -- settable property, return an 8-tuple of state
        probabilities, where the first four (HH, VV, HV, VH) are
        assuming a measurement in the HV basis and the last four (AA,
        DD, AD, DA) assume a measurement in the AD basis.
    pair_rate -- settable property, return the rate of pairs generated
        by the source. Automatically recalculate the singles_rate.
    heralding_efficiency -- settable property, return the heralding
        efficiency (ratio of pairs to singles) of the source.
        Automatically recalculate the singles_rate.
    singles_rate -- read-only property, set automatically using
        pair_rate or heralding_efficiency methods.
    """

    def __init__(self, theta=math.pi/4, phi=0,
                 pair_rate=1E5, heralding_efficiency=0.4):
        """Constructor for Source object.

        Assumes that the state in the HV basis is in the form
            sin(theta) |HH> + e^(i*phi) cos(theta) |VV>
        and hence in the AD basis is
            0.5*((sin(theta) + e^(i*phi)*cos(theta))(|DD> + |AA>)
                  + (sin(theta) - e^(i*phi)*cos(theta))(|DA> + |AD>)).
        """
        self._theta = theta
        self._phi = phi
        self.state_probs = self._calc_state_probs(theta, phi)
        # HACK: Initialise values so that property setters can be used.
        self._pair_rate = 1
        self._heralding_efficiency = 1
        self.pair_rate = pair_rate
        # singles_rate is automatically set here.
        self.heralding_efficiency = heralding_efficiency

    @staticmethod
    def _calc_state_probs(theta, phi):
        """Calculates the state probabilities from two arguments.

        Assumes that the state in the HV basis is in the form
            sin(theta) |HH> + e^(i*phi) cos(theta) |VV>
        and hence in the AD basis is
            0.5*((sin(theta) + e^(i*phi)*cos(theta))(|DD> + |AA>)
                  + (sin(theta) - e^(i*phi)*cos(theta))(|DA> + |AD>)).
        """
        s = math.sin(theta)
        c = math.cos(theta)
        e = cmath.exp(phi*1j)
        HH = math.pow(s, 2)
        VV = math.pow(c, 2)
        HV = 0
        VH = 0
        AA = math.pow(abs(0.5*(s + e*c)), 2)
        DD = AA
        AD = math.pow(abs(0.5*(s - e*c)), 2)
        DA = AD
        return al.StateProbabilities(HH, VV, HV, VH, AA, DD, AD, DA)

    @property
    def state_probs(self):
        return self._state_probs

    @state_probs.setter
    def state_probs(self, state_probs):
        for prob in state_probs:
            if not 0 <= prob <= 1:
                return ValueError("Entries must be between 0 and 1.")
        self._state_probs = al.StateProbabilities(*state_probs)
        return self._state_probs

    @property
    def pair_rate(self):
        return self._pair_rate

    @pair_rate.setter
    def pair_rate(self, pair_rate):
        """Set pair_rate and recalculate singles_rate."""
        if pair_rate < 0:
            raise ValueError("pair_rate cannot be less than 0.")
        self._pair_rate = pair_rate
        self._singles_rate = self._pair_rate / self._heralding_efficiency
        return self._pair_rate

    @property
    def heralding_efficiency(self):
        return self._heralding_efficiency

    @heralding_efficiency.setter
    def heralding_efficiency(self, heralding_efficiency):
        """Set heralding_efficiency and recalculate singles_rate"""
        if not 0 <= heralding_efficiency <= 1:
            raise ValueError("heralding_efficiency must be between 0 and 1.")
        self._heralding_efficiency = heralding_efficiency
        self._singles_rate = self._pair_rate / self._heralding_efficiency
        return self._heralding_efficiency

    # Use @property without setter to make singles_rate read-only.
    @property
    def singles_rate(self):
        """Read-only attribute, set using pair_rate or heralding_efficiency"""
        return self._singles_rate

    def splice(self, number_of_splices):
        """Return a list of sources with equal pair rates.

        Simulates the division of the Source object into several such
        objects with equally divided pair rates.
        """
        new_pair_rate = self.pair_rate/number_of_splices
        return [Source(self._theta, self._phi, new_pair_rate,
                       self.heralding_efficiency)
                for _ in range(number_of_splices)]


class Detector(al.AbstractDetector):
    """General Detector object.

    Public methods:
    paralyzable -- settable property, return True (False) if the
        detector is (not) paralyzable.
    base_efficiency -- settable property, return the base efficiency
        i.e. not accounting for dead time.
    jitter -- settable property, return detector jitter in picoseconds.
    dead_time -- settable property, return detector dead time in ps.
    dark rate -- settable property, return the rate of dark counts.
    efficiency(self, incidence_rate) --
        return the detector efficiency at a given incidence rate,
        accounting for dead time, dark rate and paralyzability.
    detection_rate(self, incidence_rate) --
        return the detection rate at a given incidence rate, accounting
        for dark counts and the detector efficiency at that incidence
        rate.
    """

    def __init__(self, paralyzable=True, base_efficiency=0.7,
                 jitter=100, dead_time=5000, dark_rate=1000):
        self.paralyzable = paralyzable
        self.base_efficiency = base_efficiency
        self.jitter = jitter
        self.dead_time = dead_time
        self.dark_rate = dark_rate

    @property
    def paralyzable(self):
        return self._paralyzable

    @paralyzable.setter
    def paralyzable(self, paralyzable):
        if paralyzable not in [True, False]:
            raise ValueError("paralyzable must be either True or False.")
        self._paralyzable = paralyzable
        return self._paralyzable

    @property
    def base_efficiency(self):
        return self._base_efficiency

    @base_efficiency.setter
    def base_efficiency(self, base_efficiency):
        if not 0 <= base_efficiency <= 1:
            raise ValueError("base_efficiency must be between 0 and 1.")
        self._base_efficiency = base_efficiency
        return self._base_efficiency

    @property
    def jitter(self):
        return self._jitter

    @jitter.setter
    def jitter(self, jitter):
        if jitter < 0:
            return ValueError("jitter cannot be less than 0.")
        self._jitter = jitter
        return self._jitter

    @property
    def dead_time(self):
        return self._dead_time

    @dead_time.setter
    def dead_time(self, dead_time):
        if dead_time < 0:
            return ValueError("dead_time cannot be less than 0.")
        self._dead_time = dead_time
        return self._dead_time

    @property
    def dark_rate(self):
        return self._dark_rate

    @dark_rate.setter
    def dark_rate(self, dark_rate):
        if dark_rate < 0:
            return ValueError("dark_rate cannot be less than 0.")
        self._dark_rate = dark_rate
        return self._dark_rate

    def _efficiency_p(self, incidence_rate):
        """Return the efficiency assuming a paralyzable detector.

        Calculate the effective efficiency from the incidence rate,
        taking into account the dead time and dark count rate.
        Assumes a paralyzable detector.
        """
        efficiency = (
            self.base_efficiency
            * math.exp(-self.dead_time*1E-12*incidence_rate))
        return efficiency

    def _efficiency_np(self, incidence_rate):
        """Return the efficiency assuming a non-paralyzable detector.

        Calculate the effective efficiency from the incidence rate,
        taking into account the dead time and dark count rate.
        Assumes a nonparalyzable detector.
        """
        efficiency = (
            self.base_efficiency
            * (1/(1 + incidence_rate*self.dead_time*1E-12)))
        return efficiency

    def efficiency(self, incidence_rate):
        """Return the detector efficiency given the incidence rate."""
        # Account for dark_rate here.
        incidence_rate += self.dark_rate
        if self.paralyzable is True:
            return self._efficiency_p(incidence_rate)
        elif self.paralyzable is False:
            return self._efficiency_np(incidence_rate)
        else:
            raise ValueError("self.paralyzable has been mis-set.")

    def detection_rate(self, incidence_rate):
        """Return the detection rate given the incidence rate.

        Accounts for dark rate and detector efficiency.
        """
        # Can't do an in-place addition because we need the raw
        # incidence_rate to plug into the efficiency method.
        total_rate = incidence_rate + self.dark_rate
        self._detection_rate = total_rate * self.efficiency(incidence_rate)
        return self._detection_rate


class Node(al.AbstractNode):
    """Node object with two combination (HA/VD) detectors.

    Public methods:
    detectors -- settable property, return the node detectors as a
        2-namedtuple (HA, VD).
    jitters -- read-only property, return the jitters of the node
        detectors as a 2-namedtuple (HA, VD).
    transmissions -- settable property, return the transmission of the
        four polarization states as a 4-namedtuple (H, A, V, D).
    noise_rate -- settable property, return the rate of noise photons
        incident on each node detector as a 2-namedtuple (HA, VD).
    basis_choice_HVtoAD -- settable property, return the probability
        of choosing to measure in the HV as opposed to AD basis (assume
        the total probability sums to 1).
    total_transmissions -- read-only property, return a 4-namedtuple
        (H, A, V, D) of functions which, when called with an incidence
        rate, return the total transmission of a photon of the given
        polarization accounting for that incidence rate on the detector
        corresponding to that photon's polarization.
    detection_rates -- read-only property, return a 2-namedtuple
        (HA, VD) of functions which, when called with an incidence
        rate, return the detection rate on that detector.
    """

    def __init__(self, detectors=(Detector(), Detector()),
                 transmissions=(0.02, 0.02, 0.02, 0.02),
                 noise_rate=(0, 0), basis_choice_HVtoAD=0.5, name=None):
        """Construtor for Node object.

        detectors parameter must be a 2-tuple of AbstractDetector
        objects.
        """
        self.detectors = detectors
        self.transmissions = transmissions
        self.noise_rate = noise_rate
        self.basis_choice_HVtoAD = basis_choice_HVtoAD
        self.name = name

    def __str__(self):
        if self.name == None:
            return self.__repr__()
        else:
            return str(self.name)

    @property
    def detectors(self):
        return self._detectors

    @detectors.setter
    def detectors(self, detectors):
        for detector in detectors:
            if not isinstance(detector, al.AbstractDetector):
                return TypeError("Entries must be AbstractDetector objects.")
        self._detectors = al.DetectorPair(*detectors)
        return self._detectors

    @property
    def jitters(self):
        jitter_HA = self.detectors.HA.jitter
        jitter_VD = self.detectors.VD.jitter
        return al.DetectorPair(jitter_HA, jitter_VD)

    @property
    def transmissions(self):
        return self._transmissions

    @transmissions.setter
    def transmissions(self, transmissions):
        for transmission in transmissions:
            if not 0 <= transmission <= 1:
                raise ValueError("transmissions must be between 0 and 1.")
        self._transmissions = al.TransmissionHAVD(*transmissions)
        return self._transmissions

    @property
    def noise_rate(self):
        return self._noise_rate

    @noise_rate.setter
    def noise_rate(self, noise_rate):
        for noise in noise_rate:
            if noise < 0:
                raise ValueError("Entries cannot be less than 0.")
        self._noise_rate = al.DetectorPair(*noise_rate)
        return self._noise_rate

    @property
    def basis_choice_HVtoAD(self):
        return self._basis_choice_HVtoAD

    @basis_choice_HVtoAD.setter
    def basis_choice_HVtoAD(self, basis_choice_HVtoAD):
        if not 0 <= basis_choice_HVtoAD <= 1:
            raise ValueError("basis_choice_HVtoAD must be between 0 and 1.")
        self._basis_choice_HVtoAD = basis_choice_HVtoAD
        return self._basis_choice_HVtoAD

    # This next bit is messy but the idea is to make total_transmission
    # return a tuple of functions which give the transmission of a HAVD
    # polarized photon (inc. detection prob) given an incidence rate.
    # Ideally these would be internal function definitions but that
    # doesn't work for some reason.
    def _trans_H(self, incidence_rate):
        # Add noise to all incidence rates
        incidence_rate += self.noise_rate.HA
        return (self.transmissions.H
                * self.detectors.HA.efficiency(incidence_rate))

    def _trans_A(self, incidence_rate):
        incidence_rate += self.noise_rate.HA
        return (self.transmissions.A
                * self.detectors.HA.efficiency(incidence_rate))

    def _trans_V(self, incidence_rate):
        incidence_rate += self.noise_rate.VD
        return (self.transmissions.V
                * self.detectors.VD.efficiency(incidence_rate))

    def _trans_D(self, incidence_rate):
        incidence_rate += self.noise_rate.VD
        return (self.transmissions.D
                * self.detectors.VD.efficiency(incidence_rate))

    @property
    def total_transmissions(self):
        """Transmissions accounting for detector efficiency."""
        return al.TransmissionHAVD(
                self._trans_H, self._trans_A, self._trans_V, self._trans_D)

    # Similar to above, this code allows self.detection_rates to create
    # a tuple of functions, one for each detector, which each calculate
    # the detection rate at that detector given an incidence rate on it.
    def _det_rate_HA(self, incidence_rate):
        # Add noise to incidence rates
        incidence_rate += self.noise_rate.HA
        return self.detectors.HA.detection_rate(incidence_rate)

    def _det_rate_VD(self, incidence_rate):
        # Add noise to incidence rates
        incidence_rate += self.noise_rate.VD
        return self.detectors.VD.detection_rate(incidence_rate)

    @property
    def detection_rates(self):
        return al.DetectorPair(self._det_rate_HA, self._det_rate_VD)


class SplittingTable(al.AbstractSplittingTable):
    """Interprets splitting tables using 'number formatting'.

    This class interprets +ve/-ve pairs of numbers (e.g. 1, -1) as
    inputs to nodes from the source indexed by that number, with the
    two signs representing the two arms of each pair. It requires as
    input a splitting table of source:share:fibre_transmission ratios
    and a separate list of appropriately indexed sources.

    !!! Do NOT input the same node twice, or the same source twice for a
    given node !!!

    Public methods:
    source_node_links(self) --
        return a dictionary encoding nodes, sources, and transmissions.
        Note that so long as the _table and _sources private attributes
        haven't been modified since initialisation, it is preferable
        to look up the attribute _source_node_links rather than call
        this method again.
    node_combinations_dictionary(self) --
        return a valueless dictionary with node pairs as keys.
    correlations(self) --
        return a dictionary of sources with node pairs as keys.
    singles_rate(self, source, node) --
        return the singles_rate from source arriving at node.
    sources -- read-only property, return the list of sources.
    pair_rate(self, source, node) --
        return the pairs from source which are sent to node.
    """

    def __init__(self, table, sources,VerifySplittingRatios = True):
        """Constructor utilising splitting table and list of sources.

        Splitting table ('table') must be of the form
        [[node1, (1, .5, .2), (2, .5, .2), (3, .5, .2)],
         [node2, (-1, .5, .2), (2, .5, .2), (-3,.5, .2)]
         [node3, (1, .5, .2), (-2, .5, .2), (-3, .5, .2)]
         [node4, (-1, .5, .2), (-2, .5, .2), (3, .5, .2)]]
        to be correctly interpreted, where the first entry of each list
        is a unique Node object, the (absolute value of the) first
        entry of each tuple represents the natural index of the Source
        object in the list of sources, the second entry of each tuple
        is the proportion of pairs sent from the source to the node,
        and the third entry is the fibre transmission to that node.
        Each source can only be assigned to a given node once
        e.g. node1 cannot have (1, .25, .2), (1, .25, .2), because the
        pair_share() method will fail to account for multiple
        assignations.
        """
        self._table = table
        self._nodes = self.find_nodes(table)
        self._sources = sources
        if VerifySplittingRatios:
            self._check_splittings_sum_to_one()
        self._correlations = self.correlations()
        self._source_node_links = self.source_node_links()

    @staticmethod
    def find_nodes(table):
        nodes = set(row[0] for row in table)
        return nodes

    def _check_splittings_sum_to_one(self):
        splittings = {}
        for index, source in enumerate(self._sources):
            index += 1
            splittings[index] = []
            splittings[-index] = []
            for row in self._table:
                for triplet in row[1:]:
                    if abs(triplet[0]) == index:
                        splittings[triplet[0]].append(triplet[1])
        sum_of_splittings = {k: sum(v) for k, v in splittings.items()}
        for k, v in sum_of_splittings.items():
            if v != 1:
                print("Source " + str(k) + "'s splittings don't sum to 1.")
        return sum_of_splittings

    def source_node_links(self):
        """Return a dictionary encoding nodes, sources, transmissions.

        Return a dictionary object with nodes as keys and dictionaries
        as values. The value-dictionaries themselves have sources as
        keys and tuples of the form (splitting, fibre_transmission) as
        values.

        The returned object can be iterated over to return
        dictionaries, which can themselves be iterated over to return
        sources
            e.g. for source in self.source_node_links()[node]
        while the tuples can be accessed using the appropriate source
            e.g. splitting = self.source_node_links()[source][0]
                 transmission = self.source_node_links()[source][1].
        """
        source_node_links = {}
        for row in self._table:
            # row[0] is a node object
            source_node_links[row[0]] = {}
            for triplet in row[1:]:
                # Abs value of first entry of triplet indexes source
                # Compensate for zero-indexing of self._sources
                index = abs(triplet[0]) - 1
                source = self._sources[index]
                transmissions = (triplet[1], triplet[2])
                source_node_links[row[0]].setdefault(source, transmissions)
        return source_node_links

    @staticmethod
    def find_nodes_linked_to_source(source_index, table):
        """Return the nodes in table linked to source_index.

        source_index is parity sensitive.
        """
        nodes = set(row[0] for row in table
                    for triplet in row[1:]
                    if triplet[0] == source_index)
        return nodes

    @classmethod
    def make_node_combinations_dictionary(cls, table):
        """Return a valueless dictionary with node pairs as keys.

        Requires a properly formatted splitting table as input.
        """
        dictionary = {}
        node_combinations = itertools.combinations_with_replacement(
                cls.find_nodes(table), 2)
        for node_combo in node_combinations:
            dictionary.setdefault(node_combo)
            dictionary.setdefault(node_combo[-1::-1])
        return dictionary

    def all_node_pairs(self):
        """Return all the possible pairs of nodes, orderless."""
        node_combinations = itertools.combinations(
                self.find_nodes(self._table), 2)
        return list(node_combinations)

    def node_combinations_dictionary(self):
        """Return a valueless dictionary with node pairs as keys.

        The dictionary can then be populated with sources (as in
        self.correlations()) or coincidence windows, etc.
        """
        return self.make_node_combinations_dictionary(self._table)

    def correlations(self):
        """Return a dictionary of sources with node pairs as keys.

        Dictionary keys are (node1, node2) tuples. Also includes the
        reversed tuples.
        Dictionary values are sets of sources which transmit signal
        photons to one node and idler photons to the other node.
        """
        # Make a valueless dictionary with the appropriate keys
        correlations = self.node_combinations_dictionary()
        # Initialise the dictionary values as sets, then add sources
        for k in correlations:
            correlations[k] = set()
        for source_index, source in enumerate(self.sources):
            source_index += 1
            # Make sure nodes are only connected by +/- index combos
            positives = self.find_nodes_linked_to_source(source_index,
                                                         self._table)
            negatives = self.find_nodes_linked_to_source(-source_index,
                                                         self._table)
            node_pairs = itertools.product(positives, negatives)
            for node_pair in node_pairs:
                correlations[node_pair].add(source)
                correlations[node_pair[-1::-1]].add(source)
        return correlations

    def singles_rate(self, source, node):
        """Return the singles_rate from source arriving at node.

        Accounts for splitting ratios and fibre transmissions.
        """
        # 'sources' is itself a dictionary with source objects as keys
        sources = self._source_node_links[node]
        if source not in sources:
            return 0
        # Multiply by the two transmissions associated with source
        return source.singles_rate * sources[source][0] * sources[source][1]

    @property
    def sources(self):
        """Return the list of sources supplied to the splitting table."""
        return self._sources

    def pair_rate(self, source, nodes):
        """Return the pair_rate from source arriving at nodes.

        Accounts for splitting ratios and fibre transmissions in both
        arms.
        """
        if source not in self._correlations[nodes]:
            return 0
        eff_node1 = (self._source_node_links[nodes[0]][source][0]
                     * self._source_node_links[nodes[0]][source][1])
        eff_node2 = (self._source_node_links[nodes[1]][source][0]
                     * self._source_node_links[nodes[1]][source][1])
        return source.pair_rate*eff_node1*eff_node2


class Network(al.AbstractNetwork):
    """Network object for two-detector nodes.

    Public methods:
    self.coincidence_windows -- dictionary with 2-tuples of node pairs
        as keys and coincidence windows between those pairs as values.
    singles_rate(self, source, node) --
        return the singles_rate from source arriving at node.
    pair_rate(self, source, node) --
        return the pairs from source which are sent to node.
    incidence_rates(self, node) --
        return the incidence rates at the two detectors in node.
    detection_rates(self, node) --
        return a DetectorPair namedtuple of detection_rates at node.
    accidentals(self, nodes) --
        return the accidental coincidence rate between nodes.
    shifted_key_rate(self, nodes) --
        return the shifted key rate genereated between nodes.
    QBER(self, nodes[, F]) --
        return the QBER between nodes. Optional argument F is the
        proportion of bits that are wrong for a miscellaneous reason
        (default F=0.01).
    secure_key_rate(self, nodes[, QBER, individual]) --
        return the secure key rate generated between nodes. QBER can be
        set manually or (by default) calculated using the method.
        individual parameter is False by default. If True, then Eve is
        assumed to be limited to individual attacks and therefore will
        have lower mutual information with Alice than otherwise.
    """

    def __init__(self, splitting_table):
        """Constructor for Network object.

        Requires a SplittingTable object as parameter.
        """
        self.splitting_table = splitting_table
        # Dictionary of all (node1, node2) pairs with corresponding
        # set of sources that link those pairs.
        # Also contains all nodes as keys with corresponding sources
        # which transmit to those nodes.
        self.sources = splitting_table.source_node_links()
        self.sources.update(splitting_table.correlations())
        # Dictionary of all (node1, node2) pairs where the values are
        # themselves dictionaries, with sources as keys and coincidence
        # windows as values
        self.coincidence_windows = (
                splitting_table.node_combinations_dictionary())
        for k in self.coincidence_windows:
            self.coincidence_windows[k] = {
                    source: self.default_coincidence_window(k)
                        for source in self.sources[k]}
        # Dictionary of all (node1, node2) pairs with value = dictionary
        # of sources as keys with value = probability of
        # coincidence between each photon polarization possibility.
        self.coincidence_probabilities = (
                splitting_table.node_combinations_dictionary())
        for k in self.coincidence_probabilities:
            self.coincidence_probabilities[k] = self.calc_cp_between_nodes(k)
        self.all_node_pairs = splitting_table.all_node_pairs()

    @staticmethod
    def default_coincidence_window(nodes):
        """Return the default coincidence window between two nodes.

        Used to initialise self.coincidence_windows dictionary.
        """
        all_jitters = nodes[0].jitters + nodes[1].jitters
        return 2*max(all_jitters)*1E-12

    def set_coincidence_windows(self, nodes, source="all", value="default"):
        """Safely change the coincidence window between nodes.

        Will check that the nodes exist and the value is valid.
        Assumes that value is in picoseconds and compensates for that.
        """
        # Check the nodes are already in the dictionary, else refuse.
        if not self.coincidence_windows.__contains__(nodes):
            raise ValueError("nodes must already be a valid key.")
        if value == "default":
            value = self.default_coincidence_window(nodes)
        if value < 0:
            raise ValueError("value must be positive.")        

        # Translate from picoseconds to seconds.
        value *= 1E-12
        reverse = nodes[-1::-1]
        cw = self.coincidence_windows[nodes]
        cwr = self.coincidence_windows[reverse]
        if source == "all":
            for source in cw:
                cw[source] = value
                cwr[source] = value
        else:
            cw[source] = value
            cwr[source] = value
        # Automatically set the coincidence probability
        self.coincidence_probabilities[nodes] = (
                self.calc_cp_between_nodes(nodes))
        self.coincidence_probabilities[reverse] = (
                self.calc_cp_between_nodes(reverse))
        return self.coincidence_windows[nodes]

    @staticmethod
    def convolved_variance(jitter1, jitter2):
        """Return the variance of a convolution of Gaussian PDFs.

        Return the variance of a convolution of two Gaussian
        probability density functions with FWHM jitter1 and jitter2.
        """
        jitter1 *= 1E-12
        jitter2 *= 1E-12
        return (jitter1**2 + jitter2**2)/(8*math.log(2))

    @staticmethod
    def gaussian(var):
        """Return a Gaussian PDF with variance var."""
        def distribution(t):
            return (1/math.sqrt(2*math.pi*var))*math.exp(-(t**2)/(2*var))
        return distribution

    @classmethod
    def coincidence_functions(cls, nodes):
        """Return an 8-tuple of coincidence functions between nodes

        Return a StateProbabilities namedtuple of functions giving the
        coincidence functions between nodes for the eight polarization
        pairs.

        Define a coincidence function C(t) for a given pair of nodes
        and photon polarizations X and Y as: the probability that the
        detector for X photons in the first node will register a count
        at precisely time t before the detector for Y photons in the
        second node registers a count, assuming that there is no noise
        and that if there were no jitter in the system then the photons
        would be detected at precisely the same time.
        """
        # Eventually, G2 = (HH, VV, HV, VH, AA, DD, AD, DA)
        cf = [0, 0, 0, 0, 0, 0, 0, 0]
        j1 = nodes[0].jitters
        j2 = nodes[1].jitters

        # HH and AA terms
        cf[0] = cf[4] = cls.gaussian(cls.convolved_variance(j1.HA, j2.HA))
        # VV and DD terms
        cf[1] = cf[5] = cls.gaussian(cls.convolved_variance(j1.VD, j2.VD))
        # HV and AD terms
        cf[2] = cf[6] = cls.gaussian(cls.convolved_variance(j1.HA, j2.VD))
        # VH and DA terms
        cf[3] = cf[7] = cls.gaussian(cls.convolved_variance(j1.VD, j2.HA))

        # Should rename StateProbabilities to something more general
        return al.StateProbabilities(*cf)

    @staticmethod
    def coincidence_probability(coincidence_function, coincidence_window):
        """Return probability of coincidence given function and window.

        Given a coincidence function and a window about the centre of
        the function, return the integral of the function over the
        window.

        Assumes the window is in seconds.
        """
        return integrate.quad(coincidence_function,
                              -coincidence_window/2,
                              coincidence_window/2)[0]

    @classmethod
    def calc_coincidence_probabilities(cls, functions, window):
        """Return an 8-tuple of coincidence probabilities.

        Return a StateProbabilities namedtuple of values giving the
        probabilities of coincidence between nodes for the eight
        polarization pairs.

        functions must be an 8-tuple.
        Assumes the window is in seconds.
        """
        cp = []
        for f in functions:
            cp.append(cls.coincidence_probability(f, window))
        return al.StateProbabilities(*cp)

    def calc_cp_between_nodes_for_a_source(self, nodes, source):
        functions = self.coincidence_functions(nodes)
        window = self.coincidence_windows[nodes][source]
        
        return self.calc_coincidence_probabilities(functions, window)

    def calc_cp_between_nodes(self, nodes):
        """Return an 8-tuple of coincidence probabilities between nodes

        Calculates the functions directly from the nodes and looks up
        the coincidence window from the object dictionary.
        """
        functions = self.coincidence_functions(nodes)
        cw = self.coincidence_windows[nodes]
        
        cp = {source: self.calc_cp_between_nodes_for_a_source(nodes, source)
                  for source in cw}
        
        return cp

    def optimize_coincidence_window(
            self, nodes, source, lbound=1, ubound="default",
            QBER=None, individual=False, nonnegative=False):
        """Sets the coincidence windows between nodes for a given source
        to the optimal value.

        optimize_coincidence_window(nodes)
            -> (optimal window, secure key rate)
        """
        if ubound == "default":
            ubound = 10*max([detector.jitter
                          for node in nodes
                          for detector in node.detectors])
        def negative_SKR_given_coincidence_window(coincidence_window):
            self.set_coincidence_windows(nodes, source, coincidence_window)
            return -self.secure_key_rate(
                    nodes, QBER, individual, nonnegative)

        optimum = optimize.minimize_scalar(
                negative_SKR_given_coincidence_window,
                bounds=(lbound, ubound), method='bounded')
        if optimum.x > 0.9*ubound:
            print("That coincidence window looks really high.")
        return (optimum.x, -optimum.fun)

    def optimize_coincidence_windows(
            self, nodes, **kwargs):
        """Sets the coincidence windows between nodes to the optimal value.

        optimize_coincidence_windows(nodes)
            -> {source: (optimal window, secure key rate)}
        """
        tuples = {source: self.optimize_coincidence_window(
                            nodes, source, **kwargs)
                  for source in self.coincidence_windows[nodes]}
        return tuples

    def singles_rate(self, source, node):
        return self.splitting_table.singles_rate(source, node)

    def total_singles_rate(self, node):
        total = 0
        for source in self.sources[node]:
            total += self.singles_rate(source, node)
        return total

    def pair_rate(self, source, nodes):
        return self.splitting_table.pair_rate(source, nodes)

    def total_pair_rate(self, nodes):
        total = 0
        for source in self.sources[nodes]:
            total += self.pair_rate(source, nodes)
        return total

    def incidence_rates(self, node):
        """Return a DetectorPair namedtuple of incidence_rates at node.

        The incidence_rates represent the rate of photons arriving at
        each detector. Accounts for all sources which transmit to node.
        Does not account for detector efficiency.
        Access the HA and VD incidence_rate at node with
            self.incidence_rates(node).HA
            self.incidence_rates(node).VD
        """
        # First entry will be HA, second will be VD
        incidence_rates = [0, 0]
        T = node.transmissions
        for source in self.sources[node]:
            singles = self.singles_rate(source, node)
            P = source.state_probs
            singles_HV = singles*node.basis_choice_HVtoAD
            singles_AD = singles*(1 - node.basis_choice_HVtoAD)
            # Rate(H) = Rate(HV)*P(H|HV measurement)*transmission(H)
            singles_H = singles_HV*(P.HH + (P.HV+P.VH)/2)*T.H
            singles_A = singles_AD*(P.AA + (P.AD+P.DA)/2)*T.A
            singles_V = singles_HV*(P.VV + (P.HV+P.VH)/2)*T.V
            singles_D = singles_AD*(P.DD + (P.AD+P.DA)/2)*T.D
            incidence_rates[0] += singles_H + singles_A
            incidence_rates[1] += singles_V + singles_D
        return al.DetectorPair(*incidence_rates)

    def detection_rates(self, node):
        """Return a DetectorPair namedtuple of detection_rates at node.

        Calls the incidence_rates(self, node) method to automatically
        determine the detection_rates. To test the detectors with a
        user-chosen incidence rate, use the Detector method directly.
        """
        det_rate_HA = node.detection_rates.HA(self.incidence_rates(node).HA)
        det_rate_VD = node.detection_rates.VD(self.incidence_rates(node).VD)
        return al.DetectorPair(det_rate_HA, det_rate_VD)

    def accidentals(self, nodes):
        """Return the accidentals rate between nodes.

        Accounts for all sources. Accounts for coincidences between
        either of one node's detectors with either of the other node's
        detectors.

        NEW IN VERSION 5: Sum over all coincidence windows for all
        sources. This accounts for the fact that different lambda
        pairs have peaks in different places, so we need to account for
        'the area under the curve' multiple times.
        """
        node1 = nodes[0]
        node2 = nodes[1]
        accidentals = 0
        for source, window in self.coincidence_windows[nodes].items():
            accidentals += (sum(self.detection_rates(node1))
                            * sum(self.detection_rates(node2))
                            * window)
        return accidentals

    def _link_shifted_key_rate(self, source, nodes):
        """Shifted key rate contribution from a source-node1-node2 link

        Rate of shifted key generation between two nodes linked by a
        source. Does not account for accidentals. Intended for use by
        self.shifted_key_rate(self, nodes) in order to calculate the
        total shifted key rate. The formula is:
        _link_shifted_key_rate =
            pair_rate * (P(both nodes choose HV basis)
                          * P(both photons detected|HV chosen))
                         + P(both nodes choose AD basis)
                            * P(both photons detected|AD chosen))
        where the pair_rate is that due to the particular source, but
        the P(both photons detected|basis chosen) depends on both the
        probability of any given pair of polarizations (characteristic
        of the source) and the different transmissions of those
        polarizations at each node, which in turn depends on the total
        incidence_rate on the individual detectors - a function of the
        dark_rate, noise_rate and singles_rate from all connected
        sources.
        """
        node1 = nodes[0]
        node2 = nodes[1]
        B_HV1 = node1.basis_choice_HVtoAD
        B_HV2 = node2.basis_choice_HVtoAD
        B_AD1 = 1 - B_HV1
        B_AD2 = 1 - B_HV2
        TT_1 = node1.total_transmissions
        TT_2 = node2.total_transmissions
        rates_1 = self.incidence_rates(node1)
        rates_2 = self.incidence_rates(node2)
        T_H1 = TT_1.H(rates_1.HA)
        T_A1 = TT_1.A(rates_1.HA)
        T_V1 = TT_1.V(rates_1.VD)
        T_D1 = TT_1.D(rates_1.VD)
        T_H2 = TT_2.H(rates_2.HA)
        T_A2 = TT_2.A(rates_2.HA)
        T_V2 = TT_2.V(rates_2.VD)
        T_D2 = TT_2.D(rates_2.VD)
        cp = self.coincidence_probabilities[nodes][source]
        # C_XY = Probability that both photons reach a detector given
        #        that both 'chose' the XY path (C = 'Coincidence')
        C_HV = (source.state_probs.HH * T_H1 * T_H2 * cp.HH
                + source.state_probs.VV * T_V1 * T_V2 * cp.VV
                + source.state_probs.HV * T_H1 * T_V2 * cp.HV
                + source.state_probs.VH * T_V1 * T_H2 * cp.VH)
        C_AD = (source.state_probs.AA * T_A1 * T_A2 * cp.AA
                + source.state_probs.DD * T_D1 * T_D2 * cp.DD
                + source.state_probs.AD * T_A1 * T_D2 * cp.AD
                + source.state_probs.DA * T_D1 * T_A2 * cp.DA)
        K = self.pair_rate(source, nodes)*(B_HV1*B_HV2*C_HV
                                           + B_AD1*B_AD2*C_AD)
        return K

    def shifted_key_rate(self, nodes):
        """Return shifted_key_rate accounting for all sources.

        Also includes all accidental coincidences.
        """
        shifted_key_rate = 0
        for source in self.sources[nodes]:
            shifted_key_rate += self._link_shifted_key_rate(source, nodes)
        # Account for accidentals here.
        shifted_key_rate += self.accidentals(nodes)
        return shifted_key_rate

    def QBER(self, nodes, F=0.0):
        """Return the Quantum Bit Error Rate between a pair of nodes.

        nodes parameter should be a 2-tuple or similar of node objects.
        F (for fudge factor) is the proportion of the shifted_key_rate
        that contributes to the QBER through inexplicable reasons.
        The formula for the QBER is
            QBER = (wrong bits) / (total bits)
        where the total is the shifted_key_rate and the wrong bits are
        a contribution from the total and additionally half of the
        accidentals (the other half will be correct by chance).
        """
        wrong = F*self.shifted_key_rate(nodes) + self.accidentals(nodes)/2
        total = self.shifted_key_rate(nodes)
        try:
            QBER = wrong/total
        except ZeroDivisionError:
            QBER = 0
        return QBER

    def secure_key_rate(self, nodes, QBER=None, individual=False,
                        nonnegative=False):
        """Return the secure key rate between a pair of nodes.

        nodes parameter should be a 2-tuple or similar of node objects.
        QBER can be set to a given value. By default it is not set and
        will call the self.QBER method.
        individual parameter, if True, restricts the space of possible
        eavesdropper attacks to 'individual attacks'. For justification
        of the formula, see e.g. Gisin (2002) Rev. Mod. Phys. By
        default it is set to False and Eve is given the maximum
        possible information.
        """
        if QBER is None:
            QBER = self.QBER(nodes)
        AliceBobMutualInfo = 1 - self.entropy(QBER)
        if individual is True:
            AliceEveMutualInfo = (
                    1 - self.entropy((0.5+math.sqrt(QBER*(1-QBER))))
                                 )
        else:
            AliceEveMutualInfo = 1 - AliceBobMutualInfo

        secure_key_rate = (self.shifted_key_rate(nodes)
                           * self.error_correction_efficiency()
                           * (AliceBobMutualInfo - AliceEveMutualInfo))

        if nonnegative == True:
            return max([0, secure_key_rate])
        elif nonnegative == False:
            return secure_key_rate
        else:
            raise ValueError("nonnegative must be True or False")

    def total_secure_key_rate(self, QBER=None, individual=False,
                              nonnegative=False):
        """Return the sum over all node pairs of secure key rates.

        
        """
        total_secure_key_rate = 0
        for nodes in self.all_node_pairs:
            total_secure_key_rate += self.secure_key_rate(
                    nodes, QBER=QBER, individual=individual, nonnegative=True)

        if nonnegative == True:
            return max([0, total_secure_key_rate])
        elif nonnegative == False:
            return total_secure_key_rate
        else:
            raise ValueError("nonnegative must be True or False")


if __name__ == "__main__":

    name = ['Alice', 'Bob', 'Chloe', 'Dave']
    source = Source(pair_rate=6E5)
    node = {name[i]: Node(name=name[i]) for i in range(4)}
    nodesA = (node['Alice'], node['Bob'])
    nodesB = (node['Alice'], node['Dave'])

    four_person_splitting1 = [
            [node['Alice'], (1, 1, .2), (2, 1, .2), (3, 1, .2)],
            [node['Bob'], (-1, 1, .2), (4, 1, .2), (5, 1, .2)],
            [node['Chloe'], (-2, 1, .2), (-4, 1, .2), (6, 1, .2)],
            [node['Dave'], (-3, 1, .2), (-5, 1, .2), (-6, 1, .2)]
                            ]

    four_person_splitting2 = [
            [node['Alice'], (1, .5, .2), (2, .5, .2), (3, .5, .2)],
            [node['Bob'], (-1, .5, .2), (2, .5, .2), (-3, .5, .2)],
            [node['Chloe'], (1, .5, .2), (-2, .5, .2), (-3, .5, .2)],
            [node['Dave'], (-1, .5, .2), (-2, .5, .2), (3, .5, .2)]
                             ]

    four_person_splitting3 = [
            [node['Alice'], (1, .5, .2), (2, .5, .2)],
            [node['Bob'], (-1, .5, .2), (2, .5, .2)],
            [node['Chloe'], (1, .5, .2), (-2, .5, .2)],
            [node['Dave'], (-1, .5, .2), (-2, .5, .2)]
                             ]

    sources1 = source.splice(6)
    sources2 = source.splice(3)
    sources3 = source.splice(2)

    s1 = SplittingTable(four_person_splitting1, sources1)
    s2 = SplittingTable(four_person_splitting2, sources2)
    s3 = SplittingTable(four_person_splitting3, sources3)
    n1 = Network(s1)
    n2 = Network(s2)
    n3 = Network(s3)

#####

    figure_number = 0
    number_of_points = 100
    QBERs = np.linspace(0, 0.15, num=number_of_points)

    figure_number += 1
    entropy = np.array([n1.entropy(QBER) for QBER in QBERs])
    plt.figure(figure_number)
    plt.plot(QBERs, entropy)
    plt.xlabel('QBER')
    plt.ylabel('Entropy')
    plt.grid(True)
    plt.show

    figure_number += 1
    SKR1 = np.array([n1.secure_key_rate(nodesA, QBER, individual=True)
                     for QBER in QBERs])
    SKR2 = np.array([n2.secure_key_rate(nodesA, QBER, individual=True)
                     for QBER in QBERs])
    SKR3A = np.array([n3.secure_key_rate(nodesA, QBER, individual=True)
                      for QBER in QBERs])
    SKR3B = np.array([n3.secure_key_rate(nodesB, QBER, individual=True)
                      for QBER in QBERs])
    plt.figure(figure_number)
    plt.plot(QBERs, SKR1, 'b', label='Six splices')
    plt.plot(QBERs, SKR2, 'r', label='Three splices')
    plt.plot(QBERs, SKR3A, 'g', label='Two splices, Alice-Bob')
    plt.plot(QBERs, SKR3B, 'c', label='Two splices, Alice-Dave')
    plt.legend()
    plt.ylim(0, .3)
    plt.xlabel('QBER')
    plt.ylabel('Secure key rate assuming individual attack(/s)')
    plt.grid(True)
    plt.show

    figure_number += 1
    SSKR1 = np.array([n1.secure_key_rate(nodesA, QBER, individual=False)
                      for QBER in QBERs])
    SSKR2 = np.array([n2.secure_key_rate(nodesA, QBER, individual=False)
                      for QBER in QBERs])
    SSKR3A = np.array([n3.secure_key_rate(nodesA, QBER, individual=False)
                       for QBER in QBERs])
    SSKR3B = np.array([n3.secure_key_rate(nodesB, QBER, individual=False)
                       for QBER in QBERs])
    plt.figure(figure_number)
    plt.plot(QBERs, SSKR1, 'b', label='Six splices')
    plt.plot(QBERs, SSKR2, 'r', label='Three splices')
    plt.plot(QBERs, SSKR3A, 'g', label='Two splices, Alice-Bob')
    plt.plot(QBERs, SSKR3B, 'c', label='Two splices, Alice-Dave')
    plt.legend()
    plt.ylim(0, .3)
    plt.xlabel('QBER')
    plt.ylabel('Strict secure key rate (/s)')
    plt.grid(True)
    plt.show

    figure_number += 1
    original_pair_rate = sources1[0].pair_rate
    pair_rates = np.logspace(1, 8, 200)
    # Secure key rate if only one source varies its pair rate
    SKRs = []
    for pair_rate in pair_rates:
        sources1[0].pair_rate = pair_rate
        SKR = n1.secure_key_rate(nodesA)
        SKRs.append(SKR)
    # Secure key rate if all sources vary pair rates
    SKRs2 = []
    for pair_rate in pair_rates:
        for source in sources1:
            source.pair_rate = pair_rate
        SKR = n1.secure_key_rate(nodesA)
        SKRs2.append(SKR)
    plt.figure(figure_number)
    plt.plot(pair_rates, SKRs, 'g', label='Only source 1 changes')
    plt.plot(pair_rates, SKRs2, 'm', label='All sources change')
    plt.legend()
    plt.xlabel('Pair rate (/s)')
    plt.ylabel('Secure key rate (/s)')
    plt.xscale('log')
    plt.yscale('log')
    plt.show
    for source in sources1:
        source.pair_rate = original_pair_rate

#    figure_number += 1
#    original_pair_rate = sources1[0].pair_rate
#    pair_rates = np.logspace(1, 8, 200)
#    # Secure key rate if only one source varies its pair rate
#    SKRs = []
#    for pair_rate in pair_rates:
#        sources1[0].pair_rate = pair_rate
#        SKR = n1.optimize_coincidence_window(nodesA)[1]
#        SKRs.append(SKR)
#    # Secure key rate if all sources vary pair rates
#    SKRs2 = []
#    for pair_rate in pair_rates:
#        for source in sources1:
#            source.pair_rate = pair_rate
#        SKR = n1.optimize_coincidence_window(nodesA)[1]
#        SKRs2.append(SKR)
#    plt.figure(figure_number)
#    plt.plot(pair_rates, SKRs, 'g', label='Only source 1 changes')
#    plt.plot(pair_rates, SKRs2, 'm', label='All sources change')
#    plt.legend()
#    plt.xlabel('Pair rate (/s)')
#    plt.ylabel('Secure key rate (/s)')
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.show
#    for source in sources1:
#        source.pair_rate = original_pair_rate

# example_xor.py
#  pyAGE - A Python implementation of the Analog Genetic Encoding
#  Copyright (C) 2010  Janosch Gräf
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from age import Population, Agent, Descriptor
from math import sqrt
from time import time
from random import uniform
import pycann
try:
    from termcolor import colored
    _HAVE_TERMCOLOR_ = True
except ImportError:
    _HAVE_TERMCOLOR_ = False

# TODO: put this in a package with pyAGE and pyCANN
#       make a function out of this for generating different boolean gates (2
#       inputs, 1 output) with a given testtable

population_size = 100
desc = Descriptor(alphabet = "ACGT",  # We use a biology-like alphabet
                  # Neurons:
                  devices = ["ACAA",  # Threshold activation function
                             "ACAC",  # Sigmoidal (exp) activation function
                             "ACAG"], # Linear activation function
                  terminal = "TT",   # Terminal marker
                  parameter = "GG",  # Parameter marker
                  elitism = 0.1,      # Elitism: 50%
                  # Mutation possibilities:
                  possibilities = {"char_delete": 0.009,
                                   "char_insert": 0.01,
                                   "char_replace": 0.01,
                                   "frag_delete": 0.01,
                                   "frag_move": 0.01,
                                   "frag_copy": 0.02,
                                   "device_insert": 0.05,
                                   "chromosome_delete": 0.00029,
                                   "chromosome_copy": 0.0003,
                                   "chromosome_crossover": 0.0003})
testtable = ((0.0, 0.0, 0.0),
             (0.0, 1.0, 1.0),
             (1.0, 0.0, 1.0),
             (1.0, 1.0, 0.0))

class Network(pycann.Network):
    """ Implements a pycann network created from a genome. """
    
    def __init__(self, genome):
        # mapping for different activation functions
        activation_functions = {"ACAA": "SIGMOID_STEP",
                                "ACAC": "SIGMOID_EXP",
                                "ACAG": "LINEAR"}
        # init network
        pycann.Network.__init__(self, 2, max((0, len(genome.devices)-3)), 1)
        # configure neurons
        for i in range(len(genome.devices)):
            di = genome.devices[i]
            self.set_activation_function(i, activation_functions[di.device])
            try:
                t = di.parameters[0][1]
            except IndexError:
                t = 1.0
            self.set_threshold(i, t)
            for j in range(len(genome.devices)):
                if (i!=j):
                    dj = genome.devices[j]
                    try:
                        w = -1.0+4.0*genome.terminal_score(di.terminals[0], dj.terminals[1])
                    except IndexError:
                        w = 1.0
                    self.set_weight(i, j, w)

    def process(self, *inputs):
        self.set_inputs(*inputs)
        self.step()
        return self.get_outputs()

def eval_func(population, agent):
    """ Evaluates an agent. """
    global testtable

    # create neural network from genome
    net = Network(agent.genome)
    agent.network = net

    # calculate adaption
    Σ = Σo = 0.0
    n = len(testtable)
    for t in testtable:
        o = net.process(t[0], t[1])[0]
        Σo += o
        Σ += (o-t[2])**2
    s = sqrt(Σ/(n-1))
    if (Σo!=0.0):
        ν = s*n/Σo
    else:
        return 0.0
    #print("---------------")
    #print("Σ(yi-Ti)² = "+str(Σ))
    #print("Σ(yi)     = "+str(Σo))
    #print("n         = "+str(n))
    #print("s         = "+str(s))
    #print("ν         = "+str(ν))

    noise = 0.005 # noise of ±0.5%
    return 1-ν#+uniform(-noise, +noise)

def print_sequence(seq):
    """ Prints a sequence of bases. If the termcolor module is available this
        functions prints bases in colors. """
    if (_HAVE_TERMCOLOR_):
        BASECOLORS = {'A': 'on_yellow',
                      'C': 'on_cyan',
                      'G': 'on_red',
                      'T': 'on_magenta'}
        cseq = ""
        for b in seq:
            cseq += colored(b, 'grey', BASECOLORS[b])
        return cseq
    else:
        return seq


print("Creating population with "+str(population_size)+" agents")
population = Population(agedesc = desc, eval_callback = eval_func)
for i in range(population_size):
    population.add(Agent(desc))

print("Evolving networks")
best = None
t_next = 0
try:
    while (best==None or best.fitness<0.99): # abort if 99% adaption reached
        population.mate()
        population.mutate()

        population.evaluate()
        best = population.get_best()

        t = time()
        if (t>t_next):
            print("[%d] Best agent: %.6f%% (%d neurons; %d chromosomes; %d bases; %d devices)"%(population.generation, best.fitness*100, best.network.size, len(best.genome.chromosomes), len(best.genome), len(best.genome.devices)))
            #for i in range(len(best.genome.chromosomes)):
            #    print("  #"+str(i)+": "+print_sequence(best.genome.chromosomes[i]))
            t_next = t+1.0
except KeyboardInterrupt:
    pass

if (best!=None):
    print()
    print("Best network:")
    print("Fitness (Error): "+str(100.0*best.fitness)+"% ("+str(100.0-100.0*best.fitness)+"%)")
    print("Chromosomes:")
    for i in range(len(best.genome.chromosomes)):
        print("#"+str(i)+": "+print_sequence(best.genome.chromosomes[i]))
    print("Size: "+str(best.network.size))
    print("Activation functions: "+repr([best.network.get_activation_function(i) for i in range(best.network.size)]))
    print("Thresholds: "+repr([best.network.get_threshold(i) for i in range(best.network.size)]))
    print("Weights:")
    for i in range(best.network.size):
        print("From "+str(i)+": "+repr([best.network.get_weight(i, j) for j in range(best.network.size)]))
    print("Output:")
    for t in testtable:
        o = best.network.process(t[0], t[1])[0]
        print("%.1f, %.1f -> %.4f"%(t[0], t[1], o))


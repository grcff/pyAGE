# age/genome.py
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


import re
import random

from .descriptor import Descriptor


class Device:
    def __init__(self, genome, device, token):
        self.genome = genome
        self.device = device
        self.token = token
        self.terminals = []
        self.parameters = []

    def __del_(self):
        pass

    def add_terminal(self, terminal):
        self.terminals.append(terminal)

    def add_parameter(self, parameter):
        self.parameters.append((parameter, self.parameter_decode(parameter, self.genome.desc.come_alpha)))

    def parameter_decode(self, parameter, alpha = 1.0):
        # use CoME
        beta = float(len(self.genome.desc.alphabet))
        return sum((self.genome.desc.alphabet.find(parameter[i])*pow(beta*alpha, -i) for i in range(len(parameter)))) / \
               ((beta-1) * sum((pow(beta*alpha, -i) for i in range(len(parameter)))))

    def __eq__(self, y):
        if (type(y)==Device):
            return self.token==y.token
        else:
            return False

    def __ne__(self, y):
        return not self.__eq__(self, y)

    def __len__(self):
        return len(self.token)

    def __str__(self):
        return self.token


class Genome:
    def __init__(self, **params):
        self.desc = params.get("desc", None)
        if (self.desc==None):
            self.desc = Descriptor(**params)
        else:
            assert type(self.desc)==Descriptor
        assert self.desc.check()

        self.chromosomes = params.get("chromosomes", [])
        assert type(self.chromosomes)==list
        for c in self.chromosomes:
            assert self.desc.check_token(c)

        self.devices = []
        self.re_find_device = re.compile("|".join(self.desc.devices))
        self.re_find_termparam = re.compile(self.desc.terminal+"|"+self.desc.parameter)

    def __del__(self):
        pass

    def get_descriptor(self):
        return self.descriptor

    def add_randomly(self, num_chromosomes = (1, 1), len_chromosomes = (10, 100)):
        if (num_chromosomes[0]==num_chromosomes[1]):
            num_chromosomes = num_chromosomes[0]
        else:
            num_chromosomes = random.randrange(num_chromosomes[0], num_chromosomes[1])
        for i in range(num_chromosomes):
            length = random.randrange(len_chromosomes[0], len_chromosomes[1])
            chromosome = ""
            for j in range(length):
                chromosome += random.choice(self.desc.alphabet)
            self.chromosomes.append(chromosome)

    def add_chromosome(self, chromosome = ""):
        assert self.desc.check_token(chromosome)
        self.chromosomes.append(chromosome)

    def remove_chromosome(self, chromosome):
        self.chromosomes.remove(chromosome)

    def get_chromosome(self, index):
        return self.chromosomes[index] if index<len(self.chromosomes) else False

    def iter_chromosomes(self):
        return self.chromosomes.__iter__()

    def get_chromosomes(self):
        return self.chromosomes

    def num_chromosomes(self):
        return len(self.chromosomes)

    def get_devices(self):
        return self.devices

    def iter_devices(self):
        return self.devices.__iter__()

    def parse(self):
        self.devices = [] # reset device list
        for c in self.chromosomes:
            self.parse_chromosome(c)

    def search_all(self, regex, s):
        matches = []
        m = None
        while (True):
            m = regex.search(s, m.end()+1 if m!=None else 0)
            if (m==None):
                break
            matches.append(m)
        return matches

    def parse_chromosome(self, chromosome):
        matches = self.search_all(self.re_find_device, chromosome)
        for i in range(len(matches)-1):
            m1 = matches[i]
            m2 = matches[i+1]
            token = chromosome[m1.start():m2.start()-1]
            self.devices.append(self.parse_device(m1.group(), token))
        if (len(matches)>0):
            token = chromosome[matches[-1].end()+1:]
            self.devices.append(self.parse_device(matches[-1].group(), token))

    def parse_device(self, device_str, token):
        device = Device(self, device_str, token)
        matches = self.search_all(self.re_find_termparam, token)
               
        for i in range(len(matches)):
            m = matches[i]
            if (i==0):
                p1 = 0
            else:
                p1 = matches[i-1].end()+1
            data = token[p1:m.start()-1]
            if (len(data)>0):
                if (m.group()==self.desc.terminal):
                    device.add_terminal(data)
                else:
                    device.add_parameter(data)
        return device

    def local_alignment_score(self, a, b):
        # using smith-waterman algorithm
        # similarity score function
        def w(a, b, scoring = None):
            if (scoring==None):
                return -1 if (a==None or b==None or a!=b) else 2
            else:
                return scoring[a][b]

        # let smaller string be a, since we dont store all lines
        if (len(b)<len(a)):
            a, b = b, a
            
        # calculate H score
        # H[0] current line
        # H[1] last line
        H = [[],[]]
        score = 0
        for i in range(len(a)):
            for j in range(len(b)):
                if (i==0 or j==0):
                    H[0].append(0)
                else:
                    h = max([0,                                           # empty suffix
                             H[1][j-1]+w(a[i], b[j], self.desc.scoring),  # math/mismatch
                             H[1][j]+w(a[i], None, self.desc.scoring),    # deletion
                             H[0][j-1]+w(None, b[j], self.desc.scoring)]) # insertion
                    if (h>score):
                        score = h
                    H[0].append(h)
            # swap lines
            H.reverse()
        return score

    def terminal_score(self, tA, tB):
        return 2.0*self.local_alignment_score(tA, tB)/(len(tA)+len(tB))

    def crossover_chromosomes(self, cA, cB):
        return cA, cB
        pA = random.randrange(len(cA))
        pB = random.randrange(len(cB))
        cA = cA[:pA]+cB[pB:]
        cB = cB[:pB]+cA[pA:]
        return (cA, cB) if random.getrandbits(1)==0 else (cB, cA)

    def crossover(gA, gB, return_genome = True):
        child_chromosomes = []
        # crossover chromosomes
        for i in range(min(len(gA.chromosomes), len(gB.chromosomes))):
            child_chromosomes.append(gA.crossover_chromosomes(gA.chromosomes[i], gB.chromosomes[i])[0])
        d = len(gA.chromosomes)-len(gB.chromosomes)
        # copy remaining chromosomes
        if (d<0):
            for i in range(len(gA.chromosomes), len(gB.chromosomes)):
                child_chromosomes.append(gB.chromosomes[i])
        elif (d>0):
            for i in range(len(gB.chromosomes), len(gA.chromosomes)):
                child_chromosomes.append(gA.chromosomes[i])
        if (return_genome):
            return Genome(chromosomes = child_chromosomes, desc = gA.desc)
        else:
            return child_chromosomes

    def mutation_occurs(self, name):
        """ Check if mutation named 'name' occurs. """
        # if genome is empty, we won't mutate anything TODO
        if (len(self.chromosomes)==0):
            return False
        # return boolean with discrete possibility
        try:
            p = self.desc.possibilities[name]
            return random.random()<p
        except KeyError:
            return False

    def mutate(self):
        """ This method applies all mutation types with specified possibilities
            to the genome. """
        
        while (self.mutation_occurs("char_delete")):
            i = random.randrange(len(self.chromosomes))
            #print("char_delete:", i, len(self.chromosomes[i]))
            if (len(self.chromosomes[i])<1):
                continue
            p = random.randrange(len(self.chromosomes[i]))
            self.chromosomes[i] = self.chromosomes[i][:p]+self.chromosomes[i][p+1:]

        while (self.mutation_occurs("char_insert")):
            i = random.randrange(len(self.chromosomes))
            #print("char_insert:", i, len(self.chromosomes[i]))
            if (len(self.chromosomes[i])<1):
                continue
            p = random.randrange(len(self.chromosomes[i]))
            c = random.choice(self.desc.alphabet)
            self.chromosomes[i] = self.chromosomes[i][:p]+c+self.chromosomes[i][p:]

        while (self.mutation_occurs("char_replace")):
            i = random.randrange(len(self.chromosomes))
            #print("char_replace:", i, len(self.chromosomes[i]))
            if (len(self.chromosomes[i])<1):
                continue
            p = random.randrange(len(self.chromosomes[i]))
            c = random.choice(self.desc.alphabet)
            self.chromosomes[i] = self.chromosomes[i][:p]+c+self.chromosomes[i][p+1:]

        while (self.mutation_occurs("frag_delete")):
            i = random.randrange(len(self.chromosomes))
            if (len(self.chromosomes[i])<2):
                continue
            p = random.randrange(len(self.chromosomes[i])-1)
            l = random.randrange(1, len(self.chromosomes[i])-p)
            self.chromosomes[i] = self.chromosomes[i][:p]+self.chromosomes[i][p+l:]

        while (self.mutation_occurs("frag_move")):
            i = random.randrange(len(self.chromosomes))
            if (len(self.chromosomes[i])<2):
                continue
            p = (random.randrange(len(self.chromosomes[i])-1),
                random.randrange(len(self.chromosomes[i])))
            l = random.randrange(1, len(self.chromosomes[i])-p[0])
            if (p[0]<p[1]):
                self.chromosomes[i] = self.chromosomes[i][:p[0]]+self.chromosomes[i][p[0]+l:p[1]]+self.chromosomes[i][p[0]:p[0]+l]+self.chromosomes[i][p[1]:]
            else:
                self.chromosomes[i] = self.chromosomes[i][:p[1]]+self.chromosomes[i][p[0]:p[0]+l]+self.chromosomes[i][p[1]:p[1]]+self.chromosomes[i][p[0]:]

        while (self.mutation_occurs("frag_copy")):
            i = random.randrange(len(self.chromosomes))
            if (len(self.chromosomes[i])<2):
                continue
            p = (random.randrange(len(self.chromosomes[i])-1),
                random.randrange(len(self.chromosomes[i])))
            l = random.randrange(1, len(self.chromosomes[i])-p[0])
            seq = self.chromosomes[i][p[0]:p[0]+l]
            self.chromosomes[i] = self.chromosomes[i][:p[1]]+seq+self.chromosomes[i][p[1]:]

        while (self.mutation_occurs("device_insert")):
            i = random.randrange(len(self.chromosomes))
            p = random.randrange(len(self.chromosomes[i]))
            l = random.randrange(max((int(0.2*len(self.chromosomes[i])), 5)))
            d = random.choice(self.desc.devices)+"".join((random.choice(self.desc.alphabet) for i in range(l)))
            self.chromosomes[i] = self.chromosomes[i][:p]+d+self.chromosomes[i][p:]

        while (self.mutation_occurs("chromosome_delete")):
            i = random.randrange(len(self.chromosomes))
            self.chromosomes.pop(i)

        while (self.mutation_occurs("chromosomes_copy")):
            i = random.randrange(len(self.chromosomes))
            self.chromosomes.append(self.chromosomes[i])

        # remove empty chromosomes
        try:
            while (True):
                self.chromosomes.remove("")
        except ValueError:
            pass
            
    def __eq__(self, y):
        if (type(y)!=Genome):
            return False
        if (len(self.chromosomes)!=len(y.chromosomes)):
            return False
        for i in range(len(self.chromosomes)):
            if (self.chromosomes[i]!=y.chromosomes[i]):
                    return False
        return True

    def __ne__(self, y):
        return not self.__eq__(y)

    def __len__(self):
        return sum(map(lambda c: len(c), self.chromosomes))
        #return len(self.chromosomes)

    def __str__(self):
        return str(self.chromosomes)

    def __getitem__(self, y):
        return self.chromosomes[y]

__all__ = ["Device", "Genome"]

import re
import random

 
class Descriptor:
    def __init__(self, **params):
        self.alphabet = params.get("alphabet", None)
        self.devices = params.get("devices", None)
        self.terminal = params.get("terminal", None)
        self.parameter = params.get("parameter", None)
        self.possibilities = params.get("possibilities", None)
        self.scoring = params.get("scoring", None)
        self.come_alpha = params.get("come_alpha", 1.0)

    def check_alphabet(self):
        if (type(self.alphabet)!=str):
            return False
        if (self.alphabet==""):
            return False
        for l in self.alphabet:
            if (self.alphabet.count(l)>1):
                return False
        return True

    def check_devices(self):
        if (type(self.devices)!=list):
            return False
        for d in self.devices:
            return self.check_token(d)

    def check_token(self, token):
        if (type(token)!=str or len(token)==0):
            return False
        for t in token:
            if (not t in self.alphabet):
                return False
        return True

    def check_possibilities(self):
        keys = ["char_delete", "char_insert", "char_replace", "frag_delete", "frag_move", "frag_copy", "device_insert", "chromosome_delete", "chromosome_copy", "chromosome_crossover"]
        if (self.possibilities==None):
            self.possibilities = {}
        elif (type(self.possibilities)!=dict):
            return False
        for p in keys:
            if (p not in self.possibilities):
                self.possibilities[p] = 0.0
            else:
                try:
                    self.possibilities[p] = float(self.possibilities[p])
                except ValueError:
                    return False
        return True

    def check_scoring_matrix(self):
        if (self.scoring==None):
            return True
        b = len(self.alphabet)
        if (type(self.scoring)!=tuple or len(self.scoring)!=b):
            return False
        for l in self.scoring:
            if (type(l)!=tuple or len(l)!=b):
                return False
            for s in l:
                if (type(s)!=float):
                    return False

    def check_come_alpha(self):
        if (type(self.come_alpha)!=float):
            try:
                self.come_alpha = float(self.come_alpha)
            except ValueError:
                return False
        return (self.come_alpha>=0.0 and self.come_alpha<=1.0)

    def check(self):
        return self.check_alphabet() \
           and self.check_devices() \
           and self.check_token(self.terminal) \
           and self.check_token(self.parameter) \
           and self.check_possibilities() \
           and self.check_scoring_matrix() \
           and self.check_come_alpha()


class Device:
    def __init__(self, genome, device, token, terminals = [], parameters = []):
        self.genome = genome
        self.device = device
        self.token = token
        self.terminals = terminals
        self.parameters = []
        for p in parameters:
            self.parameters.append((p, self.parameter_decode(p, self.genome.desc.come_alpha)))

    def __del_(self):
        pass

    def parameter_decode(self, parameter, alpha = 1.0):
        # use CoME
        beta = float(len(self.genome.desc.alphabet))
        return sum((parameter[i]*pow(beta*alpha, -i) for i in range(len(parameter)))) / \
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
        self.chromosomes = params.get("chromosomes", [])
        assert type(self.chromosomes)==list
        for c in self.chromosomes:
            assert desc.check_token(c)

        self.desc = params.get("desc", None)
        if (self.desc==None):
            self.desc = Descriptor(**params)
        else:
            assert type(self.desc)==Descriptor
        assert self.desc.check()

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

    def add_chromosomes(self, chromosome = ""):
        assert self.desc.check_token(chromosome)
        self.chromosomes.append(chromosome)

    def remove_chromosome(self, chromosome):
        self.chromosomes.remove(chromosome)

    def get_chromosomes(self, index):
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
                p1 = matches[i-1].end()
            data = token[p1:m.end()]
            if (m.group()==self.desc.terminal):
                device.terminals.append(data)
            else:
                device.parameters.append(data)
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
        for i in range(len(tA)):
            for j in range(len(tB)):
                if (i==0 or j==0):
                    H[0].append(0)
                else:
                    H[0].append(max([0,                                           # empty suffix
                                    H[1][j-1]+w(a[i], b[j], self.desc.scoring),   # math/mismatch
                                    H[1][j]+w(a[i], None, self.desc.scoring),     # deletion
                                    H[0][j-1]+w(None, b[j], self.desc.scoring)])) # insertion
            # swap lines
            H.reverse()
        return max((max(l) for l in H))

    def terminal_score(self, tA, tB):
        return 2.0*self.local_alignment_score(tA, tB)/(len(tA)+len(tB))

    def __eq__(self, y):
        if (type(y)!=Genome):
            return False
        if (self.num_chromosomes!=y.num_chromosomes):
            return False
        for i in len(self.chromosomes):
            if (self.chromosomes[i]!=y.chromosomes):
                    return False
        return True

    def __ne__(self, y):
        return not self.__eq__(y)

    def __len__(self):
        return len(self.chromosomes)

    def __str__(self):
        return str(self.chromosomes)

    def __getitem__(self, y):
        return self.chromosomes[y]

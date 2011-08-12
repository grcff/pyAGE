# age/__init__.py
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

from tarfile import TarFile, TarInfo
from os.path import splitext
from random import sample, random
from math import ceil
from .genome import *

""" This module provides classes for handling of whole population of AGE agents """

def roulette_wheel(iterable, n = 1, key = lambda x: x):
    """ Implementation of a roulette wheel with specified possibilities. """
    available = sorted(iterable, key = key, reverse = True)
    if (len(available)<n):
        raise ValueError("Not enough ("+str(n)+") elements in iterable ("+str(len(available))+")")
    S = sum(map(key, available))
    picked = []

    while (len(picked)<n):
        r = random()*S
        for e in available:
            f = key(e)
            r -= f
            if (r<=0.0):
                picked.append(e)
                available.remove(e)
                S -= f
                break

    return picked

def test_roulette_wheel():
    iterable = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    S = sum(iterable)
    n = 10000
    Fe = list(map(lambda p: p*n/S, iterable))
    Fo = [0 for i in range(len(iterable))]

    for i in range(n):
        v = roulette_wheel(iterable)[0]
        Fo[iterable.index(v)] += 1

    X_square = sum(map(lambda fe, fo: (fo-fe)**2/fe, Fe, Fo))
    # Return true if 99% sure, that roulette wheel works
    return X_square<2.088


class Agent:
    def __init__(self, agedesc, **options):
        self.agedesc = agedesc
        self.id = options.get("id")
        if ("file" in options):
            # load from file (usually a member of a TAR file)
            self.load_from_file(options["file"])
        else:
            # create new agent
            self.genome = Genome(desc = agedesc,
                                 chromosomes = options.get("chromosomes", []))
            if ("chromosomes" not in options):
                # if no chromosomes given, generate some
                self.genome.add_randomly((1, 3), (10, 200))
            self.fitness = options.get("fitness", 0)

    def load_from_file(self, f):
        # read lines
        lines = map(lambda l: l.strip(), f.readlines())
        # filter out comments
        lines = filter(lambda l: not l.startswith("#"), lines)
        # ID
        self.id = eval(lines[0], {"__builtins__": None})
        # fitness
        self.fitness = float(lines[0])
        # genome
        self.genome = age.Genome(self.agedesc, lines[1:])

    def save_to_file(self, f):
        f.write(repr(self.id)+"\n")
        f.write(str(self.fitness)+"\n")
        for c in self.genome:
            f.write(c+"\n")


class Population:
    def __init__(self, **options):
        self.eval_callback = options["eval_callback"]
        if ("file" in options):
            # load from file
            self.load_from_file(options["file"])
        else:
            # create new population
            self.agedesc = options["agedesc"]
            self.generation = options.get("generation", 0)
            self.agents = []

    def load_from_file(self, f):
        tar = TarFile(f, "r")

        # load info file
        f = tar.extractfile("info.py")
        self.agedesc, self.generation = eval(f.read(-1), {"__builtins__": None})
        f.close()

        # load agents
        for info in tar.getmembers():
            if (splitext(info.name)[1]==".agt" and info.isfile()):
                f = tar.extractfile(info)
                self.add(Agent(self.agedesc, file = f))
                f.close()

        tar.close()

    def save_to_file(self, f):
        tar = TarFile(f, "w")

        # save info file
        f = StringIO(repr((self.agedesc, self.generation)))
        info = tar.gettarinfo(None, "info.py", f)
        tar.addfile(info, f)
        f.close()

        # save agents
        for i in range(len(self.agents)):
            f = StringIO()
            self.agents[i].save_to_file(f)
            info = tar.gettarinfo(None, str(i)+".agt", f)
            tar.addfile(info, f)
            f.close()

        tar.close()                

    def add(self, *agents):
        for a in agents:
            self.agents.append(a)

    def remove(self, *agents):
        for a in agents:
            self.agents.remove(a)

    def pick(self, n = 1, rwheel = False):
        if (rwheel):
            return roulette_wheel(self.agents, n, lambda a: a.fitness)
        else:
            return sample(self.agents, n)

    def mate(self):
        n = min((2, ceil(self.agedesc.elitism*len(self.agents))))
        elite = self.pick(n, True)
        offspring = []
        while (len(offspring)-len(elite)<len(self.agents)):
            a, b = self.pick(2, True)
            c = a.genome.crossover(b.genome, False)
            offspring.append(Agent(self.agedesc, chromosomes = c))
        self.agents = elite+offspring
        self.generation += 1

    def mutate(self):
        for a in self.agents:
            a.genome.mutate()

    def evaluate(self, n = 1):
        agents = self.pick(n)
        for a in agents:
            a.genome.parse()
        F = self.eval_callback(self, *agents)
        if (type(F)!=tuple):
            F = (F,)
        for i in range(len(F)):
            agents[i].fitness = F[i]

    def get_best(self):
        return max(self.agents, key = lambda a: a.fitness)


__all__ = ["Agent", "Population"]

# age/descriptor.py
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


class Descriptor:
    def __init__(self, **params):
        self.alphabet = params.get("alphabet", None)
        self.devices = params.get("devices", None)
        self.terminal = params.get("terminal", None)
        self.parameter = params.get("parameter", None)
        self.possibilities = params.get("possibilities", None)
        self.scoring = params.get("scoring", None)
        # TODO rename scoring matrix in substitution/insert/delete matrix
        self.come_alpha = params.get("come_alpha", 1.0)
        # used in populations
        self.elitism = params.get("elitism", 0.2)

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

    def check_scoring(self):
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

    def check_elitism(self):
        return (self.elitism>0.0 and self.elitism<=1.0)

    def check(self):
        # TODO raise ValueError when check goes wrong
        return self.check_alphabet() \
           and self.check_devices() \
           and self.check_token(self.terminal) \
           and self.check_token(self.parameter) \
           and self.check_possibilities() \
           and self.check_scoring() \
           and self.check_come_alpha() \
           and self.check_elitism
    
    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "Descriptor(alphabet = "+repr(self.alphabet)+",\n" \
              +"           devices = "+repr(self.devices)+",\n" \
              +"           terminal = "+repr(self.terminal)+",\n" \
              +"           parameter = "+repr(self.parameter)+",\n" \
              +"           possibilities = "+repr(self.possibilities)+",\n" \
              +"           scoring = "+repr(self.scoring)+",\n" \
              +"           come_alpha = "+repr(self.come_alpha)+",\n" \
              +"           elitism = "+repr(self.elitism)+")"


__all__ = ["Descriptor"]

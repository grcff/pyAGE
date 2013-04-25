import age
from termcolor import colored

def print_sequence(seq):
    BASECOLORS = {'A': 'on_yellow',
                  'C': 'on_cyan',
                  'G': 'on_red',
                  'T': 'on_magenta'}
    cseq = ""
    for b in seq:
        cseq += colored(b, 'grey', BASECOLORS[b], attrs = ["bold"])
    return cseq

def print_genome(genome, name):
    print("Genome '%s' (%d):"%(name, len(genome.chromosomes)))
    for i, c in enumerate(genome.chromosomes):
        print("Chromosome %d: %s"%(i + 1, print_sequence(c)))
    for i, d in enumerate(genome.devices):
        print("Device %d: %s"%(i + 1, print_sequence(d.device)))
        for j, p in enumerate(d.parameters):
            print(" Parameter %d: %f, %s"%(j + 1, d.parameters[j][1], print_sequence(d.parameters[j][0])))
        for j in range(len(d.terminals)):
            print(" Terminal %d: %s"%(j + 1, print_sequence(d.terminals[j])))
    print()


desc = age.Descriptor(alphabet = "ACGT",
                      devices = ["ACGA"],
                      terminal = "TGC",
                      parameter = "TGA")
print(desc)

genomes = []
for i in range(2):
    g = age.Genome(desc = desc)
    genomes.append(g)
    g.add_randomly((10, 20), (20, 200))
    g.parse()
    print_genome(g, "Parent #"+str(i+1))

print("Crossover...")

child = age.Genome.crossover(genomes[0], genomes[1])
child.parse()
print_genome(child, "Child")
#child_chromosomes = age.Genome.crossover(genomes[0], genomes[1], False)
#for i in range(len(child_chromosomes)):
#    print("Chromosomes "+str(i+1)+": "+child_chromosomes[i])

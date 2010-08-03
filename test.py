import age

def print_genome(genome, name):
    print("Genome '"+name+"' ("+str(len(genome.chromosomes))+"):")
    for i in range(len(genome)):
        print("Chromosome "+str(i+1)+": "+genome[i])
    for i in range(len(genome.devices)):
        d = genome.devices[i]
        print("Device "+str(i+1)+": "+str(d))
        for j in range(len(d.parameters)):
            print(" Parameter "+str(j+1)+": "+str(d.parameters[j][1])+", "+d.parameters[j][0])
        for j in range(len(d.terminals)):
            print(" Terminal "+str(j+1)+": "+d.terminals[j])
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

import age

genome = age.Genome(alphabet = "ACGT",
                      devices = ["ACGAT"],
                      terminal = "TGC",
                      parameter = "TGA")
genome.add_randomly((10, 20), (20, 200))

for i in range(len(genome)):
    print("Chromosome "+str(i+1)+": "+genome[i])

genome.parse()

for i in range(len(genome.devices)):
    d = genome.devices[i]
    print("Device "+str(i+1)+": "+d)

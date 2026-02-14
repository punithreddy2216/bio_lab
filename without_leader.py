from Bio import Entrez, SeqIO

Entrez.email = "Punithreddy2216@gmail.com"

accession = "CAB38517.1"

# Fetch sequence
handle = Entrez.efetch(db="protein",
                       id=accession,
                       rettype="fasta",
                       retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

sequence = str(record.seq)

print("Full Protein Sequence:")
print(sequence)

# First 5 amino acids
peptide = sequence[:5]
print("\nFirst 5 Amino Acids:", peptide)

# Mass table
mass_table = {
    'G':57,'A':71,'S':87,'P':97,'V':99,
    'T':101,'C':103,'I':113,'L':113,
    'N':114,'D':115,'K':128,'Q':128,
    'E':129,'M':131,'H':137,'F':147,
    'R':156,'Y':163,'W':186
}

# Parent mass
def peptide_mass(peptide):
    return sum(mass_table[aa] for aa in peptide)

parent_mass = peptide_mass(peptide)
print("Parent Mass:", parent_mass)

# Cyclic spectrum
def cyclic_spectrum(peptide):
    prefix_mass = [0]
    for i in range(len(peptide)):
        prefix_mass.append(prefix_mass[i] + mass_table[peptide[i]])

    peptide_mass_total = prefix_mass[-1]
    spectrum = [0]

    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
            spectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < len(peptide):
                spectrum.append(peptide_mass_total - (prefix_mass[j] - prefix_mass[i]))

    return sorted(spectrum)

spectrum = cyclic_spectrum(peptide)

print("Cyclic Spectrum:")
print(spectrum)

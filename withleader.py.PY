from Bio import Entrez, SeqIO
from collections import Counter

Entrez.email = "vikkurthiyasaswi@gmail.com"

accession = "CAB38517.1"

handle = Entrez.efetch(db="protein",
                       id=accession,
                       rettype="fasta",
                       retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

sequence = str(record.seq)
peptide = sequence[:5]

print("First 5 Amino Acids:", peptide)

mass_table = {
    'G':57,'A':71,'S':87,'P':97,'V':99,
    'T':101,'C':103,'I':113,'L':113,
    'N':114,'D':115,'K':128,'Q':128,
    'E':129,'M':131,'H':137,'F':147,
    'R':156,'Y':163,'W':186
}

amino_acids = list(mass_table.keys())

def peptide_mass(peptide):
    return sum(mass_table[aa] for aa in peptide)

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

experimental_spectrum = cyclic_spectrum(peptide)
parent_mass = peptide_mass(peptide)

print("Parent Mass:", parent_mass)
print("Experimental Spectrum:", experimental_spectrum)

def score(peptide, spectrum):
    theoretical = cyclic_spectrum(peptide)
    return sum((Counter(theoretical) & Counter(spectrum)).values())

def expand(peptides):
    return [p + aa for p in peptides for aa in amino_acids]

def leaderboard_cyclopeptide(spectrum, N):
    leaderboard = [""]
    leader_peptide = ""

    while leaderboard:
        leaderboard = expand(leaderboard)
        new_leaderboard = []

        for peptide in leaderboard:
            if peptide_mass(peptide) == parent_mass:
                if score(peptide, spectrum) > score(leader_peptide, spectrum):
                    leader_peptide = peptide
            if peptide_mass(peptide) <= parent_mass:
                new_leaderboard.append(peptide)

        leaderboard = sorted(new_leaderboard,
                             key=lambda p: score(p, spectrum),
                             reverse=True)[:N]

    return leader_peptide

best_peptide = leaderboard_cyclopeptide(experimental_spectrum, N=10)

print("Best Peptide Found:", best_peptide)

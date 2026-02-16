amino_acids = ['G','A','S','P','V','T','C','I','L','N','D','K','Q','E','M','H','F','R','Y','W']

mass = {
    'G':57,'A':71,'S':87,'P':97,'V':99,'T':101,'C':103,'I':113,'L':113,
    'N':114,'D':115,'K':128,'Q':128,'E':129,'M':131,'H':137,
    'F':147,'R':156,'Y':163,'W':186
}

spectrum = [0, 57, 71, 87, 97, 128, 158, 186, 229, 300]
parent_mass = max(spectrum)

def peptide_mass(p):
    total = 0
    for x in p:
        total = total + mass[x]
    return total

def is_consistent(p):
    for i in range(len(p)):
        for j in range(i+1, len(p)+1):
            sub = p[i:j]
            if peptide_mass(sub) not in spectrum:
                return False
    return True

def branch_and_bound(k):
    peptides = [""]

    c1 = []
    c2 = []
    c3 = []

    for length in range(1, k+1):
        print("\n", length, "-mer peptides")
        new_peptides = []

        for p in peptides:
            for a in amino_acids:
                pep = p + a
                m = peptide_mass(pep)

                if m > parent_mass:
                    print(pep, "mass =", m, "INCONSISTENT")
                elif is_consistent(pep):
                    print(pep, "mass =", m, "CONSISTENT")
                    new_peptides.append(pep)

                    if length == 1:
                        c1.append(pep)
                    elif length == 2:
                        c2.append(pep)
                    elif length == 3:
                        c3.append(pep)
                else:
                    print(pep, "mass =", m, "INCONSISTENT")

        peptides = new_peptides

    print("\nConsistent 1-mers:", c1)
    print("Consistent 2-mers:", c2)
    print("Consistent 3-mers:", c3)

branch_and_bound(3)
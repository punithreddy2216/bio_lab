# Smith–Waterman Algorithm
# Local sequence alignment
def smith_waterman(seq1,seq2,match=2,mismatch=-1,gap=-1):
    n=len(seq1)
    m=len(seq2)

    score=[[0 for j in range(m+1)] for i in range(n+1)]

    max_score=0
    max_pos=None

    # Fill scoring matrix
    for i in range(1,n+1):
        for j in range(1,m+1):
            if seq1[i-1]==seq2[j-1]:
                diag=score[i-1][j-1]+match
            else:
                diag=score[i-1][j-1]+mismatch

            up=score[i-1][j]+gap
            left=score[i][j-1]+gap

            score[i][j]=max(0,diag,up,left)

            if score[i][j]>max_score:
                max_score=score[i][j]
                max_pos=(i, j)

    # Traceback
    align1=""
    align2=""
    i,j=max_pos

    while score[i][j]!=0:
        if score[i][j]==score[i-1][j-1]+(match if seq1[i-1]==seq2[j-1] else mismatch):
            align1=seq1[i-1]+align1
            align2=seq2[j-1]+align2
            i-=1
            j-=1
        elif score[i][j] == score[i-1][j]+gap:
            align1=seq1[i-1]+align1
            align2="-"+align2
            i-=1
        else:
            align1="-"+align1
            align2=seq2[j-1]+align2
            j-=1

    return align1,align2,max_score


# Example
s1 = "ACACACTA"
s2 = "AGCACACA"

a1,a2,score=smith_waterman(s1, s2)
print("Local Alignment 1:",a1)
print("Local Alignment 2:",a2)
print("Score:",score)
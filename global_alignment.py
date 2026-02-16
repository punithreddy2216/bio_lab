def needleman_wunsch(seq1,seq2,match=1,mismatch=-1,gap=-2):
    n=len(seq1)
    m=len(seq2)
    # Create score matrix
    score=[[0]*(m+1)for _ in range(n+1)]
    # Initialize first row and column
    for i in range(n+1):
        score[i][0]=i*gap
    for j in range(m+1):
        score[0][j]=j*gap
    # Fill matrix
    for i in range(1,n+1):
        for j in range(1,m+1):
            if seq1[i-1]==seq2[j-1]:
                diag=score[i-1][j-1] + match
            else:
                diag=score[i-1][j-1] + mismatch
            up=score[i-1][j]+gap
            left=score[i][j-1]+gap
            score[i][j]=max(diag,up,left)
    # Traceback
    align1 = ""
    align2 = ""
    i,j=n,m
    while i>0 and j>0:
        current=score[i][j]
        if seq1[i-1]==seq2[j-1]:
            diag=score[i-1][j-1]+match
        else:
            diag=score[i-1][j-1]+mismatch
        if current==diag:
            align1=seq1[i-1]+align1
            align2=seq2[j-1]+align2
            i-=1
            j-=1
        elif current==score[i-1][j]+gap:
            align1=seq1[i-1]+align1
            align2="-"+align2
            i-=1
        else:
            align1="-"+align1
            align2=seq2[j-1]+align2
            j-=1
    # Finish edges
    while i>0:
        align1=seq1[i-1]+align1
        align2="-"+align2
        i-=1
    while j>0:
        align1="-"+align1
        align2=seq2[j-1]+align2
        j-=1
    return score[n][m],align1,align2
# Example 
seq1="GATTACA"
seq2="GCATGCU"

score,a1,a2=needleman_wunsch(seq1,seq2)

print("Score:",score)
print(a1)
print(a2)
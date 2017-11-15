
import os

MAWFile = './patterns/filtered_hg37.maws'
originalMAWs = []
with open(MAWFile, 'r') as f:
    for motif in f:
        originalMAWs.append(motif.strip())

print 'No. MAWs originally extracted: ' + str(len(originalMAWs))

chromosomes = [str(x) for x in range(1,22)] + ['X', 'Y']
falseMAWsFile = './results/falseMAWs_%s.txt'
for chromosome in chromosomes:
    chrFile = falseMAWsFile % chromosome
    if os.path.exists(chrFile):
        with open(chrFile, 'r') as f:
            for index in f:
                originalMAWs[int(index)] = '-'
    else:
        print 'Warning: Missing file ' + chrFile

confirmedMAWs = set(originalMAWs)
confirmedMAWs.remove('-')

confirmedMAWsFile = './results/confirmedMAWs.maws'
with open(confirmedMAWsFile, 'w') as f:
    j = 0
    for maw in confirmedMAWs:
        if j != 0:
            f.write('\n')
        f.write(maw)
        j += 1

print 'No. MAWs confirmed: ' + str(j)
print 'Confirmed MAWs written to: ' + confirmedMAWsFile

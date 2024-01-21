from math import sqrt

usescipy = True
try:
    from scipy.stats import f
except ImportError:
    usescipy = False
    print('Could not import scipy.stats.f')
    print('You will need to use R to compute P value')
else:
    print('Found scipy.stats.f and will use it to compute P value')

# Dimensions
a = 2
n = __NLOCI__
p = __SAMPLESIZE__

def getDistances(fnprefix):
    d = {}
    for rep in range(n):
        d[rep] = []
        lines = open("%s%d.txt" % (fnprefix, rep+1,), "r").readlines()
        for line in lines[1:]:
            parts = line.strip().split()
            assert len(parts) == 2
            y = float(parts[1])
            d[rep].append(y)
    return d

dsmc = getDistances("smcdists")
dbeast = getDistances("beastdists")

# Calculate overall mean
M = 0.0
ncheck = 0
for i in range(n):
    ncheck += len(dsmc[i])
    M += sum(dsmc[i])
    ncheck += len(dbeast[i])
    M += sum(dbeast[i])
assert ncheck == a*n*p
M /= a*n*p

# Calculate total sum-of-squares SST
SST = 0.0
for i in range(n):
    SST += sum([pow(x - M, 2.) for x in dsmc[i]])
dfT = a*n - 1
MST = SST/dfT

# Calculate within-run sum-of-squares SSWS
SSWS = 0.0
meanrowcol = {}
for i in range(n):
    # Calculate mean for replicate i, treatment SMC
    m = sum(dsmc[i])/len(dsmc[i])
    meanrowcol[(i,0)] = m
    SSWS += sum([pow(x - m, 2.) for x in dsmc[i]])
for i in range(n):
    # Calculate mean for replicate i, treatment BEAST
    m = sum(dbeast[i])/len(dbeast[i])
    meanrowcol[(i,1)] = m
    SSWS += sum([pow(x - m, 2.) for x in dbeast[i]])
dfWS = a*n*(p-1)    
MSWS = SSWS/dfWS

# Calculate treatment sum-of-squares SSA
meancol = {}
meancol[0] = 0.0
meancol[1] = 0.0
for i in range(n):
    meancol[0] += sum(dsmc[i])    
    meancol[1] += sum(dbeast[i])
meancol[0] /= n*p
meancol[1] /= n*p
SSA = n*p*pow(meancol[0] - M, 2.) + n*p*pow(meancol[1] - M, 2.)
dfA = a - 1
MSA = SSA/dfA

# Calculate subject sum-of-squares SSS
SSS = 0.0
meanrow = {}
for i in range(n):
    meanrow[i] = sum(dsmc[i]) + sum(dbeast[i])
    meanrow[i] /= p*a
    SSS += p*a*pow(meanrow[i] - M, 2.)
dfS = n - 1
MSS = SSS/dfS

# Calculate treatment X subject sum-of-squares SSAxS
SSAxS = 0.0
for i in range(n):
    for j in range(2):
        SSAxS += p*pow(meanrowcol[i,j] - meanrow[i] - meancol[j] + M, 2.)
dfAxS = (a-1)*(n-1)
MSAxS = SSAxS/dfAxS

F = MSA/MSAxS
print('\nTest of whether SMC and BEAST differ:')
print('  MSA   = %.5f' % MSA)
print('  MSAxS = %.5f' % MSAxS)
print('  F     = %.5f' % F)
print('  dfA   = %.5f' % dfA)
print('  dfAxS = %.5f' % dfAxS)
if usescipy:
    pvalue = 1.0 - f.cdf(F, dfA, dfAxS)
    print('  P-value = %.5f' % pvalue)
    if pvalue < 0.05:
        print('  Result: significant difference between SMC and BEAST at the 5%% level')
    else:
        print('  Result: no significant difference between SMC and BEAST at the 5%% level')
else:
    print('In R, issue this command to compute P value:')
    print('1 - pf(%.5f, %g, %g)' % (F, dfA, dfAxS))

F = MSS/MSAxS
print('\nTest of whether repeated measures was needed:')
print('  MSS   = %.5f' % MSS)
print('  MSAxS = %.5f' % MSAxS)
print('  F     = %.5f' % F)
print('  dfS   = %.5f' % dfS)
print('  dfAxS = %.5f' % dfAxS)
if usescipy:
    pvalue = 1.0 - f.cdf(F, dfS, dfAxS)
    print('  P-value = %.5f' % pvalue)
    if pvalue < 0.05:
        print('  Result: repeated measures was useful')
    else:
        print('  Result: repeated measures was not necessary')
else:
    print('In R, issue this command to compute P value:')
    print('1 - pf(%.5f, %g, %g)' % (F, dfS, dfAxS))

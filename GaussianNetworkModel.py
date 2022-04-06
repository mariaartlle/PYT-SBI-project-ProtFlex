from prody import *
import matplotlib.pylab as pl

protein = parsePDB('alphafold/2yhw.pdb')
calphas = protein.select('calpha')
gnm = GNM()
gnm.buildKirchhoff(calphas)
gnm.calcModes(None)

# Fem servir el mode 0 per calcular perq és el que descriu millor el moviment

SqFlucts = calcSqFlucts(gnm[0])
NormSqFlucts = SqFlucts / (SqFlucts**2 ).sum()**0.5


pl.figure()
showCrossCorr(gnm)
pl.savefig('crosscorrelations.png')

pl.figure()
showNormedSqFlucts(gnm[0]) # plot sq fluct normalitzades
pl.savefig('normSqFlucts.png')

print(NormSqFlucts)

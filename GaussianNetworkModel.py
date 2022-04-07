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
showNormedSqFlucts(gnm[0]) # plot sq fluct normalitzades
pl.savefig('normSqFlucts.png')


print(NormSqFlucts)


def get_NormSqFluct(protein):
    calphas = protein.select('calpha')
    gnm = GNM()
    gnm.buildKirchhoff(calphas)
    gnm.calcModes(None)
    SqFlucts = calcSqFlucts(gnm[0])
    NormSqFlucts = SqFlucts / (SqFlucts**2 ).sum()**0.5
    return NormSqFlucts


## provar amb internet, no està clar 
def get_CrossCorr_figure(protein):
    pl.figure()
    showCrossCorr(gnm)
    pl.savefig('crosscorrelations.png')
    return

def get_NSF_figure(protein):

    return

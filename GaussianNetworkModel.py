import prody as pr
import matplotlib


# def GNM(protein):
#     protein = pr.parsePDB(protein)
#     calphas = protein.select('calpha')
#     gnm = pr.GNM()
#     gnm.buildKirchhoff(calphas)
#     # gnm.getKirchhoff()
#     gnm.calcModes(None)
#     return gnm
#
# print(GNM("alphafold/Q9UM73.pdb"))


protein = pr.parsePDB("alphafold/2yhw.pdb")
gnm = pr.GNM()
calphas = protein.select('calpha')
gnm.buildKirchhoff(calphas)
gnm.calcModes(50)
# hinges = pr.calcHinges(gnm)
# print(pr.calcDistFlucts(gnm))
pr.calcSqFlucts(gnm)

#pr.showCrossCorr(gnm)
pr.showSqFlucts(gnm[0])

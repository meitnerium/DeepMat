import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
import os
#import pickle
import cirpy
#import pubchempy as pcp
#import sys
#import glob
import time
import numpy as np

# get SDF
#c = pcp.Compound.from_cid(n, record_type='3d')
#pcp.download('SDF','mol.sdf',n,'cid')


def init_latex(f):
    print('OK')

def add_latex(f,pngfile,name,cid):
    #f.write(data['png']+'\\\\')
    f.write('\\begin{table}\n')
    f.write('\\begin{tabular}{ c c c }\n')
    f.write(name+' & ' +str(cid)+' & \includegraphics[scale=1.0]{'+pngfile+'}\n')
    f.write('\\end{tabular}\n')
    f.write('\\end{table}\n')

def end_latex(f):
    print('OK')

def getob(mol):
    print("test mol")
    print(mol)
    print(mol.GetNumAtoms())
    patt = Chem.MolFromSmarts('O')
    print(mol.HasSubstructMatch(patt))
    print(mol.GetSubstructMatch(patt))
    print(mol.GetAtoms())
    #Initialize count value of each atom to zero
    C = 0
    H = 0
    N = 0
    O = 0

    for atom in mol.GetAtoms():
        print(str(atom.GetAtomicNum()))
        if atom.GetAtomicNum() == 6:
            H = H +1
        if atom.GetAtomicNum() == 6:
            C = C +1
        elif atom.GetAtomicNum() == 7:
            N = N +1
        elif atom.GetAtomicNum() == 8:
            O = O + 1
    print("C = "+str(C))
    print("H = "+str(H))
    print("N = "+str(N))
    print("O = "+str(O))
    OB=-1600*(O-C-H/2)/Chem.Descriptors.ExactMolWt(mol)
    print("OB = "+str(OB))
    return OB


def getdescriptors(mol):
    # https://www.nature.com/articles/s41598-018-27344-x
    #patt = Chem.MolFromSmarts('N-N(=O)-O')
    desc = []
    descval = []
    patt = Chem.MolFromSmarts('N-N(=O)-O')
    print("test for NNO2")
    print(mol.HasSubstructMatch(patt))
    fragment = mol.GetSubstructMatches(patt)
    NNO2 = len(fragment)
    print(NNO2)
    print(fragment)
    desc.append('NNO2')
    descval.append(NNO2)

    patt = Chem.MolFromSmarts('C-N(=O)-O')
    print("test for CNO2")
    print(mol.HasSubstructMatch(patt))
    fragment = mol.GetSubstructMatches(patt)
    CNO2 = len(fragment)
    print(CNO2)
    print(fragment)
    desc.append('CNO2')
    descval.append(CNO2)

    patt = Chem.MolFromSmarts('O-N(=O)-O')
    print("test for ONO2")
    print(mol.HasSubstructMatch(patt))
    fragment = mol.GetSubstructMatches(patt)
    ONO2 = len(fragment)
    print(ONO2)
    print(fragment)
    desc.append('ONO2')
    descval.append(ONO2)

    patt = Chem.MolFromSmarts('N(=O)-O')
    print("test for NO2")
    print(mol.HasSubstructMatch(patt))
    fragment = mol.GetSubstructMatches(patt)
    # ONO2 groupment count has 2 NO2
    NO2=len(fragment)-ONO2
    print(NO2)
    print(fragment)
    desc.append('NO2')
    descval.append(NO2)
    
    patt = Chem.MolFromSmarts('C=N-F')
    print("test for CNF")
    print(mol.HasSubstructMatch(patt))
    fragment = mol.GetSubstructMatches(patt)
    CNF = len(fragment)
    print(CNF)
    print(fragment)
    desc.append('CNF')
    descval.append(CNF)

    # TODO : regler ce probleme
    # Problem with COH !!!
    #patt = Chem.MolFromSmarts('C-O-H')
    #frag='C-O-H'
    #print("test for COH")
    #fragment = getfrag(frag,mol)
    #COH = len(fragment)
    #print(mol.HasSubstructMatch(patt))
    #fragment = mol.GetSubstructMatches(patt)
    #COH = len(fragment)
    #print(COH)
    #print(fragment)
    #desc.append('COH')
    #descval.append(COH)

    print("test for NOC")
    frag='N-O-C'
    fragment = getfrag(frag,mol)
    NOC = len(fragment)
    print(NOC)
    desc.append('NOC')
    descval.append(NOC)

    print("test for CNO")
    frag='C=N-O'
    fragment = getfrag(frag,mol)
    CNO = len(fragment)
    print(CNO)
    desc.append('CNO')
    descval.append(CNO)

    print("test for CNN")
    frag='C-N=N'
    fragment = getfrag(frag,mol)
    CNN = len(fragment)
    print(CNN)
    desc.append('CNO')
    descval.append(CNO)

    # TODO : Error
    # error with CNH2
    #print("test for CNH2")
    #frag='C-N(-H)-H'
    #fragment = getfrag(frag,mol)
    #CNH2 = len(fragment)
    #print(CNH2)
    #desc.append('')
    #descval.append()

    print("test for CNOC")
    frag='C-N-O-C'
    fragment = getfrag(frag,mol)
    CNOC = len(fragment)
    print(CNOC)
    desc.append('CNOC')
    descval.append(CNOC)

    print("test for CF")
    frag='C-F'
    fragment = getfrag(frag,mol)
    CF = len(fragment)
    print(CF)
    desc.append('CF')
    descval.append(CF)

    print("test for NO")
    frag='N=O'
    fragment = getfrag(frag,mol)
    NO = len(fragment)
    print(NO)
    desc.append('NO')
    descval.append(NO)

    print("test for CO")
    frag='C=O'
    fragment = getfrag(frag,mol)
    CO = len(fragment)
    print(CO)
    desc.append('CO')
    descval.append(CO)

    #patt = Chem.MolFromSmarts('NOO')
    #print("test for nitro")
    #print(Chem.Fragments.fr_nitro(mol))
    #print(mol.HasSubstructMatch(patt))
    #print(mol.GetSubstructMatch(patt))
    # [N+](=O)[O-]
    return desc,descval

def getfrag(frag,mol):
    patt = Chem.MolFromSmarts(frag)
    return mol.GetSubstructMatches(patt)


def bagofbounds(mol2,maxentry):
    from chemml.chem import Molecule
    #caffeine_smiles = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
    #caffeine_smarts = '[#6]-[#7]1:[#6]:[#7]:[#6]2:[#6]:1:[#6](=[#8]):[#7](:[#6](=[#8]):[#7]:2-[#6])-[#6]'
    #caffeine_inchi = 'InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3'
    mol = Molecule(Chem.MolToSmiles(mol2),'smiles')
    mol.hydrogens('add')
    mol.to_xyz(optimizer='MMFF', mmffVariant='MMFF94s', maxIters=300) # 'UFF'
    print(mol)
    mol.visualize()
    #mol.visualize()



    from chemml.datasets import load_xyz_polarizability
    from chemml.chem import BagofBonds
    #coordinates, y = load_xyz_polarizability(mol)
    bob = BagofBonds(const= 1.0)
    features = bob.represent(mol)
    print(features)
    print("number of entry")
    print(features.shape[1])
    if features.shape[1] > maxentry:
        maxentry = features.shape[1]
    return features,maxentry


def sumoverbound(mol):
    #N-O, N:O, N-N, N=O, N=N, N:N, N#N, C-N, C-C, C-H, C:N, C:C, C-F, C-O, C=O, C=N, C=C, H-O, H-N, F-N.
    listbound = ['N-O', 'N:O', 'N-N', 'N=O', 'N=N', 'N:N', 'N#N', 'C-N', 'C-C', 'C-H', 'C:N', 'C:C', 
    'C-F', 'C-O', 'C=O', 'C=N', 'C=C', 'H-O', 'H-N', 'F-N']
    for bound in listbound:
        mol = rdkit.Chem.rdmolops.AddHs(mol)
        patt = Chem.MolFromSmiles(bound)
        fragment = mol.GetSubstructMatches(patt)
        print(bound+" : "+str(len(fragment)))
    #patt = Chem.MolFromSmarts('N-O')
    #fragment = mol.GetSubstructMatches(patt)
    #len(fragment+" : "+str(len(fragment) )


import csv
cids=[] # pubchem compound ID vector
names=[] # name vector
pes=[] #Explosives power vector
f=open('table.tex','w')
init_latex(f)
descriptor_file = open('descriptor_file.csv', mode='w')
descriptor_writer = csv.writer(descriptor_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
maxentry=1
bob = np.array([1,1])
m=0
with open('PE.csv') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='"')
    n=0
    for row in spamreader:
        if n != 0:
            time.sleep(5)
            cids.append(row[5])
            names.append(row[0])
            pes.append(row[3])
            print(row[0]+" "+row[3]+" "+row[5])
            #c = pcp.Compound.from_cid(n, record_type='3d')
            sdf='data/'+row[5]+'.sdf'
            if os.path.isfile(sdf) == False:
                download=True
                print("trying \"pcp.download('SDF',sdf,"+row[5]+",'cid',overwrite=True)\"")
                try:
                    pcp.download('SDF',sdf,row[5],'cid')
                except:
                    print("download not working")
                    download=False
            else:
                download=True
            if download:
                suppl = Chem.SDMolSupplier(sdf)
                for mol in suppl:
                    pngf = "data/"+row[5]+".png"
                    Draw.MolToFile(mol,pngf)
                    add_latex(f,pngf,row[0],row[5])

                    ob = getob(mol)
                    desc,descval = getdescriptors(mol)
                    desc.insert(0,'cid')
                    descval.insert(0,row[5])
                    bob2insert,maxentry = bagofbounds(mol,maxentry)
                    print("maxentry")
                    print(maxentry)
                    bob.resize([m+1,maxentry])
                    bob[m,:bob2insert.shape[1]] = bob2insert
                    m = m + 1
                    #bob2inesrt = np.array(bob2insert)
                    #bob.append(bob2insert)

                    print("Length of bob")
                    #sumbound = sumoverbound(mol)
                    #with open('eggs.csv', 'a', newline='') as csvfile:
                    #    spamwriter = csv.writer(csvfile, delimiter=',',
                    #        quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    #        spamwriter.writerow(['cid'] + [str(cid)])
                    #        spamwriter.writerow(['Spam', 'Lovely Spam', 'Wonderful Spam'])
                    if n==1:
                        descriptor_writer.writerow(desc)
                    descriptor_writer.writerow(descval)
        n=n+1
#max=0

#for i in range(len(bob[:,0])):
#    if len(bob[i,:]) > max:
#        max=len(bob[i,:])

print(bob.shape)
for i in range(len(bob[:,0])):
    if len(bob[i,:]) < maxentry:
        toinsert=zeros(maxentry-len(bob[i,:]))
        bob.append(toinsert)

end_latex(f)
f.close()

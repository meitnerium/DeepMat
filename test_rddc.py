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

    patt = Chem.MolFromSmarts('O-N(=O)-O')
    print("test for ONO2")
    print(mol.HasSubstructMatch(patt))
    fragment = mol.GetSubstructMatches(patt)
    ONO2 = len(fragment)
    print(ONO2)
    print(fragment)

    patt = Chem.MolFromSmarts('N(=O)-O')
    print("test for NO2")
    print(mol.HasSubstructMatch(patt))
    fragment = mol.GetSubstructMatches(patt)
    # ONO2 groupment count has 2 NO2
    NO2=len(fragment)-ONO2
    print(NO2)
    print(fragment)
    
    patt = Chem.MolFromSmarts('C=N-F')
    print("test for CNF")
    print(mol.HasSubstructMatch(patt))
    fragment = mol.GetSubstructMatches(patt)
    CNF = len(fragment)
    print(CNF)
    print(fragment)

    patt = Chem.MolFromSmarts('C-O-H')
    frag='C-O-H'
    print("test for COH")
    #print(mol.HasSubstructMatch(patt))
    fragment = mol.GetSubstructMatches(patt)
    COH = len(fragment)
    print(COH)
    print(fragment)

    print("test for NOC")
    frag='N-O-C'
    fragment = getfrag(frag,mol)
    NOC = len(fragment)
    print(NOC)

    print("test for CNO")
    frag='C=N-O'
    fragment = getfrag(frag,mol)
    CNO = len(fragment)
    print(CNO)

    print("test for CNN")
    frag='C-N=N'
    fragment = getfrag(frag,mol)
    CNN = len(fragment)
    print(CNN)

    print("test for CNH2")
    frag='C-N(-H)-H'
    fragment = getfrag(frag,mol)
    CNH2 = len(fragment)
    print(CNH2)

    print("test for CNOC")
    frag='C-N-O-C'
    fragment = getfrag(frag,mol)
    CNOC = len(fragment)
    print(CNOC)

    print("test for CF")
    frag='C-F'
    fragment = getfrag(frag,mol)
    CF = len(fragment)
    print(CF)

    print("test for NO")
    frag='N=O'
    fragment = getfrag(frag,mol)
    CF = len(fragment)
    print(CF)

    print("test for CO")
    frag='C=O'
    fragment = getfrag(frag,mol)
    CO = len(fragment)
    print(CO)

    #patt = Chem.MolFromSmarts('NOO')
    #print("test for nitro")
    #print(Chem.Fragments.fr_nitro(mol))
    #print(mol.HasSubstructMatch(patt))
    #print(mol.GetSubstructMatch(patt))
    # [N+](=O)[O-]

def getfrag(frag,mol):
    patt = Chem.MolFromSmarts(frag)
    return mol.GetSubstructMatches(patt)

import csv
cids=[] # pubchem compound ID vector
names=[] # name vector
pes=[] #Explosives power vector
f=open('table.tex','w')
init_latex(f)
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
                    nno2 = getdescriptors(mol)
                    #with open('eggs.csv', 'a', newline='') as csvfile:
                    #    spamwriter = csv.writer(csvfile, delimiter=',',
                    #        quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    #        spamwriter.writerow(['cid'] + [str(cid)])
                    #        spamwriter.writerow(['Spam', 'Lovely Spam', 'Wonderful Spam'])

        n=n+1


end_latex(f)
f.close()


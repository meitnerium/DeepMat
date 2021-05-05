from chemml.chem import Molecule
caffeine_smiles = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
caffeine_smarts = '[#6]-[#7]1:[#6]:[#7]:[#6]2:[#6]:1:[#6](=[#8]):[#7](:[#6](=[#8]):[#7]:2-[#6])-[#6]'
caffeine_inchi = 'InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3'
mol = Molecule(caffeine_smiles, input_type='smiles')
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


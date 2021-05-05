import chemml

from chemml.datasets import load_xyz_polarizability
from chemml.chem import BagofBonds
coordinates, y = load_xyz_polarizability()
bob = BagofBonds(const= 1.0)
features = bob.represent(coordinates)

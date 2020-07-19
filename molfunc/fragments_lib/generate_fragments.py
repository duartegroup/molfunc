from molfunc.molecules import Molecule


class Fragment:

    def __init__(self, name, smiles, aliases):

        self.name = name
        self.smiles = smiles
        self.aliases = aliases + [name.lower()]


fragments = [Fragment('NMe2', 'CN([*])C', ['dimethylamino']),
             Fragment('NO2', 'O=[N+]([*])[O-]', ['nitro']),
             Fragment('H', '[H][*]', ['hydro']),
             Fragment('Me', 'C[*]', ['ch3', 'methyl']),
             Fragment('Et', 'CC[*]', ['c2h5', 'ch2ch3', 'ethyl']),
             Fragment('CH2OH', '[*]CO', []),
             Fragment('OH', 'O[*]', ['hydroxy']),
             Fragment('OMe', '[*]OC', ['methoxy', 'och3']),
             Fragment('OPh', '[*]OC1=CC=CC=C1', ['phenoxy', 'oc6h5']),
             Fragment('OtBu', '[*]OC(C)(C)C', ['oc3h9']),
             Fragment('CO2Me', '[*]C(OC)=O', ['co2ch3']),
             Fragment('CO2Et', '[*]C(OCC)=O', ['co2ch2ch3','co2c2h5']),
             Fragment('Boc', '[*]C(OC(C)(C)C)=O', ['co2tbu', 'co2(ch3)3']),
             Fragment('Cl', 'Cl[*]', ['chloride', 'chloro']),
             Fragment('F', 'F[*]', ['fluoride', 'fluro']),
             Fragment('I', 'I[*]', ['iodide', 'iodo']),
             Fragment('Br', 'Br[*]', ['bromide', 'bromo']),
             Fragment('CN', 'N#C[*]', ['cyano', 'nitrile']),
             Fragment('iPr', 'CC([*])C', ['isopropyl', 'ch(ch3)2']),
             Fragment('tBu', 'CC([*])(C)C', ['tertbutyl', 'c(ch3)3']),
             Fragment('Ph', '[*]C1=CC=CC=C1', ['phenyl', 'c6h5']),
             Fragment('NH2', 'N[*]', ['amino']),
             Fragment('CF3', '[*]C(F)(F)F', []),
             Fragment('Ac', '[*]C(C)=O', ['acetyl', 'come', 'coch3']),
             Fragment('TMS', 'C[Si](C)([*])C',
                      ['trimethylsilyl', 'sime3', 'si(ch3)3']),
             Fragment('Bn', '[*]CC1=CC=CC=C1',
                      ['benzyl', 'ch2ph', 'ch2c6h5']),
             Fragment('Bz', 'O=C([*])C1=CC=CC=C1',
                      ['benzoyl', 'coph', 'coc6h5']),
             Fragment('Mes', 'CC1=CC(C)=CC(C)=C1[*]', ['mesityl']),
             Fragment('Ms', 'O=S(C)([*])=O',
                      ['methanesulfonyl', 'mesyl', 'so2me', 'so2ch3']),
             Fragment('Tf', 'O=S(C(F)(F)F)([*])=O', ['triflate', 'so2cf3'])]


if __name__ == '__main__':

    for fragment in fragments:

        # Generate the molecules with F rather than R atoms
        molecule = Molecule(name=fragment.name,
                            smiles=fragment.smiles.replace('*', 'F'))

        # Change back for printing
        for atom in molecule.atoms:
            if atom.label == 'F':
                atom.label = 'R'

        # Print the SMILES string and the possible aliases in the xyz files
        title_line = f'{fragment.smiles} {",".join(fragment.aliases)}'
        molecule.print_xyz_file(title_line=title_line)

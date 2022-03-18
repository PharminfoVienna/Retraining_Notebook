from rdkit import Chem
from standardiser.neutralise import run as neutralise


class Standardization:

    """
    Takes an input file and standardises the molecules obtained from the file.

    Uses Atkinson's standardiser
    """

    def __init__(self):
        """
        Create list of standardised molecules.
        :ivar standardised_mol_list: list of standardised molecules
        :param file: sdf file with the molecules
        """

    def standardise(self, mol):
        """
        Performs the standardization.
        Following functions are carried out:
            - Removing salts
            - Checking for possible components (might have been generated due to the salt removal)
                - Remove duplicate components in the same molecule
                - Check for inorganic compounds and remove molecule if yes -> cannot be predicted
                - Check for metals not forming salts in the molecule and remove molecule if yes -> cannot be predicted
                - Get the largest component (by heavy atom count) if there are still more than one component.
                - Set properties of original molecule to component -> retain all the input information
            - Check single components for inorganic compounds and remove molecule if yes -> cannot be predicted
            - Check single components for metals not forming salts in the molecule and remove molecule if yes ->
            cannot be predicted
            - neutralise molecule (performed no matter if components or single component)
        :param mol:
        :return: Tuple indicating successful standardization and molecule or error message.
        """
        mol = self.remove_salt_metals(mol)  # remove salt metals -> might produce some additional components
        # check for mixtures
        components = [c for c in Chem.GetMolFrags(mol, asMols=True,
                                                sanitizeFrags=False)]  # gets the components of the given molecule

        if len(components) > 1:  # handle mixtures
            inorganic_frag = []
            metal_frag = []
            properties = mol.GetPropsAsDict()  # save original properties to assign later on to chosen component

            try:
                # remove duplicate components possibly obtained after removing salts
                components = self.remove_duplicate_components(components)  # might fail due to neutralisation step
            except:
                return False, "Molecule could not be neutralised"

            # flag inorganics
            for comp_index in range(len(components)):  # iterate over the components
                if self.is_inorganic(components[comp_index]):  # flag inorganic components
                    inorganic_frag.append(comp_index)

            # delete all inorganic components
            if inorganic_frag:  # check if we even have inorganic parts
                inorganic_frag.reverse()  # reverse so we don't mess up the indices
                for comp_index in inorganic_frag:
                    components.pop(comp_index)

            # check if there is still more than 1 component.
            # Metals are considered as mixture, even though they are removed in the process
            # 'not used in online version -> just remove largest component'
            # if len(components) > 1:
            #     mixture.append(molecules[i])

            # flag metals
            for comp_index in range(len(components)):  # iterate over the components
                if self.contains_metal(components[comp_index]):
                    metal_frag.append(comp_index)

            # remove all metal components
            if metal_frag:  # check if we even have metals
                metal_frag.reverse()
                for comp_index in metal_frag:
                    components.pop(comp_index)

            # choose the remaining component to be our new molecule
            if len(components) > 1:  # check if there is still more than 1 component
                mol = self.get_largest_component(components)  # get largest of the remaining components
            elif len(components) > 0:  # we have more than zero, but no more than one -> exactly one component left
                mol = components[0]  # assign remaining component to molecule
            else:  # no component left -> go on with the next molecule
                return False, "Molecule contained only metals or inorganic compounds"

            for prop, value in properties.items():  # assign properties to left over component
                mol.SetProp(prop, str(value))

        else:  # do cleaning for only one component
            if self.contains_metal(mol):
                return False, "Molecule contains metal"
            if self.is_inorganic(mol):
                return False, "Molecule is inorganic"
        try:
            # neutralise molecule
            mol.UpdatePropertyCache()
            mol = neutralise(mol)  # neutralisation might throw an error due to Sanity_check
            return True, mol
        except:
            return False, "Molecule could not be neutralised"

    def is_silane(self, mol):
        """
        Checks if the molecule contains at least one carbon atom. If yes, the compound is considered to be organic (False).
        Otherwise the compound will be considered as inorganic and will be flagged as such (True).
        """
        foundCarbon = False
        foundSilicon = False
        for a in list(mol.GetAtoms()):
            if a.GetAtomicNum() == 14:  
                foundSilicon = True
            if a.GetAtomicNum() == 6: 
                foundCarbon = True
        if foundCarbon == True and foundSilicon == True:
            return True
        return False

    def is_inorganic(self, mol):
        """
        Checks if the molecule contains at least one carbon atom. If yes, the compound is considered to be organic (False).
        Otherwise the compound will be considered as inorganic and will be flagged as such (True).
        """
        for a in list(mol.GetAtoms()):
            if a.GetAtomicNum() == 6:  # returns upon the first carbon encounter -> we are not inorganic
                return False
            if a.GetAtomicNum() == 14:
                return False
        return True


    def remove_salt_metals(self, mol):
        """
        Removes the metals Li, Na, K, Mg, Ca from the molecule.
        """
        salt_metals = [3, 11, 12, 19, 20]  # li, Na, Mg, K, Ca
        to_remove = []  # maintain a list of atom indices which represent a metal

        # iterate over atoms and search for salt metals
        for a in list(mol.GetAtoms()):
            if a.GetAtomicNum() in salt_metals:
                if a.GetDegree() > 1:  # check if metal is bound to counter-ion -> ignore. Ignore only for Mg and Ca!
                    continue
                to_remove.append(a.GetIdx())  # only track the index, because molecule would change if we remove now

        # remove in a separate step, because otherwise we mess up the atom order / index
        if to_remove:
            mol = Chem.EditableMol(mol)
            to_remove.reverse()  # start from the back to not alter the indices
            for idx in to_remove:
                mol.RemoveAtom(idx)
            mol = mol.GetMol()  # get a normal molecule back
        return mol


    def contains_metal(self, mol):
        """
        Checks if the molecule contains a metal atom which we don't want.
        """
        allowed_atoms = [1, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53]
        for a in list(mol.GetAtoms()):
            if a.GetAtomicNum() not in allowed_atoms:  # returns upon encounter of the first not allowed atom
                return True
        return False


    def get_largest_component(self, components):
        """
        Keep the component with the highest number of heavy atoms.
        """
        largest_component = 0  # keep track of component index
        num_atoms = 0  # keep track of component size
        for i, comp in enumerate(components):
            if comp.GetNumHeavyAtoms() > num_atoms:  # check if molecule has more heavy atoms. If the same, keep first
                num_atoms = comp.GetNumHeavyAtoms()
                largest_component = i  # udpate index of largest component
        return components[largest_component]  # return largest component by index


    def remove_duplicate_components(self, components):
        """
        Get a hash-code of the molecules and remove duplicate components.
        """
        inchi_keys = {}
        for comp in components:
            comp.UpdatePropertyCache()
            comp = neutralise(comp)  # neutralise at this step to remove differently charged duplicates
            inchi = Chem.inchi.MolToInchi(comp)
            inchi_key = Chem.inchi.InchiToInchiKey(inchi)  # generate a hash code to filter by -> InChIKey

            if inchi_key not in inchi_keys.keys():  # checks for duplicates
                inchi_keys[inchi_key] = comp
        return [comp for comp in inchi_keys.values()]  # returns all the unique fragments
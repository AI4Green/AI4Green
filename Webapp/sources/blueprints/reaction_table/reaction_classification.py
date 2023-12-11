from rdkit import Chem
from typing import List, Dict


def define_substructs() -> Dict:
    """
    returns general SMART strings for pattern matching during reaction classification
    """
    substructs = {
        "grignard": "*[Mg][Cl,Br,I]",
        "ketone": "*[#6]C(=O)[#6]*",
        "aldehyde": "[#6X2](=O)",
        "alkene": "[#6]=[#6]*",
        "carboxylic_acid": "*C(=O)[O;!$(O(C)[#6])]",
        "amine": "*N",
        "ester": "*C(=O)O[#6]",
        "ether": "*[#6]O[#6]*",
        "boronic_acid": "B(O)O",
        "aryl_halide": "a[F,Cl,Br,I]",
        "alkyl_halide": "C[Cl,Br,I]",
        "alcohol": "*C[O;!$(O(C)[!#6])]",
        "palladium": "[Pd]",
        "ruthenium": "[Ru]",
        "baylis-hillman": "*C(O)C(=C)*",
        "amide": "*C(=O)N*",
        "biaryl": "*a-a*",
        "aryl_amine": "*Na",
        "aryl_alkene": "aC=C*",
    }

    return substructs

def define_reaction_criteria() -> Dict:
    """
    returns functional group criteria for each reaction class included in the solvent surfer
    """

    reaction_criteria = {
        "Amide bond formation": {
            "Reactants": ["carboxylic_acid", "amine"],
            "Products": ["amide"],
            "Agents": [],
        },
        "Alcohol oxidation": {
            "Reactants": ["alcohol"],
            "Products": ["ketone"],
            "Agents": [],
        },
        "Alkene metathesis": {
            "Reactants": ["alkene"],
            "Products": ["alkene"],
            "Agents": ["ruthenium"],
        },
        "Baylis-Hillman": {
            "Reactants": ["aldehyde", "alkene"],
            "Products": ["baylis-hillman"],
            "Agents": [],
        },
        "Buchwald-Hartwig coupling": {
            "Reactants": ["amine", "aryl_halide"],
            "Products": ["aryl_amine"],
            "Agents": ["palladium"],
        },
        "Ester hydrolysis": {
            "Reactants": ["ester"],
            "Products": ["carboxylic_acid", "alcohol"],
            "Agents": [],
        },
        "Grignard": {"Reactants": ["grignard"], "Products": ["alcohol"], "Agents": []},
        "Heck cross-coupling": {
            "Reactants": ["alkene", "aryl_halide"],
            "Products": ["aryl_alkene"],
            "Agents": ["palladium"],
        },
        "SNAr/SN2": {
            "Reactants": ["aryl_halide"],
            "Products": ["ether", "amine"],
            "Agents": [],
        },
        "Suzuki-Miyaura coupling": {
            "Reactants": ["boronic_acid", "aryl_halide"],
            "Products": ["biaryl"],
            "Agents": ["palladium"],
        },
    }

    return reaction_criteria

def classify_reaction(reactants: List, products: List) -> (str, List):
    """
        Classifies reactions by identifying important functional groups in reactants and products.
        Only checks for reaction classes included in the solvent surfer

        Args:
            reactants: List, List of reactants as RDKit molecules
            products: List, List of products as RDKit molecules

        Returns:
            reaction_class: Str, class of reaction identified to functional groups. Returns 'Other' if class not recognised.
            reaction_groups.keys(): List of reaction classes checked

        """

    # define functional groups
    substructs = define_substructs()

    # define reactive groups for each class
    reaction_groups = define_reaction_criteria()

    reactant_groups = []
    product_groups = []
    reaction_class = "Other"

    # match functional groups in reactants and products to dict above
    for num, mols in enumerate([reactants, products]):
        for mol in mols:
            mol.UpdatePropertyCache()
            Chem.GetSymmSSSR(mol)
            mol.GetRingInfo().NumRings()
            Chem.Kekulize(mol)
            Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_DEFAULT)
            for name, substruct in substructs.items():
                query = Chem.MolFromSmarts(substruct)
                query.UpdatePropertyCache()
                Chem.GetSymmSSSR(query)
                mol.GetRingInfo().NumRings()
                Chem.Kekulize(query)
                Chem.SetAromaticity(query)
                if mol.HasSubstructMatch(query):
                    if num == 0:
                        reactant_groups.append(name)
                    else:
                        product_groups.append(name)

    for r_class, r_dict in reaction_groups.items():
        r = r_dict["Reactants"]
        p = r_dict["Products"]
        #a = r_dict['Agents']

        if len([x for x in r if x in reactant_groups]) == len(r):
            if len([y for y in p if y in product_groups]) == len(p):
                reaction_class = r_class

    reaction_groups["Other"] = {}
    return reaction_class, reaction_groups.keys()

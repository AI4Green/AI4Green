def calculate_atom_economy(reactants_molecular_weights_sum: float, products_molecular_weights_sum: float) -> float:
    """
    Calculate the atom economy of a chemical reaction.

    Atom economy is a measure of the efficiency of a chemical reaction in terms of how well atoms are utilized. 
    It is calculated as the ratio of the molecular weight of the desired product to the total molecular weight of all reactants, multiplied by 100.

    Args:
        reactants_molecular_weights (float): The sum of molecular weights of the reactants.
        products_molecular_weights_sum (float): The sum of molecular weights of the products.
    Returns:
        float: The atom economy percentage.
    """
    
    if reactants_molecular_weights_sum == 0:
        raise ValueError("The sum of molecular weights of reactants cannot be zero.")

    atom_economy = (products_molecular_weights_sum / reactants_molecular_weights_sum) * 100
    return atom_economy


def calculate_carbon_efficiency(reactants_carbon_moles: float, products_carbon_moles_sum: float) -> float:
    """
    Calculate the carbon efficiency of a chemical reaction.

    Carbon efficiency is a measure of how effectively carbon atoms are utilized in a chemical reaction. 
    It is calculated as the ratio of the number of carbon atoms in the desired product to the total number of carbon atoms in all reactants, multiplied by 100.
    (mol product × C atoms in product) / ∑(mol reactant × C atoms in reactant) × 100 

    Args:
        reactants_carbon_moles (float): The total number of carbon moles in the reactants.
        products_carbon_moles_sum (float): The sum of carbon moles in the products.
    Returns:
        float: The carbon efficiency percentage.
    """
    
    if reactants_carbon_moles == 0:
        raise ValueError("The total number of carbon moles in reactants cannot be zero.")

    carbon_efficiency = (products_carbon_moles_sum / reactants_carbon_moles) * 100
    return carbon_efficiency

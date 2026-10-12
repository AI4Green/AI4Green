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



def calculate_atom_economy(
    reactants_molecular_weight_sum: float,
    desired_product_molecular_weight: float,
) -> float:
    """
    Calculate the atom economy (AE) of a chemical reaction.

    AE measures how many of the reactant atoms end up in the desired product,
    assuming 100% yield and exactly stoichiometric amounts (no excess).

    By-products are what is lost, so only the desired product goes on top.
        AE = MW(desired product) / sum(coefficient x MW(reactant)) x 100

    Args:
        reactants_molecular_weight_sum (float): Sum of reactant molecular
            weights (g/mol), each multiplied by its stoichiometric coefficient
            in the balanced equation.
        desired_product_molecular_weight (float): Molecular weight (g/mol) of
            the desired product only, not by-products.

    Returns:
        float: Atom economy as a percentage.

    Raises:
        ValueError: If either input is zero or negative.
    """
    if reactants_molecular_weight_sum <= 0:
        raise ValueError("Sum of reactant molecular weights must be positive.")

    if desired_product_molecular_weight <= 0:
        raise ValueError("Desired product molecular weight must be positive.")

    return (desired_product_molecular_weight / reactants_molecular_weight_sum) * 100


def calculate_carbon_efficiency(
    reactants_carbon_moles_sum: float,
    desired_product_carbon_moles: float,
) -> float:
    """
    Calculate the carbon efficiency (CE) of a chemical reaction.

    CE is the percentage of carbon in the reactants that ends up in the desired
    product. Unlike atom economy, it uses actual amounts, so it accounts for
    yield and excess reactants.
        CE = (mol product x C atoms in product) / sum(mol reactant x C atoms in reactant) x 100

    Args:
        reactants_carbon_moles_sum (float): Sum over all reactants of (moles of
            reactant x number of carbon atoms in reactant).
        desired_product_carbon_moles (float): Moles of desired product obtained x
            number of carbon atoms in the desired product.

    Returns:
        float: Carbon efficiency as a percentage.

    Raises:
        ValueError: If the reactant sum is zero or negative, or the product value
            is negative.
    """
    if reactants_carbon_moles_sum <= 0:
        raise ValueError("Sum of reactant carbon moles must be positive.")

    if desired_product_carbon_moles < 0:
        raise ValueError("Desired product carbon moles cannot be negative.")

    return (desired_product_carbon_moles / reactants_carbon_moles_sum) * 100


def calculate_reaction_mass_efficiency(
    reactants_mass_sum: float,
    desired_product_mass: float,
) -> float:
    """
    Calculate the reaction mass efficiency (RME) of a chemical reaction.

    RME is the percentage of the total mass of reactants that ends up in the
    desired product. It accounts for yield and excess reactants.
        RME = (mass of desired product) / (sum of masses of all reactants) x 100

    Args:
        reactants_mass_sum (float): Sum of masses of all reactants (g).
        desired_product_mass (float): Mass of the desired product obtained (g).

    Returns:
        float: Reaction mass efficiency as a percentage.

    Raises:
        ValueError: If the reactant sum is zero or negative, or the product mass
            is negative.
    """
    if reactants_mass_sum <= 0:
        raise ValueError("Sum of reactant masses must be positive.")

    if desired_product_mass < 0:
        raise ValueError("Desired product mass cannot be negative.")

    return (desired_product_mass / reactants_mass_sum) * 100

    
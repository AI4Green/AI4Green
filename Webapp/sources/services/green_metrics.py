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

    
def calculate_mass_intensity(
    total_mass_used: float,
    desired_product_mass: float,
) -> float:
    """
    Calculate the mass intensity (MI) of a chemical process.

    MI is the total mass of materials used per unit mass of desired product.
    Unlike RME, it includes everything that goes into the process, not just
    the reactants.
        MI = total mass used (excluding water) / mass of desired product

    Total mass used includes reactants, reagents, catalysts, solvents and
    workup/purification materials. Water is excluded.

    Args:
        total_mass_used (float): Total mass of all materials used in the
            process, excluding water (any consistent mass unit, e.g. g or kg).
        desired_product_mass (float): Mass of the desired product obtained
            (same unit as total_mass_used).

    Returns:
        float: Mass intensity as a ratio (dimensionless, e.g. kg/kg).
            An ideal process has MI = 1.

    Raises:
        ValueError: If either input is zero or negative.
    """
    if total_mass_used <= 0:
        raise ValueError("Total mass used must be positive.")

    if desired_product_mass <= 0:
        raise ValueError("Desired product mass must be positive.")

    return total_mass_used / desired_product_mass

    
def calculate_mass_productivity(
    total_mass_used: float,
    desired_product_mass: float,
) -> float:
    """
    Calculate the mass productivity (MP) of a chemical process.

    MP is the percentage of the total mass of materials used that ends up in
    the desired product. It is the reciprocal of mass intensity, expressed as
    a percentage.
        MP = mass of desired product / total mass used (excluding water) x 100
           = 100 / MI

    Args:
        total_mass_used (float): Total mass of all materials used in the
            process, excluding water (any consistent mass unit, e.g. g or kg).
        desired_product_mass (float): Mass of the desired product obtained
            (same unit as total_mass_used).

    Returns:
        float: Mass productivity as a percentage. An ideal process has
            MP = 100%.

    Raises:
        ValueError: If either input is zero or negative.
    """
    mass_intensity = calculate_mass_intensity(total_mass_used, desired_product_mass)
    return 100 / mass_intensity
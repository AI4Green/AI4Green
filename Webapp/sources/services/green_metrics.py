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


def calculate_e_factor(
    total_mass_used: float,
    desired_product_mass: float,
) -> float:
    """
    Calculate the E-factor of a chemical process.

    The E-factor is the mass of waste generated per unit mass of desired
    product. Waste is everything used in the process that does not end up
    as the desired product.
        E-factor = (total mass used - mass of desired product) / mass of desired product
                 = MI - 1

    Args:
        total_mass_used (float): Total mass of all materials used in the
            process, excluding water (any consistent mass unit, e.g. g or kg).
        desired_product_mass (float): Mass of the desired product obtained
            (same unit as total_mass_used).

    Returns:
        float: E-factor as a ratio (dimensionless, e.g. kg/kg). An ideal
            process has E-factor = 0; lower is better.

    Raises:
        ValueError: If either input is zero or negative, or if the product
            mass is greater than the total mass used.
    """
    if desired_product_mass > total_mass_used:
        raise ValueError("Desired product mass cannot exceed total mass used.")

    mass_intensity = calculate_mass_intensity(total_mass_used, desired_product_mass)
    return mass_intensity - 1


def calculate_effective_mass_yield(
    non_benign_mass_used: float,
    desired_product_mass: float,
) -> float:
    """
    Calculate the effective mass yield (EMY) of a chemical process.

    EMY is the mass of desired product as a percentage of the mass of
    non-benign materials used. Materials judged benign (e.g. water, NaCl)
    are left out of the denominator. Aocess that relies on benign
    materials scores higher.
        EMY = mass of desired product / mass of non-benign materials used x 100

    Deciding which materials count as benign is done before calling this
    function, using toxicity or hazard data. That judgement is not made here.

    Args:
        non_benign_mass_used (float): Total mass of all non-benign materials
            used in the process (any consistent mass unit, e.g. g or kg).
        desired_product_mass (float): Mass of the desired product obtained
            (same unit as non_benign_mass_used).

    Returns:
        float: Effective mass yield as a percentage. Can exceed 100% if most
            of the materials used are benign.

    Raises:
        ValueError: If the non-benign mass is zero or negative, or the product
            mass is negative.
    """
    if non_benign_mass_used <= 0:
        raise ValueError("Non-benign mass used must be positive.")

    if desired_product_mass < 0:
        raise ValueError("Desired product mass cannot be negative.")

    return (desired_product_mass / non_benign_mass_used) * 100


def calculate_optimum_efficiency(
    reaction_mass_efficiency: float,
    atom_economy: float,
) -> float:
    """
    Calculate the optimum efficiency (OE) of a chemical reaction.

    OE compares the actual mass efficiency of a reaction with its theoretical
    best. It shows how close a reaction comes to its atom economy once real
    yield and excess reactants are taken into account.
        OE = RME / AE x 100

    If yield is 100% and reactants are used in exact stoichiometric amounts,
    RME equals AE and OE is 100%. OE therefore cannot distinguish between
    reactions when only structures are known.

    Args:
        reaction_mass_efficiency (float): Reaction mass efficiency as a
            percentage, e.g. from calculate_reaction_mass_efficiency.
        atom_economy (float): Atom economy as a percentage, e.g. from
            calculate_atom_economy.

    Returns:
        float: Optimum efficiency as a percentage. An ideal reaction has
            OE = 100%.

    Raises:
        ValueError: If atom economy is zero or negative, reaction mass
            efficiency is negative, or reaction mass efficiency is greater
            than atom economy.
    """
    if atom_economy <= 0:
        raise ValueError("Atom economy must be positive.")

    if reaction_mass_efficiency < 0:
        raise ValueError("Reaction mass efficiency cannot be negative.")

    if reaction_mass_efficiency > atom_economy:
        raise ValueError("Reaction mass efficiency cannot exceed atom economy.")

    return (reaction_mass_efficiency / atom_economy) * 100
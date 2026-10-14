import math


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
        ValueError: If either input is zero or negative, or if the desired
            product molecular weight is greater than the reactant sum.
    """
    if reactants_molecular_weight_sum <= 0:
        raise ValueError("Sum of reactant molecular weights must be positive.")

    if desired_product_molecular_weight <= 0:
        raise ValueError("Desired product molecular weight must be positive.")

    if desired_product_molecular_weight > reactants_molecular_weight_sum:
        if not math.isclose(
            desired_product_molecular_weight, reactants_molecular_weight_sum
        ):
            raise ValueError(
                "Desired product molecular weight cannot exceed sum of "
                "reactant molecular weights."
            )

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
        ValueError: If the reactant sum is zero or negative, the product value
            is negative, or the product value is greater than the reactant sum.
    """
    if reactants_carbon_moles_sum <= 0:
        raise ValueError("Sum of reactant carbon moles must be positive.")

    if desired_product_carbon_moles < 0:
        raise ValueError("Desired product carbon moles cannot be negative.")

    if desired_product_carbon_moles > reactants_carbon_moles_sum:
        raise ValueError(
            "Desired product carbon moles cannot exceed sum of reactant "
            "carbon moles."
        )

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
        reactants_mass_sum (float): Sum of masses of all reactants (any
            consistent mass unit, e.g. g or kg).
        desired_product_mass (float): Mass of the desired product obtained
            (same unit as reactants_mass_sum).

    Returns:
        float: Reaction mass efficiency as a percentage.

    Raises:
        ValueError: If the reactant sum is zero or negative, the product mass
            is negative, or the product mass is greater than the reactant sum.
    """
    if reactants_mass_sum <= 0:
        raise ValueError("Sum of reactant masses must be positive.")

    if desired_product_mass < 0:
        raise ValueError("Desired product mass cannot be negative.")

    if desired_product_mass > reactants_mass_sum:
        raise ValueError("Desired product mass cannot exceed sum of reactant masses.")

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
        ValueError: If either input is zero or negative, or if the product
            mass is greater than the total mass used.
    """
    if total_mass_used <= 0:
        raise ValueError("Total mass used must be positive.")

    if desired_product_mass <= 0:
        raise ValueError("Desired product mass must be positive.")

    if desired_product_mass > total_mass_used:
        raise ValueError("Desired product mass cannot exceed total mass used.")

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
        ValueError: If either input is zero or negative, or if the product
            mass is greater than the total mass used.
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
    are left out of the denominator. A process that relies on benign
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
        if not math.isclose(reaction_mass_efficiency, atom_economy):
            raise ValueError("Reaction mass efficiency cannot exceed atom economy.")

    return (reaction_mass_efficiency / atom_economy) * 100


def calculate_process_mass_intensity(
    total_mass_used_including_water: float,
    desired_product_mass: float,
) -> float:
    """
    Calculate the process mass intensity (PMI) of a chemical process.

    PMI is the total mass of materials used per unit mass of desired product,
    including water. It has the same form as mass intensity (MI), but MI
    excludes water and PMI includes it.
        PMI = total mass used (including water) / mass of desired product

    Total mass used includes reactants, reagents, catalysts, solvents,
    workup/purification materials and water.

    Args:
        total_mass_used_including_water (float): Total mass of all materials
            used in the process, including water (any consistent mass unit,
            e.g. g or kg).
        desired_product_mass (float): Mass of the desired product obtained
            (same unit as total_mass_used_including_water).

    Returns:
        float: Process mass intensity as a ratio (dimensionless, e.g. kg/kg).
            An ideal process has PMI = 1.

    Raises:
        ValueError: If either input is zero or negative, or if the product
            mass is greater than the total mass used.
    """
    if total_mass_used_including_water <= 0:
        raise ValueError("Total mass used must be positive.")

    if desired_product_mass <= 0:
        raise ValueError("Desired product mass must be positive.")

    if desired_product_mass > total_mass_used_including_water:
        raise ValueError("Desired product mass cannot exceed total mass used.")

    return total_mass_used_including_water / desired_product_mass


def calculate_renewables_intensity(
    renewable_mass_used: float,
    desired_product_mass: float,
) -> float:
    """
    Calculate the renewables intensity (RI) of a chemical process.

    RI is the mass of renewably derivable materials used per unit mass of
    desired product. RI alone does not show what share of the inputs is
    renewable; use renewables percentage for that.
        RI = mass of renewably derivable materials used / mass of desired product

    Deciding which materials count as renewably derivable is done before
    calling this function. That judgement is not made here.

    Args:
        renewable_mass_used (float): Total mass of renewably derivable
            materials used in the process (any consistent mass unit,
            e.g. g or kg).
        desired_product_mass (float): Mass of the desired product obtained
            (same unit as renewable_mass_used).

    Returns:
        float: Renewables intensity as a ratio (dimensionless, e.g. kg/kg).

    Raises:
        ValueError: If the renewable mass is negative, or the product mass is
            zero or negative.
    """
    if renewable_mass_used < 0:
        raise ValueError("Renewable mass used cannot be negative.")

    if desired_product_mass <= 0:
        raise ValueError("Desired product mass must be positive.")

    return renewable_mass_used / desired_product_mass


def calculate_renewables_percentage(
    renewables_intensity: float,
    process_mass_intensity: float,
) -> float:
    """
    Calculate the renewables percentage (RP) of a chemical process.

    RP is the percentage of the total mass of materials used that is
    renewably derivable.
        RP = RI / PMI x 100
           = mass of renewably derivable materials / total mass used x 100

    RI and PMI must be calculated over the same set of materials (the same
    process boundaries). If PMI includes water or workup materials, the
    renewable mass behind RI must be judged over those same materials.

    Args:
        renewables_intensity (float): Renewables intensity, e.g. from
            calculate_renewables_intensity.
        process_mass_intensity (float): Process mass intensity, e.g. from
            calculate_process_mass_intensity.

    Returns:
        float: Renewables percentage. 100% means all materials used are
            renewably derivable.

    Raises:
        ValueError: If PMI is zero or negative, RI is negative, or RI is
            greater than PMI.
    """
    if process_mass_intensity <= 0:
        raise ValueError("Process mass intensity must be positive.")

    if renewables_intensity < 0:
        raise ValueError("Renewables intensity cannot be negative.")

    if renewables_intensity > process_mass_intensity:
        raise ValueError("Renewables intensity cannot exceed process mass intensity.")

    return (renewables_intensity / process_mass_intensity) * 100


def calculate_waste_intensity(
    total_mass_input: float,
    total_waste_mass: float,
) -> float:
    """
    Calculate the waste intensity (WI) of a chemical process.

    WI is the mass of waste produced per unit mass of total input. Unlike the
    E-factor, it is expressed per unit of input, not per unit of product.
        WI = total waste produced / total mass input

    Total mass input should cover the same materials as PMI, including water.

    Args:
        total_mass_input (float): Total mass of all materials used in the
            process, including water (any consistent mass unit, e.g. g or kg).
        total_waste_mass (float): Total mass of waste produced (same unit as
            total_mass_input).

    Returns:
        float: Waste intensity as a ratio (dimensionless, e.g. kg/kg). An
            ideal process has WI = 0.

    Raises:
        ValueError: If the total mass input is zero or negative, the waste
            mass is negative, or the waste mass is greater than the total
            mass input.
    """
    if total_mass_input <= 0:
        raise ValueError("Total mass input must be positive.")

    if total_waste_mass < 0:
        raise ValueError("Total waste mass cannot be negative.")

    if total_waste_mass > total_mass_input:
        raise ValueError("Total waste mass cannot exceed total mass input.")

    return total_waste_mass / total_mass_input


def calculate_waste_percentage(
    total_mass_input: float,
    total_waste_mass: float,
) -> float:
    """
    Calculate the waste percentage (WP) of a chemical process.

    WP is the percentage of the total mass input that ends up as waste.
        WP = total waste produced / total mass input x 100

    Total mass input should cover the same materials as PMI, including water,
    so that WP and PMI describe the same process boundaries.

    Args:
        total_mass_input (float): Total mass of all materials used in the
            process, including water (any consistent mass unit, e.g. g or kg).
        total_waste_mass (float): Total mass of waste produced (same unit as
            total_mass_input).

    Returns:
        float: Waste percentage. An ideal process has WP = 0%.

    Raises:
        ValueError: If the total mass input is zero or negative, the waste
            mass is negative, or the waste mass is greater than the total
            mass input.
    """
    if total_mass_input <= 0:
        raise ValueError("Total mass input must be positive.")

    if total_waste_mass < 0:
        raise ValueError("Total waste mass cannot be negative.")

    if total_waste_mass > total_mass_input:
        raise ValueError("Total waste mass cannot exceed total mass input.")

    return (total_waste_mass / total_mass_input) * 100


def calculate_turnover_number(
    catalyst_moles: float,
    product_moles: float,
) -> float:
    """
    Calculate the turnover number (TON) of a catalysed reaction.

    TON is the number of moles of product formed per mole of catalyst. It
    shows how many times each catalyst unit is used before the reaction
    stops. Higher values mean less catalyst is needed.
        TON = moles of product formed / moles of catalyst

    Args:
        catalyst_moles (float): Moles of catalyst used.
        product_moles (float): Moles of desired product formed.

    Returns:
        float: Turnover number as a ratio (dimensionless, mol/mol).

    Raises:
        ValueError: If the catalyst moles are zero or negative, or the
            product moles are negative.
    """
    if catalyst_moles <= 0:
        raise ValueError("Catalyst moles must be positive.")

    if product_moles < 0:
        raise ValueError("Product moles cannot be negative.")

    return product_moles / catalyst_moles


def calculate_space_time_yield(
    reactor_volume: float,
    reaction_time: float,
    desired_product_mass: float,
) -> float:
    """
    Calculate the space-time yield (STY) of a chemical process.

    STY is the mass of desired product made per unit reactor volume per unit
    time. It measures how productively the reactor is used.
        STY = mass of desired product / (reactor volume x reaction time)

    The result has units of (mass unit) / (volume unit x time unit), set by
    the inputs, e.g. g / (L x h) = g L^-1 h^-1. Values are only comparable
    if the same units are used.

    Args:
        reactor_volume (float): Volume of the reactor (any volume unit,
            e.g. L or m^3).
        reaction_time (float): Reaction time (any time unit, e.g. h).
        desired_product_mass (float): Mass of the desired product obtained
            (any mass unit, e.g. g or kg).

    Returns:
        float: Space-time yield, in mass / (volume x time).

    Raises:
        ValueError: If the reactor volume or reaction time is zero or
            negative, or the product mass is negative.
    """
    if reactor_volume <= 0:
        raise ValueError("Reactor volume must be positive.")

    if reaction_time <= 0:
        raise ValueError("Reaction time must be positive.")

    if desired_product_mass < 0:
        raise ValueError("Desired product mass cannot be negative.")

    return desired_product_mass / (reactor_volume * reaction_time)


def calculate_percentage_yield(
    limiting_reactant_moles: float,
    product_moles: float,
    product_per_limiting_reactant: float = 1.0,
) -> float:
    """
    Calculate the percentage yield of a chemical reaction.

    Percentage yield is the moles of desired product isolated as a percentage
    of the maximum possible from the limiting reactant.
        Yield = mol product / (mol limiting reactant x ratio) x 100

    The ratio is the stoichiometric coefficient of the product divided by
    that of the limiting reactant, e.g. 0.5 for 2 A -> B. The default of 1
    covers 1:1 reactions.

    Yields above 100% are returned, not rejected. They usually mean the
    product was wet or impure.

    Args:
        limiting_reactant_moles (float): Moles of limiting reactant used.
        product_moles (float): Moles of desired product isolated.
        product_per_limiting_reactant (float): Moles of product formed per
            mole of limiting reactant in the balanced equation. Defaults
            to 1.0.

    Returns:
        float: Percentage yield.

    Raises:
        ValueError: If the limiting reactant moles or the stoichiometric
            ratio are zero or negative, or the product moles are negative.
    """
    if limiting_reactant_moles <= 0:
        raise ValueError("Limiting reactant moles must be positive.")

    if product_per_limiting_reactant <= 0:
        raise ValueError("Stoichiometric ratio must be positive.")

    if product_moles < 0:
        raise ValueError("Product moles cannot be negative.")

    theoretical_product_moles = limiting_reactant_moles * product_per_limiting_reactant
    return (product_moles / theoretical_product_moles) * 100


def calculate_percentage_conversion(
    initial_limiting_reactant_mass: float,
    final_limiting_reactant_mass: float,
) -> float:
    """
    Calculate the percentage conversion of the limiting reactant.

    Conversion is the percentage of the limiting reactant used up in the
    reaction. It does not show whether the reactant became the desired
    product.
        Conversion = 100 - (final mass / initial mass x 100)
                   = (initial mass - final mass) / initial mass x 100

    Args:
        initial_limiting_reactant_mass (float): Mass of limiting reactant at
            the start (any consistent mass unit, e.g. g or kg).
        final_limiting_reactant_mass (float): Mass of limiting reactant left
            at the end (same unit as initial_limiting_reactant_mass).

    Returns:
        float: Percentage conversion.

    Raises:
        ValueError: If the initial mass is zero or negative, the final mass
            is negative, or the final mass is greater than the initial mass.
    """
    if initial_limiting_reactant_mass <= 0:
        raise ValueError("Initial limiting reactant mass must be positive.")

    if final_limiting_reactant_mass < 0:
        raise ValueError("Final limiting reactant mass cannot be negative.")

    if final_limiting_reactant_mass > initial_limiting_reactant_mass:
        raise ValueError("Final limiting reactant mass cannot exceed initial mass.")

    return (
        (initial_limiting_reactant_mass - final_limiting_reactant_mass)
        / initial_limiting_reactant_mass
        * 100
    )


def calculate_percentage_selectivity(
    percentage_conversion: float,
    percentage_yield: float,
) -> float:
    """
    Calculate the percentage selectivity of a chemical reaction.

    Selectivity is the percentage of the consumed limiting reactant that
    became the desired product.
        Selectivity = % yield / % conversion x 100

    This matches mol desired product / mol reactant consumed when the yield
    accounts for stoichiometry, e.g. from calculate_percentage_yield with
    the correct ratio.

    Args:
        percentage_conversion (float): Conversion as a percentage, e.g. from
            calculate_percentage_conversion.
        percentage_yield (float): Yield as a percentage, e.g. from
            calculate_percentage_yield.

    Returns:
        float: Percentage selectivity. An ideal reaction has 100%.

    Raises:
        ValueError: If conversion is zero or negative (selectivity is
            undefined at zero conversion), yield is negative, or yield is
            greater than conversion.
    """
    if percentage_conversion <= 0:
        raise ValueError("Percentage conversion must be positive.")

    if percentage_yield < 0:
        raise ValueError("Percentage yield cannot be negative.")

    if percentage_yield > percentage_conversion:
        if not math.isclose(percentage_yield, percentage_conversion):
            raise ValueError("Percentage yield cannot exceed percentage conversion.")

    return (percentage_yield / percentage_conversion) * 100

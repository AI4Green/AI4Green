import pytest
from sources import services


def test_calculate_atom_economy():
    atom_economy = services.green_metrics.calculate_atom_economy(
        reactants_molecular_weight_sum=106.12,
        desired_product_molecular_weight=88.11,
    )
    assert atom_economy == pytest.approx(83.03, rel=1e-3)


def test_calculate_carbon_efficiency():
    carbon_efficiency = services.green_metrics.calculate_carbon_efficiency(
        reactants_carbon_moles_sum=4.0,
        desired_product_carbon_moles=3.2,
    )
    assert carbon_efficiency == pytest.approx(80.0)


def test_calculate_reaction_mass_efficiency():
    reaction_mass_efficiency = (
        services.green_metrics.calculate_reaction_mass_efficiency(
            reactants_mass_sum=50.0,
            desired_product_mass=40.0,
        )
    )
    assert reaction_mass_efficiency == pytest.approx(80.0)


def test_calculate_mass_intensity():
    mass_intensity = services.green_metrics.calculate_mass_intensity(
        total_mass_used=200.0,
        desired_product_mass=40.0,
    )
    assert mass_intensity == pytest.approx(5.0)


def test_calculate_mass_productivity():
    mass_productivity = services.green_metrics.calculate_mass_productivity(
        total_mass_used=200.0,
        desired_product_mass=40.0,
    )
    assert mass_productivity == pytest.approx(20.0)


def test_calculate_e_factor():
    e_factor = services.green_metrics.calculate_e_factor(
        total_mass_used=200.0,
        desired_product_mass=40.0,
    )
    assert e_factor == pytest.approx(4.0)


def test_calculate_effective_mass_yield():
    effective_mass_yield = services.green_metrics.calculate_effective_mass_yield(
        non_benign_mass_used=80.0,
        desired_product_mass=40.0,
    )
    assert effective_mass_yield == pytest.approx(50.0)


def test_calculate_optimum_efficiency():
    optimum_efficiency = services.green_metrics.calculate_optimum_efficiency(
        reaction_mass_efficiency=40.0,
        atom_economy=80.0,
    )
    assert optimum_efficiency == pytest.approx(50.0)


def test_calculate_process_mass_intensity():
    process_mass_intensity = services.green_metrics.calculate_process_mass_intensity(
        total_mass_used_including_water=500.0,
        desired_product_mass=40.0,
    )
    assert process_mass_intensity == pytest.approx(12.5)


def test_calculate_renewables_intensity():
    renewables_intensity = services.green_metrics.calculate_renewables_intensity(
        renewable_mass_used=100.0,
        desired_product_mass=40.0,
    )
    assert renewables_intensity == pytest.approx(2.5)


def test_calculate_renewables_percentage():
    renewables_percentage = services.green_metrics.calculate_renewables_percentage(
        renewables_intensity=2.5,
        process_mass_intensity=12.5,
    )
    assert renewables_percentage == pytest.approx(20.0)


def test_calculate_waste_intensity():
    waste_intensity = services.green_metrics.calculate_waste_intensity(
        total_mass_input=200.0,
        total_waste_mass=160.0,
    )
    assert waste_intensity == pytest.approx(0.8)


def test_calculate_waste_percentage():
    waste_percentage = services.green_metrics.calculate_waste_percentage(
        total_mass_input=200.0,
        total_waste_mass=160.0,
    )
    assert waste_percentage == pytest.approx(80.0)


def test_calculate_turnover_number():
    turnover_number = services.green_metrics.calculate_turnover_number(
        catalyst_moles=0.01,
        product_moles=0.8,
    )
    assert turnover_number == pytest.approx(80.0)


def test_calculate_space_time_yield():
    space_time_yield = services.green_metrics.calculate_space_time_yield(
        reactor_volume=0.25,
        reaction_time=2.0,
        desired_product_mass=50.0,
    )
    assert space_time_yield == pytest.approx(100.0)


def test_calculate_percentage_yield():
    percentage_yield = services.green_metrics.calculate_percentage_yield(
        limiting_reactant_moles=1.0,
        product_moles=0.4,
        product_per_limiting_reactant=0.5,
    )
    assert percentage_yield == pytest.approx(80.0)


def test_calculate_percentage_conversion():
    percentage_conversion = services.green_metrics.calculate_percentage_conversion(
        initial_limiting_reactant_mass=10.0,
        final_limiting_reactant_mass=1.0,
    )
    assert percentage_conversion == pytest.approx(90.0)


def test_calculate_percentage_selectivity():
    percentage_selectivity = services.green_metrics.calculate_percentage_selectivity(
        percentage_conversion=90.0,
        percentage_yield=72.0,
    )
    assert percentage_selectivity == pytest.approx(80.0)

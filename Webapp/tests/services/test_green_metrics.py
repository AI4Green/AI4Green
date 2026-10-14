from sources import services


def test_calculate_atom_economy():
    atom_economy = services.green_metrics.calculate_atom_economy(reactants_molecular_weight_sum=46.0, desired_product_molecular_weight=46.0)

    assert atom_economy == 100.0
    assert 5==4, 

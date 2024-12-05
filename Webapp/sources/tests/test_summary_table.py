import flask_testing
from sources import app
import unittest


class MyTestCase(flask_testing.TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        app.config['LOGIN_DISABLED'] = True
        return app

    def setUp(self):
        pass

    def send_to_summary(self, amountunit, volumeunit, massunit, solventVolumeUnit, productMassUnit,
                        reactants, reactantMolecularWeights,
                        reactantDensities, reactantConcentrations, reactantEquivalents,
                        reactantAmounts, roundedReactantAmounts, reactantVolumes,
                        roundedReactantVolumes, reactantMasses, roundedReactantMasses,
                        reactantHazards, reactantPhysicalForms,
                        reactantMassSum,reactantMolecularWeightSum,
                        numberOfReagents,
                        reagents, reagentTableNumbers, reagentMolecularWeights,
                        reagentDensities, reagentConcentrations, reagentEquivalents,
                        reagentAmounts, roundedReagentAmounts, reagentVolumes,
                        roundedReagentVolumes, reagentMasses, roundedReagentMasses,
                        reagentHazards, reagentPhysicalForms,
                        reagentMassSum, reagentMolecularWeightSum,
                        numberOfSolvents, solvents, solventColors, solventTableNumbers,
                        roundedSolventConcentrations, solventConcentrations,
                        solventVolumes, solventPhysicalForms, solventHazards,
                        products, productTableNumbers,
                        productMasses, roundedProductMasses,
                        productMolecularWeights,
                        productHazards, productPhysicalForms,
                        mainProductTableNumber, reactantPrimaryKeys,
                        reagentPrimaryKeys, solventPrimaryKeys, productPrimaryKeys, js_summary_table_data,
                        reactionSmiles):
        return self.client.post(
            '/_summary',
            data=dict(amountUnit=amountunit, volumeUnit=volumeunit, massUnit=massunit,
                      solventVolumeUnit=solventVolumeUnit, productMassUnit=productMassUnit,
                      reactants=reactants,
                      reactantMolecularWeights=reactantMolecularWeights,
                      reactantDensities=reactantDensities,
                      reactantConcentrations=reactantConcentrations,
                      reactantEquivalents=reactantEquivalents,
                      reactantAmounts=reactantAmounts,
                      roundedReactantAmounts=roundedReactantAmounts,
                      reactantVolumes=reactantVolumes,
                      roundedReactantVolumes=roundedReactantVolumes,
                      reactantMasses=reactantMasses,
                      roundedReactantMasses=roundedReactantMasses,
                      reactantHazards=reactantHazards,
                      reactantPhysicalForms=reactantPhysicalForms,
                      reactantMassSum=reactantMassSum, reactantMolecularWeightSum=reactantMolecularWeightSum,
                      numberOfReagents=numberOfReagents,
                      reagents=reagents,
                      reagentTableNumbers=reagentTableNumbers,
                      reagentMolecularWeights=reagentMolecularWeights,
                      reagentDensities=reagentDensities,
                      reagentConcentrations=reagentConcentrations,
                      reagentEquivalents=reagentEquivalents,
                      reagentAmounts=reagentAmounts,
                      roundedReagentAmounts=roundedReagentAmounts,
                      reagentVolumes=reagentVolumes,
                      roundedReagentVolumes=roundedReagentVolumes,
                      reagentMasses=reagentMasses,
                      roundedReagentMasses=roundedReagentMasses,
                      reagentHazards=reagentHazards,
                      reagentPhysicalForms=reagentPhysicalForms,
                      reagentMassSum=reagentMassSum, reagentMolecularWeightSum=reagentMolecularWeightSum,
                      numberOfSolvents=numberOfSolvents,
                      solvents=solvents,
                      solventColors=solventColors,
                      solventTableNumbers=solventTableNumbers,
                      roundedSolventConcentrations=roundedSolventConcentrations,
                      solventConcentrations=solventConcentrations,
                      solventVolumes=solventVolumes,
                      solventPhysicalForms=solventPhysicalForms,
                      solventHazards=solventHazards,
                      products=products,
                      productTableNumbers=productTableNumbers,
                      productMasses=productMasses,
                      roundedProductMasses=roundedProductMasses,
                      productMolecularWeights=productMolecularWeights,
                      productHazards=productHazards,
                      productPhysicalForms=productPhysicalForms,
                      mainProductTableNumber=mainProductTableNumber,
                      reactantPrimaryKeys=reactantPrimaryKeys,
                      reagentPrimaryKeys=reagentPrimaryKeys,
                      solventPrimaryKeys=solventPrimaryKeys,
                      productPrimaryKeys=productPrimaryKeys,
                      js_summary_table_data=js_summary_table_data,
                      reactionSmiles=reactionSmiles
                      ),
            follow_redirects=True
        )

    def test_summary(self):
        """Here we test that the data for is passed correctly"""

        response = self.send_to_summary('mol', 'mL', 'mg', 'ml', 'mg',
                                        'reactant',
                                        0.3, 0.45, 0.33, 1.2,
                                        0.34, 0.45, 0.66, 0.6, 0.45, 0.45,
                                        'H330', 'Gas',
                                        0.766, 0.44,
                                        1, 'reagent', 1, 0.3, 0.45, 0.33, 1.2,
                                        0.34, 0.45, 0.66, 0.6, 0.45, 0.45,
                                        'H334', 'Dense solid',
                                        0.755, 0.788,
                                        1, 'Ethanol', 'green', 1, 0.3, 0.45, 0.33,
                                        'Volatile liquid', 'H315',
                                        10, 1, 0.44, 0.33, 0.54,
                                        'H240', 'Dusty Solid',
                                        1, [2], [6],
                                        [3], [1], "no data", 'CCO>>CCN')
        print('test', response.data)
        self.assertIn(b'H330 Fatal if inhaled', response.data)


if __name__ == '__main__':
    unittest.main()

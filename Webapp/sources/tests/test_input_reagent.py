from unittest import mock, main
from sources import app, db
import flask_testing
from database_setup import test_database_create
from pony.orm import db_session


@mock.patch('sources.input_reagent.routes.current_user')
class MyTestCase(flask_testing.TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def setUp(self):
        pass

    def tearDown(self):
        restore_db()

    def input_reagent(self, commonName, iupac, hPhrase, cas, density, boilingPoint,
                      molWeight, concentration, amountUnit, volumeUnit, massUnit,
                      reactantPrimaryKeys, limitingReactantTableNumber,
                      numberOfReactants, reactants, reactantMolecularWeights,
                      reactantDensities, reactantConcentrations, reactantEquivalents,
                      reactantAmounts, roundedReactantAmounts, reactantVolumes,
                      roundedReactantVolumes, reactantMasses, roundedReactantMasses,
                      reactantHazards, reactantPhysicalForms, numberOfReagents,
                      reagents, reagentTableNumbers, reagentMolecularWeights,
                      reagentDensities, reagentConcentrations, reagentEquivalents,
                      reagentAmounts, roundedReagentAmounts, reagentVolumes,
                      roundedReagentVolumes, reagentMasses, roundedReagentMasses,
                      reagentHazards, reagentPhysicalForms, solventVolumeUnit,
                      numberOfSolvents, solvents, solventColors, solventTableNumbers,
                      roundedSolventConcentrations, solventConcentrations,
                      solventVolumes, solventPhysicalForms, solventHazards,
                      numberOfProducts, products, productTableNumbers, productMolecularWeights,
                      productAmountUnit, productMassUnit, productAmounts,
                      roundedProductAmounts, productMasses, roundedProductMasses,
                      productPhysicalForms, productHazards, productPrimaryKeys):
        return self.client.post(
            '/_input_reagent',
            data=dict(
                commonName=commonName,
                iupac=iupac,
                hPhrase=hPhrase,
                cas=cas,
                density=density,
                boilingPoint=boilingPoint,
                molWeight=molWeight,
                concentration=concentration,
                amountUnit=amountUnit,
                volumeUnit=volumeUnit,
                massUnit=massUnit,
                limitingReactantTableNumber=limitingReactantTableNumber,
                numberOfReactants=numberOfReactants,
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
                reactantPrimaryKeys=reactantPrimaryKeys,

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

                solventVolumeUnit=solventVolumeUnit,
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
                numberOfProducts=numberOfProducts,
                productTableNumbers=productTableNumbers,
                productMolecularWeights=productMolecularWeights,
                productAmountUnit=productAmountUnit,
                productMassUnit=productMassUnit,
                productAmounts=productAmounts,
                roundedProductAmounts=roundedProductAmounts,
                productMasses=productMasses,
                roundedProductMasses=roundedProductMasses,
                productPhysicalForms=productPhysicalForms,
                productHazards=productHazards,
                productPrimaryKeys=productPrimaryKeys,
            ),
            follow_redirects=True
        )

    def test_successful_input_reagent(self, mock_user):
        """First, we input a new reagent using the tested route"""
        mock_user.email = 'PI@test.com'
        response = self.input_reagent('this', 'this', 'H111', '111', '1', '100',
                                      '18', '1', 'mol', 'ml', 'g', '2, 3', 1,
                                      '1', 'oxygen', '16',
                                      '100', '200', '1',
                                      '10', '10', '0.1',
                                      '0.1', '2', '2',
                                      'H100', 'gas', '0',
                                      '', '', '',
                                      '', '', '',
                                      '', '', '',
                                      '', '', '',
                                      '', '', '',
                                      '0', '', '', '',
                                      '', '', '',
                                      '', '', '',
                                      'that', '2', '23',
                                      'mol', 'mg', '1.012',
                                      '1', '2.0123', '2',
                                      'solid', 'H331', '9')
        """Then, we delete the new reagent in the database to make the tests independent from each other"""
        #self.rg_con.execute("DELETE from reagent_table WHERE iupac = 'this';")
        #self.rg_con.commit()
        print(response.data)
        self.assertIn(b'your reagent is added to the database', response.data)  # checks the response message

    @db_session
    def test_failed_input_reagent(self, mock_user):
        """First, we manually input a new reagent"""
        mock_user.email = 'PI@test.com'
        new_compound = db.NovelCompound(name='this', workbook=1)
        # self.rg_con.execute("INSERT INTO reagent_table VALUES ('this', 'this', 'H111', '111', '1', '100', '18', '1')")
        # self.rg_con.commit()
        """Then, we try to input the new reagent using the tested route"""
        response = self.input_reagent('this', 'this', 'H111', '111', '1', '100',
                                      '18', '1', 'mol', 'ml', 'g',
                                      '1', 'oxygen', '16',
                                      '100', '200', '1',
                                      '10', '10', '0.1',
                                      '0.1', '2', '2',
                                      'H100', 'gas', '0',
                                      '', '', '',
                                      '', '', '',
                                      '', '', '',
                                      '', '', '',
                                      '', '', '',
                                      '0', '', '', '',
                                      '', '', '', ''
                                      '', '', '',
                                      'that', '2', '23',
                                      'mol', 'mg', '5', '5', '1',
                                      '1', '2', '2',
                                      'solid', 'H331', '9')
        """Here we delete the new reagent in the database to make the tests independent from each other"""
        self.assertIn(b'your reagent is already in the database', response.data)  # checks the response message


if __name__ == '__main__':
    main()

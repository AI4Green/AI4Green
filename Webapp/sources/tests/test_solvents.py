import flask_testing
from sources import app, db
from pony.orm import select, db_session
import unittest


class TestSolvents(flask_testing.TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        app.config['LOGIN_DISABLED'] = True
        return app

    def setUp(self):
        pass

    def send_to_solvents(self, solvent, number):
        # stops redirect loop error being raised when running from test loader
        return self.client.post(
            '/_solvents',
            data=dict(solvent=solvent, number=number, workgroup='Test-Workgroup', workbook='Test-Workbook'))

    def test_sever(self):
        """Here we test the response works"""
        response = self.send_to_solvents('Ethanol', 1)
        self.assertEqual(response.status_code, 200)

    @db_session
    def test_data(self):
        """Here we test that the data for ethanol is found in the response."""
        response = self.send_to_solvents('Ethanol', 1)
        ethanol_pk = str(select(x.id for x in db.Compound if x.name == 'Ethanol').first()).encode('utf-8')
        expected_response = b'{"alert_message":"","flag":"hazard-acceptable","hazards":"H225-H302-H319-H371","new_solvent":false,"num":"1",' \
                            b'"primary_key":' + ethanol_pk + b',"solvent":"Ethanol"}'
        self.assertIn(expected_response, response.data)

    def test_solvent_doesnt_exist(self):
        """Here we test an error is raised when a solvent not present is added"""
        response = self.send_to_solvents('imaginary_solvent', 1)
        expected_response = b''

        self.assertIn(expected_response, response.data)

    def test_no_input(self):
        """Here we test the values are set to 0 if the default select value is passed"""
        response = self.send_to_solvents('-select-', 1)
        expected_response = b''

        self.assertIn(expected_response, response.data)

    @db_session
    def test_different_number(self):
        """Here we test the response is correct using a different number"""
        response = self.send_to_solvents('Ethanol', 4)
        ethanol_pk = str(select(x.id for x in db.Compound if x.name == 'Ethanol').first()).encode('utf-8')
        expected_response = b'{"alert_message":"","flag":"hazard-acceptable","hazards":"H225-H302-H319-H371","new_solvent":false,"num":"4",' \
                            b'"primary_key":' + ethanol_pk + b',"solvent":"Ethanol"}'
        self.assertIn(expected_response, response.data)

    @db_session
    def test_different_solvent(self):
        """Here we test the response is correct using a different solvent"""
        response = self.send_to_solvents('Water', 1)
        water_pk = str(select(x.id for x in db.Compound if x.name == 'Water').first()).encode('utf-8')
        expected_response = b'{"alert_message":"","flag":"hazard-acceptable","hazards":"Not Hazardous","new_solvent":false,"num":"1",' \
                            b'"primary_key":'+ water_pk + b',"solvent":"Water"}'
        self.assertIn(expected_response, response.data)


if __name__ == '__main__':
    unittest.main()

from sources import app, db
from pony.orm import db_session, select
import unittest


class TestSolvents(unittest.TestCase):
    def setUp(self):
        app.config.from_object('config.TestConfig')
        app.config['LOGIN_DISABLED'] = True
        self.app = app.test_client()

    def send_to_solvents(self, solvent, number, workgroup='Test-Workgroup', workbook='Test-Workbook'):
        # stops redirect loop error being raised when running from test loader
        return self.app.post(
            '/_solvents',
            data=dict(solvent=solvent, number=number, workgroup=workgroup, workbook=workbook))

    def test_server(self):
        """Here we test the response works"""
        response = self.send_to_solvents('64-17-5', 1)
        self.assertEqual(response.status_code, 200)

    @db_session
    def test_data(self):
        """Here we test that the data for ethanol is found in the response."""
        ethanol_pk = str(select(x.id for x in db.Compound if x.cas == '64-17-5').first()).encode('utf-8')
        response = self.send_to_solvents('64-17-5', 1)
        expected_response = \
            b'{"alert_message":"","flag":"hazard-acceptable","hazards":"H225","new_solvent":false,"num":"1","primary_key":'\
            + ethanol_pk + b',"solvent":"Ethanol"}\n'
        self.assertIn(expected_response, response.data)

    def test_solvent_doesnt_exist(self):
        """Here we test an error is raised when a solvent not present is added"""
        response = self.send_to_solvents('12-34-5', 1)
        expected_response = b'{"alert_message":"The solvent is not found in the available databases. You can add it as ' \
                            b'a new compound.","flag":"white","hazards":"","new_solvent":true,"num":"1",' \
                            b'"primary_key":"","solvent":"12-34-5"}\n'
        self.assertIn(expected_response, response.data)

    @db_session
    def test_different_solvent(self):
        """Here we test the response is correct using a different solvent"""
        water_pk = str(select(x.id for x in db.Compound if x.cas == '7732-18-5').first()).encode('utf-8')
        response = self.send_to_solvents('7732-18-5', 1)
        expected_response = \
            b'{"alert_message":"","flag":"hazard-acceptable","hazards":"No hazard codes found","new_solvent":false,"num":"1",' \
            b'"primary_key":' + water_pk + b',"solvent":"Water"}\n'

        self.assertIn(expected_response, response.data)


if __name__ == '__main__':
    unittest.main()

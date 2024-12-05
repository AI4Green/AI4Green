import flask_testing
from sources import app
import unittest


class TestReagents(flask_testing.TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        app.config['LOGIN_DISABLED'] = True
        return app

    def setUp(self):
        pass

    def send_to_reagents(self, reagent, number, workbook='Test-Workbook', workgroup='Test-Workgroup'):
        # stops redirect loop error being raised when running from test loader
        return self.client.post(
            '/_reagents',
            data=dict(reagent=reagent, number=number, workbook=workbook, workgroup=workgroup))

    def test_sever(self):
        """Here we test the response works"""
        response = self.send_to_reagents('Aniline', 1)
        self.assertEqual(response.status_code, 200)

    def test_data(self):
        """Here we test that the data for aniline is found in the response."""
        response = self.send_to_reagents('Aniline', 1)
        aniline_data = b'{"concentration":null,"density":1.022,"hazards":' \
                       b'"H301-H311-H317-H318-H331-H341-H351-H372-H400",' \
                       b'"molWeight":93.13,"name":"Aniline","number":"1","primary_key":'
        self.assertIn(aniline_data, response.data)

    def test_reagent_doesnt_exist(self):
        """Here we test it is recognised if a reagent is not found"""
        response = self.send_to_reagents('imaginary_reagent', 1)
        self.assertIn(b'{"identifiers":[],"reagent":"imaginary_reagent","reagent_not_found":true}', response.data)

    def test_different_number(self):
        """Here we test the response is correct using a different number"""
        response = self.send_to_reagents('Aniline', 4)
        expected_response = b'{"concentration":null,"density":1.022,"hazards":' \
                            b'"H301-H311-H317-H318-H331-H341-H351-H372-H400",' \
                            b'"molWeight":93.13,"name":"Aniline","number":"4","primary_key":'
        self.assertIn(expected_response, response.data)

    def test_different_reagent(self):
        """Here we test the response is correct using a different reagent"""
        response = self.send_to_reagents('Methanol', 1)
        expected_response = b'{"concentration":null,"density":0.792,"hazards":"H225-H301-H311-H331-H370",' \
                            b'"molWeight":32.042,"name":"Methanol","number":"1","primary_key":'
        self.assertIn(expected_response, response.data)


if __name__ == '__main__':
    unittest.main()

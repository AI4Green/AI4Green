import sources.auxiliary as aux
from sources import app
import unittest


class Test_GetData(unittest.TestCase):
    def setUp(self):
        self.app = app.test_client()

    def tearDown(self):
        pass

    def test_sharp_replace(self):
        SMILE = "CsharpN"
        SMILE = aux.replace_symbols(SMILE)
        SMILEnew = "C%23N"
        self.assertEqual(SMILE, SMILEnew)

    def test_plus_replace(self):
        SMILE = "CplusN"
        SMILE = aux.replace_symbols(SMILE)
        SMILEnew = "C%2BN"
        self.assertEqual(SMILE, SMILEnew)

    def test_minus_replace(self):
        SMILE = "CminusN"
        SMILE = aux.replace_symbols(SMILE)
        SMILEnew = "C%2DN"
        self.assertEqual(SMILE, SMILEnew)

    def test_hash_replace(self):
        SMILE = "C#N"
        SMILE = aux.code_symbols(SMILE)
        SMILEnew = "C%23N"
        self.assertEqual(SMILE, SMILEnew)

    def test_plus_replace1(self):
        SMILE = "C+N"
        SMILE = aux.code_symbols(SMILE)
        SMILEnew = "C%2BN"
        self.assertEqual(SMILE, SMILEnew)

    def test_minus_replace1(self):
        SMILE = "C-N"
        SMILE = aux.code_symbols(SMILE)
        SMILEnew = "C%2DN"
        self.assertEqual(SMILE, SMILEnew)

    def test_green(self):
        metric = 90.1
        flag = aux.metric_flag(metric)
        self.assertEqual(flag, "hazard-acceptable")

    def test_yellow_upper(self):
        metric = 90
        flag = aux.metric_flag(metric)
        self.assertEqual(flag, "hazard-warning")

    def test_yellow_lower(self):
        metric = 70
        flag = aux.metric_flag(metric)
        self.assertEqual(flag, "hazard-warning")

    def test_red(self):
        metric = 69.9
        flag = aux.metric_flag(metric)
        self.assertEqual(flag, "hazard-hazardous")


if __name__ == '__main__':
    unittest.main()

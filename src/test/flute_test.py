import unittest
from ..store.flute import StoreFlute

class TestStoreFlute(unittest.TestCase):
    testFlute = StoreFlute
    def testEvaluteWL(self):
        mockXcoordinates = [2, 0, 3, 1, 4, 6]
        mockYcoordinates = [0, 1, 2, 3, 5, 3]
        wl = self.testFlute.evaluteWL(mockXcoordinates, mockYcoordinates, len(mockXcoordinates), 3)
        self.assertEqual(wl, 13)

if __name__ == '__main__':
    unittest.main()

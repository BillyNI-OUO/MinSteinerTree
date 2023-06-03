from .. import repo
from ..repo.flute import Flute
from typing import List, Union
from ..store.logger import Logger

class storeFlute:
    def createMST(self, xCoordintaes:List[int], yCoordinates:List[int], degree:int, accuracy:int = 3):
        try:
            if len(xCoordintaes) != len(yCoordinates):    
                raise ValueError("The number of x_coordinate should be same as y_cordinate.")
        except ValueError:
            Logger.error("The number of x_coordinate should be same as y_cordinate.")
            return None
        try:
            if degree > len(xCoordintaes):
                raise ValueError("The value of degree shouldn't larger than the number of coordinates.")
        except ValueError:
            Logger.error("The value of degree shouldn't larger than the number of coordinates.")
            return None
        try:
            if degree < 1 or len(xCoordintaes) < 1:
                raise ValueError("The value of degree shouldn't less than 0.")
        except ValueError:
            Logger.error("The value of degree shouldn't less than 0.")
            return None
        try:
            tree = Flute.runFlute(xCoordintaes, yCoordinates, degree, accuracy)
        except Exception as e:
            Logger.error(str(e))
            return None
        try:
            if tree == None:
                Logger.error("StoreFlute run Flute failed.")
            res = Flute.treeToPairArray(tree)
        except Exception as e:
            Logger.error(str(e))
            return None
        
        return res

StoreFlute = storeFlute()
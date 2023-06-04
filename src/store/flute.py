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
    
    def evaluteWL(self, xCoordintaes:List[int], yCoordinates:List[int], degree:int, accuracy:int = 3):
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
            Logger.info("Start evalute MST wire length")
            res = Flute.runFluteWl(xCoordintaes, yCoordinates, degree, accuracy)
            Logger.info("Finsh evalute MST wire length.")
            if res != None:
                Logger.info(f"WireLength: {res}")
        except Exception as e:
            Logger.error(str(e))
            return None
        return res
    
    def createMSTAndWL(self, xCoordintaes:List[int], yCoordinates:List[int], degree:int, accuracy:int = 3):
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
            Logger.info("Start create MST")
            tree = Flute.runFlute(xCoordintaes, yCoordinates, degree, accuracy)
            Logger.info("Finish create MST")
        except Exception as e:
            Logger.error(str(e))
            return None
        try:
            if tree == None:
                Logger.error("StoreFlute run Flute failed.")
            Logger.info("Start transform tree to coordinates")
            res = Flute.treeToPairArray(tree)
            Logger.info("Finish transform tree to coordinates")
        except Exception as e:
            Logger.error(str(e))
            return None
        try:
            if tree == None:
                Logger.error("StoreFlute run Flute failed.")
            Logger.info("Start caculate wireLength from tree")
            resWL = Flute.caculateTreeWireLength(tree)
            Logger.info("Finsh caculate wireLength from tree")
            if resWL != None:
                Logger.info(f"Wirelength: {resWL}")
        except Exception as e:
            Logger.error(str(e))
            return None
        
        return res, resWL

StoreFlute = storeFlute()
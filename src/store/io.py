from .parser import StoreParser
from .logger import Logger
from .plot import StorePlot
import sys
from typing import List, Tuple

class storeIO:
    def __init__(self) -> None:
        pass
    def input(self) -> Tuple[List[int], List[int]]:
        if StoreParser.args.input == None:
            Logger.info("Parsing data from stdin.")
            res = self.parseFromStdin()
            Logger.info("End parsing data from stdin.")
            if res == None:
                Logger.error("Parsing fail from stdin.")
            return res
        else:
            Logger.info(f"Parsing data from {StoreParser.args.input}.")
            res = self.parseFromFile(StoreParser.args.input)
            Logger.info(f"End parsing data from {StoreParser.args.input}.")
            if res == None:
                Logger.error(f"Parsing fail from {StoreParser.args.input}.")
            return res 
    
    def parseFromStdin(self) -> Tuple[List[int], List[int]]:
        Logger.info("Input x y coordinate seperate by space and line, end with EOF")
        res = ([], [])
        inputLines = sys.stdin.readlines()
        for line in inputLines:
            coordinates = line.split(" ")
            if(len(coordinates) != 2):
                Logger.error("Input form invalid.")
                return None
            try:
                res[0].append(int(coordinates[0]))
                res[1].append(int(coordinates[1]))
            except Exception as e:
                Logger.error(str(e))
                return None
        return res
    
    def parseFromFile(self, fileName:str)-> Tuple[List[int], List[int]]:
        Logger.info("Input x y coordinate seperate by space and line, end with EOF")
        res = ([], [])
        with open(fileName) as fp:
            inputLines = fp.readlines()
            for line in inputLines:
                coordinates = line.split(" ")
                if(len(coordinates) != 2):
                    Logger.error("Input form invalid.")
                    return None
                try:
                    res[0].append(int(coordinates[0]))
                    res[1].append(int(coordinates[1]))
                except Exception as e:
                    Logger.error(str(e))
                    return None
            return res
    
    def output(self)->None:
        if StoreParser.args.output != None:
            try:
                StorePlot.savePicture(StoreParser.args.output)
            except Exception as e:
                Logger.error(str(e))
        else:
            Logger.info("End with no output")
        if StoreParser.args.visible:
            try:
                StorePlot.show()
            except Exception as e:
                Logger.error(str(e))
        else:
            Logger.info("End with no show")
            


IO = storeIO()
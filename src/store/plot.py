import matplotlib.pyplot as plt
from ..constants import Constants
from .logger import Logger
from typing import List, Tuple
class storePlot:
    def __init__(self) -> None:
        pass

    def setup(self, xSize:float=Constants.PLOT_X_SIZE, ySize:float=Constants.PLOT_Y_SIZE) -> None:
        plt.figure(figsize=(xSize, ySize))


    def drawMST(self, xCoordinates:List[int], yCoordinate:List[int], edges:Tuple[List[int], List[int]]) -> None:
        plt.scatter(xCoordinates, yCoordinate, s=10, color='r')
        for i in range(int(len(edges[0])/2)):
            plt.plot(edges[0][2*i:2*i+2], edges[1][2*i:2*i+2], color='k')    

        plt.xlim(float(min(xCoordinates)-2), float(max(xCoordinates)+2))
        plt.ylim(float(min(yCoordinate)-2), float(max(yCoordinate)+2))
        plt.xlabel(r'$X$', size=16)
        plt.ylabel(r'$Y$', size=16)
        plt.grid(b=True, which='major', color = 'r', linestyle= '--', linewidth = '0.5')
        plt.tight_layout()
    
    def savePicture(self, fileName:str)->None:
        Logger.info(f"Output to {fileName}")
        plt.savefig(fileName)

    def show(self)->None:
        Logger.info("Showing the graph")
        plt.show()
        Logger.info("Terminate successful")

StorePlot=storePlot()
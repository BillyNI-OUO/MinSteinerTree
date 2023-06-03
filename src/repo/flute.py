from . import _flute
from ._flute import Tree
from ..store.logger import Logger
from typing import List, Union

class flute(object):
	
	def __init__(self):
		_flute.readLUT()

	def runFlute(self, x_coordintaes:List[int], y_coordintates:List[int], degree:int, accuracy:int = 3)->Tree:
		try:
			res = _flute.flute(degree, x_coordintaes, y_coordintates, accuracy)
		except Exception as e:
			Logger.error("ReporunFlute failed.\nx_coordintaes: {x_coordintaes}\ny_coordintates: {y_coordintates}")
			Logger.error(str(e))	
			return None
		return res
	
	def treeToPairArray(self, tree: Tree)->Union[List[int], List[int]]: 
		try:
			res = _flute.treeToPairArray(tree)
		except Exception as e:
			Logger.error("RepotreeToPairArray failed.")
			Logger.error(str(e))	
			return None
		return res.x, res.y

Flute=flute()
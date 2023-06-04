import src
from src import store
from src.store import StoreFlute, StorePlot
from src.store import Logger, IO, StoreParser

x, y = IO.input()
res, wl = StoreFlute.createMSTAndWL(x, y, len(x), 3)
StorePlot.drawMST(x, y, res)
IO.output()
#res = StoreFlute.createMST(x, y, len(x), 3)
#StorePlot.drawMST(x, y, res)
#IO.output()


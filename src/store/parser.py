import argparse
from typing import Any
class storeParser:
    parser: argparse.ArgumentParser
    def __init__(self) -> None:
        self.setup()

    def setup(self) -> None:
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("-o", "--output", help="the output png filename, default = None", default=None)
        self.parser.add_argument("-i", "--input", help="input file name, default = None", default=None)
        self.parser.add_argument("-v", "--visible", help="visible", default=True, action="store_true")
    @property
    def args(self)->Any:
        return self.parser.parse_args()


StoreParser = storeParser()

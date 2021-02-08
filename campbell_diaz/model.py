import sys, os
from pcse.engine import Engine

this_dir = os.path.dirname(__file__)


class CampbellDiazModel(Engine):
    """Implementation of the Campbell-Diaz water productivity model

    """
    config = os.path.join(this_dir, "campbell.conf")

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)
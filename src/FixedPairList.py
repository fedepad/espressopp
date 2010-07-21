from espresso import pmi
import _espresso 
import espresso
from espresso.esutil import cxxinit

class FixedPairListLocal(_espresso.FixedPairList):
    'The (local) fixed pair list.'

    def __init__(self, storage):
        'Local construction of a fixed pair list'
        cxxinit(self, _espresso.FixedPairList, storage)

    def add(self, pid1, pid2):
        'add pair to fixed pair list'
        return self.cxxclass.add(self, pid1, pid2)

if pmi.isController:
    class FixedPairList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.FixedPairListLocal',
            pmicall = [ "add" ]
            )
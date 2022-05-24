import os
import sys
from .data_update import update
from .iri import IRI, MissingDataError

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

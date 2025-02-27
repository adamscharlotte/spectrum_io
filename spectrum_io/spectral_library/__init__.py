"""Initialize spectral library."""
import logging

from . import digest
from .dlib import DLib
from .msp import MSP
from .spectronaut import Spectronaut

logger = logging.getLogger(__name__)

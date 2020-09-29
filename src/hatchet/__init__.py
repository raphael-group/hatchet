__version__ = '0.1.1'

import os.path
from importlib_resources import path
import hatchet
from hatchet.utils.config import Config


with path(hatchet, 'hatchet.ini') as ini_file:
    filenames = [ini_file]
    # If a hatchet.ini file exists where we were imported from, use it after the pre-packaged .ini file
    # This allows overriding of values without having to alter the pre-packaged .ini file, which may
    # be buried deep inside the site-packages folder.s
    if os.path.exists('hatchet.ini'):
        filenames.append('hatchet.ini')
    config = Config('hatchet', filenames)

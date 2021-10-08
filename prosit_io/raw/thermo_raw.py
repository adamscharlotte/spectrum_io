from typing import Optional
from .msraw import MSRaw
import os
import pathlib
from sys import platform

import logging
logger = logging.getLogger(__name__)

class ThermoRaw(MSRaw):

    @staticmethod
    def convert_raw_mzml(input_path: str, output_path: Optional[str] = None, gzip = False):
        """
        Converts a ThermoRaw file to mzML

        :param input_path: File path of the Thermo Rawfile
        :param output_path: File path of the mzML path
        """
        if output_path is None:
            output_path = os.path.splitext(input_path)[0] + ".mzML"
        
        if os.path.isfile(output_path):
            logger.info(f"Found converted file at {output_path}, skipping conversion")
            return output_path
        
        if gzip:
            gzip = "-g"
        else:
            gzip = ""

        if "linux" in platform:
            mono = "mono"
        elif "win" in platform:
            mono = ""
        
        exec_path = pathlib.Path(__file__).parent.absolute() # get path of parent directory of current file
        exec_command = f"{mono} {exec_path}/utils/ThermoRawFileParser/ThermoRawFileParser.exe {gzip} -i {input_path} -b {output_path}.tmp"
        logger.info(f"Converting thermo rawfile to mzml with the command: '{exec_command}'")
        os.system(exec_command)
        
        # only rename the file now, so that we don't have a partially converted file if something fails
        os.rename(f"{output_path}.tmp", output_path)
        
        return output_path

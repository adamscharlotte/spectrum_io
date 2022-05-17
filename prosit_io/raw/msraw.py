import logging
import os
import warnings
from abc import abstractmethod
from typing import Union, List, Optional, Dict

import pandas as pd

from fundamentals.constants import MZML_DATA_COLUMNS

logger = logging.getLogger(__name__)


class MSRaw:
    path: str
    output_path: str

    def __init__(self, path: Optional[str] = None, output_path: Optional[str] = None):
        self.path = path
        self.output_path = output_path

    @abstractmethod
    def convert_raw_mzml(self, input_path, output_path):
        """
        Use https://github.com/compomics/ThermoRawFileParser for conversion
        """
        raise NotImplementedError

    @staticmethod
    def read_mzml(
        source: Union[str, List[str]],
        ext: str = 'mzml',
        package: str = 'pyteomics',
        scanidx: Optional[List] = None,
        *args,
        **kwargs
    ):
        """
        Reads mzml and generates a dataframe containing intensities and m/z values.

        :param source: A directory containing mzml files, a list of files or a single file.
        :param ext: File extension for searching a specified directory.
        :param package: Package for parsing the mzml file. Can eiter be "pymzml" or "pyteomics"

        :return: Pandas DataFrame
        """

        if isinstance(source, str):
            file_list = []
            if os.path.isdir(source):
                # if string is provided and is a directory, search all mzml files with provided extension
                for file in os.listdir(source):
                    if file.lower().endswith(ext.lower()):
                        file_list.append(file)

            else:
                file_list = [source]
            source = file_list
        data = {}
        print(package)
        if package == 'pymzml':
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=ImportWarning)
                for file_path in source:
                    logger.info(f"Reading mzML file: {file_path}")
                    MSRaw._get_scans_pymzml(file_path, data, scanidx, *args, **kwargs)
        elif package == 'pyteomics':
            from pyteomics import mzml
            print(source)
            for file_path in source[:1]:
                logger.info(f"Reading mzML file: {file_path}")
                data_iter = mzml.read(source=file_path, *args, **kwargs)
                file_name = os.path.splitext(os.path.basename(file_path))[0]
                for spec in data_iter:
                    id = spec['id'].split('scan=')[-1]
                    key = f"{file_name}_{id}"
                    data[key] = [file_name, id, spec['intensity array'], spec['m/z array'], int(spec['precursorList']['precursor'][0]['activation']['collision energy'])]
                data_iter.close()
        else:
            assert False, "Choose either 'pymzml' or 'pyteomics'"

        data = pd.DataFrame.from_dict(data, orient='index', columns=MZML_DATA_COLUMNS)
        data['SCAN_NUMBER'] = pd.to_numeric(data['SCAN_NUMBER'])
        return data
    
    @staticmethod
    def _get_scans_pymzml(
        file_path: str,
        data: Dict,
        scanidx: Optional[List] = None,
        *args,
        **kwargs
    ) -> None:
        """
        Reads mzml and generates a dataframe containing intensities and m/z values.

        :param file_path: path to a single mzml file.
        :param data: dictionary to be added to by this function
        :param scanidx: optional list of scan numbers to extract. if not specified, all scans will be extracted
        """
        import pymzml
        data_iter = pymzml.run.Reader(file_path, args=args, kwargs=kwargs)
        file_name = os.path.splitext(os.path.basename(file_path))[0]
        if scanidx is None:
            for spec in data_iter:
                key = f"{file_name}_{spec.ID}"
                data[key] = [file_name, spec.ID, spec.i, spec.mz]
        else:
            for idx in scanidx:
                spec = data_iter[idx]
                key = f"{file_name}_{spec.ID}"
                data[key] = [file_name, spec.ID, spec.i, spec.mz]
        data_iter.close()


if __name__ == "__main__":
    from sys import argv
    if len(argv) == 3:
        data = {}
        scanidx = list(map(int, argv[2].split(",")))
        MSRaw()._get_scans_pymzml(argv[1], data, scanidx)
        for _, spec_id, intensities, mzs in data.values():
            print(spec_id)
            for m, i in zip(mzs, intensities):
                print(m, i, sep='\t')
    else:
        print("Please specify a mzml file and a comma separated list of scan numbers to extract")

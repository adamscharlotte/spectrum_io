import logging
import os

import sqlite3
import pandas as pd
from timspy.df import TimsPyDF
from mgf_filter.masterSpectrum import MasterSpectrum

logger = logging.getLogger(__name__)

def binning(inp, ignoreCharges):
    ms = MasterSpectrum()
    ms.load_from_tims(inp, ignoreCharges)
    ms.export_to_csv(temp_path)
    comb_ms = pd.read_csv(temp_path)
    scan = inp["Scannumber"]
    comb_ms["Scannumber"] = scan
    comb_ms = comb_ms.drop(columns=["counts", "left border", "right border", "start_mz", "ms1_charge", "rel_intensity_ratio", "counts_ratio"])
    return comb_ms

class TimsD:
    """Main to read .d folder and generate dataframe containing intensities and m/z values."""

    d_path: Optional[str]
    output_path: Optional[str]
    maxquant_path: Optional[str]

    def __init__(self, d_path: Optional[str] = None, output_path: Optional[str] = None, maxquant_path: Optional[str] = None):
        """
        :param path: path to .d folder
        :param output_path: path to save the output file
        """
        self.d_path = d_path
        self.output_path = output_path
        self.maxquant_path = maxquant_path

    def extract_pasef(self):
        tdf_path = self.d_path + "/analysis.tdf"
        con = sqlite3.connect(tdf_path)
        cur = con.cursor()
        pasef_df = pd.read_sql_query("SELECT * from PasefFrameMsMsInfo", con)
        con.close()

    def extract_d_psm(self, pasef_df):
        try:
            import opentims_bruker_bridge
            all_columns = ('frame','scan','tof','intensity','mz','inv_ion_mobility','retention_time')
        except ModuleNotFoundError:
            print("Without Bruker proprietary code we cannot yet perform tof-mz and scan-dt transformations.")
            print("Download 'opentims_bruker_bridge' if you are on Linux or Windows.")
            print("Otherwise, you will be able to use only these columns:")
            all_columns = ('frame','scan','tof','intensity','retention_time')
            self.d_path
            path_precursors = self.maxquant_path + "accumulatedMsmsScans.txt"
            df_precursors = pd.read_csv(path_precursors, sep = "\t")
            df_precursors.columns = df_precursors.columns.str.replace(" ", "")
            df_precursors_filtered = df_precursors[df_precursors.Scannumber.isin(SCAN_NUMBERS)]
            df_precursors_filtered.PASEFprecursorIDs = df_precursors_filtered.PASEFprecursorIDs.str.split(";").apply(lambda s: [int(x) for x in s])
            df_precursors_filtered = df_precursors_filtered.explode("PASEFprecursorIDs")
            path_pasef = self.maxquant_path + "pasefMsmsScans.txt"
            df_pasef = pd.read_csv(path_pasef, sep = "\t")
            df_pasef.columns = df_pasef.columns.str.replace(" ", "")
            df_pasef_filtered = df_pasef[df_pasef.Precursor.isin(df_precursors_filtered["PASEFprecursorIDs"])]
            frames = df_pasef_filtered.Frame.to_numpy()
            D = TimsPyDF(self.d_path)
            df_tims = D.query(frames=frames, columns=all_columns)
            df_raw_pasef_all_scans = df_tims.merge(df_pasef_filtered, left_on='frame', right_on='Frame')
            df_raw_pasef = df_raw_pasef_all_scans[df_raw_pasef_all_scans.ScanNumBegin <= df_raw_pasef_all_scans.scan]
            df_raw_pasef = df_raw_pasef[df_raw_pasef.scan <= df_raw_pasef.ScanNumEnd]
            df_raw_pasef_mapped = df_raw_pasef.drop(columns=["scan"]).merge(df_precursors_filtered, left_on='Precursor', right_on='PASEFprecursorIDs')
            df_intensity = df_raw_pasef_mapped.groupby(['Scannumber'])['intensity'].apply(list).reset_index(name='combined_INTENSITIES')
            df_tof = df_raw_pasef_mapped.groupby(['Scannumber'])['tof'].apply(list).reset_index(name='combined_MZ').rename(columns = {"Charge": "CHARGE"}, inplace=True)
            df_ms_values = df_raw_pasef_mapped[['Scannumber', 'Sequence', 'Modifiedsequence', 'CHARGE']].drop_duplicates().merge(df_tof, left_on=['Scannumber'], right_on=['Scannumber']).merge(df_intensity, left_on=['Scannumber'], right_on=['Scannumber'])
            # Combine Scans
            bin_result_df = pd.DataFrame()
            for index, line in df_ms_values.iterrows():
                bin_result = binning(line, True)
                bin_result_df = bin_result_df.append(bin_result)
            bin_result_df_collapsed = bin_result_df.groupby("Scannumber").agg(list)
            un_annot_df_combined = pd.merge(df_ms_values, bin_result_df_collapsed, on="Scannumber")




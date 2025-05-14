import click
import csv
import pandas as pd


@click.command()
@click.option('--sicstats', type=click.Path(exists=True), required=True, help='Path to SICstats file from MASIC.')
@click.option('--syn', type=click.Path(exists=True), required=True, help='Path to syn file from MS-GF+.')
@click.option('--output', type=click.Path(writable=True), required=True, help='Path to write the combined output CSV.')
def merge_files(sicstats, syn, output):
    """
    Merges SICstats and syn files side-by-side and writes the combined DataFrame to a TSV output file.
    """

    df_sic = pd.read_csv(sicstats, sep='\t')
    df_syn = pd.read_csv(syn, sep='\t')

    # Syn file columns of interest
    columns_syn = [
        'ResultID', 
        'Scan',
        'FragMethod',
        'SpecIndex',
        'Charge',
        'PrecursorMZ',
        'DelM',
        'DelM_PPM',
        'MH',
        'Peptide',
        'Protein',
        'NTT',
        'DeNovoScore',
        'MSGFScore',
        'MSGFDB_SpecEValue',
        'Rank_MSGFDB_SpecEValue',
        'EValue',
        'QValue',
        'PepQValue',
        'IsotopeError'
    ]

    # SICstats file columns of interest
    columns_sic = [
        'OptimalPeakApexScanNumber',
        'PeakMaxIntensity',
        'PeakSignalToNoiseRatio',
        'FWHMInScans',
        'PeakArea',
        'ParentIonIntensity',
        'MZ',
        'StatMomentsArea',
        'PeakScanStart',
        'PeakScanEnd',
        'FragScanNumber'
    ]

    df_sic = df_sic[columns_sic]
    df_syn  = df_syn[columns_syn]
    
    # Add empty columns to df_syn to preserve original column order
    df_syn['ElutionTime'] = ''
    df_syn['ScanType'] = ''
    df_syn['TotalIonIntensity'] = ''
    df_syn['BasePeakIntensity'] = ''
    df_syn['BasePeakMZ'] = ''

    df_sic = df_sic.rename(columns={'MZ': 'ParentIonMZ', 'OptimalPeakApexScanNumber': 'Optimal_Scan_Number'})

    merged_df = pd.merge(
        df_syn,
        df_sic,
        left_on='Scan',
        right_on='FragScanNumber',
        how='inner'
    )

    # Scan and FragScanNumber are the same, so we can drop one of them
    merged_df.drop(columns='FragScanNumber', inplace=True)

    # Added here to preserve the original column order
    merged_df['PeakWidthMinutes'] = 0

    merged_df = merged_df.sort_values(by='ResultID')

    merged_df.to_csv(output, index=False, sep='\t')

if __name__ == '__main__':
    merge_files()


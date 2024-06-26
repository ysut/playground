# Standard modules
from dataclasses import dataclass
from logging import getLogger, config
import os

# Third party modules
import gffutils
import pybedtools
import pandas as pd
import numpy as np
import yaml
from pathlib2 import Path
from tqdm import tqdm
from typing import NamedTuple

# Local modules
from .libs.args import parser_setting
from .libs.utils import load_config, dfs_to_excel, OutputSettings
from .libs.preprocess import PreProcessExomeSummary
from .libs.modesamples import ModeSamples
from .libs.annolib.anno import Anno
from .libs.annolib.genebased import GeneBasedAnno
from .libs.filter.maffilter import MafFilter
from .libs.filter.typefilter import TypeFilter
from .libs.filter.gtfilter import GtFilter
from .libs.filter.qcfilter import QcFilter
from .libs.filter.hardfilter import HardFilter
from .libs.filter.counter import counter

# Settings
tqdm.pandas()


def main():
    #----- STEP 0. Logging settings
    parent_directory = os.path.dirname(os.path.dirname(__file__))
    config_path: str = os.path.join(parent_directory, 'config/logging.yaml')
    with open(config_path, 'r') as f:
        config.dictConfig(yaml.safe_load(f))
    logger = getLogger(__name__)

    #----- STEP 1. Argument settings
    logger.info('STEP 1. Argument settings')
    args = parser_setting()

    #----- STEP 2. Load exome_summary file and config file
    logger.info('STEP 2. Load exome_summary file and config file')
    df: pd.DataFrame = pd.read_table(args['input'], header=0, dtype=str)
    configs: dict = load_config(args['config'])

    #----- STEP 3. Get Mode and Samples information
    logger.info('STEP 3. Get Mode and Samples information')
    modesamples: NamedTuple = ModeSamples(df=df, args=args)
    mode_samples_info = modesamples.get_mode_samples_info()
    logger.info(f'Analysis mode: {mode_samples_info.mode}')
    logger.info(f'Analysis samples: {modesamples.analysis_samples}')

    #----- STEP 4. Output settings
    logger.info('STEP 4. Output settings')
    output_settings = OutputSettings(
        args=args, mode_samples_info=mode_samples_info
        )
    output_file_path: str = output_settings.get_saving_file_path()

    #----- STEP 5. Pre-processing
    logger.info('STEP 5. Pre-processing')
    preprocessing = PreProcessExomeSummary(
        df=df, args=args, mode_samples_info=mode_samples_info
        )
    df = preprocessing.all_pre_processing()

    #-----   STEP 6. Annotation
    logger.info('STEP 6. Annotation')
    logger.info('STEP 6-1. Gene-based annotation')
    gba = GeneBasedAnno(args)
    df = gba.anno_hgmd(df=df)
    df = gba.dm_filter(df=df)
    df = gba.anno_dcpr(df=df)
    df['MOI'] = df.apply(gba.summarize_moi, axis=1)

    # logger.info('STEP 6-2. Variant-based annotation')

    #-----   STEP 7. Filtering
    qcfilter = QcFilter(configs=configs)
    df = qcfilter.exclude_low_quality(df)

    maffilter = MafFilter(
        df=df, mode_samples_info=mode_samples_info, configs=configs)
    df = maffilter.all_filtering()

    typefilter = TypeFilter(df=df)
    df = typefilter.exclude_hlamuc_and_exonicsyno()

    gtfilter = GtFilter(
        df=df, mode_samples_info=mode_samples_info)
    dfs = gtfilter.genotypeing_filter()


    #-----   STEP 8. Count variants of filtering process
    logger.info('STEP 8. Count variants of filtering process')
    countsummery_file = str(Path(output_file_path).parent) + '/CountSummary.xlsx'
    filtered_dfs = counter(dfs=dfs, output_excel=countsummery_file)
    

    #----  STEP 10. Hard filtering
    logger.info('STEP 10. Hard filtering')
    hard = HardFilter(dfs=filtered_dfs)
    filtered_dfs = hard.hard_filtering()


    #-----   STEP 9. Output as an Excel file
    dfs_to_excel(filtered_dfs, f"{output_file_path}.xlsx")
    logger.info(f"Output file: {output_file_path}.xlsx has been saved.")

    

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

# Local modules
from .libs.args import parser_setting
from .libs.utils import load_config, OutputSettings
from .libs.preprocess import PreProcessExomeSummary
from .libs.modesamples import ModeSamples
from .libs.annolibs.anno import Anno
from .libs.filter.maffilter import MafFilter
from .libs.filter.typefilter import TypeFilter
from .libs.filter.gtfilter import GtFilter
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
    modesamples = ModeSamples(df=df, args=args)
    mode_samples_info = modesamples.get_mode_samples_info()
    logger.info(f'Analyze mode: {mode_samples_info["mode"]}')
    logger.info(f'Samples: {mode_samples_info["samples"]}')

    #----- STEP 4. Output settings
    output_settings = OutputSettings(
        args=args, mode_samples_info=mode_samples_info
        )
    output_file_path: str = output_settings.get_saving_file_path()

    #----- STEP 5. Pre-processing
    preprocessing = PreProcessExomeSummary(
        df=df, mode_samples_info=mode_samples_info
        )
    df = preprocessing.all_pre_processing()




    

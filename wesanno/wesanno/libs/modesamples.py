import sys
import pandas as pd
from typing import NamedTuple, Optional


class ModeSamplesInfo(NamedTuple):
    mode: Optional[str] = None
    proband_id: Optional[str] = None
    father_id: Optional[str] = None
    mother_id: Optional[str] = None
    sibling_id: Optional[str] = None
    parent_id: Optional[str] = None


class ModeSamples:
    def __init__(self, df: pd.DataFrame, args: dict):
        self.df: pd.DataFrame = df
        self.input_mode: str = args['mode']
        self.input_samples: tuple = args['samples']

    def __varidate_mode(self) -> str:
        format_index: int = self.df.columns.get_loc('FORMAT')
        self.samples: list = self.df.columns[format_index+1:].tolist()
        if len(self.samples) == 1:
            return 'proband'
        elif len(self.samples) == 2:
            return 'duo'
        elif len(self.samples) == 3:
            return 'trio'
        elif len(self.samples) == 4:
            # When 'auto' mode, return 'quad_affected'
            return 'quad_affected'
        elif len(self.samples) == 0:
            print("Error: This file does not have sample columns")
            sys.exit(1)
        else:
            # Have not been implemented.
            return 'quad_affected'

    def get_mode_samples_info(self) -> ModeSamplesInfo:
        # Determine analysis mode
        if self.input_mode == 'auto':
            mode = self.__varidate_mode()
        else:
            mode = self.input_mode
            
        # Get Sample info
        if self.input_samples == 'auto':
            format_index: int = self.df.columns.get_loc('FORMAT')
            samples: list = self.df.columns[format_index+1:].tolist()
            sample_generator = (sample for sample in samples)
            if mode == 'proband':
                return ModeSamplesInfo(mode=mode,
                                   proband_id=next(sample_generator))
            elif mode == 'duo':
                return ModeSamplesInfo(mode=mode,
                                   proband_id=next(sample_generator),
                                   parent_id=next(sample_generator))
            elif mode == 'trio':
                return ModeSamplesInfo(mode=mode,
                                   proband_id=next(sample_generator),
                                   father_id=next(sample_generator),
                                   mother_id=next(sample_generator))
            elif ((mode == 'quad_affected') 
                  | (mode == 'quad_unaffected')):
                return ModeSamplesInfo(mode=mode,
                                   proband_id = next(sample_generator),
                                   father_id = next(sample_generator),
                                   mother_id = next(sample_generator),
                                   sibling_id = next(sample_generator))

        elif ((len(self.input_samples) == 1) 
              | (len(self.input_samples) == 2)
              | (len(self.input_samples) == 3)
              | (len(self.input_samples) == 4)):
            # mode = self.__varidate_mode()
            if mode == 'proband':
                return ModeSamplesInfo(mode=mode, 
                                   proband_id=self.input_samples[0]
                                   )
            elif mode == 'duo':
                return ModeSamplesInfo(mode=mode,
                                   proband_id=self.input_samples[0],
                                   parent_id=self.input_samples[1]
                                   )
            elif mode == 'trio':
                return ModeSamplesInfo(mode=mode,
                                proband_id=self.input_samples[0],
                                father_id=self.input_samples[1],
                                mother_id=self.input_samples[2]
                                )
            elif ((mode == 'quad_affected') 
                  | (mode == 'quad_unaffected')):
                    return ModeSamplesInfo(mode=self.input_mode,
                                       proband_id=self.input_samples[0],
                                       father_id=self.input_samples[1],
                                       mother_id=self.input_samples[2],
                                       sibling_id=self.input_samples[3]
                                       )
            
        else:
            print('Error: It has not been implemented ' 
                  'in more than 5 sample cases.')
            sys.exit()

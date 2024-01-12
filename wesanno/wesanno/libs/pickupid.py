import sys

class ModeID:
    def __init__(self, dataframe, mode):
        print(f"\nInput analysis mode     : {mode}")
        self.columns = dataframe.columns
        self.dnaIDs = [column for column in self.columns if column.startswith('Sample_')]
        if mode == 'auto':
            num_samples = len(self.dnaIDs)
            if num_samples == 1:
                self.analysis_mode = 'proband'
            elif num_samples == 2:
                self.analysis_mode = 'duo'
            elif num_samples == 3:
                self.analysis_mode = 'trio'
            elif num_samples == 4:
                self.analysis_mode = 'quad'
            elif num_samples == 0:
                print('Sample infomation is not included in this file.')
                sys.exit()
            else:
                print('Too many sammple IDs.')
                sys.exit()            
        else:
            self.analysis_mode = mode
        print(f"Determined analysis mode: {self.analysis_mode}\n")


    def __dnaIdGenerator(self):
        for dnaID in self.dnaIDs:
            yield dnaID


    def pickup_dnaID(self):
        self.id = self.__dnaIdGenerator()
        if self.analysis_mode == 'proband':
            proband_id = next(self.id)
            return proband_id, # return tuple!!

        elif self.analysis_mode == 'trio':
            proband_id = next(self.id)
            father_id = next(self.id)
            mother_id = next(self.id)
            return proband_id, father_id, mother_id
        
        elif self.analysis_mode == 'quad':
            proband_id = next(self.id)
            father_id = next(self.id)
            mother_id = next(self.id)
            sibling_id = next(self.id)
            return proband_id, father_id, mother_id, sibling_id

        else:
            proband_id = next(self.id)
            parent_id = next(self.id)
            return proband_id, parent_id

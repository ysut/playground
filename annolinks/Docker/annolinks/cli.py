import os
import sys
from logging import getLogger, config

import yaml

sys.path.append(os.path.join(os.path.dirname(__file__), 'libs'))
import hyperlink
from args import parser_setting, analyze_args

parent_directory = os.path.dirname(os.path.dirname(__file__))
config_path = os.path.join(parent_directory, 'config/logging_config.yaml')


def main():
    #------- logging setting -------#
    with open(config_path, 'r') as f:
        config.dictConfig(yaml.safe_load(f))
        
    logger = getLogger(__name__)

    #------- args setting -------#
    raw_args = parser_setting()
    args = analyze_args(raw_args)
    
    #------- Step 1. Insert hyperlink columns to each dataframe -------#
    logger.info('Insert hyperlink columns to each dataframe')
    
    #1-1. Load excel file as df(dataclass) and get sheet names
    logger.info(f'Load {args["input"]}')
    dfs = hyperlink.load_excel_as_dataclass(args['input'])
    
    #1-2. Generate a sheet names list for annotation
    anno_sheets: list = hyperlink.generate_anno_sheets_list(
        input_sheets=list(dfs.sheets.keys()), 
        skip_sheets=args['skip_sheets']
        )
    
    #1-3. Insert hyperlink columns to each dataframe
    hl = hyperlink.Hyperlink(
        gene_symbol_col=args['gene_col'], 
        alt_col=args['alt_col'],
        franklin_page=args['franklin_page'],
        assembly=args['assembly'],
        skip_sites=args['skip_sites'],
        pos19=args['pos19'],
        pos38=args['pos38'],
        ucsc_width=args['ucsc_width'],
        spliceai_maskFlag=args['spliceai_maskFlag'],
        spliceai_dist=args['spliceai_dist'],
        windowsFlag=args['windowsFlag']
        )
    
    for sheet in anno_sheets:
        logger.info(f'Insert hyperlink columns to {sheet}')
        dfs.sheets[sheet] = hl.insert_hyperlinks(dfs.sheets[sheet])
    
    #------- Step 2. Write to temporary excel file -------#
    tmp_file_path: str = args['input'].stem + '_tmp.xlsx'
    logger.info(f'Write to {tmp_file_path} for temporary use')
    hyperlink.df_to_excel(dfs, tmp_file_path)

    #------- Step 3. Convert URL to Hyperlink using Openpyxl -------#
    logger.info('Convert URL to Hyperlink using Openpyxl')
    wb = hyperlink.convert_url_to_hyperlink(
        input_excel=tmp_file_path, 
        anno_sheets=anno_sheets, 
        skip_sites=args['skip_sites'])
    
    #------- Step 4. Save the result as a new excel file -------#
    logger.info('Save the result as a new excel file')
    wb.save(args['output'])
    wb.close()
    os.remove(tmp_file_path)






    

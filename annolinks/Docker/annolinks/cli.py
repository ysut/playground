import os
import sys
from logging import getLogger, config

import yaml

sys.path.append(os.path.join(os.path.dirname(__file__), 'libs'))
import preprocess
import errcheck
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
    
    #------- Step 1. Load an input file and Check the file format -------#
    #1-1. Load an excel file and get sheet names
    logger.info(f"Load {args['input']}")
    dfs = hyperlink.load_excel_as_dataclass(args['input'])
    
    #1-2. Generate a sheet names list for annotation
    anno_sheets: list = hyperlink.generate_anno_sheets_list(
        input_sheets=list(dfs.sheets.keys()), 
        skip_sheets=args['skip_sheets']
        )
    
    #1-3. Error check
    ec = errcheck.ErrorCheck(dfs=dfs, anno_sheets=anno_sheets, 
                             skip_sites=args['skip_sites'], args=args)
    
    logger.info('Check specified columns (HGMD, DECIPHER, ...)')    
    if ec.has_specified_columns():
        sys.exit(1)
    
    logger.info("Check required columns "
                "(CHROM, POS, REF, and ALT (or sepcified ALT column name))")
    if not ec.has_required_columns():
        sys.exit(1)
    
    logger.info('Check hyperlink limitation')
    if not ec.is_within_limit():
        sys.exit(1)

    #------- Step 2. Preprocess -------#
    #2-1. Insert liftover columns to each dataframe
    logger.info('Preprocess (liftover)')
    dfs = preprocess.liftover_to_hg38(
        dfs=dfs, args=args, anno_sheets=anno_sheets)
    
    #2-2. Split ALT column
    logger.info('Preprocess (split ALT)')
    dfs = preprocess.split_alt_col(
        dfs=dfs, args=args, anno_sheets=anno_sheets)


    #------- Step 3. Insert URLs -------#
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
    
    # Insert URLs to each dataframe
    for sheet in anno_sheets:
        logger.info(f'Insert URL columns to {sheet}')
        dfs.sheets[sheet] = hl.insert_urls(dfs.sheets[sheet])
    
    #------- Step 4. Write to temporary excel file -------#
    tmp_file_path: str = args['input'].stem + '_tmp.xlsx'
    logger.info(f'Write to {tmp_file_path} for temporary use')
    hyperlink.df_to_excel(dfs, tmp_file_path)

    #------- Step 5. Convert URL to Hyperlink using Openpyxl -------#
    logger.info('Convert URLs to Hyperlinks using Openpyxl')
    wb = hyperlink.convert_url_to_hyperlink(
        input_excel=tmp_file_path, 
        anno_sheets=anno_sheets, 
        skip_sites=args['skip_sites'])
    
    #------- Step 6. Save the result as a new excel file -------#
    logger.info('Saving the result as a new excel file......')
    wb.save(args['output'])
    wb.close()
    os.remove(tmp_file_path)
    logger.info(f"Saved as {args['output']}")
    logger.info('Done')





    

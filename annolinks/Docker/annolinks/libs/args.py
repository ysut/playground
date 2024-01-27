import argparse
import pathlib2
from logging import getLogger

logger = getLogger(__name__)

def parser_setting() -> dict:
    parser = argparse.ArgumentParser(description = 'WES')
    
    #1 Input
    parser.add_argument('--input', '-I', required=True, type=pathlib2.Path, 
                        help='path to a excel file')
    
    #2 Gene symbol column name
    parser.add_argument('--gene-col', '-G', required=True, type=str, 
                        help='The column name of gene symbol')
    
    #3 Using ALT column name
    parser.add_argument('--alt-col', '-L', required=False, type=str,
                        default='ALT', help='The column name of ALT')
    
    #4 Output
    parser.add_argument('--output', '-O', required=False, type=pathlib2.Path, 
                        help='path to the output file')

    #5 Assembly
    parser.add_argument('--assembly', '-A', required=False, type=str, 
                        choices=['hg19', 'GRCh37', 'hg38', 'GRCh38'], 
                        default='hg19', 
                        help='Genome assembly')  

    #6 Skip excel sheets
    parser.add_argument('--skip-sheets', required=False, type=str,
                        help='Skip annotation in selected excel sheets')

    #7 Skip sites
    parser.add_argument('--skip-sites', required=False, type=str,
                        help=('Skip annotation to selected site links.'
                              'Select from "HGMD", "Franklin", "DECIPHER", '
                              '"SpliceAI", and, "UCSC".'
                              ' (e.g. Franklin,UCSC,HGMD etc.)'
                              ))
    
    #8 USCS genome browser display width
    parser.add_argument('--ucsc-width', required=False, type=int, 
                        help='The width of displaying the UCSC Genome Browser')
    
    #9 POS column name for hg19
    parser.add_argument('--pos19', required=False, type=str, 
                        help='The column name of POS of hg19')
    
    #10 POS column name for hg38
    parser.add_argument('--pos38', required=False, type=str, 
                        help='The column name of POS of hg38')
    
    #11 Franklin page to access
    parser.add_argument('--franklin-page', required=False, type=str,
                        choices=[
                              'assessment-tools', 'variant-interpretation', 
                              'publications', 'gene-assessment', 
                              'conditions', 'clinical-evidence', 
                              'community-classification-demo-app'
                              ], 
                        default='assessment-tools', 
                        help='Franklin page to access') 
    
    #12 SpliceAI lookup option - Masked or raw scores (default: masked)
    parser.add_argument('--spliceai-raw', required=False,
                        action='store_false', default=True, dest='maskFlag',
                        help='Use raw scores')
    
    #13 SpliceAI lookup option - Max distance (default: 10000)
    parser.add_argument('--spliceai-dist', required=False, type=int,
                        default=10000, help='Max distance for SpliceAI lookup')
    

    #14 Liftover option
    parser.add_argument('--liftover', required=False,
                        action='store_true', default=False, dest='liftoverFlag',
                        help='Add liftover process (hg19 -> hg38)')
    
    #15 Windows option
    parser.add_argument('--windows', required=False,
                        action='store_true', default=False, dest='windowsFlag',
                        help='For using on Windows OS')
                        
    return vars(parser.parse_args())


def analyze_args(args: dict) -> dict:
    parsed_args ={}

    parsed_args['input'] = args['input']
    parsed_args['gene_col'] = args['gene_col']
    parsed_args['alt_col'] = args['alt_col']
    parsed_args['assembly'] = args['assembly']
    parsed_args['franklin_page'] = args['franklin_page']
    parsed_args['spliceai_maskFlag'] = args['maskFlag']
    parsed_args['spliceai_dist'] = args['spliceai_dist']
    parsed_args['liftover'] = args['liftoverFlag']
    parsed_args['windowsFlag'] = args['windowsFlag']

    # Output    
    if args['output'] is None:
        input_str = str(parsed_args['input'])
        parsed_args['output'] = input_str.replace('.xlsx', '_hyperlinked.xlsx')
    else:
        parsed_args['output'] = args['output']

    # Skip sheets
    if args['skip_sheets'] is None:
        parsed_args['skip_sheets'] = []
    else:
        parsed_args['skip_sheets'] = args['skip_sheets'].split(',')

    # Skip sites
    if args['skip_sites'] is None:
        parsed_args['skip_sites'] = []
    else:
        parsed_args['skip_sites'] = args['skip_sites'].split(',')

    # UCSC genome browser displaying width
    if args['ucsc_width'] is None:
        parsed_args['ucsc_width'] = 45
    else:
        parsed_args['ucsc_width'] = args['ucsc_width']

    # Position columns
    if ((args['assembly'] == 'hg19') | (args['assembly'] == 'GRCh37')):
        if args['pos19'] is None:
            parsed_args['pos19'] = 'POS'
        else:
            parsed_args['pos19'] = args['pos19']
        if args['pos38'] is None:
            parsed_args['pos38'] = 'POS_hg38'
        else:
            parsed_args['pos38'] = args['pos38']
            
    elif ((args['assembly'] == 'hg38') | (args['assembly'] == 'GRCh38')):
        if args['pos19'] is None:
            parsed_args['pos19'] = 'POS_hg19'
        else:
            parsed_args['pos19'] = args['pos19']
        if args['pos38'] is None:
            parsed_args['pos38'] = 'POS'
        else:
            parsed_args['pos38'] = args['pos38']

    else:
        logger.error(f"Invalid assembly: {args['assembly']}")
        raise ValueError('Invalid assembly')
    
    # Showing arguments in log
    if parsed_args['spliceai_maskFlag']:
        spliceai_score = 'masked'
    else:
        spliceai_score = 'raw'

    logger.info(
        f"""

                              ## Arguments ##
                              Input         : {parsed_args['input']}
                              Output        : {parsed_args['output']}
                              Assembly      : {parsed_args['assembly']}
                              Gene column   : {parsed_args['gene_col']}
                              Skip sheets   : {parsed_args['skip_sheets']}
                              Skip sites    : {parsed_args['skip_sites']}
                              UCSC width    : {parsed_args['ucsc_width']}
                              POS19         : {parsed_args['pos19']}
                              POS38         : {parsed_args['pos38']}
                              SpliceAI score: {spliceai_score}
                              SpliceAI dist : {parsed_args['spliceai_dist']}
                              Windows URL   : {parsed_args['windowsFlag']}
        """
        )

    return parsed_args

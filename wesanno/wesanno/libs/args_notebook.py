import argparse
import pathlib2
import os

from logging import getLogger

logger = getLogger(__name__)

params = [
    '--input', '/Volumes/vol/work/Github/TestData/proband/31070/annovar/exome_summary.20220607_145914.txt',
    '--xhmm', '/work/Github/TestData/proband/xhmm/data.segdup.strvar.haplo.deciph.omim.xcnv.gene.uniq',
    # '--resources', '/Volumes/resources'
    '--mode', 'quad_unaffected'
    ]


def print_args(args: dict):
    logger.info('Arguments are as follows:\n')
    for k, v in args.items():
        print(f'     {k}: {v}')
    print('\n')


def parser_setting() -> dict:
    module_path = pathlib2.Path(os.path.dirname(__file__))

    parser = argparse.ArgumentParser(description = 'WES')
    parser.add_argument('--input', '-i', required=True, 
                        type=pathlib2.Path, help='path to exome_summary file')
    parser.add_argument('--output', '-o', required=False,
                        type=pathlib2.Path, help='path to root directory for output files')
    parser.add_argument('--xhmm', '-x', required=False, 
                        type=pathlib2.Path, help='path to a directry including XHMM files')
    # parser.add_argument('--vcf', '-v', required=False, 
    #                     type=pathlib2.Path, help='path to vcf file')
    parser.add_argument('--phenotype', '-p', required=False, 
                        type=pathlib2.Path, help='path to a file including phenotypes')
    parser.add_argument('--mode', '-m', required=False, default='auto', 
                        choices=[
                            'proband', 'trio', 'duo', 
                            'quad_affected', 'quad_unaffected', 'auto'
                            ], 
                        help='analyze mode')                        
    parser.add_argument('--samples', '-s', required=False, default='auto',
                        nargs='*', help='order of samples')      

    parser.add_argument('--assembly', '-a', required=False, default='hg19',
                        choices=['hg19', 'hg38', 'GRCh37', 'GRCh38'], 
                        help='genome assembly')        
    
    parser.add_argument('--config', '-c', required=False, 
                        default=f'{module_path.parents[1]}/config/config.toml',
                        type=pathlib2.Path,  help='path to config file')
    parser.add_argument('--resources', '-r', required=False, 
                        default=f'{module_path.parents[1]}/resources',
                        type=pathlib2.Path,  help='path to root directory for databases')
    
    #Options setting for additional annotations
    parser.add_argument('--no-gnomad', action='store_true', required=False, 
                        help='Not annotate gnomAD constraint scores')
    parser.add_argument('--no-hgmd', action='store_true', required=False, 
                        help='Not annotate HGMD infomation such as DM count') 
    parser.add_argument('--no-decipher', action='store_true', required=False, 
                        help='Not annotate DECIPHER data')
    parser.add_argument('--no-ddg2p', action='store_true', required=False, 
                        help='Not annotate DDG2P data')
    parser.add_argument('--no-jarvis', action='store_true', required=False, 
                        help='Not annotate JARVIS data')
    parser.add_argument('--no-spliceai', action='store_true', required=False, 
                        help='Not annotate SpliceAI pre-computed scores') 
    parser.add_argument('--no-syno', required=False, 
                        help='Not create synonymous analysis sheet')
    parser.add_argument('--no-alphamissense', action='store_true', required=False, 
                    help='Not annotate SpliceAI pre-computed scores') 
    parser.add_argument('--no-revel', action='store_true', required=False, 
                    help='Not annotate SpliceAI pre-computed scores') 

    parser.add_argument('--no-trap', action='store_true', required=False, 
                    help='Not annotate SpliceAI pre-computed scores') 
    
    # Option setting for Excel formatings
    parser.add_argument('--excel-formating', required=False, 
                        default=True, help='Not decorate Excel')
    

    # args = vars(parser.parse_args())
    args = vars(parser.parse_args(params))

    print_args(args)

    return args
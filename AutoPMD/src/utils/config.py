import os

class config: 

    S_GLOBAL_ARGS = {
            'PRIMER_OPT_SIZE': 29,   
            'PRIMER_MIN_SIZE': 27,
            'PRIMER_MAX_SIZE': 36,
            'PRIMER_OPT_TM': 65.0,
            'PRIMER_MIN_TM': 60.0,
            'PRIMER_MAX_TM': 75.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
    }


    D_GLOBAL_ARGS = {
                'PRIMER_OPT_SIZE': 29,
                'PRIMER_MIN_SIZE': 27,
                'PRIMER_MAX_SIZE': 36,
                'PRIMER_OPT_TM': 65.0,
                'PRIMER_MIN_TM': 60.0,
                'PRIMER_MAX_TM': 75.0,
                'PRIMER_MIN_GC': 20.0,
                'PRIMER_MAX_GC': 80.0,
        }

    #单双点设计引物的全局参数  
    GLOBAL_ARGS = {
                'PRIMER_PICK_ANYWAY':0,
                'PRIMER_PRODUCT_SIZE_RANGE': 0,
                'PRIMER_NUM_RETURN':10
        }  

    #测序引物设计的全局参数
    Q_ARGS = {
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_MIN_SIZE': 18,   
                'PRIMER_MAX_SIZE': 25,
                'PRIMER_OPT_TM': 65.0,
                'PRIMER_MIN_TM': 55.0,
                'PRIMER_MAX_TM': 75.0,   
                'PRIMER_MIN_GC': 20.0,
                'PRIMER_MAX_GC': 80.0,
    }

    Q_GLOBAL_ARGS = {   
                'PRIMER_PICK_ANYWAY':0,  
                'PRIMER_TASK': 'pick_primer_list',
        }
      
    #提取target_gene_seq和marker_seq需要的常量
    AMPLICONIC_GENE_TARGET_SEQ_LENGTH = 20
    
    #提取marker时需要的常量
    AMPLICONIC_MARKER_SEQ_LENGTH = 100
    AMPLICONIC_MARKER_SEQ_START = 0

    #
    TARGETGENE_AFTER_BEFORE_SEQ_N = 20 

    INPUT_FILE_PATH = ''   
    OUTPUT_FILE_PATH =''    

    PLASMID_FILE_NAME = ''    
    INPUT_FILE_NAME = ''  
    NEW_OUTPUT_FILE_PATH=''
    
    OUTPUT_SINGLE_UP_NAME = 'all_single_primer_up_success.xlsx' 
    OUTPUT_SINGLE_DOWN_NAME = 'all_single_primer_down_success.xlsx'   
    OUTPUT_DOUBLE_UP_NAME = 'all_double_primer_up_success.xlsx'   
    OUTPUT_DOUBLE_DOWN_NAME = 'all_double_primer_down_success.xlsx' 
    SEQUENCING_PRIMER_SUCCESS = 'sequencing_primer_success.xlsx'
    SEQUENCING_PRIMER_FAILTRUE = 'sequencing_primer_failtrue.xlsx'
    PCR_PRIMER_FAILTRUE = 'primer_failtrue.xlsx'
    SUMMARY_PTIMER = 'primer_information_summary.xlsx'
    OUTPUT_ORDERS_NAME = 'order.xlsx'       
    IMG_ORDERS_NAME = 'order1.png'

    OUTPUT_FILE_NAME_PLASMID_MUTATION = "plasmid_mutation.gb" 
    DATA_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")) + '/data'

      
import sys,os
from utils.config import config
from utils import util
import module.from_gbFile_to_seq as gb_seq
import module.organize_plasmid_primers as opp
import module.sequencing_primer as sp         
import module.single_double as sd
import pandas as pd  
from Bio import SeqIO
import json
import time
import argparse
from utils.config import config
from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath


def batch_plasmid(many_plasmid_mute_result,primer_output_path):
    '''
        根据突变质粒的结果字典, 批量输出突变质粒的GB文件
        params:
            many_plasmid_mute_result ([dict]) -- 突变质粒的结果字典,k--质粒模板的名字,v--质粒模板下的突变质粒字典:k--突变位点,v--突变质粒的GB文件
            primer_output_path ([String]) -- 
        returns:
            无
    '''
    for k,v in many_plasmid_mute_result.items():
        #判断是否创建了一个文件夹
        if os.path.exists(primer_output_path+k):
            print(f'{k}文件夹已存在')
        else:   
            os.makedirs(primer_output_path + k)  
            print(f'{k}文件夹创建')
        for site,plasmid in v.items():   
            # 输出gb文件   
            SeqIO.write(plasmid,primer_output_path + k + '/' +site, "genbank")  

#聚集引物输出      
def gather_primer(many_pcr_primer_result, many_sequencing_primer_result,num, outputdir = config.OUTPUT_FILE_PATH):
    '''
        批量输出PCR引物、PCR备用引物、测序引物
        params:
            输出路径
        returns: 
                
    '''
    #所有质粒模板的测序引物的设计成功、失败
    suc_sequencing_primer_df,sequencing_fail_primer_df = sequencing_primer(many_sequencing_primer_result,outputdir=config.OUTPUT_FILE_PATH)
    #所有质粒模板的pcr引物的设计成功、失败
    suc_pcr,fail_pcr = pcr_primer(many_pcr_primer_result,outputdir=config.OUTPUT_FILE_PATH)
    plasmid_primer_df,spare_plasmid_primer_df = suc_pcr
    pcr_s_u_fail_primer_df, pcr_s_d_fail_primer_df,pcr_d_u_fail_primer_df,pcr_d_d_fail_primer_df = fail_pcr

    fai_sequencing_primer_df = sequencing_fail_primer_df

    suc_sequencing_primer_df.reset_index(drop=True, inplace=True)
    spare_plasmid_primer_df.reset_index(drop=True, inplace=True)
    suc_sequencing_primer_df.reset_index(drop=True, inplace=True)
    fai_sequencing_primer_df.reset_index(drop=True, inplace=True)


    if not os.path.exists(outputdir):  
        os.makedirs(outputdir)  

    with pd.ExcelWriter(outputdir+'primer_information_summary.xlsx') as writer:  
        plasmid_primer_df.to_excel(writer,sheet_name = 'plasmid_primer',index_label='No.')      
        spare_plasmid_primer_df.to_excel(writer,sheet_name = 'spare_plasmid_primer',index_label='No.') 
        suc_sequencing_primer_df.to_excel(writer,sheet_name='sequencing_primer',index_label='No.') 
        fai_sequencing_primer_df.to_excel(writer,sheet_name='fail_sequencing_primer',index_label='No.')

    if num ==1:   
        #输出失败
        with pd.ExcelWriter(outputdir + 'primer_failtrue.xlsx') as writer:
            pcr_s_u_fail_primer_df = pcr_s_u_fail_primer_df.to_excel(writer,sheet_name = 'pcr_s_u_fail_primer',index_label='No.')
            pcr_s_d_fail_primer_df = pcr_s_d_fail_primer_df.to_excel(writer,sheet_name = 'pcr_s_d_fail_primer',index_label='No.')
            pcr_d_u_fail_primer_df = pcr_d_u_fail_primer_df.to_excel(writer,sheet_name = 'pcr_d_u_fail_primer',index_label='No.')
            pcr_d_d_fail_primer_df = pcr_d_d_fail_primer_df.to_excel(writer,sheet_name = 'pcr_d_d_fail_primer',index_label='No.')
            sequencing_fail_primer_df = sequencing_fail_primer_df.to_excel(writer,sheet_name = 'sequencing_fail_primer',index_label='No.')
          
def normalize_primer(primer_df,input_path_name,name='suc'):  
    #读取输入文件
    input_df = pd.read_excel(input_path_name)
    input_df.columns = input_df.columns.str.lower()
    input_df = input_df[['sample','template','mutagenesis','design']]
    


    if len(primer_df) >0:
        all_primer_df=input_df.merge(primer_df,left_on=['template','sample'],right_on=['template','id'])
        all_primer_df = all_primer_df.drop(columns=['id'])
        all_primer_df.fillna('null',inplace=True)

        #规范化处理每一个gb文件对应的引物集合
        if 'suc' in name:
            all_primer_df = all_primer_df.rename(columns={'sample':'Sample',
                                    'template':'Template',
                                    'mutagenesis':'Mutation',
                                    'design':'Design',
                                    'success_or_failtrue':'Success_or_Failtrue',
                                    "primer_up_f_seq_(5'-3')":"Primer-U-F_Sequence（5'-3')",
                                    "primer_up_r_seq_(5'-3')":"Primer-U-R_Sequence（5'-3')",
                                    "primer_up_f_Tm":"Primer-U-R_F_Tm",
                                    "primer_up_r_Tm":"Primer-U-R_R_Tm",
                                    "up_product_sequence_length":"Product_Length-U",
                                    "up_product_sequence":"Product_Seq-U",
                                    "primer_down_f_seq_(5'-3')":"Primer-D-F_Sequence（5'-3')",
                                    "primer_down_r_seq_(5'-3')":"Primer-D-R_Sequence（5'-3')",
                                    "primer_down_f_Tm":"Primer-D-F_Tm",
                                    "primer_down_r_Tm":"Primer-D-R_Tm",
                                    "down_product_sequence_length":"Product_Length-D",
                                    "down_product_sequence":"Product_Seq-D",
                            })
        elif 'u_fail' in name:
            all_primer_df = all_primer_df.rename(columns={  'sample':'Sample',
                                                            'template':'Template',
                                                            'mutagenesis':'Mutation',
                                                            'design':'Design',
                                                            'PRIMER_LEFT_EXPLAIN':'Pimer_U_F_Explain',
                                                            'PRIMER_RIGHT_EXPLAIN':'Primer_U_R_Explain',
                                                            'PRIMER_PAIR_EXPLAIN':'Primer_U_Pair_Explain'
                                                         })
        elif 'd_fail' in name:
            all_primer_df = all_primer_df.rename(columns={  'sample':'Sample',
                                                            'template':'Template',
                                                            'mutagenesis':'Mutation',
                                                            'design':'Design',
                                                            'PRIMER_LEFT_EXPLAIN':'Primer_D_F_Explain',
                                                            'PRIMER_RIGHT_EXPLAIN':'Primer_D_R_Explain',
                                                            'PRIMER_PAIR_EXPLAIN':'Primer_D_Pair_Explain'
                                                         })    
        
    else:
        all_primer_df = primer_df       

    return all_primer_df

def pcr_primer(many_pcr_primer_result,outputdir = config.OUTPUT_FILE_PATH):

    pcr_suc_primer_df = pd.DataFrame()
    pcr_suc_spare_primer_df = pd.DataFrame()
    pcr_s_u_fail_primer_df = pd.DataFrame()
    pcr_s_d_fail_primer_df = pd.DataFrame()
    pcr_d_u_fail_primer_df = pd.DataFrame()
    pcr_d_d_fail_primer_df = pd.DataFrame()

    

    for k,v in many_pcr_primer_result.items():   
        suc_tuple,fail_dict = v
             
        #成功的
        suc_primer_df, suc_spare_primer_df = suc_tuple
        pcr_suc_primer_df = pcr_suc_primer_df.append(suc_primer_df)
        pcr_suc_spare_primer_df = pcr_suc_spare_primer_df.append(suc_spare_primer_df)

        #失败的  
        if 'single_dict_up_primers_failtrue' in fail_dict.keys():
            pcr_s_u_fail_primer_df = pcr_s_u_fail_primer_df.append(fail_dict['single_dict_up_primers_failtrue'])
        if 'single_dict_down_primers_failtrue' in fail_dict.keys():
            pcr_s_d_fail_primer_df = pcr_s_d_fail_primer_df.append(fail_dict['single_dict_down_primers_failtrue'])
        if 'double_dict_up_primers_failtrue' in fail_dict.keys():
            pcr_d_u_fail_primer_df = pcr_d_u_fail_primer_df.append(fail_dict['double_dict_up_primers_failtrue'])
        if 'double_dict_down_primers_failtrue' in fail_dict.keys():
            pcr_d_d_fail_primer_df = pcr_d_d_fail_primer_df.append(fail_dict['double_dict_down_primers_failtrue'])   
       
    #根据索引整合input文件
    input_path_name = config.INPUT_FILE_PATH +config.INPUT_FILE_NAME  

    pcr_suc_primer_df = normalize_primer(pcr_suc_primer_df, input_path_name,'pcr_suc')   
    pcr_suc_spare_primer_df = normalize_primer(pcr_suc_spare_primer_df,input_path_name,'pcr_suc_spare')
    pcr_s_u_fail_primer_df = normalize_primer(pcr_s_u_fail_primer_df,input_path_name,'pcr_s_u_fail')
    pcr_s_d_fail_primer_df = normalize_primer(pcr_s_d_fail_primer_df,input_path_name,'pcr_s_d_fail')
    pcr_d_u_fail_primer_df = normalize_primer(pcr_d_u_fail_primer_df,input_path_name,'pcr_d_u_fail')
    pcr_d_d_fail_primer_df = normalize_primer(pcr_d_d_fail_primer_df,input_path_name,'pcr_d_d_fail')
  
    suc = pcr_suc_primer_df,pcr_suc_spare_primer_df
    fail = pcr_s_u_fail_primer_df, pcr_s_d_fail_primer_df,pcr_d_u_fail_primer_df,pcr_d_d_fail_primer_df
   
    return suc,fail   

#测序引物
def sequencing_primer(many_sequencing_primer_result,outputdir):

    '''
        
    '''  
    # v[0]:succcess,    v[1]:failtrue
    sequencing_primer = pd.DataFrame()
    failtrue_primer = pd.DataFrame()

    #给测序引物按照读入gb文件的顺序输出
    file_path = os.path.join(config.INPUT_FILE_PATH, config.INPUT_FILE_NAME)
    df = pd.read_excel(file_path)
    temp_df = df.drop_duplicates(subset='template')
    template = list(temp_df['template'])


    # for k,v in many_sequencing_primer_result.items():
    for i in template:
        k = i + '_sequencing_primer_result'
        v=many_sequencing_primer_result.get(k)
        if v != None:
            if len(v[0])!=0:
                sequencing_primer = sequencing_primer.append(v[0])
            if len(v[1])!=0:
                failtrue_primer = failtrue_primer.append(v[1])
    


    #规范化s
    if len(sequencing_primer) > 0:
        for i in range(1,len(sequencing_primer.columns)):  
            cols = list(sequencing_primer.columns)   
            num = cols[i].split('_')[-1]
            
            if int(num) == 2:   
                sequencing_primer[f'Sequencing_primer_{num}_name'] = sequencing_primer['site_id']+f'_seq{num}_R'
            else:
                sequencing_primer[f'Sequencing_primer_{num}_name'] = sequencing_primer['site_id']+f'_seq{num}_F' 

            if 'SEQUENCING_PRIMER' in cols[i]:  
                sequencing_primer = sequencing_primer.rename(columns={
                    cols[i]:f"Sequencing_primer_{num}_Sequence(5'-3')"
                })

        sequencing_primer = sequencing_primer.rename(columns={
                'site_id':'Sample'
        })


        #定义顺序
        columns = list(sequencing_primer.columns)
        aa = sorted(columns)
        sequencing_primer = sequencing_primer[aa]
    #输出
    # sequencing_primer.to_excel(outputdir+config.SEQUENCING_PRIMER_SUCCESS,index=False)
    # if len(failtrue_primer)!=0:
    #     failtrue_primer.to_excel(outputdir+config.SEQUENCING_PRIMER_FAILTRUE,index=False)

   
    return sequencing_primer, failtrue_primer    

def position_check(input_mute_tuple_df, before_processed_seq_dict):

    single_mute_df, double_mute_df = input_mute_tuple_df
    target_gene_seq = before_processed_seq_dict['target_gene_seq']
    single_failture_df =  pd.DataFrame()
    double_failture_df = pd.DataFrame()

    if len(single_mute_df) > 0:
        for i,v in single_mute_df.iterrows():
            position = int(v['position1']) -1
            before = v['before1']
            seq = target_gene_seq[position:position+3]
            if before != seq:
                single_failture_df = single_failture_df.append(v)

    if len(double_mute_df) > 0:
        for i,v in double_mute_df.iterrows():
            position1 = int(v['position1']) - 1
            before1 = v['before1']  
            seq1 =  target_gene_seq[position1:position1+3]

            position2 = int(v['position2']) - 1
            before2 =  v['before2']
            seq2 =  target_gene_seq[position2:position2+3]

            if before1 != seq1 or before2 != seq2:
                double_failture_df = double_failture_df.append(v)
    
    return single_failture_df, double_failture_df
  
def everyone_plasmid_process(num, input_file_path=config.INPUT_FILE_PATH + config.INPUT_FILE_NAME, primer_output_path= config.OUTPUT_FILE_PATH):
    '''
        根据引物的设计结果生成突变质粒  
    '''

    input_file_path = config.INPUT_FILE_PATH + config.INPUT_FILE_NAME
    primer_output_path = config.OUTPUT_FILE_PATH 

    #批量生成引物
    # batch_primer(input_file_path,primer_output_path) 

    #提取多个质粒信息  
    many_gb_dict = util.extract_many_gb(input_file_path = config.INPUT_FILE_PATH + config.INPUT_FILE_NAME)
    many_plasmid_mute_result = {}
    many_sequencing_primer_result = {}
    many_pcr_primer_result = {}


    #位置问题
    failture_position_df = pd.DataFrame()

    for k,v in many_gb_dict.items():  
         #提取gb文件信息
        before_processed_seq_dict, after_processed_seq_dict = gb_seq.get_data_from_genebank(   
                                                                    infile = f'{k}.gb',
                                                                    marker=v['marker'],
                                                                    target_gene=v['target'])
        #读取每个质粒的点突变输入数据,且转换氨基酸成核苷酸
        input_mute_tuple_df = util.read_input_file_by_gb(input_file_path, gb_name = v['gb_name'])
        single_mute_df, double_mute_df = input_mute_tuple_df

        single_failture_df, double_failture_df = position_check(input_mute_tuple_df, before_processed_seq_dict)
        single_double_failture_df = double_failture_df.append(single_failture_df)
        failture_position_df = failture_position_df.append(single_double_failture_df)
        failture_position_df['template'] = v['gb_name']

        
        #取出位置错误的df 
        temp = pd.concat([single_mute_df, single_failture_df])
        single_mute_df = temp.drop_duplicates(['sample','target'],keep=False)
        temp= pd.concat([double_mute_df, double_failture_df])
        double_mute_df = temp.drop_duplicates(['sample','target'],keep=False)
        input_mute_tuple_df = single_mute_df,double_mute_df


        #设计引物 
        if len(single_mute_df) >0 or len(double_mute_df) > 0:  
            result_output_success_failtrue_dict = sd.design_process(input_mute_tuple_df,after_processed_seq_dict, before_processed_seq_dict)
            primer_success_dict_df, primer_failtrue_dict_df = util.one_plasmid_primer_result_to_df(v['gb_name'],result_output_success_failtrue_dict) 
            primer_df,spare_primer_df = util.merge_primer(primer_success_dict_df,gb_name = v['gb_name'],type='one')

            #重组质粒
            mute_after_plasmid_seq_dict = opp.concat_mute_plasmid_seq(input_mute_tuple_df,before_processed_seq_dict) 

            #针对每个质粒生成点突变的gb文件
            primer_gb_file_dict = opp.write_gb_file(primer_df, mute_after_plasmid_seq_dict, f'{k}.gb')  
            spare_primer_gb_file_dict = opp.write_gb_file(spare_primer_df, mute_after_plasmid_seq_dict,f'{k}.gb')
    
            #合并gb文件，生成点突变质粒
            site_mute_plasmid = opp.merge_gb_output_file(primer_gb_file_dict, spare_primer_gb_file_dict) 

            #设计测序引物
            suc_primer_df, fail_primer_df = sp.high_sequencing_primer(mute_after_plasmid_seq_dict,input_mute_tuple_df)
            if len(suc_primer_df)>0:
                suc_primer_df = suc_primer_df.set_index(['site_id',"primer_id"])['primer'].unstack().rename_axis(columns=None).reset_index()

            #生成一个质粒的设计结果
            many_plasmid_mute_result[v['gb_name']+'_site_mute'] = site_mute_plasmid
            many_sequencing_primer_result[v['gb_name']+'_sequencing_primer_result'] = (suc_primer_df,fail_primer_df)

            #引物设计结果
            suc_primer_df, suc_spare_primer_df = util.merge_primer(primer_success_dict_df, gb_name = v['gb_name'], type='two')
            many_pcr_primer_result[v['gb_name']+'_pcr_primer_result'] = (suc_primer_df, suc_spare_primer_df), primer_failtrue_dict_df
        

    #批量生成gb文件
    batch_plasmid(many_plasmid_mute_result, primer_output_path)
    #生成汇总引物信息
    gather_primer(many_pcr_primer_result, many_sequencing_primer_result, num, outputdir = config.OUTPUT_FILE_PATH)  
    return many_sequencing_primer_result, many_plasmid_mute_result, many_pcr_primer_result ,failture_position_df

  
def main(input_file_path):
    with open(input_file_path,'r',encoding='utf8') as fp:
        json_params = json.load(fp)
        config.OUTPUT_FILE_PATH = config.DATA_ROOT + json_params['outputdir']
        config.INPUT_FILE_PATH = config.DATA_ROOT + json_params['inputdir']  
        config.INPUT_FILE_NAME = json_params['input_mute_name']
        config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH = json_params['global_params']['AMPLICONIC_GENE_TARGET_SEQ_LENGTH']
        config.AMPLICONIC_MARKER_SEQ_LENGTH = json_params['global_params']['AMPLICONIC_MARKER_SEQ_START_LENGTH'][1]
        config.AMPLICONIC_MARKER_SEQ_START = json_params['global_params']['AMPLICONIC_MARKER_SEQ_START_LENGTH'][0]  
        config.S_GLOBAL_ARGS = json_params['pcr_single_primer_params']    
        config.D_GLOBAL_ARGS = json_params['pcr_double_primer_params']
        config.Q_ARGS = json_params['seq_primer_params']
        config.TARGETGENE_AFTER_BEFORE_SEQ_N = json_params['targetGene_after_before_seq_n']
        
        # TARGETGENE_AFTER_BEFORE_SEQ_N = 20 
        # config.NEW_OUTPUT_FILE_PATH = config.DATA_ROOT + json_params['new_outputdir']  

        #第一次调用主程序配置
        config.GLOBAL_ARGS['PRIMER_PICK_ANYWAY'] = 0
        config.Q_GLOBAL_ARGS['PRIMER_PICK_ANYWAY'] = 0   
        many_sequencing_primer_result, many_plasmid_mute_result, many_pcr_primer_result, failture_position_df = everyone_plasmid_process(num=1)
        #第二次调用主程序配置  
        config.GLOBAL_ARGS['PRIMER_PICK_ANYWAY'] = 1
        config.Q_GLOBAL_ARGS['PRIMER_PICK_ANYWAY'] = 1  
        num=2   
        many_sequencing_primer_result, many_plasmid_mute_result, many_pcr_primer_result, failture_position_df = everyone_plasmid_process(num)   

        #输出新的汇总   
        if len(failture_position_df) > 0:
            failture_position_df['failture_reason'] = 'The position of mutation is wrong' 
            all_input_df = pd.read_excel(config.INPUT_FILE_PATH + config.INPUT_FILE_NAME)
            failture_position_df.reset_index(drop=True,inplace=True)
            all_input_df.reset_index(drop=True,inplace=True)   
            failture_position_df = pd.merge(all_input_df,failture_position_df, on=['sample','template'], how='inner', suffixes=('', '_'))
            failture_position_df = failture_position_df[['sample', 'template', 'mutagenesis','failture_reason']]      
            failture_position_df.rename(columns={'sample':'Sample','template':"Template",'mutagenesis':'Mutagenesis'},inplace=True)
           
        util.read_failtrue_df_and_merge(sum_info_file = config.SUMMARY_PTIMER, failtrue_file = config.PCR_PRIMER_FAILTRUE, outputdir = config.OUTPUT_FILE_PATH, failture_position_df = failture_position_df)  


if  __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', help='input params file', required=True)     
    arguments = parser.parse_args()
    input_file_path = arguments.input
    
    time1=time.time()    
    main(input_file_path)
    time2=time.time()
    print('高通量引物设计已完成用时：'+str(time2-time1)+'s')     
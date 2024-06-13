import primer3
import pandas as pd
import sys,os
sys.path.append('../')
from utils.config import config   
from utils import util as ut
import warnings   
warnings.filterwarnings('ignore')

import re
    
#读取输入文件 
def read_input_file(input_file_path=config.INPUT_FILE_PATH+config.INPUT_FILE_NAME):
    input_df = pd.read_excel(input_file_path)  
    input_df.columns = input_df.columns.str.lower()
    try :
        double_mute_df = input_df[input_df['mutagenesis'].str.contains(';')]   
        single_mute_df = input_df[~input_df['mutagenesis'].str.contains(';')]    

        if len(single_mute_df) > 0 :
            single_mute_df['mutagenesis']= single_mute_df.mutagenesis.str.split('-')
            single_mute_df = single_mute_df[['sample','marker','target','mutagenesis']]
            new_single_df = pd.DataFrame(single_mute_df.mutagenesis.to_list(),columns=['before','position','after'])
            #重置索引
            single_mute_df.reset_index(drop=True,inplace=True)
            new_single_df.reset_index(drop=True,inplace=True)
            single_mute_df = pd.concat([single_mute_df,new_single_df],axis=1)
 
        def work(x):
            arr = []
            for item in x:
                arr.extend(item.split('-'))
            return arr
        double_mute_df['mutagenesis']=double_mute_df.mutagenesis.str.split(';').apply(lambda x: work(x))
        double_mute_df = double_mute_df[['sample','marker','target','mutagenesis']]
        double_mute_df.reset_index(drop=True,inplace=True)
        new_double_df = pd.DataFrame(double_mute_df.mutagenesis.to_list(),columns=['before1','position1','after1','before2','position2','after2'])
        double_mute_df = pd.concat([double_mute_df,new_double_df],axis=1)
        return single_mute_df, double_mute_df
    
    except Exception as e:
        print(e)   
        single_mute_df = input_df
        single_mute_df['mutagenesis']= single_mute_df.mutagenesis.str.split('-')
        single_mute_df = single_mute_df[['sample','marker','target','mutagenesis']]
        new_single_df = pd.DataFrame(single_mute_df.mutagenesis.to_list(),columns=['before','position','after'])
        single_mute_df = pd.concat([single_mute_df,new_single_df],axis=1)
        return single_mute_df
   
#取单位点突变的引物模板
def single_input_to_dict_temp(  single_mute_df,
                                dict_plasmid_seq,
                                dict_left_primers_succ=None  
                                ):     
    dict_temp={}
    marker_seq=dict_plasmid_seq['marker_seq']
    target_gene_up_seq=dict_plasmid_seq['target_gene_up_seq']
    target_gene_down_seq=dict_plasmid_seq['target_gene_down_seq']     
    target_gene_seq=dict_plasmid_seq['target_gene_seq']
   

    temp_num =  + config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH
      
    for index,row in single_mute_df.iterrows():    
            site_id = row['sample']
            position = int(row['position1'])  
            seq_after_mut = row['after1']
            print(single_mute_df)   

            if dict_left_primers_succ: #确定target下游引物设计模板，需用target_up的左引物在marker上定位
                print(dict_left_primers_succ.keys())   
               
                if row['sample'] in list(dict_left_primers_succ.keys()):
                    up_left_primer=dict_left_primers_succ[site_id]['PRIMER_LEFT_0_SEQUENCE']
                    down_temp_right_point=marker_seq.find(up_left_primer)+ 20      #加入重叠区域20bp
                #    seq_temp = target_gene_seq[position+10:position+19]+seq_after_mut+target_gene_seq[position+22:]+target_gene_down_seq+marker_seq[:down_temp_right_point]
                    seq_temp = target_gene_seq[position + config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH - 10 : position + config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH - 1 ] + seq_after_mut + target_gene_seq[position+ config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH + 2:] + target_gene_down_seq + marker_seq[:down_temp_right_point]

                    dict_temp[site_id]=seq_temp   
                
            else:   #确定target上游引物设计模板
                temp = marker_seq+target_gene_up_seq
                
                # seq_temp=temp + target_gene_seq[:position+19] + seq_after_mut + target_gene_seq[position+22:position+31] 
                seq_temp=temp + target_gene_seq[:position + config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH - 1 ] + seq_after_mut + target_gene_seq[position + config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH + 2 :position + config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH + 2 + 9  ] 

                dict_temp[site_id]=seq_temp   

        
    return dict_temp   
  
#取双位点突变引物模板
def double_input_to_dict_temp(  double_mute_df,  
                                dict_plasmid_seq,
                                stype='target',
                                dict_left_primers_succ=None
                                ):

    
    dict_temp={}
    marker_seq=dict_plasmid_seq['marker_seq']
    target_gene_up_seq=dict_plasmid_seq['target_gene_up_seq']
    target_gene_down_seq=dict_plasmid_seq['target_gene_down_seq']     
    target_gene_seq=dict_plasmid_seq['target_gene_seq']
    plasmid_seq=target_gene_down_seq+marker_seq+target_gene_up_seq
   
    for index,row in double_mute_df.iterrows():  
        site_id = row['sample']   

        if site_id == 'site20':
            print("jkgfdsvjkl")
            pass


        position1 = int(row['position1'])   
        seq_after_mut1=row['after1']
        position2=int(row['position2'])   
        seq_after_mut2=row['after2']
        distance_bp = abs(int(position2) - int(position1+2)) -1    #突变位点之间的碱基个数

        if distance_bp <= 120:   
            overlap = distance_bp + 6
            max_overlap = 36
            surplus_bp = max_overlap - overlap
            if surplus_bp % 2 != 0:
                surplus_bp = surplus_bp + 1
            one_surplus_bp = int(surplus_bp / 2)

            if 21 <= distance_bp and distance_bp <= 120:               #双点突变的第一种情况
                if dict_left_primers_succ:   #确定下游引物
                    if row['sample'] in list(dict_left_primers_succ.keys()):
                        up_left_primer=dict_left_primers_succ[site_id]['PRIMER_LEFT_0_SEQUENCE']
                        down_temp_right_point = marker_seq.find(up_left_primer) + 20
                        seq_temp = seq_after_mut2 + target_gene_seq[  config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position2 + 3  :] + target_gene_down_seq + marker_seq[:down_temp_right_point]
                        dict_temp[site_id]=seq_temp
                else: #确定上游引物模板 
                        seq_temp = marker_seq + target_gene_up_seq + target_gene_seq[:  config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position1 ] + seq_after_mut1  
                        dict_temp[site_id]=seq_temp
            elif 0 <=distance_bp and distance_bp <21:                   #双点突变的第二种情况  
                if dict_left_primers_succ:   #确定下游引物
                    if row['sample'] in list(dict_left_primers_succ.keys()):
                        up_left_primer=dict_left_primers_succ[site_id]['PRIMER_LEFT_0_SEQUENCE']
                        down_temp_right_point = marker_seq.find(up_left_primer) + 20      #重叠区域20bp
    
                        temp1 = target_gene_seq[  config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position1 - one_surplus_bp :   config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position1  ] + seq_after_mut1 
                        distance_seq =  target_gene_seq[config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position1 + 3 : config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position1 + 3 +  distance_bp ]
                        temp2 =  seq_after_mut2  + target_gene_seq[ config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position2 + 3 :  ]  +  target_gene_down_seq + marker_seq[:down_temp_right_point] 

                        dict_temp[site_id] = temp1 + distance_seq + temp2  
                else: #确定上游引物模板 
                        temp1 = marker_seq + target_gene_up_seq + target_gene_seq[:  config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position1 ] + seq_after_mut1 
                        distance_seq = target_gene_seq[config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position1 + 3 : config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position1 + 3 +  distance_bp ]
                        temp2 = seq_after_mut2 + target_gene_seq[ config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position2 +3 : config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position2 +3 + one_surplus_bp  ]
                        seq_temp = temp1 + distance_seq + temp2
                        dict_temp[site_id]=seq_temp
        else:                                                           #双点突变的第三种情况
            if stype=='target': 
                        temp1 = target_gene_seq[ config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position1 - 9  : config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position1 ] + seq_after_mut1
                        distance_seq =  target_gene_seq[config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position1 + 3 : config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position1 + 3 +  distance_bp ]
                        temp2 = seq_after_mut2 + target_gene_seq[ config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position2 + 3  : config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position2 + 3 +9 ]
                        temp_seq = temp1 + distance_seq + temp2
                        # temp_seq = target_gene_seq[position1+10:position1+19] + seq_after_mut1 + distance_seq + seq_after_mut2 + target_gene_seq[position2+22:position2+31]
            elif stype=='plasmid':
                        temp1 = target_gene_seq[ config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position2 - 9 : config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position2 ] + seq_after_mut2 +   target_gene_seq[ config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position2  + 3 : ]
                        temp2 = plasmid_seq + target_gene_seq[:config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position1: ] + seq_after_mut1 + target_gene_seq[  config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position1 + 3 : config.AMPLICONIC_GENE_TARGET_SEQ_LENGTH-1 + position1 + 3 + 9 ]
                        temp_seq = temp1 + temp2
                        # temp_seq = target_gene_seq[position2+10:position2+19] + seq_after_mut2 + target_gene_seq[position2+22:] + plasmid_seq + target_gene_seq[:position1+19] + seq_after_mut1 + target_gene_seq[position1+22:position1+31]
            dict_temp[site_id]=temp_seq
   
    return dict_temp 


def filter_primer_design_result(primer3_result):

    primer_desgin = primer3_result
    numbers=[''.join(re.findall(r'\d+', k))  for k in primer_desgin.keys()]  
    numbers = set(numbers)

    numbers = [item for item in numbers if item != '']

    left_primer = ''
    right_primer = ''

    li = []
    for i,item in enumerate(numbers):
        
        temp = {  f'PRIMER_LEFT_{i}_SEQUENCE': primer_desgin[ f'PRIMER_LEFT_{i}_SEQUENCE' ],
                f'PRIMER_RIGHT_{i}_SEQUENCE': primer_desgin[ f'PRIMER_RIGHT_{i}_SEQUENCE' ],
                f'PRIMER_LEFT_{i}_TM': primer_desgin[ f'PRIMER_LEFT_{i}_TM'],
                f'PRIMER_RIGHT_{i}_TM': primer_desgin[ f'PRIMER_RIGHT_{i}_TM'],
                f'PRIMER_LEFT_{i}': primer_desgin[ f'PRIMER_LEFT_{i}' ],
                f'PRIMER_RIGHT_{i}': primer_desgin[ f'PRIMER_RIGHT_{i}'],
                f'PRIMER_PAIR_{i}_PRODUCT_SIZE' : primer_desgin[ f'PRIMER_PAIR_{i}_PRODUCT_SIZE' ]
            }
        
        if i == 0:
            left_primer = temp[ f'PRIMER_LEFT_{i}_SEQUENCE']
            right_primer = temp[ f'PRIMER_RIGHT_{i}_SEQUENCE']
            li.append(temp)
        else:
            if left_primer != temp[ f'PRIMER_LEFT_{i}_SEQUENCE'] and right_primer != temp[ f'PRIMER_RIGHT_{i}_SEQUENCE']:
                li.append(ut.replace_digits_in_keys(temp, '1'))
                break
            elif i == len(numbers)-1:
                li.append(ut.replace_digits_in_keys(temp, '1'))
    
    li[0].update(li[1])
    primer3_result = li[0]

    return primer3_result

#设计引物
def primer_design(seqId,
                  seqTemplate,
                  stype,
                  mute_type='single',
                  global_args=config.GLOBAL_ARGS):   
    if mute_type == 'single':
        config.GLOBAL_ARGS.update(config.S_GLOBAL_ARGS)
        global_args=config.GLOBAL_ARGS    
    elif mute_type =='double':
        config.GLOBAL_ARGS.update(config.D_GLOBAL_ARGS)
        global_args=config.GLOBAL_ARGS
    elif mute_type == 'sequencing': 
        config.Q_GLOBAL_ARGS.update(config.Q_ARGS)
        global_args=config.Q_GLOBAL_ARGS  
    
    #序列参数
    seqlength = len(seqTemplate)
    seq_args = {
                'SEQUENCE_ID': seqId,
                'SEQUENCE_TEMPLATE': seqTemplate,   
                'SEQUENCE_FORCE_LEFT_START':-1,
                'SEQUENCE_FORCE_RIGHT_START': seqlength-1,
                'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [[0,0,0,0]]
        }
    #选择类型，设定序列参数
    if mute_type == 'single': #单点突变
        if stype == "left":   #target上游引物设计
           
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [[0,100,-1,-1]]
            size_range = [seqlength-73,seqlength]
        elif stype == "right":   #target下游引物设计
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [[-1,-1,seqlength-36,36]]  
            seq_args['SEQUENCE_FORCE_LEFT_START'] = 0 
            size_range = [seqlength-2,seqlength]     
        primer_num = 10      
    elif mute_type == 'double3':  #双点突变
        seq_args['SEQUENCE_FORCE_LEFT_START'] = 0  
        seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [[0,35,-1,-1]]
        size_range = [seqlength,seqlength]     
        primer_num = 10
    elif mute_type == 'sequencing':  #测序引物  
        seq_args['SEQUENCE_ID'] = seqId
        seq_args['SEQUENCE_TEMPLATE'] = seqTemplate
        size_range = [25,seqlength]   
        primer_num = 1
    elif mute_type == 'double1_2':  #双点突变中的第一、二情况
        if stype == 'left':
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST']= [[0,100,seqlength-36,36]]
            size_range = [seqlength-73,seqlength]    #设置引物的最小长度27
        else:
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [[-1,-1,seqlength-36,36]]
            seq_args['SEQUENCE_FORCE_LEFT_START'] = 0   
            size_range = [seqlength,seqlength]
        primer_num = 10
  
    #设定全局参数       
    global_args['PRIMER_PRODUCT_SIZE_RANGE']=size_range
    global_args['PRIMER_NUM_RETURN']= primer_num
    #调用工具
    primer3_result = primer3.bindings.designPrimers(seq_args, global_args) 

    if len(primer3_result)>10:
        primer3_result = filter_primer_design_result(primer3_result)
   
    return primer3_result  
      
#取引物属性  
def get_primer_attribute(dict_primers_succ,dict_temp,double_df=pd.DataFrame(),primer_num_return = 2,type='double2',stype='up'):
    if type == 'double2' and len(double_df) > 0 :
        double_df = double_df.set_index('sample')
        double_df['distance'].astype(int)
  
    primers_dict={}
    for key in dict_primers_succ:
        if type == 'double2' and  len(double_df) > 0 :
            distance_seq = double_df.loc[key,'distance_seq']
            distance = double_df.loc[key,'distance']  
            overlap_seq= ''
            if int(distance) % 2 == 0:
                overlap_seq = distance_seq[int(distance/2-10):int(distance/2+10)]
            else:
                distance = distance + 1
                overlap_seq = distance_seq[int(distance/2-10):int(distance/2+10)]


        li = []
        for i in range(primer_num_return):
        
            #取左右引物tm
            left_primer_tm= dict_primers_succ[key][f'PRIMER_LEFT_{i}_TM']
            right_primer_tm = dict_primers_succ[key][f'PRIMER_RIGHT_{i}_TM']
            #取左右引物seq
            left_primer_seq=dict_primers_succ[key][f'PRIMER_LEFT_{i}_SEQUENCE']
            right_primer_seq=dict_primers_succ[key][f'PRIMER_RIGHT_{i}_SEQUENCE']
            #取互补链
            right_primer_rev_seq=ut.revComp(right_primer_seq)
            temp_seq=dict_temp[key]
            #取引物的起始-----终止坐标
            product_left_point=temp_seq.find(left_primer_seq)
            product_right_point=temp_seq.find(right_primer_rev_seq)+len(right_primer_seq)-1
            #取产物序列
            product_seq=temp_seq[product_left_point:product_right_point+1]

            if type == 'double2' and stype == 'up' and len(double_df) > 0:
                real_right_primer_rev_seq =right_primer_rev_seq + distance_seq[:distance_seq.find(overlap_seq)+len(overlap_seq)]
                right_primer_seq = ut.revComp(real_right_primer_rev_seq)
                product_seq = product_seq + distance_seq[:distance_seq.find(overlap_seq)+len(overlap_seq)]
            elif type == 'double2' and stype == 'down' and len(double_df) > 0:
                left_primer_seq = distance_seq[distance_seq.find(overlap_seq):]+left_primer_seq
                product_seq = distance_seq[distance_seq.find(overlap_seq):] +product_seq
            
            li.append({'product_seq':product_seq,
                       'product_seq_len':str(len(product_seq)),
                       'left_primer_seq':left_primer_seq,
                       'right_primer_seq':right_primer_seq,
                       'left_primer_tm':left_primer_tm,
                       'right_primer_tm':right_primer_tm,  
                      })

        primers_dict[key]=li
    return primers_dict

#输出设计失败的
def result_output_failtrue(primers_dict_failture,outputdir):
    df = pd.DataFrame()
    for k,v in primers_dict_failture.items():
        temp = pd.DataFrame([v])
        temp.insert(0,'id',k)
        df = df.append(temp)
    df.to_excel(outputdir)
        
#高通量生成引物  
def create_primers(dict_temp_template,type='left',mute_type='single'):
    dict_temp_primer_success={}
    dict_temp_primer_failure={}
    for key in dict_temp_template:
        template_seq = dict_temp_template[key]
       
        primers_dict=primer_design(key,template_seq,type,mute_type)

        if len(primers_dict)<10:
            #设计失败
            dict_temp_primer_failure[key] = primers_dict   
        else:   
            #设计成功
                           
            dict_temp_primer_success[key] = primers_dict

          
    return dict_temp_primer_success,dict_temp_primer_failure

#判断设计引物类型         
def judge_primer_type_by_input(mute_tuple_df, before_processed_seq_dict):
    single_mute_df,double_mute_df = mute_tuple_df
    if len(single_mute_df) == 0:
        mute_type = 'double'
    elif len(double_mute_df) == 0:
        mute_type = 'single'
    else:
        mute_type = 'single_double'
    print(mute_type)   
    if mute_type == 'double' or mute_type == 'single_double':   
        double1, double2, double3 = double_mute_class(double_mute_df, before_processed_seq_dict)
        double = double1, double2, double3
    else:
        double = 'NONE'   
    
    return single_mute_df,double_mute_df,double,mute_type


#判断double的类型，并按照双点距离之间碱基的个数进行分类
def double_mute_class(double,before_processed_seq_dict):
    def work(x1,x2):
        if x2 > x1:
            distance_seq = before_processed_seq_dict['target_gene_seq'][int(x1)+2:int(x2)-1]
        else:
            distance_seq = before_processed_seq_dict['target_gene_seq'][int(x2)+2:int(x1)-1]
        return distance_seq
   
    double['distance_seq'] = double.apply(lambda x: work(x['position1'],x['position2']),axis=1)
    double['distance'] =  abs(double['position2'].astype(int) - double['position1'].astype(int))-3
    double1 = double[(0<=double['distance']) & (double['distance']<=21)]
    double2 = double[(21<=double['distance']) & (double['distance']<=120)]
    double3 = double[double['distance']>120]
    return double1, double2, double3

#流程
def design_process(input_mute_tuple_df,dict_plasmid_seq,before_processed_seq_dict):  
    #判断设计类型
    single_mute_df,double_mute_df,double,mute_type = judge_primer_type_by_input(input_mute_tuple_df, before_processed_seq_dict)
    if double != 'NONE':
        double1, double2, double3 = double
    
    result_output_success_failtrue_dict = {}

    single_dict_left_primers_failtrue = 'NONE'  
    single_dict_left_primers_attribute = 'NONE'
    single_dict_right_primers_failtrue = 'NONE'
    single_dict_right_primers_attribute = 'NONE'

    double_dict_left_primers_failtrue = 'NONE'
    double_dict_left_primers_attribute = 'NONE'
    double_dict_right_primers_failtrue = 'NONE'
    double_dict_right_primers_attribute = 'NONE'



    if mute_type == 'single' or mute_type == 'single_double':
        #设计左引物模板
        dict_left_primer_template = single_input_to_dict_temp(single_mute_df,dict_plasmid_seq)
        #生成左引物
        dict_left_primers = create_primers( dict_temp_template = dict_left_primer_template, type = 'left', mute_type='single' )
        #取出左引物失败的情况
        if len(dict_left_primers[1]) > 0:
            single_dict_left_primers_failtrue = dict_left_primers[1]     
        else:
            single_dict_left_primers_failtrue = 'NONE'
            
        #在左引物设计成功的前提下设计右引物   
        if len(dict_left_primers[0]) > 0:   
            dict_left_primers_success = dict_left_primers[0]    
            #取出左引物属性
            single_dict_left_primers_attribute = get_primer_attribute(dict_left_primers_success, dict_left_primer_template)      
            #
            
            #设计右引物模板
            dict_right_primer_template=single_input_to_dict_temp(   single_mute_df,    
                                                                    dict_plasmid_seq, 
                                                                    dict_left_primers_succ = dict_left_primers_success,
                                                                )
            #生成右引物
            dict_right_primers = create_primers( dict_temp_template = dict_right_primer_template,
                                                 type = 'right',
                                                 mute_type='single'   
                                               ) 
            #判断右引物是否设计成功
            if len(dict_right_primers[1]) > 0:
                single_dict_right_primers_failtrue = dict_right_primers[1]
            else:
                single_dict_right_primers_failtrue = 'NONE'
                
     
            if len(dict_right_primers[0]) > 0:   
                dict_right_primers_success = dict_right_primers[0]  
                #取出右引物属性
                single_dict_right_primers_attribute = get_primer_attribute(dict_right_primers_success, dict_right_primer_template)
            else:
                single_dict_right_primers_attribute = 'NONE'
        else:
            single_dict_left_primers_attribute = 'NONE'

    if mute_type == 'double' or mute_type == 'single_double':

        double1_dict_left_primers_failtrue = 'NONE'
        double1_dict_right_primers_failtrue = 'NONE'
        double2_dict_left_primers_failtrue = 'NONE'
        double2_dict_right_primers_failtrue = 'NONE'
        double3_dict_left_primers_failtrue = 'NONE'
        double3_dict_right_primers_failtrue = 'NONE'

        double1_dict_left_primers_attribute = 'NONE'
        double1_dict_right_primers_attribute = 'NONE'
        double2_dict_left_primers_attribute = 'NONE'
        double2_dict_right_primers_attribute = 'NONE'
        double3_dict_left_primers_attribute = 'NONE'
        double3_dict_right_primers_attribute = 'NONE'




        if len(double1)>0:
            #设计上游引物模板
            dict_left_primer_template = double_input_to_dict_temp(double1,dict_plasmid_seq)  
            #生成引物
            dict_left_primers = create_primers(dict_left_primer_template,type='left',mute_type='double1_2')


            #找出设计失败的，设计成功的
            if len(dict_left_primers[1]) > 0:
                double1_dict_left_primers_failtrue = dict_left_primers[1]
                #输出
            else:
                double1_dict_left_primers_failtrue = "NONE"
            if len(dict_left_primers[0]) > 0:
                dict_left_primers_success = dict_left_primers[0]
                  
                #取出左(target)引物属性 (上游引物)
                double1_dict_left_primers_attribute = get_primer_attribute( dict_left_primers_success, dict_left_primer_template)
                #输出
                #在上游引物设计成功的基础上设计下游引物模板
                
                dict_right_primer_template = double_input_to_dict_temp( double1, dict_plasmid_seq, dict_left_primers_succ = dict_left_primers_success)
                #生成右(plasmid)引物
                dict_right_primers = create_primers( dict_right_primer_template, type = 'right', mute_type='double1_2')
                #判断
                if len(dict_right_primers[1])>0:
                    double1_dict_right_primers_failtrue = dict_right_primers[1]
                #输出
                else:
                    double1_dict_right_primers_failtrue = 'NONE'
                if len(dict_right_primers[0]) > 0:  
                    dict_right_primers_success = dict_right_primers[0]
                    #取出右(plasmid)引物属性
                    double1_dict_right_primers_attribute = get_primer_attribute(dict_right_primers_success,dict_right_primer_template)
                   
                    #输出   
                else:
                    double1_dict_right_primers_attribute = 'NONE'
            else:     
                double1_dict_left_primers_attribute = 'NONE'
            
           
        if len(double2)>0:
            #设计上游引物模板
            dict_left_primer_template = double_input_to_dict_temp(double2,dict_plasmid_seq)
            #生成引物
            dict_left_primers = create_primers(dict_left_primer_template,type='left',mute_type='double1_2')
            #找出设计失败的，设计成功的
            if len(dict_left_primers[1]) > 0:
                double2_dict_left_primers_failtrue = dict_left_primers[1]
                #输出
            else:
                double2_dict_left_primers_failtrue = "NONE"
            if len(dict_left_primers[0]) > 0:
                dict_left_primers_success = dict_left_primers[0]
                #取出左(target)引物属性 (上游引物)
                double2_dict_left_primers_attribute = get_primer_attribute( dict_left_primers_success, dict_left_primer_template,double2,stype='up')
                #输出
                #在上游引物设计成功的基础上设计下游引物模板
                dict_right_primer_template = double_input_to_dict_temp( double2, dict_plasmid_seq, dict_left_primers_succ = dict_left_primers_success)
                #生成右(plasmid)引物
                dict_right_primers = create_primers( dict_right_primer_template, type = 'right', mute_type='double1_2')
                #判断
                if len(dict_right_primers[1])>0:
                    double2_dict_right_primers_failtrue = dict_right_primers[1]
                #输出
                else:
                    double2_dict_right_primers_failtrue = 'NONE'
                if len(dict_right_primers[0]) > 0:
                    dict_right_primers_success = dict_right_primers[0]
                    #取出右(plasmid)引物属性
                    double2_dict_right_primers_attribute = get_primer_attribute(dict_right_primers_success,dict_right_primer_template,double2,stype='down')
                    #输出 
                else:
                    double2_dict_right_primers_attribute = 'NONE'
            else:
                double2_dict_left_primers_attribute = 'NONE'


        if len(double3)>0:
        
            #设计左(target)引物模板
            dict_left_primer_template = double_input_to_dict_temp(double3,dict_plasmid_seq,stype='plasmid')
            #生成左(target)引物
            dict_left_primers = create_primers(dict_temp_template = dict_left_primer_template, type = 'plasmid', mute_type='double3')
            #找出设计失败的，设计成功的
            if len(dict_left_primers[1]) > 0:
                double3_dict_left_primers_failtrue = dict_left_primers[1]
                #输出
            else:
                double3_dict_left_primers_failtrue = "NONE"
            if len(dict_left_primers[0]) > 0:
                dict_left_primers_success = dict_left_primers[0]
                #取出左(target)引物属性
                double3_dict_left_primers_attribute = get_primer_attribute( dict_left_primers_success, dict_left_primer_template)
                #输出 
            else:
                double3_dict_left_primers_attribute = 'NONE'
            #设计右(plasmid)引物模板
            dict_right_primer_template=double_input_to_dict_temp( double3, dict_plasmid_seq, stype='target')
            #生成右(plasmid)引物
            dict_right_primers = create_primers( dict_temp_template = dict_right_primer_template, type = 'target', mute_type='double3')
            if len(dict_right_primers[1])>0:
                double3_dict_right_primers_failtrue = dict_right_primers[1]
                #输出
            else:
                double3_dict_right_primers_failtrue = 'NONE'
            if len(dict_right_primers[0]) > 0:
                dict_right_primers_success = dict_right_primers[0]
                #取出右(plasmid)引物属性
                double3_dict_right_primers_attribute = get_primer_attribute(dict_right_primers_success,dict_right_primer_template)
                #输出 
            else:
                double3_dict_right_primers_attribute = 'NONE'

        #输出双点突变引物设计结果
        #汇总double
        #失败
        if double1_dict_left_primers_failtrue=='NONE' and double2_dict_left_primers_failtrue == 'NONE' and double3_dict_left_primers_failtrue == 'NONE':
            double_dict_left_primers_failtrue = 'NONE'
        else:
            double_dict_left_primers_failtrue={}
            for i in [double1_dict_left_primers_failtrue,double2_dict_left_primers_failtrue,double3_dict_left_primers_failtrue]:
                if i != 'NONE':
                    double_dict_left_primers_failtrue.update(i)

        if double1_dict_right_primers_failtrue=='NONE' and double2_dict_right_primers_failtrue == 'NONE' and double3_dict_right_primers_failtrue == 'NONE':
            double_dict_right_primers_failtrue = 'NONE'
        else:
            double_dict_right_primers_failtrue={}
            for i in [double1_dict_right_primers_failtrue,double2_dict_right_primers_failtrue,double3_dict_right_primers_failtrue]:
                if i != 'NONE':
                    double_dict_right_primers_failtrue.update(i)


        #成功        
        if double1_dict_left_primers_attribute=='NONE' and double2_dict_left_primers_attribute == 'NONE' and double3_dict_left_primers_attribute == 'NONE':
            double_dict_left_primers_attribute = 'NONE'
        else:
            double_dict_left_primers_attribute={}
            for i in [double1_dict_left_primers_attribute,double2_dict_left_primers_attribute,double3_dict_left_primers_attribute]:
                if i != 'NONE':
                    double_dict_left_primers_attribute.update(i)

        if double1_dict_right_primers_attribute=='NONE' and double2_dict_right_primers_attribute == 'NONE' and double3_dict_right_primers_attribute == 'NONE':
            double_dict_right_primers_attribute = 'NONE'
        else:
            double_dict_right_primers_attribute={}
            for i in [double1_dict_right_primers_attribute,double2_dict_right_primers_attribute,double3_dict_right_primers_attribute]:
                if i != 'NONE':
                    double_dict_right_primers_attribute.update(i)


    #输出单点突变引物设计结果           
    result_output_success_failtrue_dict['single_dict_up_primers_failtrue'] = single_dict_left_primers_failtrue
    result_output_success_failtrue_dict['single_dict_up_primers_success'] = {'dict_primers_attribute':single_dict_left_primers_attribute,'primer_num_return':2,'type':'up'}
    result_output_success_failtrue_dict['single_dict_down_primers_failtrue'] = single_dict_right_primers_failtrue
    result_output_success_failtrue_dict['single_dict_down_primers_success'] = {'dict_primers_attribute':single_dict_right_primers_attribute,'primer_num_return':2,'type':'down'}     

    result_output_success_failtrue_dict['double_dict_up_primers_failtrue'] = double_dict_left_primers_failtrue
    result_output_success_failtrue_dict['double_dict_up_primers_success']={'dict_primers_attribute':double_dict_left_primers_attribute,'primer_num_return':2,'type':'up'}
    result_output_success_failtrue_dict['double_dict_down_primers_failtrue'] = double_dict_right_primers_failtrue
    result_output_success_failtrue_dict['double_dict_down_primers_success'] = {'dict_primers_attribute':double_dict_right_primers_attribute,'primer_num_return':2,'type':'down'}

       
            
    return result_output_success_failtrue_dict   

      

          

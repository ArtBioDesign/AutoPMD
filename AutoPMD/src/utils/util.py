# import utils.config as config  
import pandas as pd  
from Bio import SeqIO 
import primer3
import sys,os
sys.path.append('../')
import json   
import module.organize_plasmid_primers as opp
import shutil   
from utils.config import config    
import zipfile   
from zipfile import ZipFile
import os
import zipfile
import module.from_gbFile_to_seq as gb_seq
import itertools


#将字典中的某个key的值替换成另一个值
def replace_digits_in_keys(dictionary, replacement_char):
    new_dict = {}
    for key in dictionary:
        new_key = ''.join(replacement_char if char.isdigit() else char for char in key)
        new_dict[new_key] = dictionary[key]
    return new_dict



#取互补链
def revComp(seq):
    complementSeq=seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))
    revcompSeq = complementSeq[::-1]
    return revcompSeq

#重新排列df
def reorgnize(df,columns =[]):
    df = df[columns]  
    return df   
  
#换名字
def replace_primer3Name_to_peopleReadName(df,type='up'):
    names = df.columns.to_list()
    print(names)      
    df =df.rename(columns={
                            names[2]:f"primer_{type}_f_seq_(5'-3')",
                            names[3]:f"primer_{type}_r_seq_(5'-3')",
                            names[4]:f"primer_{type}_f_Tm",
                            names[5]:f"primer_{type}_r_Tm",
                            names[1]:f"{type}_product_sequence_length",
                            names[0]:f"{type}_product_sequence"
                        } )
                        
    return df   
    
#坐标转序列  
def coordinate_2_seq(coordinate,seq):
    return str(seq[coordinate[0]:coordinate[1]])
  
#以坐标轴为1的坐标，转换为坐标轴为0的坐标  
def coordinate_2_coordinate(start,length=3):
    return int(start)-1,(int(start)-1)+length

#提取多个plasmid template file
def extract_many_gb(input_file_path = config.INPUT_FILE_PATH + config.INPUT_FILE_NAME):
   
    input_file_path = config.INPUT_FILE_PATH + config.INPUT_FILE_NAME
  
    #提取input_excel信息
    input_df = pd.read_excel(input_file_path)
    input_df.columns = input_df.columns.str.lower()
   
  
    #分离gb文件的相关信息
      
    gb_files_name = set(input_df['template'].values.tolist())
   
    gb_dict = {}

    for gb_name in gb_files_name: 

        gb_df = input_df[input_df['template']==gb_name]
      
        marker = list(gb_df['marker'].values)[0]
        target = list(gb_df['target'].values)[0]
        gb_path = config.INPUT_FILE_PATH  
        gb_dict[gb_path+gb_name]={  'marker':marker,  
                            'target':target,
                            'gb_name':gb_name,
                            'input_mute_df':gb_df
                        }
    return gb_dict    
  

def convert_nucleotide_to_aminoacid(df):
    def work(mutagenesis,material_type):
 
        if ';' in mutagenesis:
            mute1,mute2 = mutagenesis.split(';')
            mute1_arr1,mute1_arr2,mute1_arr3  =  mute1.split('-')
            mute2_arr1,mute2_arr2,mute2_arr3  =  mute2.split('-')

            if material_type == 'amino acid':
                mute1_arr2 = int(mute1_arr2) * 3 - 3 +1
                mute2_arr2 = int(mute2_arr2) * 3 - 3 +1

                if mute1_arr2 <= mute2_arr2:
                    mutagenesis = mute1_arr1 + '-' + str(mute1_arr2) + '-'+ str(mute1_arr3) + ';' + mute2_arr1 + '-' + str(mute2_arr2) + '-'+ mute2_arr3    
                else:
                    mutagenesis = mute2_arr1 + '-' + str(mute2_arr2) + '-'+ mute2_arr3 + ';' + mute1_arr1 + '-' + str(mute1_arr2) + '-'+ str(mute1_arr3)
            return mutagenesis



        else:
            arr1 = mutagenesis.split('-')[0]
            arr2 = int(mutagenesis.split('-')[1])
            arr3 = mutagenesis.split('-')[2]
            if material_type == 'amino acid':
                arr2 = arr2 * 3 - 3+ 1
                mutagenesis = arr1 + '-' + str(arr2) + '-'+arr3
            return mutagenesis
        
    df = lambda2cols(df,lambdaf=work, in_coln=['mutagenesis', 'type'], to_colns=['Mutagenesis'])
    df.drop(columns='mutagenesis',inplace=True)
    df.rename(columns={'Mutagenesis':'mutagenesis'},inplace=True)
    return df



#读取输入文件 
def read_input_file_by_gb(input_file_path=config.INPUT_FILE_PATH + config.INPUT_FILE_NAME, gb_name = 'icd-28a-new'):

    '''
        根据突变种类分离输入文件中每一个质粒模板的GOI
        params: 
            input_file_path
            gb_name
        returns:
            single_mute_df, double_mute_df    
    
    '''
      
    input_file_path = config.INPUT_FILE_PATH + config.INPUT_FILE_NAME

    input_df = pd.read_excel(input_file_path)  
    input_df.columns = input_df.columns.str.lower()

    #根据名字取出gb文件的内容
    input_df = input_df[input_df['template']==gb_name]

    #分离
    double_mute_df = input_df[input_df['mutagenesis'].str.contains(';')]
    single_mute_df = input_df[~input_df['mutagenesis'].str.contains(';')]
     

    if len(single_mute_df) > 0 :
            single_mute_df = convert_nucleotide_to_aminoacid(single_mute_df)
            single_mute_df['mutagenesis']= single_mute_df.mutagenesis.str.split('-')
            single_mute_df = single_mute_df[['sample','marker','target','mutagenesis']]
            new_single_df = pd.DataFrame(single_mute_df.mutagenesis.to_list(),columns=['before1','position1','after1'])
            #重置索引
            
            single_mute_df.reset_index(drop=True,inplace=True)
            new_single_df.reset_index(drop=True,inplace=True)
            single_mute_df = pd.concat([single_mute_df,new_single_df],axis=1)   
            print(single_mute_df)  
    if len(double_mute_df) > 0:
        def work(x):
                arr = []
                for item in x:   
                    arr.extend(item.split('-'))
                return arr
        double_mute_df = convert_nucleotide_to_aminoacid(double_mute_df)
        double_mute_df['mutagenesis']=double_mute_df.mutagenesis.str.split(';').apply(lambda x: work(x))
        double_mute_df = double_mute_df[['sample','marker','target','mutagenesis']]
        double_mute_df.reset_index(drop=True,inplace=True)
        new_double_df = pd.DataFrame(double_mute_df.mutagenesis.to_list(),columns=['before1','position1','after1','before2','position2','after2'])
        double_mute_df = pd.concat([double_mute_df,new_double_df],axis=1)
    
    return single_mute_df, double_mute_df
        


#针对多个gb在一个excel时，读取输入文件，提取多个gb文件
def read_input_file_to_many_gb(input_file_path=config.INPUT_FILE_PATH+config.INPUT_FILE_NAME):
    
    input_df = pd.read_excel(input_file_path)
    input_df.columns = input_df.columns.str.lower()
    #分离gb文件的相关信息
    gb_files_name = set(input_df['template'].values.tolist())
    gb_dict = {}
    for gb_name in gb_files_name:
        gb_df = input_df[input_df['template']==gb_name]
        marker = list(gb_df['marker'].values)[0]
        target = list(gb_df['target'].values)[0]
        gb_dict[gb_name]={'marker':marker,'target':target,'gb':gb_df}
    return gb_dict 


#设计引物
def primer_design(seqId,
                  seqTemplate,
                  stype,
                  mute_type='single',
                  global_args=config.GLOBAL_ARGS):

    if mute_type == 'single':
        config.GLOBAL_ARGS.update(config.S_GLOBAL_ARGS)
        global_args = config.GLOBAL_ARGS 
    elif mute_type =='double':
        config.GLOBAL_ARGS.update(config.D_GLOBAL_ARGS)
        global_args = config.GLOBAL_ARGS
    elif mute_type == 'sequencing': 
        config.Q_GLOBAL_ARGS.update(config.Q_ARGS)
        global_args =config.Q_GLOBAL_ARGS 

              
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
            size_range = [seqlength-66,seqlength-7]
        elif stype == "right":   #target下游引物设计
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [[-1,-1,seqlength-36,36]]
            seq_args['SEQUENCE_FORCE_LEFT_START'] = 0
            size_range = [seqlength-2,seqlength]
        primer_num = 2    
    elif mute_type == 'double':  #双点突变                                   
        if stype == "target":   #target引物设计
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [[0,35,-1,-1]]
        elif stype == "plasmid":   #plasmid引物设计
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [[-1,-1,seqlength-36,36]]     
        size_range = [seqlength,seqlength]
        primer_num = 2
    elif mute_type == 'sequencing':  #测序引物
        seq_args['SEQUENCE_ID'] = seqId
        seq_args['SEQUENCE_TEMPLATE'] = seqTemplate  
        size_range = [25,seqlength]
        primer_num = 1

    #设定全局参数   
    global_args['PRIMER_PRODUCT_SIZE_RANGE']=size_range
    global_args['PRIMER_NUM_RETURN']= primer_num
    #调用工具
    primer3_result = primer3.bindings.designPrimers(seq_args, global_args)    
    return primer3_result

#输出设计引物失败  
def result_output_failtrue_df(plasmid_name,primers_dict_failture):

    df = pd.DataFrame()
    for k,v in primers_dict_failture.items():
        temp = pd.DataFrame([v])
        temp.insert(0,'id',k)
        df = df.append(temp)
    df['template']=plasmid_name
    return df  


#输出引物设计成功的
def result_output_success_df(plasmid_name,primers_dict,type='down'):
    all_df = pd.DataFrame()
    
    for key in primers_dict:
        df =pd.DataFrame(primers_dict[key])    
        df['id'] = key   
        all_df=all_df.append(df)  
    #换名子
    all_df = replace_primer3Name_to_peopleReadName(all_df,type)  
    columns_name = all_df.columns.to_list()
    # #列的顺序
    all_df = all_df[[columns_name[-1],columns_name[2], columns_name[3],columns_name[4],columns_name[5],columns_name[0],columns_name[1]]]
    all_df['template']=plasmid_name
    return all_df  
  

#合并引物设计成功的情况
def merge_primer(primer_dict,gb_name,type='one'):
     
  
    up_primer_df = pd.DataFrame()   
    down_primer_df = pd.DataFrame()
    primer_df = pd.DataFrame()
    for k,v in primer_dict.items():
        if type == 'one':
            if 'up' in k:
                primer_df_everyone = reorgnize(v,columns=['id',"primer_up_f_seq_(5'-3')","primer_up_r_seq_(5'-3')"])
                up_primer_df = up_primer_df.append(primer_df_everyone)
            elif 'down' in k:
                primer_df_everyone = reorgnize(v,columns=['id',"primer_down_f_seq_(5'-3')","primer_down_r_seq_(5'-3')"])
                down_primer_df = down_primer_df.append(primer_df_everyone)
        else:
            if 'up' in k:
                v.drop(columns='template',inplace=True)
                up_primer_df = up_primer_df.append(v)
            elif 'down' in k:
                v.drop(columns='template',inplace=True)
                down_primer_df = down_primer_df.append(v)

    #合并 
    # print(up_primer_df.columns,down_primer_df.columns)
    if len(up_primer_df)>0 and len(down_primer_df)>0:
        primer_df = pd.merge(up_primer_df,down_primer_df,on='id',how='outer')     
    elif len(up_primer_df) == 0 and len(down_primer_df)>0:
        primer_df = down_primer_df
    elif len(up_primer_df) > 0 and len(down_primer_df)==0:
        primer_df = up_primer_df

              
     #填充
    primer_df = primer_df.fillna('null')
       
    primer_df['template'] = gb_name        
    
    plasmid_primer_df = primer_df[primer_df.index%2 == 0]   
    spare_plasmid_primer = primer_df[primer_df.index%2 == 1]

    #去重
    plasmid_primer_df.drop_duplicates(subset='id', inplace=True, ignore_index=False)  
    spare_plasmid_primer.drop_duplicates(subset='id', inplace=True, ignore_index=False)  
    return plasmid_primer_df,spare_plasmid_primer    
 

def one_plasmid_primer_result_to_df(gb_name,result_output_success_failtrue_dict):
        plasmid_name = gb_name      

        success_result_dict = {}
        failtrue_result_dict = {}     
   

        for k,v in result_output_success_failtrue_dict.items():
           
            
            if 'failtrue' in k and v != 'NONE':
                temp= result_output_failtrue_df(plasmid_name,v)
                failtrue_result_dict[k] = temp   
            elif 'failtrue' in k and v == 'NONE':
                failtrue_result_dict[k] = pd.DataFrame()
            elif 'success' in k and v['dict_primers_attribute'] != 'NONE':
                
                temp = result_output_success_df(plasmid_name, v['dict_primers_attribute'], v['type'])
                success_result_dict[k] = temp
            elif 'success' in k and v['dict_primers_attribute'] != 'NONE':
                success_result_dict[k] = pd.DataFrame()
       

        return success_result_dict, failtrue_result_dict       


# #打包文件为zip
def zip_ya(startdir, file_news,num=1):
    z = zipfile.ZipFile(file_news, 'a', zipfile.ZIP_DEFLATED)
    for dirpath, dirnames, filenames in os.walk(startdir):
        fpath = dirpath.replace(startdir, '')  # 这一句很重要，不replace的话，就从根目录开始复制
        fpath = fpath and fpath + os.sep or ''
        for filename in filenames:  
            if num == 1:   
                if filename == 'failtrue.xlsx':
                    z.write(os.path.join(dirpath, filename), fpath + filename,zipfile.ZIP_DEFLATED)
            elif num == 2:
                if filename =='success.xlsx':
                    z.write(os.path.join(dirpath, filename), fpath + filename,zipfile.ZIP_DEFLATED)
            else:
                fpath=startdir.split('/')[-2]+'/'
                z.write(os.path.join(dirpath, filename), fpath + filename,zipfile.ZIP_DEFLATED)

# #汇总信息
def read_failtrue_df_and_merge(sum_info_file='',failtrue_file ='',outputdir = config.OUTPUT_FILE_PATH, failture_position_df = pd.DataFrame()):

    #读取合并表格与失败表格
    plasmid_primer = pd.read_excel(outputdir+sum_info_file,sheet_name='plasmid_primer')
    spare_plasmid_primer = pd.read_excel(outputdir+sum_info_file,sheet_name='spare_plasmid_primer')
    sequencing_primer = pd.read_excel(outputdir+sum_info_file,sheet_name='sequencing_primer')
    sequencing_fail_primer_df = pd.read_excel(outputdir+sum_info_file,sheet_name='fail_sequencing_primer')

    xl = pd.ExcelFile(outputdir + failtrue_file)
    sheet_names = xl.sheet_names   
    
    #合并失败引物
    u_df = pd.DataFrame(columns=['Sample','Template','Pimer_U_F_Explain','Primer_U_R_Explain','Primer_U_Pair_Explain'])
    d_df = pd.DataFrame(columns=['Sample','Template','Primer_D_F_Explain','Primer_D_R_Explain','Primer_D_Pair_Explain'])
    for name in sheet_names:
        if 'u_fail' in name and len(xl.parse(name)!=0):
            df1 = xl.parse(name)[['Sample','Template','Pimer_U_F_Explain','Primer_U_R_Explain','Primer_U_Pair_Explain']]
            u_df = u_df.append(df1)
        elif 'd_fail' in name and len(xl.parse(name)!=0):    
            df2 = xl.parse(name)[['Sample','Template','Primer_D_F_Explain','Primer_D_R_Explain','Primer_D_Pair_Explain']]
            d_df = d_df.append(df2)
    
    u_d_df = pd.merge(u_df,d_df,on=['Sample','Template'],how='outer')


   
    #失败与所有引物合并
    all_df=pd.merge(plasmid_primer,u_d_df,how='outer',on=['Sample','Template']) 
    all_df = all_df.fillna('null')
    #标记出哪些是失败的
    def work(x1,x2,x3,x4,x5,x6):
        if x1=='null' and x4=='null':
            return 'NO'   
        elif x1=='null' and x4!='null':
            return 'D_YES'
        elif x1!='null' and x4 == 'null':  
            return 'U_YES'
        else:
            return 'YES'  
    all_df['Warning'] = all_df.apply(lambda x:work(x['Pimer_U_F_Explain'],x['Primer_U_R_Explain'],x['Primer_U_Pair_Explain'],x['Primer_D_F_Explain'],x['Primer_D_R_Explain'],x['Primer_D_Pair_Explain']),axis=1)

    #新的汇总表----success.xlsx
    plasmid_primer = all_df[['No.', 'Sample', 'Template', 'Mutation', 'Design',
                        "Primer-U-F_Sequence（5'-3')", "Primer-U-R_Sequence（5'-3')",
                       'Product_Length-U', "Primer-D-F_Sequence（5'-3')",
                        "Primer-D-R_Sequence（5'-3')", 'Product_Length-D', 'Warning']]
     #处理备用引物  
    temp = plasmid_primer[['Sample','Template','Warning']] 
   
    spare_plasmid_primer = pd.merge(spare_plasmid_primer,temp,on=['Sample','Template'])    

    #新汇总表------failtrue.xlsx
    fail_u_d_df = pd.merge(plasmid_primer,u_d_df,how='inner')   
    
    #增加target_gene的先后nbp
    add_seq_to_success_df = add_seq_to_success(n=config.TARGETGENE_AFTER_BEFORE_SEQ_N, input_file_path = config.INPUT_FILE_PATH + config.INPUT_FILE_NAME)
    plasmid_primer = plasmid_primer.merge(add_seq_to_success_df,left_on=['Sample','Template'],right_on=['sample','template'],how='inner')
    plasmid_primer = plasmid_primer.drop(columns=['template','sample','mutagenesis'])
    spare_plasmid_primer = spare_plasmid_primer.merge(add_seq_to_success_df,left_on=['Sample','Template'],right_on=['sample','template'],how='inner')
    spare_plasmid_primer = spare_plasmid_primer.drop(columns=['template','sample','mutagenesis'])
    sequencing_alignment_sequence_df = plasmid_primer[['No.','Sample','targetGene_after_before_seq_n']]



    plasmid_primer['No.'] = plasmid_primer['No.'] + 1
    spare_plasmid_primer['No.'] = spare_plasmid_primer['No.'] + 1
    sequencing_primer['No.'] = sequencing_primer['No.'] + 1
    sequencing_alignment_sequence_df['No.'] = sequencing_alignment_sequence_df['No.']  + 1



    #输出
    with pd.ExcelWriter(outputdir+'success.xlsx') as writer:
        plasmid_primer.to_excel(writer,sheet_name = 'plasmid_primer', index=False)
        spare_plasmid_primer.to_excel(writer,sheet_name = 'spare_plasmid_primer', index=False)
        sequencing_primer.to_excel(writer,sheet_name = 'sequencing_primer', index=False)
        sequencing_alignment_sequence_df.to_excel(writer,sheet_name='sequencing_alignment_sequence', index=False) 

    
    fail_u_d_df = fail_u_d_df.append(failture_position_df)

    fail_u_d_df['No.'] = range(1,len(fail_u_d_df)+1)
    sequencing_fail_primer_df['No.'] = sequencing_fail_primer_df['No.'] +1
    with pd.ExcelWriter(outputdir+'failtrue.xlsx') as writer:
        fail_u_d_df.to_excel(writer,sheet_name = 'plasmid_primer', index=False) 
        sequencing_fail_primer_df.to_excel(writer,sheet_name='sequencing_primer',index=False)


    
    #寻找错误的gb文件，且输出到文件夹中
    #获取所有设计失败的订单

    #获取所有设计失败的gb名字
    gb_name_set = set(u_d_df['Template'])
    #遍历每个gb
    for name in list(gb_name_set):
        #在所有文件中取一个gb文件下的所有site
        all_fail_set = set(u_d_df[u_d_df['Template']==name]['Sample']+'_plasmid_mutation.gb')
        gb_name = name + '_site_mute'
        if os.path.exists(outputdir+gb_name):
        
            #读取文件列表
            all_file_list=os.listdir(outputdir+gb_name)
            if len(set(all_file_list)) >= len(all_fail_set):
                
                for gb_fail_name in all_fail_set:
                
                    #源dir名称
                    s_name = gb_name
                    #新dir
                    d_name = 'failtrue_'+gb_name
                    if os.path.exists(outputdir+d_name+'/'):
                        shutil.copy(outputdir+s_name+'/'+gb_fail_name,outputdir+d_name+'/')
                        
                    else:  
                        os.mkdir(outputdir+d_name+'/')   
                        shutil.copy(outputdir+s_name+'/'+gb_fail_name,outputdir+d_name+'/')    

    #打包输出
    #处理失败的zip
    for template in set(fail_u_d_df['Template']):
        gb_name ='failtrue_'+template+'_site_mute'
        zip_ya(outputdir+gb_name+'/', outputdir+'failure_site_mute.zip',num=0)

    zip_ya(outputdir, outputdir+'failure_site_mute.zip',num=1)

    #处理成功的zip
    #处理失败的zip  
    for template in set(all_df['Template']):
        gb_name =template+'_site_mute'
        zip_ya(outputdir+gb_name+'/', outputdir+'success_site_mute.zip',num=0)

    zip_ya(outputdir, outputdir+'success_site_mute.zip',num=2)


def add_seq_to_success(n=config.TARGETGENE_AFTER_BEFORE_SEQ_N, input_file_path = config.INPUT_FILE_PATH + config.INPUT_FILE_NAME):


    many_gb_dict = extract_many_gb(input_file_path)
    add_seq_to_success_df =  pd.DataFrame() 
    for k,v in many_gb_dict.items():
        before_processed_seq_dict, after_processed_seq_dict =gb_seq.get_data_from_genebank(      
                                                                infile = f'{k}.gb',
                                                                marker=v['marker'],
                                                                target_gene=v['target']) 
        target_gene_seq = before_processed_seq_dict['target_gene_seq']
        target_gene_up_seq =  before_processed_seq_dict['target_gene_up_seq']
        target_gene_down_seq = before_processed_seq_dict['target_gene_down_seq']   
        if len(target_gene_up_seq) > n : 
            target_gene_up_seq_n = target_gene_up_seq[-n:]
        else:
            target_gene_up_seq_n = target_gene_up_seq
            target_gene_up_seq_n = target_gene_down_seq[-(n-len(target_gene_up_seq)):] + target_gene_up_seq_n
        if len(target_gene_down_seq) > n:
            target_gene_down_seq_n = target_gene_down_seq[:n] 
        else:
            target_gene_down_seq_n = target_gene_down_seq
            target_gene_down_seq_n = target_gene_down_seq_n + target_gene_up_seq[:(n-len(target_gene_down_seq))]
        input_mute_df = v['input_mute_df'][['sample','template','mutagenesis']]
          
        def work(x):   
            gene_seq = '' 
            if ';' in x:
                arr = x.split(';')
                for item in arr:  
                    ar = item.split('-')
                    gene_seq = target_gene_seq[:int(ar[1])-1] + ar[2] + target_gene_seq[int(ar[1])+3-1:]
                return target_gene_up_seq_n + gene_seq +  target_gene_down_seq_n
            else:
                ar = x.split('-')
                gene_seq = target_gene_seq[:int(ar[1])-1] + ar[2] + target_gene_seq[int(ar[1])+3-1:]
                return target_gene_up_seq_n + gene_seq +  target_gene_down_seq_n
              
        input_mute_df['targetGene_after_before_seq_n'] = input_mute_df.mutagenesis.apply(lambda x: work(x))
        add_seq_to_success_df = add_seq_to_success_df.append(input_mute_df)  
    return add_seq_to_success_df
    

def lambda2cols(df,lambdaf,in_coln,to_colns):         #apply函数的助手

    if len(in_coln) == 2:
        df_=df.apply(lambda x: lambdaf(x[in_coln[0]], x[in_coln[1]]),
                     axis=1).apply(pd.Series)
    if len(in_coln) == 3:
        df_=df.apply(lambda x: lambdaf(x[in_coln[0]],x[in_coln[1]],x[in_coln[2]]),
                     axis=1).apply(pd.Series)
    elif len(in_coln) == 4:
        df_=df.apply(lambda x: lambdaf(x[in_coln[0]],x[in_coln[1]],x[in_coln[2]],x[in_coln[3]]),
                     axis=1).apply(pd.Series)
    elif len(in_coln) == 5:
         df_=df.apply(lambda x: lambdaf(x[in_coln[0]],x[in_coln[1]],x[in_coln[2]],x[in_coln[3]],x[in_coln[4]]),
                     axis=1).apply(pd.Series)
    df_.columns = to_colns
    df = df.join(df_)        
    return df


  



    



  




 


    
    
    

   
    
    







    



  


  
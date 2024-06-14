from asyncore import write
import sys,os
from symbol import single_input
sys.path.append('../')
from utils.config import config   
from utils import util  
import pandas as pd  
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import math     
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC

#重新组合突变质粒片段  
def concat_mute_plasmid_seq(input_mute_df,before_processed_seq_dict):

    double_mute_df =pd.DataFrame(columns=['sample', 'marker', 'target', 'mutagenesis', 'before1', 'position1','after1', 'before2', 'position2', 'after2'])
    single_mute_df = pd.DataFrame(columns=['sample', 'marker', 'target', 'mutagenesis', 'before', 'position','after'])
    mute_before_target_gene_seq = before_processed_seq_dict['target_gene_seq']

    #区分单双点突变
    single_mute_df = input_mute_df[0]
    double_mute_df = input_mute_df[1]   
        #gene_seq坐标
    mute_before_target_gene_seq_start = 0
    mute_before_target_gene_seq_end = len(mute_before_target_gene_seq)   
           
    #合并单点、双点突变

    
    if len(single_mute_df) > 0:
        mute1_df = single_mute_df[['before1', 'position1', 'after1']]
        mute1_df = pd.concat([single_mute_df.iloc[:,0:1],mute1_df],axis=1)
    else:
        mute1_df = pd.DataFrame(columns=['before1', 'position1', 'after1'])

    
    if len(double_mute_df) > 0:
        mute2_df = pd.concat([double_mute_df.iloc[:,0:1], double_mute_df[['before1', 'position1', 'after1', 'before2', 'position2', 'after2']]] ,axis=1)
    else:
        mute2_df = pd.DataFrame(columns=['before1', 'position1', 'after1', 'before2', 'position2', 'after2'] )  


    print(mute1_df, '\n' ,mute2_df)
   
    mute_df = pd.concat([mute1_df,mute2_df])  
    mute_df.reset_index(drop=True,inplace=True)  


    mute_after_plasmid_seq_dict = {}
   
    
    #单、双点突变组装回原target_gene
    for index,row in mute_df.iterrows():
            #突变位点（密码子）在文件中的坐标
        
        position1 = util.coordinate_2_coordinate(row['position1'],length=3) #length:突变长度

        mute1_seq = row['after1']
       
        mute1_up_seq = util.coordinate_2_seq((mute_before_target_gene_seq_start,position1[0]),seq=mute_before_target_gene_seq) 
       
        if pd.isna(mute_df.iloc[index,4]):
            #只有一个突变位点
            mute1_down_seq = util.coordinate_2_seq((position1[1],mute_before_target_gene_seq_end),seq=mute_before_target_gene_seq)
            mute_after_target_gene_seq = mute1_up_seq + mute1_seq + mute1_down_seq
        else:
            #有两个突变位点
            position2 = util.coordinate_2_coordinate(row['position2'],length=3) #length:突变长度
            mute2_seq = row['after2']
            mute1_down_mute2_up_seq = util.coordinate_2_seq( (position1[1],position2[0]), seq = mute_before_target_gene_seq)
            mute2_down_seq = util.coordinate_2_seq((position2[1],mute_before_target_gene_seq_end), seq = mute_before_target_gene_seq)
            #拼接序列
            mute_after_target_gene_seq = mute1_up_seq + mute1_seq + mute1_down_mute2_up_seq + mute2_seq + mute2_down_seq

        #拼接突变后的序列   
        target_gene_up_seq = before_processed_seq_dict['target_gene_up_seq']
        target_gene_down_seq = before_processed_seq_dict['target_gene_down_seq']
        mute_after_plasmid_seq = target_gene_up_seq + mute_after_target_gene_seq + target_gene_down_seq
        temp_dict = {
                                        'target_gene_up_seq':target_gene_up_seq, 
                                        'mute_after_target_gene_seq':mute_after_target_gene_seq,   
                                        'target_gene_down_seq':target_gene_down_seq,
                                        'mute_after_plasmid_seq':mute_after_plasmid_seq       
                                    }
                                    

        mute_after_plasmid_seq_dict[row['sample']]= temp_dict      

    return mute_after_plasmid_seq_dict
          
#读引物输出文件
def read_primer_by_site_from_dir(
                              s_d_site,
                              sheet_name='spare_plasmid_primer',
                              dir_path=config.OUTPUT_FILE_PATH,
                              fileName=[config.OUTPUT_SINGLE_UP_NAME,  
                                        config.OUTPUT_SINGLE_DOWN_NAME,  
                                        config.OUTPUT_DOUBLE_DOWN_NAME,
                                        config.OUTPUT_DOUBLE_UP_NAME
                                       ]):  
    files= os.listdir(dir_path) #得到文件夹下的所有文件名称
    #挑选文件       input_from_primer_output=[ i for i in files if 'single'in i or 'double' in i and 'failtrue' not in i]
    up_primer_df = pd.DataFrame()
    down_primer_df = pd.DataFrame()
    primer_df = pd.DataFrame()
    for name in fileName:
        if name in files:
            primer_df_everyone = pd.read_excel(config.OUTPUT_FILE_PATH+name,sheet_name=sheet_name)
            #取每个gb文件对应的site
            if len(s_d_site) > 0:  
                s_d_site_df = pd.DataFrame(data=list(s_d_site),columns=['sample'])
                primer_df_everyone = primer_df_everyone.merge(s_d_site_df,how='inner',left_on='id',right_on='sample')     
                if 'up' in name:   
                    #重排位置
                    primer_df_everyone = util.reorgnize(primer_df_everyone,columns=['id',"primer_up_f_seq_(5'-3')","primer_up_r_seq_(5'-3')"])
                    up_primer_df = up_primer_df.append(primer_df_everyone)
                elif 'down' in name:
                    #重排位置
                    primer_df_everyone = util.reorgnize(primer_df_everyone,columns=['id',"primer_down_f_seq_(5'-3')","primer_down_r_seq_(5'-3')"])
                    down_primer_df = down_primer_df.append(primer_df_everyone)
            else:
                if 'up' in name:   
                    #重排位置
                    primer_df_everyone = util.reorgnize(primer_df_everyone,columns=['id',"primer_up_f_seq_(5'-3')","primer_up_r_seq_(5'-3')",'template','product_sequence_length'])
                    up_primer_df = up_primer_df.append(primer_df_everyone)
                elif 'down' in name:
                    #重排位置
                    primer_df_everyone = util.reorgnize(primer_df_everyone,columns=['id',"primer_down_f_seq_(5'-3')","primer_down_r_seq_(5'-3')",'template','product_sequence_length'])
                    down_primer_df = down_primer_df.append(primer_df_everyone)
            
    #合并
    primer_df = pd.merge(up_primer_df,down_primer_df,on='id',how='outer')     
    #填充
    primer_df = primer_df.fillna('null')
    return primer_df

def write_gb_file(primer_df,mute_after_plasmid_seq_dict,gb_path_name):
    output_gb_file_dict={}
     #写出GB文件

    columns =  primer_df.columns
       
    for i,row in primer_df.iterrows():
        #获取重组质粒序列
        mute_after_plasmid_seq = mute_after_plasmid_seq_dict[row['id']]['mute_after_plasmid_seq']
  
        #质粒文件地址
        plasmid_file_path = gb_path_name
        plasmid = SeqIO.read(plasmid_file_path, "genbank")
      
        features = []
        for ele in plasmid.features:
            if ele.type != "primer_bind":
                features.append(ele)
        plasmid.features = features
        plasmid.seq=Seq(mute_after_plasmid_seq, IUPAC.unambiguous_dna)  
        plasmid.id = f'site{i}'+config.OUTPUT_FILE_NAME_PLASMID_MUTATION
        plasmid.name = f'site{i}'+config.OUTPUT_FILE_NAME_PLASMID_MUTATION  

        temp_features = []
         
        if "primer_up_f_seq_(5'-3')" in columns and "primer_up_r_seq_(5'-3')" in columns and row[1]!='null' and row[2] != 'null':   
            #获取位置信息
            primer_up_f_seq_coordinate = (mute_after_plasmid_seq.find(row[1]), mute_after_plasmid_seq.find(row[1]) + len(row[1]) )
            primer_up_r_seq_coordinate = (mute_after_plasmid_seq.find( util.revComp(row[2]) ), mute_after_plasmid_seq.find( util.revComp(row[2])) + len(row[2]) )
            #构建特征
            primer1_features = SeqFeature(FeatureLocation(primer_up_f_seq_coordinate[0], primer_up_f_seq_coordinate[1]), type="primer_bind", strand=1, id='primer1', qualifiers={'label':['primer_1'] } )
            primer2_features = SeqFeature(FeatureLocation(primer_up_r_seq_coordinate[0], primer_up_r_seq_coordinate[1]), type="primer_bind", strand=-1, id = 'primer2', qualifiers={'label':['Primer_2']} )  
            temp_features.extend([primer1_features,primer2_features])

        if "primer_down_f_seq_(5'-3')" in columns and "primer_down_r_seq_(5'-3')" in columns and row[3]!='null' and row[4] != 'null':  
            primer_down_f_seq_coordinate = (mute_after_plasmid_seq.find(row[3]), mute_after_plasmid_seq.find(row[3]) + len(row[3]) )
            primer_down_r_seq_coordinate = (mute_after_plasmid_seq.find( util.revComp(row[4]) ), mute_after_plasmid_seq.find( util.revComp(row[4]) ) + len(row[4]) )  
            primer3_features = SeqFeature(FeatureLocation(primer_down_f_seq_coordinate[0], primer_down_f_seq_coordinate[1]), type="primer_bind", strand=1, id = 'primer3', qualifiers={'label':['Primer_3']} )
            primer4_features = SeqFeature(FeatureLocation(primer_down_r_seq_coordinate[0], primer_down_r_seq_coordinate[1]), type="primer_bind", strand=-1, id= 'primer4', qualifiers={'label':['Primer_4']} ) 
            temp_features.extend([primer3_features,primer4_features])

        plasmid.features.extend(temp_features)   
        #存储gb文件
        output_gb_file_dict[row[0]]=plasmid

    return output_gb_file_dict  
         
#primer_df引物、genebank_db：原质粒
def merge_gb_output_file(primer_gb_file_dict, spare_primer_gb_file_dict):

    site_mute_plasmid = {}
    for k,v in spare_primer_gb_file_dict.items():
        primer_bind_features=[]
        for item in v.features:
            if item.type == 'primer_bind':
                primer_bind_features.append(item)
        #换名字
        for item in primer_bind_features:
            item.id='spare_'+item.id
            item.qualifiers={'label':[item.id]}
        #加入备用引物    
        primer_gb_file_dict[k].features.extend(primer_bind_features) 
        site_mute_plasmid[f'{k}_'+config.OUTPUT_FILE_NAME_PLASMID_MUTATION]= primer_gb_file_dict[k]
        # #输出gb文件   
        # SeqIO.write(primer_gb_file_dict[k],
        #             outputdir+f'{k}_'+config.OUTPUT_FILE_NAME_PLASMID_MUTATION,
        #             "genbank")
    return site_mute_plasmid
     
#整理读入的引物
def read_primer_from_primer_df(primer_df):
    #读取引物 
    f_list = list(primer_df.iloc[:,1:2].values.flatten())
    f_list.extend(list(primer_df.iloc[:,3:4].values.flatten()))
    f_df = pd.DataFrame(columns=["碱基序列(5'to3')"],data=f_list)
    f_df = f_df[f_df["碱基序列(5'to3')"]!='null'].drop_duplicates()
    f_df.insert(0,'板号','f')
    #读取引物
    r_list = list(primer_df.iloc[:,2:3].values.flatten())
    r_list.extend(list(primer_df.iloc[:,3:4].values.flatten()))
    r_df = pd.DataFrame(columns=["碱基序列(5'to3')"],data=r_list)
    r_df = r_df[r_df["碱基序列(5'to3')"]!='null'].drop_duplicates()
    r_df.insert(0,'板号','r')
    #合并
    primer_orders_df=pd.concat([f_df,r_df]).reset_index(drop=True)
    return primer_orders_df

#引物重命名
def order_primer_rename_by_ways(primer_orders_df,label=['A','B','C','D','E','F','G','H']):    
    primer_name_df = pd.DataFrame()
    cell = math.ceil(len(primer_orders_df) / len(label))
    cell_remainder = len(primer_orders_df) % len(label)
    big_cell = math.ceil( cell/12 )
    big_cell_remainder = cell % 12

    for i in range(1,big_cell+1):
        for j in range(1,13):
            primer_name_df =primer_name_df.append(pd.DataFrame(label)+str(j))
            
    #重置索引        
    primer_name_df.reset_index(drop=True,inplace=True)

    #减去余数    
    if big_cell_remainder !=0:
        primer_name_df=primer_name_df.iloc[:-(big_cell_remainder*12*8),:]
        for j in range(1,big_cell_remainder+1):
            primer_name_df =primer_name_df.append(pd.DataFrame(label)+str(j))

    if cell_remainder != 0:
        primer_name_df=primer_name_df.iloc[:-(len(label)-cell_remainder),:]
        
    name = primer_name_df[0].values.tolist()
#     primer_name_df.rename(columns={0:'引物名称'},inplace=True)
    return name
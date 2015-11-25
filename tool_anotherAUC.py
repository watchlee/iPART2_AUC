#!/usr/bin/env python
# coding=utf-8
# variable detial
#   SARA_TEST_DATA_PATH = SARA當初測試用的資料
#   TtoR_list = 存放取得的TtoR pdb資料
#   TtoR_document_path = 已經alignment好的所有TtoR資料所放的檔案錄記
#   pdbpath = 放置原始的pdb檔案
#   oneDseq_path = 放置原始檔案的1D序列
#   Complete_path = 完整TtoR檔案路徑
"""使用###進行註解的地方是用於debug使用"""
import subprocess
#######################################################################
#   當初d<=2的測試用的資料 是錯的但沒辦法還原 所以使用此方法代替
read_file = open('../processed_eachto418_d2','r')
bigfamily_compare = read_file.read()
read_file.close()
#######################################################################
class PDB_DATA:
    big_family = '' # 祖先
    family = ''     # d = 0的family
    pdb_name = ''
    def __init__(self):
        pass
    def set_data(self,big_family_data,family_data,pdb_name_data):
        self.big_family = big_family_data
        self.family = family_data
        self.pdb_name = pdb_name_data
    def get_family(self):
        return self.family
    def get_bigfamily(self):
        return self.big_family
    def get_pdb(self):
        return self.pdb_name
    def test(self):
        return self.big_family+' '+self.family+' '+self.pdb_name
#######################################################################
class Compare_pdb:
    RMSD= 0
    num_gap1 = 0
    num_gap2 = 0
    length_align =0
    def __init__(self,RMSD,num_gap1,num_gap2,length_align):
        self.RMSD = RMSD
        self.num_gap1 = num_gap1
        self.num_gap2 = num_gap2
        self.length_align = length_align
    def getRMSD(self):
        return self.RMSD
    def get_gap1(self):
        return self.num_gap1
    def get_gap2(self):
        return self.num_gap2
    def get_align(self):
        return self.length_align
#######################################################################
#    Profit計算 
class Profit_pdb:
    # index表示是目前比對的第幾筆結果
    index = 0
    reference = ''
    mobile = ''
    readalignment=''
    path = ''
    def __init__(self,seq1_name,seq2_name,path,index):
        self.reference = 'reference /home/watchlee/Research_Programming/RMSD/pdb/'+seq1_name+'.pdb'
        self.mobile = 'mobile /home/watchlee/Research_Programming/RMSD/pdb/'+seq2_name+'.pdb'
        self.readalignment = 'readalignment '+path+'/ori_ali_seq.pir'+str(index)
        self.index = index
        self.path = path
        ###print 'profitW_log'+str(self.index)
    def write_file(self):
        with open(self.path+'/pf_script'+str(self.index),'w') as file:
            file.write(self.reference+'\n')
            file.write(self.mobile+'\n')
            file.write(self.readalignment+'\n')
            file.write('atoms *\n')
            file.write('ignoremissing\n')
            file.write('fit\n')
            file.write('quit\n')
        subprocess.call('./profit2.5.3 < '+self.path+'/pf_script'+str(self.index)+' > '+self.path+'/profitW_log'+str(self.index),shell=True)
    def read_RMSD(self):
        try:
            #預設是先找proftW_log的RMSD若檔案找不到則跳到except或是下一個profit_log
            with open(self.path+'/profit_log'+str(self.index),'r') as file:
                for each_line in file:
                    if(each_line.find('RMS')!=-1):
                        return each_line.replace('RMS: ','')#.replace('\n','')
            #有可能找不到RMSD的值則使用profit_log
            with open(self.path+'/profitW_log'+str(self.index),'r') as file:
                for each_line in file:
                    if(each_line.find('RMS')!=-1):
                        return each_line.replace('RMS:','')
            return 'Error! '+self.path+'/profitW_log'+str(self.index)
        except:
            with open(self.path+'/profitW_log'+str(self.index),'r') as file:
                for each_line in file:
                    if(each_line.find('RMS')!=-1):
                        return each_line.replace('RMS:','')
            return 'Error! '+self.path+'/profit_log'+str(self.index)
#######################################################################
def read_file(SARA_TEST_DATA_PATH):
    f_list = []
    with open(SARA_TEST_DATA_PATH,'r') as file:
        for each_line in file:
            pdb = PDB_DATA()
            bigfamily = each_line.split("|")[0].split(":")[0]
            family = each_line.split("|")[0].split(":")[1]
            pdb_name = each_line.split("|")[0].split(":")[2]+'_'+each_line.split("|")[1]
            pdb.set_data(bigfamily,family,pdb_name)
            f_list.append(pdb)
    return f_list
#######################################################################
def Raw_FSCOR_Process(FSCOR_list,FSCOR_document_path,outfile):
    ###count = 0
    sort_list = []
    sortbig_list= []
    for index in range(len(FSCOR_list)):
        for inner_index in range(index+1,len(FSCOR_list)):
            if(inner_index!=index):
                context_list = []
                context_length = 0
                file_document_name = FSCOR_list[index].get_pdb()+'_to_'+FSCOR_list[inner_index].get_pdb()
                file_document_name2 = FSCOR_list[inner_index].get_pdb()+'_to_'+FSCOR_list[index].get_pdb()
                temp_max = 0
                try:
                    with open(FSCOR_document_path+file_document_name+'/semiG_result.php','r') as  file:
                        for each_line in file:
                            context_list.append(each_line)
                            context_length+=1
                    times = context_length/ 7
                    temp_max = float(context_list[6]) 
                    for loop in range(1,times):
                        if(temp_max < float(context_list[loop*7+6])):
                            temp_max = float(context_list[loop*7+6])
                except:
                    with open(FSCOR_document_path+file_document_name2+'/semiG_result.php','r') as file:
                        for each_line in file:
                            context_list.append(each_line)
                            context_length+=1
                    times = context_length/ 7
                    temp_max = float(context_list[6]) 
                    for loop in range(1,times):
                        if(temp_max < float(context_list[loop*7+6])):
                            temp_max = float(context_list[loop*7+6])
                if(FSCOR_list[index].get_family()==FSCOR_list[inner_index].get_family()):
                    sort_list.append('p,p,'+str(temp_max))
                else:
                    sort_list.append('n,p,'+str(temp_max))
                pdb_name = FSCOR_list[index].get_pdb()+'_to_'+FSCOR_list[inner_index].get_pdb()
                pdb_name2=FSCOR_list[inner_index].get_pdb()+'_to_'+FSCOR_list[index].get_pdb()
                result = search_family(pdb_name,pdb_name2)
                sortbig_list.append(result+str(temp_max))
    ###WRITE_FILE('big_log',big_log)
    WRITE_FILE(outfile+'_0_RAW_another',sort_list)
    WRITE_FILE(outfile+'_2_RAW_another',sortbig_list)
            ###print FSCOR_list[index].get_family()+' '+FSCOR_list[inner_index].get_family()
            ###print FSCOR_list[index].get_pdb()+' '+FSCOR_list[inner_index].get_pdb()
            ###count+=1
    return 
#######################################################################
def Raw_TtoR_Process(TtoR_list,TtoR_document_path,outfile):
    count= 0
    sort_list = []
    sortbig_list = []
    for T_count in range(0,227):
        for R_count in range(227,419):
            count+=1
            #   檔案名稱
            file_document_name = TtoR_list[T_count].get_pdb()+'_to_'+TtoR_list[R_count].get_pdb()
            file_document_name2 = TtoR_list[R_count].get_pdb()+'_to_'+TtoR_list[T_count].get_pdb()
            temp_max = 0
            context_list=[]
            context_length = 0
            try :
                with open(TtoR_document_path+file_document_name+'/semiG_result.php','r') as file:
                    for each_line in file:
                        context_list.append(each_line)
                        context_length+=1
                times = context_length/7
                temp_max = float(context_list[6])
                print 'test'+file_document_name
                for loop in range(1,times):
                    if(temp_max < float(context_list[loop*7+6])):
                        temp_max = float(context_list[loop*7 + 6])
            except:
                with open(TtoR_document_path+file_document_name2+'/semiG_result.php','r') as file:
                    for each_line in file:
                        context_list.append(each_line)
                        context_length+=1
                times = context_length/7
                temp_max = float(context_list[6])
                for loop in range(1,times):
                    if(temp_max < float(context_list[loop*7+6])):
                        temp_max = float(context_list[loop*7+6])
            if(TtoR_list[T_count].get_family()==TtoR_list[R_count].get_family()):
                sort_list.append('p,p,'+str(temp_max))
            else:
                sort_list.append('n,p,'+str(temp_max))
            pdb_name = TtoR_list[T_count].get_pdb()+'_to_'+TtoR_list[R_count].get_pdb()
            pdb_name2 = TtoR_list[R_count].get_pdb()+'_to_'+TtoR_list[T_count].get_pdb()
            result = search_family(pdb_name,pdb_name2)
            sortbig_list.append(result+str(temp_max))
    WRITE_FILE(outfile+'_0_RAW_another',sort_list)
    WRITE_FILE(outfile+'_2_RAW_another',sortbig_list)
    ###print count 
    return 
#######################################################################
def FSCOR_Process(FSCOR_list,FSCOR_document_path,outfile):
    #count = 0
    d_list = []
    dSAS_list=[]
    dSI_list = []
    dMI_list = []
    d2SAS_list = []
    d2MI_list = []
    d2SI_list = []
    for index in range(len(FSCOR_list)):
        for inner_index in range(index+1,len(FSCOR_list)):
            if(index!=inner_index):
                file_document_name = FSCOR_list[index].get_pdb()+'_to_'+FSCOR_list[inner_index].get_pdb()
                file_document_name2 = FSCOR_list[inner_index].get_pdb()+'_to_'+FSCOR_list[index].get_pdb()
                context_length = 0
                min = 0
                align_length = 0
                gap_num_seq1 = 0
                gap_num_seq2 = 0
                compare_list=[]
                seq1=''
                seq2=''
                try:
                    with open(FSCOR_document_path+file_document_name+'/semiG_result.php','r') as file:
                        for each_line in file:
                            context_length+=1
                    times = context_length / 7
                    for loop in range(times):
                        temp_RMSD = 0
                        temp_gap_seq1 = 0
                        temp_gap_seq2 = 0
                        temp_length = 0
                        with open(FSCOR_document_path+file_document_name+'/profit_log'+str(loop),'r') as file:
                            for each_line in file:
                                if(each_line.find('RMS')!=-1):
                                    temp_RMSD = float(each_line.replace('RMS:',''))
                        context_list = []
                        with open(FSCOR_document_path+file_document_name+'/ori_ali_seq.pir'+str(loop),'r') as file:
                            for each_line in file:
                                context_list.append(each_line)
                        seq1 = context_list[2]
                        seq2 = context_list[5]
                        for i in range(0,len(seq1)-2):
                            if(seq1[i]!='-'):
                                temp_gap_seq1+=1
                            if(seq2[i]!='-'):
                                temp_gap_seq2+=1
                            if(seq1[i]!='-' and seq2[i]!='-'):
                                temp_length+=1
                        ###print file_document_name+' length='+str(temp_length)+' gap1='+str(temp_gap_seq1)+' gap2='+str(temp_gap_seq2)
                        pdb = Compare_pdb(temp_RMSD,temp_gap_seq1,temp_gap_seq2,temp_length)
                        compare_list.append(pdb)
                except:
                    with open(FSCOR_document_path+file_document_name2+'/semiG_result.php','r') as file:
                        for each_line in file:
                            context_list.append(each_line)
                            context_length+=1
                        times = context_length / 7
                    for loop in range(times):
                        temp_RMSD = 0
                        temp_gap_seq1 = 0
                        temp_gap_seq2 = 0
                        temp_length = 0
                        with open(FSCOR_document_path+file_document_name2+'/profit_log'+str(loop),'r') as file:
                            for each_line in file:
                                if(each_line.find('RMS')!=-1):
                                    temp_RMSD = float(each_line.replace('RMS:',''))
                        context_list = []
                        with open(FSCOR_document_path+file_document_name2+'/ori_ali_seq.pir'+str(loop),'r') as file:
                            for each_line in file:
                                context_list.append(each_line)
                        seq1 = context_list[2]
                        seq2 = context_list[5]
                        for i in range(0,len(seq1)-2):
                            if(seq1[i]!='-'):
                                temp_gap_seq1+=1
                            if(seq2[i]!='-'):
                                temp_gap_seq2+=1
                            if(seq1[i]!='-' and seq2[i]!='-'):
                                temp_length+=1
                        ###print file_document_name2+' length = '+str(temp_length)+' gap1= '+str(temp_gap_seq1)+' gap2 ='+str(temp_gap_seq2)
                        pdb = Compare_pdb(temp_RMSD,temp_gap_seq1,temp_gap_seq2,temp_length)
                        compare_list.append(pdb)
                if(len(seq1)!=len(seq2)):
                    print file_document_name+' not equal'+str(len(seq1))+' '+str(len(seq2))
                else:
                    print file_document_name+' equal'
                min = compare_list[0].getRMSD()
                gap_num_seq1 = compare_list[0].get_gap1()
                gap_num_seq2 = compare_list[0].get_gap2()
                align_length = compare_list[0].get_align()
                for i in range(1,len(compare_list)):
                    if(min > compare_list[i].getRMSD()):
                        min = compare_list[i].getRMSD()
                        gap_num_seq1 = compare_list[i].get_gap1()
                        gap_num_seq2 = compare_list[i].get_gap2()
                        align_length = compare_list[i].get_align()
                try:
                    SAS = min*100 / align_length
                    SAS = 0 - SAS
                    SI = (min * MIN(gap_num_seq1,gap_num_seq2))/align_length
                    SI = 0 - SI
                    MI = 1 - ((1+align_length)/((1+(min/1.5))*(1+MIN(gap_num_seq1,gap_num_seq2))))
                    MI = 0 - MI
                except:
                    SI = 0
                    MI = 0
                    SAS = 0
                if(FSCOR_list[index].get_family()==FSCOR_list[inner_index].get_family()):
                    dSAS_list.append('p,p,'+str(SAS))
                    dSI_list.append('p,p,'+str(SI))
                    dMI_list.append('p,p,'+str(MI))
                else:
                    dSAS_list.append('n,p,'+str(SAS))
                    dSI_list.append('n,p,'+str(SI))
                    dMI_list.append('n,p,'+str(MI))
                pdb_name = FSCOR_list[index].get_pdb()+'_to_'+FSCOR_list[inner_index].get_pdb()
                pdb_name2 = FSCOR_list[inner_index].get_pdb()+'_to_'+FSCOR_list[index].get_pdb()
                d_list.append(pdb_name+' pdb_min = '+str(min)+' align= '+str(align_length)+' gpa1='+str(gap_num_seq1)+' gap2='+str(gap_num_seq2)) 
                result = search_family(pdb_name,pdb_name2)
                d2SAS_list.append(result+str(SAS))
                d2SI_list.append(result+str(SI))
                d2MI_list.append(result+str(MI))
    ###WRITE_FILE('center46_FSCOR_0_log',dlist)
    ###WRITE_FILE('center46_FSCOR_2_log',d2list)
    WRITE_FILE(outfile+'_log_another',d_list)
    WRITE_FILE(outfile+'_0_SAS_another',dSAS_list)
    WRITE_FILE(outfile+'_0_SI_another',dSI_list)
    WRITE_FILE(outfile+'_0_MI_another',dMI_list)
    WRITE_FILE(outfile+'_2_SI_another',d2SI_list)
    WRITE_FILE(outfile+'_2_MI_another',d2MI_list)
    WRITE_FILE(outfile+'_2_SAS_another',d2SAS_list)
#######################################################################
def TtoR_Process(TtoR_list,TtoR_document_path,outfile):
    count = 0
    d_list = []
    SAS_list= []
    MI_list = []
    SI_list = []
    SAS2_list= []
    MI2_list = []
    SI2_list = []
### test
    for T_count in range(0,227):
        SAS = 0
        SI = 0
        MI = 0
        for R_count in range(227,419):
            count+=1
            #   檔案名稱
            file_document_name = TtoR_list[T_count].get_pdb()+'_to_'+TtoR_list[R_count].get_pdb()
            file_document_name2 = TtoR_list[R_count].get_pdb()+'_to_'+TtoR_list[T_count].get_pdb()
            ###print file_document_name
            context_length = 0
            min = 0
            align_length = 0
            gap_num_seq1 = 0
            gap_num_seq2 = 0
            compare_list = []
            try:
                with open(TtoR_document_path+file_document_name+'/semiG_result.php','r') as file:
                    for each_line in file:
                        context_length+=1
                    times = context_length / 7
                for loop in range(0,times):
                    temp_RMSD = 0
                    temp_gap_seq1=0
                    temp_gap_seq2= 0
                    temp_length = 0
                    context_list=[]
                    with open(TtoR_document_path+file_document_name+'/ori_ali_seq.pir'+str(loop),'r') as file2:
                        for each_line in file2:
                            context_list.append(each_line)
                    seq1 = context_list[2]
                    seq2 = context_list[5]
                    for index in range(len(seq1)-2):
                        if(seq1[index]!='-'):
                            temp_gap_seq1+=1
                        if(seq2[index]!='-'):
                            temp_gap_seq2+=1
                        if(seq1[index]!='-' and seq2[index]!='-'):
                            temp_length+=1
                    with open(TtoR_document_path+file_document_name+'/profit_log'+str(loop),'r') as file3:
                        for each_line in file3:
                            if(each_line.find('RMS')!=-1):
                                temp_RMSD = float(each_line.replace('RMS:',''))
                    pdb = Compare_pdb(temp_RMSD,temp_gap_seq1,temp_gap_seq2,temp_length)
                    compare_list.append(pdb) 
            except:
                with open(TtoR_document_path+file_document_name2+'/semiG_result.php','r') as file:
                    for each_line in file:
                        context_length+=1
                    times = context_length / 7
                for loop in range(0,times):
                    temp_RMSD = 0
                    temp_gap_seq1=0
                    temp_gap_seq2= 0
                    temp_length = 0
                    context_list=[]
                    with open(TtoR_document_path+file_document_name2+'/ori_ali_seq.pir'+str(loop),'r') as file2:
                        for each_line in file2:
                            context_list.append(each_line)
                    seq1 = context_list[2]
                    seq2 = context_list[5]
                    for index in range(len(seq1)-2):
                        if(seq1[index]!='-'):
                            temp_gap_seq1+=1
                        if(seq2[index]!='-'):
                            temp_gap_seq2+=1
                        if(seq1[index]!='-' and seq2[index]!='-'):
                            temp_length+=1
                    with open(TtoR_document_path+file_document_name2+'/profit_log'+str(loop),'r') as file3:
                        for each_line in file3:
                            if(each_line.find('RMS')!=-1):
                                temp_RMSD = float(each_line.replace('RMS:',''))
                    pdb = Compare_pdb(temp_RMSD,temp_gap_seq1,temp_gap_seq2,temp_length)
                    compare_list.append(pdb) 
            min = compare_list[0].getRMSD()
            gap_num_seq1 = compare_list[0].get_gap1()
            gap_num_seq2 = compare_list[0].get_gap2()
            align_length = compare_list[0].get_align()
            for index in range(1,len(compare_list)):
                if(min > compare_list[index].getRMSD()):
                    min = compare_list[index].getRMSD()
                    gap_num_seq1 = compare_list[index].get_gap1()
                    gap_num_seq2 = compare_list[index].get_gap2()
                    align_length = compare_list[index].get_align()
            try:
                SAS = (min*100)/align_length
                SAS = 0 - SAS
                SI = (min * MIN(gap_num_seq1,gap_num_seq2))/ align_length
                SI = 0 - SI
                MI = 1-((1+align_length)/((1+(min/1.5))*(1+MIN(gap_num_seq1,gap_num_seq2))))
                MI = 0 - MI
            except:
                SI = 0 
                MI = 0
                SAS = 0
            pdb_name = TtoR_list[T_count].get_pdb()+'_to_'+TtoR_list[R_count].get_pdb()
            pdb_name2 = TtoR_list[R_count].get_pdb()+'_to_'+TtoR_list[T_count].get_pdb()
            d_list.append(pdb_name+' pdb_min = '+str(min)+' align= '+str(align_length)+' gpa1='+str(gap_num_seq1)+' gap2='+str(gap_num_seq2)) 
   ###         print pdb_name+" pdb_min = "+str(min)+' align='+str(align_length)+' gap1 = '+str(gap_num_seq1)+' gap2='+str(gap_num_seq2)
            result = search_family(pdb_name,pdb_name2)
        ### test
        ###
            SAS2_list.append(result+str(SAS))
            SI2_list.append(result+str(SI))
            MI2_list.append(result+str(MI))
            if(TtoR_list[T_count].get_family()==TtoR_list[R_count].get_family()):
                SAS_list.append('p,p,'+str(SAS))
                SI_list.append('p,p,'+str(SI))
                MI_list.append('p,p,'+str(MI))
            else:
                SAS_list.append('n,p,'+str(SAS))
                SI_list.append('n,p,'+str(SI))
                MI_list.append('n,p,'+str(MI))
    WRITE_FILE(outfile+'_log_another',d_list)
    WRITE_FILE(outfile+'_0_SAS_another',SAS_list)
    WRITE_FILE(outfile+'_0_SI_another',SI_list)
    WRITE_FILE(outfile+'_0_MI_another',MI_list)
    WRITE_FILE(outfile+'_2_SAS_another',SAS2_list)
    WRITE_FILE(outfile+'_2_SI_another',SI2_list)
    WRITE_FILE(outfile+'_2_MI_another',MI2_list)
    return 
#######################################################################
def search_family(pdb_name,pdb_name2):
    result = ''
    if(bigfamily_compare.find(pdb_name)!=-1):
        count = bigfamily_compare.find(pdb_name)
        for index in range(count+17,count+21):
            result+=bigfamily_compare[index]
    elif(bigfamily_compare.find(pdb_name2)!=-1):
        count = bigfamily_compare.find(pdb_name2)
        for index in range(count+17,count+21):
            result+=bigfamily_compare[index]
    return result
#######################################################################
def MIN(numberA,numberB):
    if(numberA<numberB):
        return numberA
    else:
        return numberB
#######################################################################
###  做測試
def TEST(list):
    for count in range(0,len(list)):
        print list[count].test()
#######################################################################
### 輸出檔案
def WRITE_FILE(outname,list):
    with open(outname,'w') as file:
        for index in range(len(list)):
            file.write(list[index]+'\n')
#######################################################################
if __name__ =='__main__':
    FSCOR_list = read_file("/home/watchlee/Research_Programming/research_data/inputDataset/SARA_FSCOR.sa")
    ###TEST(FSCOR_list)
    TtoR_list = read_file("/home/watchlee/Research_Programming/research_data/inputDataset/SARA_TtoR-FSCOR.sa")
    ###TEST(TtoR_list)
###FSCOR setting input file
    FSCOR_file = ['/home/watchlee/Research_Programming/RMSD/alignment_main/23-4L_matrix-O6E1.5-SARA_FSCOR_23C_4L_result-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/4L_matrix-O15E1-SARA_FSCOR_4L_result-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/matrix-O4E3.5-FSCOR-semiG.job/'
                  ,'/home/bingts/iPARTS2_training/alignment_main/final_matrix.txt-O5E0.5-SARA_FSCOR_over1k_23c-semiG.job/'
                  ,'/home/bingts/iPARTS2_training/alignment_main/matrix.txt-O12E1-SARA_FSCOR_over1k_69c-semiG.job/'
                  ,'/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/iPARTS_BLOSUM-like_SM-O6E1-new_encoded_iPARTS-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/5dims_c46_K10_matrix-O8E0.5-46C_K10-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS_BLOSUM-like_SM-O6E1-SARA_FSCOR.sa-semiG.job/'
                  ,'/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/5dims_c46_K30_matrix-O9E1-46c_30k-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/5dims_c46_K60_matrix-O7E2-iPARTS_5D_46C_K60_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/5dims_c46_K90_matrix-O6E2.5-iPARTS_5D_46C_K90_SARA_FSCOR-semiG.job/'
                  
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter02_matrix-O8E1.5-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter03_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter04_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter05_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter02_matrix-O13E4-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter06_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter06_matrix-O13E3.5-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter07_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/']

    FSCOR_output_file = ['23C_4L_FSCOR','4L_23C_FSCOR','46C_FSCOR','23C_FSCOR','69C_FSCOR','iPARTS_FSCOR','5K_46C_K10_FSCOR','iPARTS_FSCOR_old','5K_46C_K30_FSCOR','5K_46C_K60_FSCOR','5K_46C_K90_FSCOR','iter02_23C_FSCOR','iter03_23C_FSCOR','iter04_23C_FSCOR','iter05_23C_FSCOR','iter02_MI_FSCOR','iter06_23C_FSCOR','iter06_MI_FSCOR','iter07_23_FSCOR']
###TtoR setting input file 
    TtoR_output_file = ['23C_4L_TtoR','4L_23C_TtoR','46C_TtoR','23C_TtoR','69C_TtoR','iPARTS_TtoR','5K_46C_K10_TtoR','5K_46C_K30_TtoR','new_23C_TtoR','5K_46C_K60_TtoR','5K_46C_K90_TtoR','iter02_23C_TtoR','iter03_23C_TtoR','iter04_23C_TtoR','iter05_23C_TtoR','iter02_MI_TtoR','iter07_23_TtoR']
    TtoR_file = ["/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/23-4L_matrix-O6E1.5-SARA_FSCOR_23C_4L_result-semiG.job/"
                 ,"/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/4L_matrix-O15E1-SARA_FSCOR_4L_result-semiG.job/"
                 ,"/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/matrix-O4E3.5-FSCOR-semiG.job/"
                 ,"/home/bingts/iPARTS2_training/alignment_main/final_matrix.txt-O5E0.5-SARA_FSCOR_over1k_23c-semiG.job/"
                 ,"/home/bingts/iPARTS2_training/alignment_main/matrix.txt-O12E1-SARA_FSCOR_over1k_69c-semiG.job/"
                 ,'/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/iPARTS_BLOSUM-like_SM-O6E1-new_encoded_iPARTS-semiG.job/'
                 ,'/home/watchlee/Research_Programming/RMSD/alignment_main/5dims_c46_K10_matrix-O8E0.5-46C_K10-semiG.job/'
                 ,'/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/5dims_c46_K30_matrix-O9E1-46c_30k-semiG.job/'
                 ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPART2_23c_matrix.txt-O5E0.5-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                 ,'/home/watchlee/Research_Programming/RMSD/alignment_main/5dims_c46_K60_matrix-O7E2-iPARTS_5D_46C_K60_SARA_FSCOR-semiG.job/'
                 ,'/home/watchlee/Research_Programming/RMSD/alignment_main/5dims_c46_K90_matrix-O6E2.5-iPARTS_5D_46C_K90_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter02_matrix-O8E1.5-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter03_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter04_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter05_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter02_matrix-O13E4-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter07_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/']

    F_index =18
    T_index =10  
    TtoR_document_path = TtoR_file[T_index]
    FSCOR_document_path = FSCOR_file[F_index]
    pdbpath = '../pdb/'
    oneDseq_path = '../1Dseq/'
    #TtoR_Process(TtoR_list,TtoR_file[T_index],TtoR_output_file[T_index])
    FSCOR_Process(FSCOR_list,FSCOR_file[F_index],FSCOR_output_file[F_index])
    #Raw_TtoR_Process(TtoR_list,TtoR_file[T_index],TtoR_output_file[T_index])
    Raw_FSCOR_Process(FSCOR_list,FSCOR_file[F_index],FSCOR_output_file[F_index])
   # for index in range(len(FSCOR_output_file)):
    #for index in range(0,6):
    #    FSCOR_Process(FSCOR_list,FSCOR_file[index],FSCOR_output_file[index])
    #    TtoR_Process(TtoR_list,TtoR_file[index],TtoR_output_file[index])
    #    Raw_TtoR_Process(TtoR_list,TtoR_file[index],TtoR_output_file[index])
    #    Raw_FSCOR_Process(FSCOR_list,FSCOR_file[index],FSCOR_output_file[index])

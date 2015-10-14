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
    big_log=[] 
    sortbig_list= []
    for index in range(0,len(FSCOR_list)):
        pdb_max = 0
        family1= 0
        family2 =0 
        for inner_index in range(0,len(FSCOR_list)):
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
                if(pdb_max < temp_max):
                    pdb_max = temp_max
                    family1 = index
                    family2 = inner_index
        ###print str(family1)+' '+str(family2)+' '+str(pdb_max)    
        if(FSCOR_list[family1].get_family()==FSCOR_list[family2].get_family()):
            sort_list.append('p,p,'+str(pdb_max))
        else:
            sort_list.append('n,p,'+str(pdb_max))
        pdb_name= FSCOR_list[family1].get_pdb()+'_to_'+FSCOR_list[family2].get_pdb()
        pdb_name2= FSCOR_list[family2].get_pdb()+'_to_'+FSCOR_list[family1].get_pdb()
        result = search_family(pdb_name,pdb_name2)
        sortbig_list.append(result+str(pdb_max))
        '''
        if(FSCOR_list[family2].get_bigfamily()==FSCOR_list[family1].get_bigfamily()):
            i+=1
            sortbig_list.append('p,p,'+str(pdb_max))
            big_log.append(FSCOR_list[family1].get_pdb()+' '+FSCOR_list[family2].get_pdb())
            big_log.append(FSCOR_list[family1].get_bigfamily()+'\n'+FSCOR_list[family2].get_bigfamily())
            big_log.append('p,p,'+str(pdb_max)+' '+str(i)+'\n')
        else:
            j+=1
            sortbig_list.append('n,p,'+str(pdb_max))
            big_log.append(FSCOR_list[family1].get_pdb()+' '+FSCOR_list[family2].get_pdb())
            big_log.append(FSCOR_list[family1].get_bigfamily()+'\n'+FSCOR_list[family2].get_bigfamily())
            big_log.append('n,p,'+str(pdb_max)+' '+str(j)+'\n')
        '''
    ###WRITE_FILE('big_log',big_log)
    
    WRITE_FILE(outfile+'_0_RAW',sort_list)
    WRITE_FILE(outfile+'_2_RAW',sortbig_list)
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
        pdb_max = 0
        family1 = 0
        family2 = 0
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
                temp_max = int(context_list[6])
                for loop in range(1,times):
                    if(temp_max < int(context_list[loop*7+6])):
                        temp_max = int(context_list[loop*7 + 6])
            except:
                with open(TtoR_document_path+file_document_name2+'/semiG_result.php','r') as file:
                    for each_line in file:
                        context_list.append(each_line)
                        context_length+=1
                times = context_length/7
                temp_max = int(context_list[6])
                for loop in range(1,times):
                    if(temp_max < int(context_list[loop*7+6])):
                        temp_max = int(context_list[loop*7+6])
            if(pdb_max < temp_max):
                pdb_max = temp_max
                family1 = R_count
                family2 = T_count
        if(TtoR_list[family1].get_family()==TtoR_list[family2].get_family()):
            sort_list.append('p,p,'+str(pdb_max))
        else:
            sort_list.append('n,p,'+str(pdb_max))
        pdb_name = TtoR_list[family1].get_pdb()+'_to_'+TtoR_list[family2].get_pdb()
        pdb_name2 = TtoR_list[family2].get_pdb()+'_to_'+TtoR_list[family1].get_pdb()
        result = search_family(pdb_name,pdb_name2)
        sortbig_list.append(result+str(pdb_max))

    WRITE_FILE(outfile+'_0_RAW',sort_list)
    WRITE_FILE(outfile+'_2_RAW',sortbig_list)
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
        pdb_min = 999999
        family1 = 0
        family2 = 0
        gap_pdb1 = 0
        gap_pdb2 = 0
        align = 0
        for inner_index in range(len(FSCOR_list)):
            if(index!=inner_index):
                file_document_name = FSCOR_list[index].get_pdb()+'_to_'+FSCOR_list[inner_index].get_pdb()

                file_document_name2 = FSCOR_list[inner_index].get_pdb()+'_to_'+FSCOR_list[index].get_pdb()
                context_length = 0
                min = 0
                align_length = 0
                gap_num_seq1 = 0
                gap_num_seq2 = 0
                compare_list=[]
                
                try:
                    context_list = []
                    with open(FSCOR_document_path+file_document_name+'/semiG_result.php','r') as file:
                        for each_line in file:
                            context_list.append(each_line)
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
                        seq1 = context_list[0]
                        seq2 = context_list[1]
                        print str(len(seq1))+' '+str(len(seq2))
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
                    context_list = []
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
                        seq1 = context_list[0]
                        seq2 = context_list[1]
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
                if(pdb_min >min):
                    pdb_min = min
                    gap_pdb1 = gap_num_seq1
                    gap_pdb2 = gap_num_seq2
                    align = align_length
                    family1 = index
                    family2 = inner_index
        try:
            SAS = pdb_min * 100 / align
            SAS = 0 - SAS
            SI = (pdb_min * MIN(gap_pdb1,gap_pdb2))/ align
            SI = 0 - SI
            MI = 1 - ((1+align)/((1+(pdb_min/1.5))*(1+MIN(gap_pdb1,gap_pdb2))))
            MI = 0 - MI
        except:
            SI = 0
            MI = 0
            SAS = 0
        if(FSCOR_list[family1].get_family()==FSCOR_list[family2].get_family()):
            dSAS_list.append('p,p,'+str(SAS))
            dSI_list.append('p,p,'+str(SI))
            dMI_list.append('p,p,'+str(MI))
        else:
            dSAS_list.append('n,p,'+str(SAS))
            dSI_list.append('n,p,'+str(SI))
            dMI_list.append('n,p,'+str(MI))
        pdb_name = FSCOR_list[family1].get_pdb()+'_to_'+FSCOR_list[family2].get_pdb()
        pdb_name2 = FSCOR_list[family2].get_pdb()+'_to_'+FSCOR_list[family1].get_pdb()
        d_list.append(pdb_name+' pdb_min = '+str(pdb_min)+' align= '+str(align)+' gpa1='+str(gap_pdb1)+' gap2='+str(gap_pdb2)) 
        result = search_family(pdb_name,pdb_name2)
        d2SAS_list.append(result+str(SAS))
        d2SI_list.append(result+str(SI))
        d2MI_list.append(result+str(MI))
    ###WRITE_FILE('center46_FSCOR_0_log',dlist)
    ###WRITE_FILE('center46_FSCOR_2_log',d2list)
    WRITE_FILE(outfile+'_log',d_list)
    WRITE_FILE(outfile+'_0_SAS',dSAS_list)
    WRITE_FILE(outfile+'_0_SI',dSI_list)
    WRITE_FILE(outfile+'_0_MI',dMI_list)
    WRITE_FILE(outfile+'_2_SI',d2SI_list)
    WRITE_FILE(outfile+'_2_MI',d2MI_list)
    WRITE_FILE(outfile+'_2_SAS',d2SAS_list)
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
    SAS_test_list = []
    SI_test_list = []
    for T_count in range(0,227):
        pdb_min = 999999
        family1 = 0
        family2 = 0
        gap_pdb1 = 0
        gap_pdb2 = 0
        align = 0
        SAS = 0
        SI = 0
        MI = 0
        SAS2 = 0
        SI2 = 0
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
                        
            if(pdb_min > min):
                pdb_min  =min
                family1 = T_count
                family2 = R_count
                gap_pdb1 = gap_num_seq1
                gap_pdb2 = gap_num_seq2
                align = align_length
        try:
            SAS = (pdb_min * 100) / align
            SAS2 = pdb_min * 100 / align
            SAS = 0 - SAS
            SAS2 = 0 - SAS2
            SI = (pdb_min * MIN(gap_pdb1,gap_pdb2))/ align
            SI2 = pdb_min * MIN(gap_pdb1,gap_pdb2)/ align
            SI = 0 - SI
            SI2 = 0 - SI2
            MI = 1 - ((1+align)/((1+pdb_min/1.5)*(1+MIN(gap_pdb1,gap_pdb2))))
            MI = 0 - MI
        except:
            SAS = 0
            SAS2 = 0
            SI2 = 0
            SI = 0
            MI = 0
        pdb_name = TtoR_list[family1].get_pdb()+'_to_'+TtoR_list[family2].get_pdb()
        pdb_name2 = TtoR_list[family2].get_pdb()+'_to_'+TtoR_list[family1].get_pdb()
        d_list.append(pdb_name+' pdb_min = '+str(pdb_min)+' align= '+str(align)+' gpa1='+str(gap_pdb1)+' gap2='+str(gap_pdb2)) 
        print pdb_name+" pdb_min = "+str(pdb_min)+' align='+str(align)+' gap1 = '+str(gap_pdb1)+' gap2='+str(gap_pdb2)
        result = search_family(pdb_name,pdb_name2)
        ### test
        SAS_test_list.append(result+str(SAS2))
        SI_test_list.append(result+str(SI2))
        ###
        SAS2_list.append(result+str(SAS))
        SI2_list.append(result+str(SI))
        MI2_list.append(result+str(MI))
        if(TtoR_list[family1].get_family()==TtoR_list[family2].get_family()):
            SAS_list.append('p,p,'+str(SAS))
            SI_list.append('p,p,'+str(SI))
            MI_list.append('p,p,'+str(MI))
        else:
            SAS_list.append('n,p,'+str(SAS))
            SI_list.append('n,p,'+str(SI))
            MI_list.append('n,p,'+str(MI))
    WRITE_FILE(outfile+'_log',d_list)
    WRITE_FILE(outfile+'_2_test_SAS',SAS_test_list)
    WRITE_FILE(outfile+'_2_test_SI',SI_test_list)
    WRITE_FILE(outfile+'_0_SAS',SAS_list)
    WRITE_FILE(outfile+'_0_SI',SI_list)
    WRITE_FILE(outfile+'_0_MI',MI_list)
    WRITE_FILE(outfile+'_2_SAS',SAS2_list)
    WRITE_FILE(outfile+'_2_SI',SI2_list)
    WRITE_FILE(outfile+'_2_MI',MI2_list)
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
    FSCOR_file = ['/home/watchlee/Research_Programming/RMSD/alignment_main/23-4L_matrix-O6E1.5-SARA_FSCOR_23C_4L_result-semiG.job/','/home/watchlee/Research_Programming/RMSD/alignment_main/4L_matrix-O15E1-SARA_FSCOR_4L_result-semiG.job/','/home/watchlee/Research_Programming/RMSD/alignment_main/matrix-O4E3.5-FSCOR-semiG.job/','/home/bingts/iPARTS2_training/alignment_main/final_matrix.txt-O5E0.5-SARA_FSCOR_over1k_23c-semiG.job/','/home/bingts/iPARTS2_training/alignment_main/matrix.txt-O12E1-SARA_FSCOR_over1k_69c-semiG.job/','/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/iPARTS_BLOSUM-like_SM-O6E1-SARA_FSCOR.sa-semiG.job/']
    FSCOR_output_file = ['23C_4L_FSCOR','4L_23C_FSCOR','46C_FSCOR','23C_FSCOR','69C_FSCOR','iPARTS_FSCOR']

###TtoR setting input file 
    TtoR_output_file = ['23C_4L_TtoR','4L_23C_TtoR','46C_TtoR','23C_TtoR','69C_TtoR','iPARTS_TtoR']
    TtoR_file = ["/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/23-4L_matrix-O6E1.5-SARA_FSCOR_23C_4L_result-semiG.job/","/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/4L_matrix-O15E1-SARA_FSCOR_4L_result-semiG.job/","/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/matrix-O4E3.5-FSCOR-semiG.job/","/home/bingts/iPARTS2_training/alignment_main/final_matrix.txt-O5E0.5-SARA_FSCOR_over1k_23c-semiG.job/","/home/bingts/iPARTS2_training/alignment_main/matrix.txt-O12E1-SARA_FSCOR_over1k_69c-semiG.job/",'/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/iPARTS_BLOSUM-like_SM-O6E1-SARA_FSCOR.sa-semiG.job/']

    F_index = 5 
    T_index = 5
    TtoR_document_path = TtoR_file[T_index]
    FSCOR_document_path = FSCOR_file[F_index]

    pdbpath = '../pdb/'
    oneDseq_path = '../1Dseq/'
   # TtoR_Process(TtoR_list,TtoR_file[T_index],TtoR_output_file[T_index])
    FSCOR_Process(FSCOR_list,FSCOR_file[F_index],FSCOR_output_file[F_index])
    #Raw_TtoR_Process(TtoR_list,TtoR_file[T_index],TtoR_output_file[T_index])
  #  Raw_FSCOR_Process(FSCOR_list,FSCOR_file[F_index],FSCOR_output_file[F_index])
   # for index in range(len(FSCOR_output_file)):
    #for index in range(0,5):
    #    FSCOR_Process(FSCOR_list,FSCOR_file[index],FSCOR_output_file[index])
        #TtoR_Process(TtoR_list,TtoR_file[index],TtoR_output_file[index])
        #Raw_TtoR_Process(TtoR_list,TtoR_file[index],TtoR_output_file[index])
    #    Raw_FSCOR_Process(FSCOR_list,FSCOR_file[index],FSCOR_output_file[index])
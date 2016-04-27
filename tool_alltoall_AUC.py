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
#   發現processed_eachto418_d2的資料有誤，改用新的
read_file = open('../new_processed_eachto418_d2','r')
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
#-----------------------Converting string into floating or integer number-------------------------------#
def num(value):
    try:
        return int(value)
    except ValueError:
        return float(value)
#######################################################################
class Compare_pdb:
    RMSD= 0
    num_gap1 = 0
    num_gap2 = 0
    length_align =0
    match = 0
    def __init__(self,RMSD,num_gap1,num_gap2,length_align,match):
        self.RMSD = RMSD
        self.num_gap1 = num_gap1
        self.num_gap2 = num_gap2
        self.length_align = length_align
        self.match = match
    def getMatch(self):
        return self.match
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
def other_PSI_Process(FSCOR_list,FSCOR_document_path,outfile):
    #count = 0
    dPSI_list=[]
    d2PSI_list = []
    dPSI_analysis_list = []
    d2PSI_analysis_list = []
    Positive_PSI_analysis_list=[]
    Negative_PSI_analysis_list=[]
    Positive2_PSI_analysis_list=[]
    Negative2_PSI_analysis_list=[]
    sameParent_list = []
    for index in range(len(FSCOR_list)):
        for inner_index in range(index+1,len(FSCOR_list)):
            if(index!=inner_index):
                #--------------------profit SARA_FSCOR
                file_document_name = FSCOR_list[index].get_pdb()+'_to_'+FSCOR_list[inner_index].get_pdb()
                file_document_name2 = FSCOR_list[inner_index].get_pdb()+'_to_'+FSCOR_list[index].get_pdb()
                #----------------SARA_FSCOR_253
                #file_document_name = FSCOR_list[index].get_pdb()+'-'+FSCOR_list[inner_index].get_pdb()
                #file_document_name2 = FSCOR_list[inner_index].get_pdb()+'-'+FSCOR_list[index].get_pdb()
                max = 0
                #---------------------計算ＰＳＩ用的value
                PSI_times=0
                psi_direct= 0.0
                #---------------因為不知道A_to_B or B_to_A才這樣寫---------------------#
                #----------------------------讀取setter_result.php中 seq長度------------------#
                #----------------------------讀取profit_log------------------#
                #----------------------------跑出result------------------#
                try:
                    '''
                    with open(FSCOR_document_path+file_document_name+'/log','r') as file:
                        for each_line in file:
                            if(each_line.find('PSI')!=-1):
                                psi_direct=float(each_line.replace('<PSI>','').replace('</PSI>',''))
                                print each_line.replace('\n','')+' '+str(psi_direct)
                    In order to handle the reference problem, we have to use another way to calculate AUC, which calculating by PSI
                    with open(FSCOR_document_path+file_document_name+'/profit_log_P'+str(loop),'r') as file:
                    with open(FSCOR_document_path+file_document_name+'/profit_log_C3'+str(loop),'r') as file:
                    with open(FSCOR_document_path+file_document_name+'/profit_log_S'+str(loop),'r') as file:
                    with open(FSCOR_document_path+file_document_name+'/profit_log','r') as file:
                        #要找每一個核甘酸的RMSD 整體的不要找
                    '''
                    with open(FSCOR_document_path+file_document_name+'/profit_log','r') as file:
                        for each_line in file:
                            try:
                                if(each_line.find('RMS')!=-1):
                                    if(num(each_line.split(':')[2].replace(' ',''))<=4.0):
                                        PSI_times+=1
                            except:
                                #-----第一個RMSD忽略掉
                                pass
                            #找出RMS判斷是否<=4
                            #若為0則全部設為0
                        ###print file_document_name+' length='+str(temp_length)+' gap1='+str(temp_gap_seq1)+' gap2='+str(temp_gap_seq2)
                except:
                    '''
                    with open(FSCOR_document_path+file_document_name2+'/log','r') as file:
                        for each_line in file:
                            if(each_line.find('PSI')!=-1):
                                psi_direct=float(each_line.replace('<PSI>','').replace('</PSI>',''))
                                print each_line.replace('\n','')+' '+str(psi_direct)
                    In order to handle the reference problem, we have to use another way to calculate AUC, which calculating by PSI
                    with open(FSCOR_document_path+file_document_name+'/profit_log_P'+str(loop),'r') as file:
                    with open(FSCOR_document_path+file_document_name+'/profit_log_C3'+str(loop),'r') as file:
                    with open(FSCOR_document_path+file_document_name+'/profit_log_S'+str(loop),'r') as file:
                    '''
                    with open(FSCOR_document_path+file_document_name2+'/profit_log','r') as file:
                        #要找每一個核甘酸的RMSD 整體的不要找
                        for each_line in file:
                            try:
                                if(each_line.find('RMS')!=-1):
                                    if(num(each_line.split(':')[2].replace(' ',''))<=4.0):
                                        PSI_times+=1
                            except:
                                #-----第一個RMSD忽略掉
                                pass
                            #找出RMS判斷是否<=4
                            #若為0則全部設為0
                psi_direct=float(psi_direct*100)
                print FSCOR_document_path+file_document_name2+' '+str(psi_direct)                            
                        ###print file_document_name2+' length = '+str(temp_length)+' gap1= '+str(temp_gap_seq1)+' gap2 ='+str(temp_gap_seq2)
                        #------------------------處理PSI情況----------------------------#
                #----get PSI times
                max = PSI_times 
                #----------------PSI = percentage of surperposed nucleotides or base pair within a given distance cut-off.
                #----------------PSI = 100* number of aligned nucleotides within a threshold of 4.0 A / the length of the shorter of the tow RAN structures
                #----------------------------讀取setter_result.php中 seq長度------------------#
                gap_seq=''
                with open('../1Dseq/'+FSCOR_list[index].get_pdb()+'.seq','r') as file:
                    gap_seq=file.read()
                gap_seq2=''
                with open('../1Dseq/'+FSCOR_list[inner_index].get_pdb()+'.seq','r') as file:
                    gap_seq2=file.read()
                try:
                    min_seq = MIN(len(gap_seq),len(gap_seq2))
                    PSI = float(max*100)/float(min_seq) 
                    print str(float(max*100)/float(min_seq))+' '+' PSI='+str(max)+' MIN_SEQ='+str(min_seq)
                except:
                    PSI=0
                temp_deter=''
                if(FSCOR_list[index].get_family()==FSCOR_list[inner_index].get_family()):
                #    if(min<3):
                    #special_list.append(file_document_name+' pdb_min:'+str(min)+' align:'+str(align_length)+' gpa1:'+str(gap_num_seq1)+' gap2:'+str(gap_num_seq2)+' '+str(PSI)) 
                    dPSI_list.append('p,p '+str(PSI))
                    #----改成直截拿他的來用
                    #dPSI_list.append('p,p '+str(psi_direct))
                #分析用
                    dPSI_analysis_list.append('p,p '+str(PSI)+' '+file_document_name)
                    Positive_PSI_analysis_list.append('p,p '+str(psi_direct)+' '+file_document_name)
                    sameParent_list.append(file_document_name+' p,p')
                    temp_deter=' d=0: p,p '
                else:
                #分析用
                    dPSI_analysis_list.append('n,p '+str(PSI)+' '+file_document_name)
                    Negative_PSI_analysis_list.append('n,p '+str(psi_direct)+' '+file_document_name)
                    dPSI_list.append('n,p '+str(PSI))
                    #----改成直截拿他的來用
                    #dPSI_list.append('n,p '+str(psi_direct))
                    temp_deter=' d=0: n,p '
                pdb_name = FSCOR_list[index].get_pdb()+'_to_'+FSCOR_list[inner_index].get_pdb()
                pdb_name2 = FSCOR_list[inner_index].get_pdb()+'_to_'+FSCOR_list[index].get_pdb()
                result = search_family(pdb_name,pdb_name2)
                temp_deter=temp_deter+'d<=2: '+result
                #----改成直截拿他的來用
                #d2PSI_list.append(result+str(psi_direct))
                d2PSI_list.append(result+str(PSI))
                d2PSI_analysis_list.append(result+str(PSI)+' '+file_document_name)
                #分析用
    atom='S'
    WRITE_FILE(outfile+'_0_PSI_'+atom+'_another_analysis',dPSI_analysis_list)
    WRITE_FILE(outfile+'_0_PSI_'+atom+'_Positive_another_analysis',Positive_PSI_analysis_list)
    WRITE_FILE(outfile+'_0_PSI_'+atom+'_Negative_another_analysis',Negative_PSI_analysis_list)
    WRITE_FILE(outfile+'_2_PSI_'+atom+'_another_analysis',d2PSI_analysis_list)
    WRITE_FILE(outfile+'_0_PSI_'+atom+'_another',dPSI_list)
    WRITE_FILE(outfile+'_2_PSI_'+atom+'_another',d2PSI_list)
def PSI_Process(FSCOR_list,FSCOR_document_path,outfile,atom):
    #count = 0
    d_list = []
    dPSI_list=[]
    d2PSI_list = []
    d2PSI_analysis_list = []
    dPSI_analysis_list = []
    special_list = []
    sameParent_list = []
    d0TPSI_list=[]
    d2TPSI_list=[]
    translate_result_d2=9
    translate_result=9
    for index in range(len(FSCOR_list)):
        for inner_index in range(index+1,len(FSCOR_list)):
            if(index!=inner_index):
                file_document_name = FSCOR_list[index].get_pdb()+'_to_'+FSCOR_list[inner_index].get_pdb()
                file_document_name2 = FSCOR_list[inner_index].get_pdb()+'_to_'+FSCOR_list[index].get_pdb()
                context_length = 0
                max = 0
                align_length = 0
                gap_num_seq1 = 0
                gap_num_seq2 = 0
                compare_list=[]
                seq1=''
                seq2=''
                #---------------------計算ＰＳＩ用的value
                PSI_times=0
                #---------------因為不知道A_to_B or B_to_A才這樣寫---------------------#
                try:
                    with open(FSCOR_document_path+file_document_name+'/semiG_result.php','r') as file:
                        for each_line in file:
                            context_length+=1
                    times = context_length / 7
                    for loop in range(times):
                        PSI_times=0
                        temp_gap_seq1 = 0
                        temp_gap_seq2 = 0
                        temp_length = 0
                        temp_match = 0
                        '''
                        In order to handle the reference problem, we have to use another way to calculate AUC, which calculating by PSI
                        with open(FSCOR_document_path+file_document_name+'/profit_log_P'+str(loop),'r') as file:
                        with open(FSCOR_document_path+file_document_name+'/profit_log_C3'+str(loop),'r') as file:
                        with open(FSCOR_document_path+file_document_name+'/profit_log_S'+str(loop),'r') as file:
                        '''
                        with open(FSCOR_document_path+file_document_name+'/profit_log_'+atom+str(loop),'r') as file:
                            #要找每一個核甘酸的RMSD 整體的不要找
                            for each_line in file:
                                try:
                                    if(each_line.find('RMS')!=-1):
                                        if(num(each_line.split(':')[2].replace(' ',''))<=4.0):
                                            PSI_times+=1
                                except:
                                    #-----第一個RMSD忽略掉
                                    pass
                                #找出RMS判斷是否<=4
                                #若為0則全部設為0
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
                            if(seq1[i]!='-' and seq2[i]!='-' and seq1[i]==seq2[i]):
                                temp_match+=1
                        ###print file_document_name+' length='+str(temp_length)+' gap1='+str(temp_gap_seq1)+' gap2='+str(temp_gap_seq2)
                        pdb = Compare_pdb(PSI_times,temp_gap_seq1,temp_gap_seq2,temp_length,temp_match)
                        compare_list.append(pdb)
                except:
                    with open(FSCOR_document_path+file_document_name2+'/semiG_result.php','r') as file:
                        for each_line in file:
                            context_list.append(each_line)
                            context_length+=1
                        times = context_length / 7
                    for loop in range(times):
                        PSI_times=0
                        temp_gap_seq1 = 0
                        temp_gap_seq2 = 0
                        temp_length = 0
                        temp_match = 0
                        '''
                        In order to handle the reference problem, we have to use another way to calculate AUC, which calculating by PSI
                        with open(FSCOR_document_path+file_document_name+'/profit_log_P'+str(loop),'r') as file:
                        with open(FSCOR_document_path+file_document_name+'/profit_log_C3'+str(loop),'r') as file:
                        with open(FSCOR_document_path+file_document_name+'/profit_log_S'+str(loop),'r') as file:
                        '''
                        with open(FSCOR_document_path+file_document_name2+'/profit_log_'+atom+str(loop),'r') as file:
                            #要找每一個核甘酸的RMSD 整體的不要找
                            for each_line in file:
                                try:
                                    if(each_line.find('RMS')!=-1):
                                        if(num(each_line.split(':')[2].replace(' ',''))<=4.0):
                                            PSI_times+=1
                                except:
                                    #-----第一個RMSD忽略掉
                                    pass
                                #找出RMS判斷是否<=4
                                #若為0則全部設為0
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
                            if(seq1[i]!='-' and seq2[i]!='-' and seq1[i]==seq2[i]):
                                temp_match+=1
                        ###print file_document_name2+' length = '+str(temp_length)+' gap1= '+str(temp_gap_seq1)+' gap2 ='+str(temp_gap_seq2)
                        #------------------------處理PSI情況----------------------------#
                        pdb = Compare_pdb(PSI_times,temp_gap_seq1,temp_gap_seq2,temp_length,temp_match)
                        compare_list.append(pdb)
                if(len(seq1.replace('\n',''))!=len(seq2.replace('\n',''))):
                    print file_document_name+' not equal'+str(len(seq1))+' '+str(len(seq2))
                else:
                    print file_document_name+' equal'
                #----get PSI times
                max = compare_list[0].getRMSD()
                gap_num_seq1 = compare_list[0].get_gap1()
                gap_num_seq2 = compare_list[0].get_gap2()
                align_length = compare_list[0].get_align()
                match = compare_list[0].getMatch() 
                for i in range(1,len(compare_list)):
                    if(max < compare_list[i].getRMSD()):
                        max = compare_list[i].getRMSD()
                        gap_num_seq1 = compare_list[i].get_gap1()
                        gap_num_seq2 = compare_list[i].get_gap2()
                        align_length = compare_list[i].get_align()
                        match = compare_list[i].getMatch() 
                #----------------PSI = percentage of surperposed nucleotides or base pair within a given distance cut-off.
                #----------------PSI = 100* number of aligned nucleotides within a threshold of 4.0 A / the length of the shorter of the tow RAN structures
                gap_seq=''
                with open('../1Dseq/'+FSCOR_list[index].get_pdb()+'.seq','r') as file:
                    gap_seq=file.read()
                gap_seq2=''
                with open('../1Dseq/'+FSCOR_list[inner_index].get_pdb()+'.seq','r') as file:
                    gap_seq2=file.read()
                min_seq = MIN(len(gap_seq),len(gap_seq2))
                try:
                    PSI = float(max*100)/float(min_seq) 
                    print str(float(max*100)/float(min_seq))+' '+' PSI='+str(max)+' MIN_SEQ='+str(min_seq)
                except:
                    PSI=0
                    print 'ERROR'
                temp_deter=''
            
                if(FSCOR_list[index].get_family()==FSCOR_list[inner_index].get_family()):
                #    if(min<3):
                    #special_list.append(file_document_name+' pdb_min:'+str(min)+' align:'+str(align_length)+' gpa1:'+str(gap_num_seq1)+' gap2:'+str(gap_num_seq2)+' '+str(PSI)) 
                    dPSI_list.append('p,p '+str(PSI))
                #分析用
                    dPSI_analysis_list.append('p,p '+str(PSI)+' '+file_document_name)
                    sameParent_list.append(file_document_name+' p,p')
                    d0TPSI_list.append('1 '+str(PSI)) 
                    temp_deter=' d=0: p,p '
                else:
                #分析用
                    dPSI_analysis_list.append('n,p '+str(PSI)+' '+file_document_name)
                    dPSI_list.append('n,p '+str(PSI))
                    d0TPSI_list.append('0 '+str(PSI)) 
                    temp_deter=' d=0: n,p '

                    
                pdb_name = FSCOR_list[index].get_pdb()+'_to_'+FSCOR_list[inner_index].get_pdb()
                pdb_name2 = FSCOR_list[inner_index].get_pdb()+'_to_'+FSCOR_list[index].get_pdb()
                d_list.append(pdb_name+' PSI_max = '+str(max)+' align= '+str(align_length)+' gpa1='+str(gap_num_seq1)+' gap2='+str(gap_num_seq2)) 
                result = search_family(pdb_name,pdb_name2)
                temp_deter=temp_deter+'d<=2: '+result
                d2PSI_list.append(result+str(PSI))
                if(result.find('n,p')!=-1):
                    d2TPSI_list.append('0 '+str(PSI))
                elif(result.find('p,p')!=-1):
                    d2TPSI_list.append('1 '+str(PSI))
                #分析用
                d2PSI_analysis_list.append(result+str(PSI)+' '+file_document_name)
                special_list.append(file_document_name+' PSI_times:'+str(max)+' align:'+str(align_length)+' first_seq:'+str(gap_num_seq1)+' second_seq:'+str(gap_num_seq2)+' PSI:'+str(PSI)+temp_deter+' match:'+str(match)) 
    WRITE_FILE(outfile+'_0_PSI_'+atom+'_another_analysis',dPSI_analysis_list)
    WRITE_FILE(outfile+'_2_PSI_'+atom+'_another_analysis',d2PSI_analysis_list)
    WRITE_FILE(outfile+'_sameParent',sameParent_list)
    WRITE_FILE(outfile+'_PSI_'+atom+'_special_another',special_list)
    WRITE_FILE(outfile+'_log_another',d_list)
    WRITE_FILE(outfile+'_0_PSI_'+atom+'_another',dPSI_list)
    WRITE_FILE(outfile+'_2_PSI_'+atom+'_another',d2PSI_list)
    WRITE_FILE(outfile+'_0_PSI_'+atom+'_pyformat',d0TPSI_list)
    WRITE_FILE(outfile+'_2_PSI_'+atom+'_pyformat',d2TPSI_list)
#######################################################################
#-------------------------應付comment---------------------日後可能會用到--------------------#
def FORCOMMENT_TFSCOR_Process(FSCOR_list,FSCOR_document_path,outfile):
    #count = 0
    d_list = []
    dSAS_list=[]
    dSI_list = []
    dMI_list = []
    d2SAS_list = []
    d2SAS_analysis_list = []
    dSAS_analysis_list = []
    d2MI_list = []
    d2SI_list = []
    special_list = []
    sameParent_list = []
    SAS_family=[]
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
                #---------------因為不知道A_to_B or B_to_A才這樣寫---------------------#
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
                        temp_match = 0
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
                            if(seq1[i]!='-' and seq2[i]!='-' and seq1[i]==seq2[i]):
                                temp_match+=1
                        ###print file_document_name+' length='+str(temp_length)+' gap1='+str(temp_gap_seq1)+' gap2='+str(temp_gap_seq2)
                        pdb = Compare_pdb(temp_RMSD,temp_gap_seq1,temp_gap_seq2,temp_length,temp_match)
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
                        temp_match = 0
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
                            if(seq1[i]!='-' and seq2[i]!='-' and seq1[i]==seq2[i]):
                                temp_match+=1
                        ###print file_document_name2+' length = '+str(temp_length)+' gap1= '+str(temp_gap_seq1)+' gap2 ='+str(temp_gap_seq2)
                        #------------------------處理一般RMSD情況----------------------#
                        pdb = Compare_pdb(temp_RMSD,temp_gap_seq1,temp_gap_seq2,temp_length,temp_match)
                        compare_list.append(pdb)
                if(len(seq1.replace('\n',''))!=len(seq2.replace('\n',''))):
                    print file_document_name+' not equal'+str(len(seq1))+' '+str(len(seq2))
                else:
                    print file_document_name+' equal'
                min = compare_list[0].getRMSD()
                gap_num_seq1 = compare_list[0].get_gap1()
                gap_num_seq2 = compare_list[0].get_gap2()
                align_length = compare_list[0].get_align()
                match = compare_list[0].getMatch() 
                for i in range(1,len(compare_list)):
                    if(min > compare_list[i].getRMSD()):
                        min = compare_list[i].getRMSD()
                        gap_num_seq1 = compare_list[i].get_gap1()
                        gap_num_seq2 = compare_list[i].get_gap2()
                        align_length = compare_list[i].get_align()
                        match = compare_list[i].getMatch() 
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
                temp_deter=''
                if(FSCOR_list[index].get_family()==FSCOR_list[inner_index].get_family()):
                #    if(min<3):
                    #special_list.append(file_document_name+' pdb_min:'+str(min)+' align:'+str(align_length)+' gpa1:'+str(gap_num_seq1)+' gap2:'+str(gap_num_seq2)+' '+str(SAS)) 
                    dSAS_list.append('p,p '+str(SAS))
                #分析用
                    dSAS_analysis_list.append('p,p '+str(SAS)+' '+file_document_name)
                    dSI_list.append('p,p,'+str(SI))
                    dMI_list.append('p,p,'+str(MI))
                    sameParent_list.append(file_document_name+' p,p')
                    temp_deter=' d=0: p,p '
                else:
                #分析用
                    dSAS_analysis_list.append('n,p '+str(SAS)+' '+file_document_name)
                    dSAS_list.append('n,p '+str(SAS))
                    dSI_list.append('n,p,'+str(SI))
                    dMI_list.append('n,p,'+str(MI))
                    temp_deter=' d=0: n,p '
                pdb_name = FSCOR_list[index].get_pdb()+'_to_'+FSCOR_list[inner_index].get_pdb()
                pdb_name2 = FSCOR_list[inner_index].get_pdb()+'_to_'+FSCOR_list[index].get_pdb()
                d_list.append(pdb_name+' pdb_min = '+str(min)+' align= '+str(align_length)+' gpa1='+str(gap_num_seq1)+' gap2='+str(gap_num_seq2)) 
                result = search_family(pdb_name,pdb_name2)
                temp_deter=temp_deter+'d<=2: '+result
                SAS_family.append(FSCOR_list[index].get_pdb()+'|'+FSCOR_list[inner_index].get_pdb()+' '+str(SAS))
                d2SAS_list.append(result+str(SAS))
                #分析用
                d2SAS_analysis_list.append(result+str(SAS)+' '+file_document_name)
                d2SI_list.append(result+str(SI))
                d2MI_list.append(result+str(MI))
                special_list.append(file_document_name+' RMSD:'+str(min)+' align:'+str(align_length)+' first_seq:'+str(gap_num_seq1)+' second_seq:'+str(gap_num_seq2)+' SAS:'+str(SAS)+temp_deter+' match:'+str(match)) 
    ###WRITE_FILE('center46_FSCOR_0_log',dlist)
    ###WRITE_FILE('center46_FSCOR_2_log',d2list)
    WRITE_FILE(outfile+'_0_SAS_another_analysis',dSAS_analysis_list)
    WRITE_FILE(outfile+'_2_SAS_another_analysis',d2SAS_analysis_list)
    WRITE_FILE(outfile+'_sameParent',sameParent_list)
    WRITE_FILE(outfile+'_SAS_special_another',special_list)
    WRITE_FILE(outfile+'_log_another',d_list)
    WRITE_FILE(outfile+'_0_SAS_another',dSAS_list)
    WRITE_FILE(outfile+'_0_SI_another',dSI_list)
    WRITE_FILE(outfile+'_0_MI_another',dMI_list)
    WRITE_FILE(outfile+'_2_SI_another',d2SI_list)
    WRITE_FILE(outfile+'_2_MI_another',d2MI_list)
    WRITE_FILE(outfile+'_2_SAS_another',d2SAS_list)
    WRITE_FILE(outfile+'_SAS_result',SAS_family)
#######################################################################
def FSCOR_Process(FSCOR_list,FSCOR_document_path,outfile):
    #count = 0
    d_list = []
    dSAS_list=[]
    dSI_list = []
    dMI_list = []
    d2SAS_list = []
    d2SAS_analysis_list = []
    dSAS_analysis_list = []
    d2MI_list = []
    d2SI_list = []
    special_list = []
    sameParent_list = []
    SAS_family=[]
    compare_family_list=[]
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
                #---------------------計算ＰＳＩ用的value
                #PSI_times=0
                #---------------因為不知道A_to_B or B_to_A才這樣寫---------------------#
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
                        temp_match = 0
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
                            if(seq1[i]!='-' and seq2[i]!='-' and seq1[i]==seq2[i]):
                                temp_match+=1
                        ###print file_document_name+' length='+str(temp_length)+' gap1='+str(temp_gap_seq1)+' gap2='+str(temp_gap_seq2)
                        pdb = Compare_pdb(temp_RMSD,temp_gap_seq1,temp_gap_seq2,temp_length,temp_match)
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
                        temp_match = 0
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
                            if(seq1[i]!='-' and seq2[i]!='-' and seq1[i]==seq2[i]):
                                temp_match+=1
                        ###print file_document_name2+' length = '+str(temp_length)+' gap1= '+str(temp_gap_seq1)+' gap2 ='+str(temp_gap_seq2)
                        #------------------------處理一般RMSD情況----------------------#
                        pdb = Compare_pdb(temp_RMSD,temp_gap_seq1,temp_gap_seq2,temp_length,temp_match)
                        compare_list.append(pdb)
                if(len(seq1.replace('\n',''))!=len(seq2.replace('\n',''))):
                    print file_document_name+' not equal'+str(len(seq1))+' '+str(len(seq2))
                else:
                    print file_document_name+' equal'
                min = compare_list[0].getRMSD()
                gap_num_seq1 = compare_list[0].get_gap1()
                gap_num_seq2 = compare_list[0].get_gap2()
                align_length = compare_list[0].get_align()
                match = compare_list[0].getMatch() 
                for i in range(1,len(compare_list)):
                    if(min > compare_list[i].getRMSD()):
                        min = compare_list[i].getRMSD()
                        gap_num_seq1 = compare_list[i].get_gap1()
                        gap_num_seq2 = compare_list[i].get_gap2()
                        align_length = compare_list[i].get_align()
                        match = compare_list[i].getMatch() 
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
                temp_deter=''
                pdb_name = FSCOR_list[index].get_pdb()+'_to_'+FSCOR_list[inner_index].get_pdb()
                pdb_name2 = FSCOR_list[inner_index].get_pdb()+'_to_'+FSCOR_list[index].get_pdb()
                temp_class = pdb_name+' '
                if(FSCOR_list[index].get_family()==FSCOR_list[inner_index].get_family()):
                #    if(min<3):
                    #special_list.append(file_document_name+' pdb_min:'+str(min)+' align:'+str(align_length)+' gpa1:'+str(gap_num_seq1)+' gap2:'+str(gap_num_seq2)+' '+str(SAS)) 
                    dSAS_list.append('p,p '+str(SAS))
                #分析用
                    dSAS_analysis_list.append('p,p '+str(SAS)+' '+file_document_name)
                    dSI_list.append('p,p,'+str(SI))
                    dMI_list.append('p,p,'+str(MI))
                    temp_class+='p,p '
                    sameParent_list.append(file_document_name+' p,p')
                    temp_deter=' d=0: p,p '
                else:
                #分析用
                    dSAS_analysis_list.append('n,p '+str(SAS)+' '+file_document_name)
                    dSAS_list.append('n,p '+str(SAS))
                    dSI_list.append('n,p,'+str(SI))
                    dMI_list.append('n,p,'+str(MI))
                    temp_deter=' d=0: n,p '
                    temp_class+='n,p '
                d_list.append(pdb_name+' pdb_min = '+str(min)+' align= '+str(align_length)+' gpa1='+str(gap_num_seq1)+' gap2='+str(gap_num_seq2)) 
                result = search_family(pdb_name,pdb_name2)
                temp_deter=temp_deter+'d<=2: '+result
                temp_class+=result
                compare_family_list.append(temp_class)
                SAS_family.append(FSCOR_list[index].get_pdb()+'|'+FSCOR_list[inner_index].get_pdb()+' '+str(SAS))
                d2SAS_list.append(result+str(SAS))
                #分析用
                d2SAS_analysis_list.append(result+str(SAS)+' '+file_document_name)
                d2SI_list.append(result+str(SI))
                d2MI_list.append(result+str(MI))
                special_list.append(file_document_name+' RMSD:'+str(min)+' align:'+str(align_length)+' first_seq:'+str(gap_num_seq1)+' second_seq:'+str(gap_num_seq2)+' SAS:'+str(SAS)+temp_deter+' match:'+str(match)) 
    ###WRITE_FILE('center46_FSCOR_0_log',dlist)
    ###WRITE_FILE('center46_FSCOR_2_log',d2list)
    WRITE_FILE(outfile+'_all_comparsion',compare_family_list)
    WRITE_FILE(outfile+'_0_SAS_another_analysis',dSAS_analysis_list)
    WRITE_FILE(outfile+'_2_SAS_another_analysis',d2SAS_analysis_list)
    WRITE_FILE(outfile+'_sameParent',sameParent_list)
    WRITE_FILE(outfile+'_SAS_special_another',special_list)
    WRITE_FILE(outfile+'_log_another',d_list)
    WRITE_FILE(outfile+'_0_SAS_another',dSAS_list)
    #WRITE_FILE(outfile+'_0_SI_another',dSI_list)
    #WRITE_FILE(outfile+'_0_MI_another',dMI_list)
    #WRITE_FILE(outfile+'_2_SI_another',d2SI_list)
    #WRITE_FILE(outfile+'_2_MI_another',d2MI_list)
    WRITE_FILE(outfile+'_2_SAS_another',d2SAS_list)
    WRITE_FILE(outfile+'_SAS_result',SAS_family)
#######################################################################
#計算六萬筆的TtoR PSI
def lost_PSI_TtoR_Process(TtoR_document_path,outfile,atom):
    count = 0
    d_list = []
    PSI_list= []
    PSI2_list= []
    analysis=[]
### test
    lost_TtoR_list=[]
    with open('./lost_TtoR_list','r') as file:
        for line in file:
            lost_TtoR_list.append(line)
    for index in range(len(lost_TtoR_list)):
        count+=1
        #   檔案名稱
        current_comparsion=lost_TtoR_list[index].split(' ')
        file_document_name = current_comparsion[0]
        file_document_name2 = current_comparsion[0].split('_to_')[1]+'_to_'+current_comparsion[0].split('_to_')[0]
        print file_document_name
        print file_document_name2
        parent_result = current_comparsion[1]+' '
        class_result = current_comparsion[2]+' '
        ###print file_document_name
        context_length = 0
        max= 0
        align_length = 0
        gap_num_seq1 = 0
        gap_num_seq2 = 0
        compare_list = []
        temp_match = 0
        PSI_times=0
        try:
            with open(TtoR_document_path+file_document_name+'/semiG_result.php','r') as file:
                for each_line in file:
                    context_length+=1
                times = context_length / 7
            for loop in range(0,times):
                PSI_times= 0
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
                        temp_match+=1
                with open(TtoR_document_path+file_document_name+'/profit_log_'+atom+str(loop),'r') as file3:
                    for each_line in file3:
                        try:
                            if(each_line.find('RMS')!=-1):
                                if(num(each_line.split(':')[2].replace(' ',''))<=4.0):
                                    PSI_times+=1
                        except:
                            pass
                print 'hhhhh'
                pdb = Compare_pdb(PSI_times,temp_gap_seq1,temp_gap_seq2,temp_length,temp_match)
                compare_list.append(pdb) 
        except:
            with open(TtoR_document_path+file_document_name2+'/semiG_result.php','r') as file:
                for each_line in file:
                    context_length+=1
                times = context_length / 7
            for loop in range(0,times):
                PSI_times= 0
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
                        temp_match+=1
                with open(TtoR_document_path+file_document_name2+'/profit_log_'+atom+str(loop),'r') as file3:
                    for each_line in file3:
                        try:
                            if(each_line.find('RMS')!=-1):
                                if(num(each_line.split(':')[2].replace(' ',''))<=4.0):
                                    PSI_times+=1
                        except:
                            pass
                pdb = Compare_pdb(PSI_times,temp_gap_seq1,temp_gap_seq2,temp_length,temp_match)
                compare_list.append(pdb) 
        max= compare_list[0].getRMSD()
        gap_num_seq1 = compare_list[0].get_gap1()
        gap_num_seq2 = compare_list[0].get_gap2()
        align_length = compare_list[0].get_align()
        for index in range(1,len(compare_list)):
            if(max< compare_list[index].getRMSD()):
                max= compare_list[index].getRMSD()
                gap_num_seq1 = compare_list[index].get_gap1()
                gap_num_seq2 = compare_list[index].get_gap2()
                align_length = compare_list[index].get_align()
        gap_seq=''
        with open('../1Dseq/'+file_document_name.split('_to_')[0]+'.seq','r') as file:
            gap_seq=file.read()
        gap_seq2=''
        with open('../1Dseq/'+file_document_name.split('_to_')[1]+'.seq','r') as file:
            gap_seq2=file.read()
        try:
            min_seq = MIN(len(gap_seq),len(gap_seq2))
            PSI = float(PSI_times*100)/float(min_seq) 
        except:
            PSI=0
        d_list.append(file_document_name+' pdb_min = '+str(min)+' align= '+str(align_length)+' gpa1='+str(gap_num_seq1)+' gap2='+str(gap_num_seq2)) 
###         print pdb_name+" pdb_min = "+str(min)+' align='+str(align_length)+' gap1 = '+str(gap_num_seq1)+' gap2='+str(gap_num_seq2)
    ### test
    ###
        PSI2_list.append(class_result+str(PSI))
        PSI_list.append(parent_result+str(PSI))
        analysis.append(class_result+str(PSI)+' '+file_document_name)
    #WRITE_FILE(outfile+'_PSI_'+atom+'_log_another_lost',d_list)
    WRITE_FILE(outfile+'_PSI_'+atom+'_another_lost_analysis',analysis)
    WRITE_FILE(outfile+'_0_PSI_'+atom+'_another_lost',PSI_list)
    WRITE_FILE(outfile+'_2_PSI_'+atom+'_another_lost',PSI2_list)
#######################################################################
#計算6萬筆資料用的TtoRcode
def lost_TtoR_Process(TtoR_document_path,outfile):
    count = 0
    d_list = []
    SAS_list= []
    MI_list = []
    SI_list = []
    SAS2_list= []
    MI2_list = []
    SI2_list = []
### test
    lost_TtoR_list=[]
    with open('./lost_TtoR_list','r') as file:
        for line in file:
            lost_TtoR_list.append(line)
    for index in range(len(lost_TtoR_list)):
        count+=1
        #   檔案名稱
        current_comparsion=lost_TtoR_list[index].split(' ')
        file_document_name = current_comparsion[0]
        file_document_name2 = current_comparsion[0].split('_to_')[1]+'_to_'+current_comparsion[0].split('_to_')[0]
        print file_document_name
        print file_document_name2
        parent_result = current_comparsion[1]+' '
        class_result = current_comparsion[2]+' '
        ###print file_document_name
        context_length = 0
        min = 0
        align_length = 0
        gap_num_seq1 = 0
        gap_num_seq2 = 0
        compare_list = []
        temp_match = 0
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
                        temp_match+=1
                with open(TtoR_document_path+file_document_name+'/profit_log'+str(loop),'r') as file3:
                    for each_line in file3:
                        if(each_line.find('RMS')!=-1):
                            temp_RMSD = float(each_line.replace('RMS:',''))
                print str(temp_RMSD)+' 1'
                print 'wtf'
                pdb = Compare_pdb(temp_RMSD,temp_gap_seq1,temp_gap_seq2,temp_length,temp_match)
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
                        temp_match+=1
                with open(TtoR_document_path+file_document_name2+'/profit_log'+str(loop),'r') as file3:
                    for each_line in file3:
                        if(each_line.find('RMS')!=-1):
                            temp_RMSD = float(each_line.replace('RMS:',''))
                print str(temp_RMSD)+' 2'
                pdb = Compare_pdb(temp_RMSD,temp_gap_seq1,temp_gap_seq2,temp_length,temp_match)
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
        d_list.append(file_document_name+' pdb_min = '+str(min)+' align= '+str(align_length)+' gpa1='+str(gap_num_seq1)+' gap2='+str(gap_num_seq2)) 
###         print pdb_name+" pdb_min = "+str(min)+' align='+str(align_length)+' gap1 = '+str(gap_num_seq1)+' gap2='+str(gap_num_seq2)
    ### test
    ###
        SAS2_list.append(class_result+str(SAS))
        SI2_list.append(class_result+str(SI))
        MI2_list.append(class_result+str(MI))
        SAS_list.append(parent_result+str(SAS))
        SI_list.append(parent_result+str(SI))
        MI_list.append(parent_result+str(MI))
    WRITE_FILE(outfile+'_log_another_lost',d_list)
    WRITE_FILE(outfile+'_0_SAS_another_lost',SAS_list)
    WRITE_FILE(outfile+'_2_SAS_another_lost',SAS2_list)
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

    compare_TtoR_list=[]
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
                    pdb = Compare_pdb(temp_RMSD,temp_gap_seq1,temp_gap_seq2,temp_length,0)
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
                    pdb = Compare_pdb(temp_RMSD,temp_gap_seq1,temp_gap_seq2,temp_length,0)
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
            temp_class=pdb_name+' ' 
            SAS2_list.append(result+str(SAS))
            SI2_list.append(result+str(SI))
            MI2_list.append(result+str(MI))
            if(TtoR_list[T_count].get_family()==TtoR_list[R_count].get_family()):
                SAS_list.append('p,p,'+str(SAS))
                SI_list.append('p,p,'+str(SI))
                MI_list.append('p,p,'+str(MI))
                temp_class+='p,p '
            else:
                SAS_list.append('n,p,'+str(SAS))
                SI_list.append('n,p,'+str(SI))
                MI_list.append('n,p,'+str(MI))
                temp_class+='n,p '
            temp_class+=result
            compare_TtoR_list.append(temp_class)
    WRITE_FILE(outfile+'_all_comparsion',compare_TtoR_list)
    WRITE_FILE(outfile+'_log_another',d_list)
    WRITE_FILE(outfile+'_0_SAS_another',SAS_list)
    WRITE_FILE(outfile+'_0_SI_another',SI_list)
    WRITE_FILE(outfile+'_0_MI_another',MI_list)
    WRITE_FILE(outfile+'_2_SAS_another',SAS2_list)
    WRITE_FILE(outfile+'_2_SI_another',SI2_list)
    WRITE_FILE(outfile+'_2_MI_another',MI2_list)
    return 
#######################################################################
#計算PSI TtoR
def PSI_TtoR_Process(TtoR_list,TtoR_document_path,outfile,atom):
    count = 0
    d_list = []
    PSI_list= []
    PSI2_list= []
    d0TPSI_list=[]
    d2TPSI_list=[]
    PSI_times=0
    translate_result_d2=9
### test
    for T_count in range(0,227):
        PSI = 0
        for R_count in range(227,419):
            count+=1
            #   檔案名稱
            file_document_name = TtoR_list[T_count].get_pdb()+'_to_'+TtoR_list[R_count].get_pdb()
            file_document_name2 = TtoR_list[R_count].get_pdb()+'_to_'+TtoR_list[T_count].get_pdb()
            ###print file_document_name
            context_length = 0
            max = 0
            align_length = 0
            gap_num_seq1 = 0
            gap_num_seq2 = 0
            compare_list = []
            PSI_times=0
            try:
                with open(TtoR_document_path+file_document_name+'/semiG_result.php','r') as file:
                    for each_line in file:
                        context_length+=1
                    times = context_length / 7
                for loop in range(0,times):
                    PSI_times = 0
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
                    with open(TtoR_document_path+file_document_name+'/profit_log_'+atom+str(loop),'r') as file3:
                        for each_line in file3:
                            try:
                                if(each_line.find('RMS')!=-1):
                                    if(num(each_line.split(':')[2].replace(' ',''))<=4.0):
                                        PSI_times+=1
                            except:
                                pass
                                
                    pdb = Compare_pdb(PSI_times,temp_gap_seq1,temp_gap_seq2,temp_length,0)
                    compare_list.append(pdb) 
            except:
                with open(TtoR_document_path+file_document_name2+'/semiG_result.php','r') as file:
                    for each_line in file:
                        context_length+=1
                    times = context_length / 7
                for loop in range(0,times):
                    PSI_times= 0
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
                    with open(TtoR_document_path+file_document_name2+'/profit_log_'+atom+str(loop),'r') as file3:
                        for each_line in file3:
                            try:
                                if(each_line.find('RMS')!=-1):
                                    if(num(each_line.split(':')[2].replace(' ',''))<=4.0):
                                        PSI_times+=1
                            except:
                                pass
                    pdb = Compare_pdb(PSI_times,temp_gap_seq1,temp_gap_seq2,temp_length,0)
                    compare_list.append(pdb) 
            max = compare_list[0].getRMSD()
            gap_num_seq1 = compare_list[0].get_gap1()
            gap_num_seq2 = compare_list[0].get_gap2()
            align_length = compare_list[0].get_align()
            for index in range(1,len(compare_list)):
                if(max < compare_list[index].getRMSD()):
                    max = compare_list[index].getRMSD()
                    gap_num_seq1 = compare_list[index].get_gap1()
                    gap_num_seq2 = compare_list[index].get_gap2()
                    align_length = compare_list[index].get_align()

            file_document_name = TtoR_list[T_count].get_pdb()+'_to_'+TtoR_list[R_count].get_pdb()
            gap_seq=''
            with open('../1Dseq/'+TtoR_list[T_count].get_pdb()+'.seq','r') as file:
                gap_seq=file.read()
            gap_seq2=''
            with open('../1Dseq/'+TtoR_list[R_count].get_pdb()+'.seq','r') as file:
                gap_seq2=file.read()
            min_seq = MIN(len(gap_seq),len(gap_seq2))

            try:
                PSI=float(max*100)/float(min_seq)
            except:
                PSI=0

            pdb_name = TtoR_list[T_count].get_pdb()+'_to_'+TtoR_list[R_count].get_pdb()
            pdb_name2 = TtoR_list[R_count].get_pdb()+'_to_'+TtoR_list[T_count].get_pdb()
   ###         print pdb_name+" pdb_min = "+str(min)+' align='+str(align_length)+' gap1 = '+str(gap_num_seq1)+' gap2='+str(gap_num_seq2)
            result = search_family(pdb_name,pdb_name2)
        ### test
        ###
            PSI2_list.append(result+str(PSI))
            if(TtoR_list[T_count].get_family()==TtoR_list[R_count].get_family()):
                d0TPSI_list.append('1 '+str(PSI))
                PSI_list.append('p,p,'+str(PSI))
            else:
                d0TPSI_list.append('0 '+str(PSI))
                PSI_list.append('n,p,'+str(PSI))
            if(result.find('n,p')!=-1):
                translate_result_d2=0
            elif(result.find('p,p')!=-1):
                translate_result_d2=1
            d2TPSI_list.append(str(translate_result_d2)+' '+str(PSI))

    WRITE_FILE(outfile+'_0_PSI_'+atom+'_pyformat',d0TPSI_list)
    WRITE_FILE(outfile+'_2_PSI_'+atom+'_pyformat',d2TPSI_list)
    WRITE_FILE(outfile+'_0_PSI_'+atom+'_another',PSI_list)
    WRITE_FILE(outfile+'_2_PSI_'+atom+'_another',PSI2_list)
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
    result=result.replace('\n','')+' '
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
#計算RASS那部分的FSCOR and TtoR
def RASS_PSI_Process(FSCOR_document_path,outfile,atom,mode):
    #count = 0
    #--------------新增轉換格式-----------------#
    d0TPSI_list=[]
    d2TPSI_list=[]


    dSAS_list=[]
    d2SAS_list = []
    test_PSI_list=[]
    lost_FSCOR_AUC_list = []
    if(mode=='TtoR'):
        path = './lost_TtoR_list'
    elif(mode=='FSCOR'):
        path = './lost_FSCOR_list'
    with open(path,'r') as file: 
    #with open('./lost_TtoR_list','r') as file: 
    #with open('./Empty_Error','r') as file:
        for line in file:
            lost_FSCOR_AUC_list.append(line)
    translate_result=9
    translate_result_d2=9
    for line in lost_FSCOR_AUC_list:
        translate_result=9
        translate_result_d2=9
        compare_result= line.split(' ')[1].replace('\n','')
        compare_pdb_name = line.split(' ')[0].replace('\n','')
        compare_pdb_name2= line.split(' ')[0].split('_to_')[1]+'_to_'+line.split(' ')[0].split('_to_')[0]
        #the result of compare family
        PSI_times=0
        print FSCOR_document_path+compare_pdb_name+'/profit_log'
        try:
            with open(FSCOR_document_path+compare_pdb_name+'/profit_log','r') as file:
                for each_line in file:
                    try:
                        if(each_line.find('RMS')!=-1):
                            if(num(each_line.split(':')[2].replace(' ',''))<=4.0):
                                PSI_times+=1
                    except:
                        pass
        except:
            with open(FSCOR_document_path+compare_pdb_name2+'/profit_log','r') as file:
                for each_line in file:
                    try:
                        if(each_line.find('RMS')!=-1):
                            if(num(each_line.split(':')[2].replace(' ',''))<=4.0):
                                PSI_times+=1
                    except:
                        pass
        gap_seq=''
        with open('../1Dseq/'+compare_pdb_name.split('_to_')[0]+'.seq','r') as file:
            gap_seq=file.read()
        gap_seq2=''
        with open('../1Dseq/'+compare_pdb_name.split('_to_')[1]+'.seq','r') as file:
            gap_seq2=file.read()
        try:
            min_seq=MIN(len(gap_seq),len(gap_seq2))
            PSI = float(PSI_times*100)/float(min_seq) 
        except:
            PSI=0.0
            print 'ERROR'
        #test_PSI_list.append('number of PSI='+str(PSI_times)+' min_seq='+str(min_seq)+' seq1='+str(len(gap_seq))+' seq2='+str(len(gap_seq2))+' PSI='+str(PSI)) 

        #-----轉換為Python ROC 可以跑的格式
        if(compare_result=='n,p'):
            translate_result = 0
        elif(compare_result=='p,p'):
            translate_result = 1 
        d0TPSI_list.append(str(translate_result)+' '+str(PSI))

        dSAS_list.append(compare_result+' '+str(PSI))
        test_PSI_list.append(compare_result+' '+str(PSI)+compare_pdb_name)
        ###d_list.append(pdb_name+' pdb_min = '+str(min)+' align= '+str(align_length)+' gpa1='+str(gap_num_seq1)+' gap2='+str(gap_num_seq2)) 
        result = search_family(compare_pdb_name,compare_pdb_name2)
        #------轉換為d<=2的python roc format
        if(result.find('n,p')!=-1):
            translate_result_d2=0
        elif(result.find('p,p')!=-1):
            translate_result_d2=1
        d2TPSI_list.append(str(translate_result_d2)+' '+str(PSI))

        d2SAS_list.append(result+str(PSI))
    WRITE_FILE(outfile+'_PSI_'+atom+'_another_lost_analysis',test_PSI_list)
    WRITE_FILE(outfile+'_0_PSI_'+atom+'_another_lost',dSAS_list)
    WRITE_FILE(outfile+'_2_PSI_'+atom+'_another_lost',d2SAS_list)
    WRITE_FILE(outfile+'_0_PSI_'+atom+'_pyformat_lost',d0TPSI_list)
    WRITE_FILE(outfile+'_2_PSI_'+atom+'_pyformat_lost',d2TPSI_list)
#######################################################################
def lost_FSCOR_PSI_Process(FSCOR_document_path,outfile,atom):
    #count = 0
    d_list = []
    dPSI_list=[]
    d2PSI_list = []
    d0TPSI_list=[]
    d2TPSI_list=[]
    lost_FSCOR_AUC_list = []
    translate_result=9
    translate_result_d2=9
    analysis=[]
    #someone give me the fail file!
    with open('./lost_FSCOR_list','r') as file: 
    #with open('./Empty_Error','r') as file:
        for line in file:
            lost_FSCOR_AUC_list.append(line)
    for line in lost_FSCOR_AUC_list:
        #the result of compare family
        #-----python format variable
        translate_result=9
        translate_result_d2=9
        compare_result= line.split(' ')[1].replace('\n','')
        compare_pdb_name = line.split(' ')[0].replace('\n','')
        context_length = 0
        max = 0
        align_length = 0
        gap_num_seq1 = 0
        gap_num_seq2 = 0
        compare_list=[]
        seq1=''
        seq2=''
        PSI_times=0
        with open(FSCOR_document_path+compare_pdb_name+'/semiG_result.php','r') as file:
            for each_line in file:
                context_length+=1
            times = context_length / 7
            for loop in range(times):
                PSI_times= 0
                temp_gap_seq1 = 0
                temp_gap_seq2 = 0
                temp_length = 0
                with open(FSCOR_document_path+compare_pdb_name+'/profit_log_'+atom+str(loop),'r') as file:
                    for each_line in file:
                        try:
                            if(each_line.find('RMS')!=-1):
                                if(num(each_line.split(':')[2].replace(' ',''))<=4.0):
                                    PSI_times+=1
                        except:
                            pass

                context_list = []
                with open(FSCOR_document_path+compare_pdb_name+'/ori_ali_seq.pir'+str(loop),'r') as file:
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
                pdb = Compare_pdb(PSI_times,temp_gap_seq1,temp_gap_seq2,temp_length,0)
                compare_list.append(pdb)
        if(len(seq1.replace('\n',''))!=len(seq2.replace('\n',''))):
            print compare_pdb_name+' not equal'+str(len(seq1))+' '+str(len(seq2))
        else:
            print compare_pdb_name+' equal'
        max = compare_list[0].getRMSD()
        gap_num_seq1 = compare_list[0].get_gap1()
        gap_num_seq2 = compare_list[0].get_gap2()
        align_length = compare_list[0].get_align()
        for i in range(1,len(compare_list)):
            if(max < compare_list[i].getRMSD()):
                max = compare_list[i].getRMSD()
                gap_num_seq1 = compare_list[i].get_gap1()
                gap_num_seq2 = compare_list[i].get_gap2()
                align_length = compare_list[i].get_align()

        gap_seq=''
        with open('../1Dseq/'+compare_pdb_name.split('_to_')[0]+'.seq','r') as file:
            gap_seq=file.read()
        gap_seq2=''
        with open('../1Dseq/'+compare_pdb_name.split('_to_')[1]+'.seq','r') as file:
            gap_seq2=file.read()

        try:
            min_seq = MIN(len(gap_seq),len(gap_seq2))                
            PSI = float(max*100)/float(min_seq)
        except:
            PSI = 0

        if(compare_result=='n,p'): 
            translate_result=0
        elif(compare_result=='p,p'):
            translate_result=1
        d0TPSI_list.append(str(translate_result)+' '+str(PSI))
        dPSI_list.append(compare_result+' '+str(PSI))
        analysis.append(compare_result+' '+str(PSI)+' '+compare_pdb_name)

        ###d_list.append(pdb_name+' pdb_min = '+str(min)+' align= '+str(align_length)+' gpa1='+str(gap_num_seq1)+' gap2='+str(gap_num_seq2)) 
        

        result = search_family(compare_pdb_name,compare_pdb_name)
        if(result.find('n,p')!=-1):
            translate_result_d2=0
        elif(result.find('p,p')!=-1):
            translate_result_d2=1
        d2TPSI_list.append(str(translate_result_d2)+' '+str(PSI))

        d_list.append(compare_pdb_name+' '+compare_result+' '+result+' '+str(min)+' '+str(PSI)+' '+str(align_length)+' '+str(gap_num_seq1)+' '+str(gap_num_seq2))
        d2PSI_list.append(result+str(PSI))
    ###WRITE_FILE('center46_FSCOR_0_log',dlist)
    ###WRITE_FILE('center46_FSCOR_2_log',d2list)
    #WRITE_FILE(outfile+'_log_another_lost',d_list)
    WRITE_FILE(outfile+'_PSI_'+atom+'_another_lost_analysis',analysis)
    WRITE_FILE(outfile+'_0_PSI_'+atom+'_another_lost',dPSI_list)
    WRITE_FILE(outfile+'_2_PSI_'+atom+'_another_lost',d2PSI_list)
    WRITE_FILE(outfile+'_0_PSI_'+atom+'_pyformat_lost',dPSI_list)
    WRITE_FILE(outfile+'_2_PSI_'+atom+'_pyformat_lost',d2PSI_list)
#######################################################################
def lost_FSCOR_Process(FSCOR_document_path,outfile):
    #count = 0
    d_list = []
    dSAS_list=[]
    dSI_list = []
    dMI_list = []
    d2SAS_list = []
    d2MI_list = []
    d2SI_list = []
    lost_FSCOR_AUC_list = []
    #someone give me the fail file!
    with open('./lost_FSCOR_list','r') as file: 
    #with open('./Empty_Error','r') as file:
        for line in file:
            lost_FSCOR_AUC_list.append(line)
    for line in lost_FSCOR_AUC_list:
        #the result of compare family
        compare_result= line.split(' ')[1].replace('\n','')
        compare_pdb_name = line.split(' ')[0].replace('\n','')
        context_length = 0
        min = 0
        align_length = 0
        gap_num_seq1 = 0
        gap_num_seq2 = 0
        compare_list=[]
        seq1=''
        seq2=''
        with open(FSCOR_document_path+compare_pdb_name+'/semiG_result.php','r') as file:
            for each_line in file:
                context_length+=1
            times = context_length / 7
            for loop in range(times):
                temp_RMSD = 0
                temp_gap_seq1 = 0
                temp_gap_seq2 = 0
                temp_length = 0
                print FSCOR_document_path+compare_pdb_name+'/profit_log'+str(loop)
                with open(FSCOR_document_path+compare_pdb_name+'/profit_log'+str(loop),'r') as file:
                    for each_line in file:
                        if(each_line.find('RMS')!=-1):
                            print each_line.replace('RMS:','')
                            temp_RMSD = float(each_line.replace('RMS:',''))
                context_list = []
                with open(FSCOR_document_path+compare_pdb_name+'/ori_ali_seq.pir'+str(loop),'r') as file:
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
                pdb = Compare_pdb(temp_RMSD,temp_gap_seq1,temp_gap_seq2,temp_length,0)
                compare_list.append(pdb)
        if(len(seq1.replace('\n',''))!=len(seq2.replace('\n',''))):
            print compare_pdb_name+' not equal'+str(len(seq1))+' '+str(len(seq2))
        else:
            print compare_pdb_name+' equal'
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
        dSAS_list.append(compare_result+' '+str(SAS))
        dSI_list.append(compare_result+' '+str(SI))
        dMI_list.append(compare_result+' '+str(MI))
        ###d_list.append(pdb_name+' pdb_min = '+str(min)+' align= '+str(align_length)+' gpa1='+str(gap_num_seq1)+' gap2='+str(gap_num_seq2)) 
        result = search_family(compare_pdb_name,compare_pdb_name)
        d_list.append(compare_pdb_name+' '+compare_result+' '+result+' '+str(min)+' '+str(SAS)+' '+str(align_length)+' '+str(gap_num_seq1)+' '+str(gap_num_seq2))
        d2SAS_list.append(result+str(SAS))
        d2SI_list.append(result+str(SI))
        d2MI_list.append(result+str(MI))
    ###WRITE_FILE('center46_FSCOR_0_log',dlist)
    ###WRITE_FILE('center46_FSCOR_2_log',d2list)
    WRITE_FILE(outfile+'_log_another_lost',d_list)
    WRITE_FILE(outfile+'_0_SAS_another_lost',dSAS_list)
    WRITE_FILE(outfile+'_0_SI_another_lost',dSI_list)
    WRITE_FILE(outfile+'_0_MI_another_lost',dMI_list)
    WRITE_FILE(outfile+'_2_SI_another_lost',d2SI_list)
    WRITE_FILE(outfile+'_2_MI_another_lost',d2MI_list)
    WRITE_FILE(outfile+'_2_SAS_another_lost',d2SAS_list)
###############################################################
if __name__ =='__main__':
    FSCOR_list = read_file("/home/watchlee/Research_Programming/research_data/inputDataset/SARA_FSCOR.sa")
    TFSCOR_list = read_file("/home/watchlee/Research_Programming/research_data/inputDataset/SARA_T")
    #FSCOR_list = read_file("/home/watchlee/Research_Programming/research_data/inputDataset/SARA_FSCOR_watchlee_under1K")
    ###TEST(FSCOR_list)
    TtoR_list = read_file("/home/watchlee/Research_Programming/research_data/inputDataset/SARA_TtoR-FSCOR.sa")
    ###TEST(TtoR_list)
###FSCOR setting input file
    FSCOR_file = ['/home/watchlee/Research_Programming/RMSD/alignment_main/23-4L_matrix-O6E1.5-SARA_FSCOR_23C_4L_result-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/4L_matrix-O15E1-SARA_FSCOR_4L_result-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/matrix-O4E3.5-FSCOR-semiG.job/'
                 ,"/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/final_matrix.txt-O5E0.5-SARA_FSCOR_over1k_23c-semiG.job/"
                  ,'/home/bingts/harry_before/iPARTS2_training/alignment_main/matrix.txt-O12E1-SARA_FSCOR_over1k_69c-semiG.job/'
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
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter07_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter08_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter09_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter09_matrix-O14E2-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter10_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter10_matrix-O15E2-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter11_matrix-O1E1.5-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter11_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter11_matrix-O1E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter09_matrix-O1E1.5-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter12_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_true_iter09_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_true_iter09_matrix-O1E1.5-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_true_iter10_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_true_iter11_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_true_iter12_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter13_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter14_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_new_23C_4L_matrix-O15E0.5-SARA_FSCOR_new_23C_4L-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/23-4L_matrix-O9E1-SARA_FSCOR_23C_4L_result-semiG.job/']
    FSCOR_output_file = ['23C_4L_FSCOR','4L_23C_FSCOR','46C_FSCOR','23C_FSCOR','69C_FSCOR','iPARTS_FSCOR','5K_46C_K10_FSCOR','iPARTS_FSCOR_old','5K_46C_K30_FSCOR','5K_46C_K60_FSCOR','5K_46C_K90_FSCOR','iter02_23C_FSCOR','iter03_23C_FSCOR','iter04_23C_FSCOR','iter05_23C_FSCOR','iter02_MI_FSCOR','iter06_23C_FSCOR','iter06_MI_FSCOR','iter07_23_FSCOR','iter08_23_FSCOR','iter09_23_FSCOR','iter09_23_MI_FSCOR','iter10_23_FSCOR','iter10_23_MI_FSCOR','iter11_SARA_FSCOR','iter11_23_FSCOR','iter11_SARA2_FSCOR','iter09_SARA_FSCOR','iter12_23_FSCOR','iter09_true_FSCOR','iter09_true_SARA_FSCOR','iter10_true_FSCOR','iter11_true_FSCOR','iter12_true_FSCOR','iter13_23_FSCOR','iter14_23_FSCOR','new_23C_4L_FSCOR','final_23C_4L_FSCOR']
###TtoR setting input file 
    TtoR_output_file = ['23C_4L_TtoR','4L_23C_TtoR','46C_TtoR','23C_TtoR','69C_TtoR','iPARTS_TtoR','5K_46C_K10_TtoR','5K_46C_K30_TtoR','new_23C_TtoR','5K_46C_K60_TtoR','5K_46C_K90_TtoR','iter02_23C_TtoR','iter03_23C_TtoR','iter04_23C_TtoR','iter05_23C_TtoR','iter02_MI_TtoR','iter07_23_TtoR','iter08_23_TtoR','iter09_23_TtoR','iter09_23_MI_TtoR','iter10_23_TtoR','iter10_23_MI_TtoR','iter11_SARA_TtoR','iter11_23_TtoR','iter11_SARA2_TtoR','iter09_SARA_TtoR','iter12_23_TtoR','iter09_true_TtoR','iter09_true_SARA_TtoR','iter10_true_TtoR','iter11_true_TtoR','iter12_true_TtoR','iter13_23_TtoR','iter14_23_TtoR','final_23C_4L_TtoR']
    TtoR_file = ["/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/23-4L_matrix-O6E1.5-SARA_FSCOR_23C_4L_result-semiG.job/"
                 ,"/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/4L_matrix-O15E1-SARA_FSCOR_4L_result-semiG.job/"
                 ,"/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/matrix-O4E3.5-FSCOR-semiG.job/"
                 ,"/home/watchlee/Research_Programming/iPARTS2_training/alignment_main/final_matrix.txt-O5E0.5-SARA_FSCOR_over1k_23c-semiG.job/"
                 ,"/home/bingts/harry_before/iPARTS2_training/alignment_main/matrix.txt-O12E1-SARA_FSCOR_over1k_69c-semiG.job/"
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
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter07_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter08_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter09_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter09_matrix-O14E2-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter10_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter10_matrix-O15E2-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter11_matrix-O1E1.5-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter11_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter11_matrix-O1E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter09_matrix-O1E1.5-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter12_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_true_iter09_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_true_iter09_matrix-O1E1.5-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_true_iter10_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_true_iter11_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_true_iter12_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter13_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/iPARTS2_23_iter14_matrix-O8E1-iPARTS2_23C_SARA_FSCOR-semiG.job/'
                  ,'/home/watchlee/Research_Programming/RMSD/alignment_main/23-4L_matrix-O9E1-SARA_FSCOR_23C_4L_result-semiG.job/']
    F_index =0
    T_index =17  
    TtoR_document_path = TtoR_file[T_index]
    FSCOR_document_path = FSCOR_file[F_index]
    pdbpath = '../pdb/'
    oneDseq_path = '../1Dseq/'
    #TtoR_Process(TtoR_list,TtoR_file[T_index],TtoR_output_file[T_index])
    #FSCOR_Process(FSCOR_list,FSCOR_file[F_index],FSCOR_output_file[F_index])
    #Raw_TtoR_Process(TtoR_list,TtoR_file[T_index],TtoR_output_file[T_index])
    #Raw_FSCOR_Process(FSCOR_list,FSCOR_file[F_index],FSCOR_output_file[F_index])
   # for index in range(len(FSCOR_output_file)):
    #for index in range(0,6):
    #    FSCOR_Process(FSCOR_list,FSCOR_file[index],FSCOR_output_file[index])
    #    TtoR_Process(TtoR_list,TtoR_file[index],TtoR_output_file[index])
    #    Raw_TtoR_Process(TtoR_list,TtoR_file[index],TtoR_output_file[index])
    #    Raw_FSCOR_Process(FSCOR_list,FSCOR_file[index],FSCOR_output_file[index])
    #lost_FSCOR_Process(FSCOR_file[0],FSCOR_output_file[0])
    #lost_FSCOR_Process(FSCOR_file[3],FSCOR_output_file[3])
    #lost_FSCOR_Process(FSCOR_file[5],FSCOR_output_file[5])
    #lost_FSCOR_Process(FSCOR_file[-1],FSCOR_output_file[-1])
    #FSCOR_Process(FSCOR_list,FSCOR_file[-1],FSCOR_output_file[-1])
    #FSCOR_Process(FSCOR_list,FSCOR_file[0],FSCOR_output_file[0])
    #FSCOR_Process(FSCOR_list,FSCOR_file[3],FSCOR_output_file[3])
    #FSCOR_Process(FSCOR_list,FSCOR_file[5],FSCOR_output_file[5])
    #Raw_FSCOR_Process(FSCOR_list,FSCOR_file[0],FSCOR_output_file[0])
    #Raw_FSCOR_Process(FSCOR_list,FSCOR_file[5],FSCOR_output_file[5])
    #FSCOR_Process(FSCOR_list,FSCOR_file[3],FSCOR_output_file[3])
    #Raw_FSCOR_Process(FSCOR_list,FSCOR_file[3],FSCOR_output_file[3])
    #TtoR_Process(TtoR_list,TtoR_file[-1],TtoR_output_file[-1])
    #TtoR_Process(TtoR_list,TtoR_file[0],TtoR_output_file[0])
    #TtoR_Process(TtoR_list,TtoR_file[3],TtoR_output_file[3])
    #TtoR_Process(TtoR_list,TtoR_file[5],TtoR_output_file[5])
    #Raw_TtoR_Process(TtoR_list,TtoR_file[0],TtoR_output_file[0])
    #Raw_TtoR_Process(TtoR_list,TtoR_file[3],TtoR_output_file[3])
    #Raw_TtoR_Process(TtoR_list,TtoR_file[5],TtoR_output_file[5])
    

    #-------------someone give me that shit data, which he don't want to process by himself. So I fucked it up.
    #other_PSI_Process(FSCOR_list,'/home/millard/iPARTS2/alignment/setter/SARA_FSCOR_PSI_ALLATOM/','SETTER_FSCOR')
    #other_PSI_Process(FSCOR_list,'/home/millard/iPARTS2/alignment/setter/SARA_FSCOR_253/','SETTER_FSCOR')
    #-------------fuck fuck fuck fuck fuck  ah....this code is used to run RASS data
    RASS_PSI_Process('/home/millard/iPARTS2/alignment/RASS/SARA_FSCOR_PSI_ALL/','RASS_TtoR','S','TtoR')
    RASS_PSI_Process('/home/millard/iPARTS2/alignment/RASS/SARA_FSCOR_PSI_C3/','RASS_TtoR','C3','TtoR')
    RASS_PSI_Process('/home/millard/iPARTS2/alignment/RASS/SARA_FSCOR_PSI_P/','RASS_TtoR','P','TtoR')
    RASS_PSI_Process('/home/millard/iPARTS2/alignment/RASS/SARA_FSCOR_PSI_ALL/','RASS_FSCOR','S','FSCOR')
    RASS_PSI_Process('/home/millard/iPARTS2/alignment/RASS/SARA_FSCOR_PSI_C3/','RASS_FSCOR','C3','FSCOR')
    RASS_PSI_Process('/home/millard/iPARTS2/alignment/RASS/SARA_FSCOR_PSI_P/','RASS_FSCOR','P','FSCOR')
    RASS_PSI_Process('/home/millard/iPARTS2/alignment/setter/SARA_FSCOR_PSI_ALLATOM/','SETTER_TtoR','S','TtoR')
    RASS_PSI_Process('/home/millard/iPARTS2/alignment/setter/SARA_FSCOR_PSI_ATOM_C3/','SETTER_TtoR','C3','TtoR')
    RASS_PSI_Process('/home/millard/iPARTS2/alignment/setter/SARA_FSCOR_PSI_ATOM_P/','SETTER_TtoR','P','TtoR')
    RASS_PSI_Process('/home/millard/iPARTS2/alignment/setter/SARA_FSCOR_PSI_ALLATOM/','SETTER_FSCOR','S','FSCOR')
    RASS_PSI_Process('/home/millard/iPARTS2/alignment/setter/SARA_FSCOR_PSI_ATOM_C3/','SETTER_FSCOR','C3','FSCOR')
    RASS_PSI_Process('/home/millard/iPARTS2/alignment/setter/SARA_FSCOR_PSI_ATOM_P/','SETTER_FSCOR','P','FSCOR')

    #FORCOMMENT_TFSCOR_Process(TFSCOR_list,'/home/watchlee/Research_Programming/RMSD/alignment_main/comment_23-4L-O15E5-SARA_T_23-4L-semiG.job/','23-4L_TFSCOR')
    #FORCOMMENT_TFSCOR_Process(TFSCOR_list,'/home/watchlee/Research_Programming/RMSD/alignment_main/comment_iPARTS-O15E5-SARA_T_iPARTS-semiG.job/','iPARTS_TFSCOR')
    #lost_TtoR_Process(TtoR_file[-1],TtoR_output_file[-1])
    #PSI_Process(FSCOR_list,FSCOR_file[5],FSCOR_output_file[5],'C3')
    #PSI_Process(FSCOR_list,FSCOR_file[5],FSCOR_output_file[5],'P')
    #PSI_Process(FSCOR_list,FSCOR_file[5],FSCOR_output_file[5],'S')
    '''
    PSI_TtoR_Process(TtoR_list,TtoR_file[3],TtoR_output_file[3],'C3')
    PSI_TtoR_Process(TtoR_list,TtoR_file[3],TtoR_output_file[3],'P')
    PSI_TtoR_Process(TtoR_list,TtoR_file[3],TtoR_output_file[3],'S')
    '''
    lost_PSI_TtoR_Process(TtoR_file[5],TtoR_output_file[5],'C3')
    lost_PSI_TtoR_Process(TtoR_file[5],TtoR_output_file[5],'P')
    lost_PSI_TtoR_Process(TtoR_file[5],TtoR_output_file[5],'S')
    lost_FSCOR_PSI_Process(FSCOR_file[5],FSCOR_output_file[5],'C3')
    lost_FSCOR_PSI_Process(FSCOR_file[5],FSCOR_output_file[5],'S')
    lost_FSCOR_PSI_Process(FSCOR_file[5],FSCOR_output_file[5],'P')
    lost_PSI_TtoR_Process(TtoR_file[-1],TtoR_output_file[-1],'C3')
    lost_PSI_TtoR_Process(TtoR_file[-1],TtoR_output_file[-1],'P')
    lost_PSI_TtoR_Process(TtoR_file[-1],TtoR_output_file[-1],'S')
    lost_FSCOR_PSI_Process(FSCOR_file[-1],FSCOR_output_file[-1],'C3')
    lost_FSCOR_PSI_Process(FSCOR_file[-1],FSCOR_output_file[-1],'S')
    lost_FSCOR_PSI_Process(FSCOR_file[-1],FSCOR_output_file[-1],'P')
    '''
    #lost_TtoR_Process(TtoR_file[5],TtoR_output_file[5])
    #lost_TtoR_Process(TtoR_file[-1],TtoR_output_file[-1])

    '''

    '''
    處理完整跟不完整(6萬筆,跟廢物RASS比較用的)資料有兩種方法可以用 
    1.直接讀取每一筆資料夾的資料去處理 <----------要一直讀取檔案 浪費時間
    2.從已經產生好的(完整8萬筆)結果中去parsing 來處理 <------------一瞬間就好了 廢話只是處理字串不必一直做I/O處理
    '''

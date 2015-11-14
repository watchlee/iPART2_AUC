#!/usr/bin/env python
# coding=utf-8
import sys
import subprocess

def calculate_AUC():
    list = ['./23C'
            ,'./23C_4L'
            ,'./46C'
            ,'./69C'
            ,'./iPARTS'
            ,'./4L_23C'
            ,'./5K_46C_K10'
            ,'./5K_46C_K30'
            ,'./5K_46C_K60'
            ,'./5K_46C_K90'
            ,'./iter02_23C']
    Research_information = '_FSCOR_'
    distanceZERO = '0_'
    distanceTWO = '2_'
    another = '_another'
    output_file = ' > FSCOR_result'
    output_path = 'FSCOR_result'
    Cal_list = ['SAS','SI','MI','RAW']
    result_list = []
    for line in Cal_list:
        for cal_line in list:
            cmd = 'ROC.sh '+cal_line+Research_information+distanceZERO+line+output_file
            subprocess.call(cmd,shell=True)
            result_list.append(cmd)
            with open(output_path,'r') as file:
                for line2 in file:
                    if(line2.find('AUC')!=-1 or line2.find('ROC p')!=-1):
                        result_list.append(line2)
                    
    for line in Cal_list:
        for cal_line in list:
            cmd = 'ROC.sh '+cal_line+Research_information+distanceZERO+line+another+output_file
            subprocess.call(cmd,shell=True)
            result_list.append(cmd)
            with open(output_path,'r') as file:
                for line2 in file:
                    if(line2.find('AUC')!=-1 or line2.find('ROC p')!=-1):
                        result_list.append(line2)
    for line in Cal_list:
        for cal_line in list:
            cmd ='ROC.sh '+cal_line+Research_information+distanceTWO+line+output_file
            subprocess.call(cmd,shell=True)
            result_list.append(cmd)
            with open(output_path,'r') as file:
                for line2 in file:
                    if(line2.find('AUC')!=-1 or line2.find('ROC p')!=-1):
                        result_list.append(line2)
    for line in Cal_list:
        for cal_line in list:
            cmd='ROC.sh '+cal_line+Research_information+distanceTWO+line+another+output_file
            subprocess.call(cmd,shell=True)
            result_list.append(cmd)
            with open(output_path,'r') as file:
                for line2 in file:
                    if(line2.find('AUC')!=-1 or line2.find('ROC p')!=-1):
                        result_list.append(line2)

    with open(output_path,'w') as file:
        for line in result_list:
            file.write(line.replace('\n','')+'\n')

if __name__=='__main__':
    try:
        list = []
        value=['_0','_2']
        type = ['_SAS','_SAS_another','_SI','_SI_another','_MI','_MI_another','_RAW','_RAW_another']
        for count in value:
            for count2 in type:
                list.append(str(sys.argv[1])+count+count2)
                cmd = 'ROC.sh '+str(sys.argv[1])+count+count2+' > temp_AUC'
                print 'input file is '+str(sys.argv[1])+count+count2
                subprocess.call(cmd,shell=True)
                try:
                    print 'output file is '+str(sys.argv[2])
                    with open('temp_AUC','r') as file:
                        for line in file:
                            if(line.find('AUC')!=-1 or line.find('ROC p')!=-1):
                                list.append(line)
                except:
                    with open('temp_AUC','r') as file:
                        for line in file:
                            if(line.find('AUC')!=-1 or line.find('ROC p')!=-1):
                                print line.replace('\n','')
        try : 
            with open(str(sys.argv[2]),'w') as file:
                for line in list:
                    if(line.find('AUC')!=-1):
                        file.write(line.replace('\n','')+'\n\n')
                    else:
                        file.write(line.replace('\n','')+'\n')
        except:
            pass
        subprocess.call('rm temp_AUC',shell=True) 


    except:
        print 'defualt mode'
        calculate_AUC()

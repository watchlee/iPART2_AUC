#!/usr/bin/env python
# coding=utf-8

import subprocess

list = ['./23C','./23C_4L','./46C','./69C','./iPARTS','./4L_23C','./5K_46C_K10','./5K_46C_K30','./5K_46C_K60','./5K_46C_K90']

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

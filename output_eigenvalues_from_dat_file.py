# -*- coding: cp936 -*-
##��ȡ������������ֵ##


#����Abaqus�е�job��
datname = raw_input('''����Abaqus�е�job������������׺����
��Ĭ��·��Ϊ��E:\Study\ABAQUS\Temp\����''')


#�򿪽��.dat�ļ�
datfile = open('E:\Study\ABAQUS\Temp\%s.dat'%datname)


#����.dat�ļ�
for line in datfile:
    #�жϼ���������ֵ�����
    if line == ' MODE NO      EIGENVALUE\n':
        i = 1
        print '\n\n MODE NO      EIGENVALUE'
        #���ǰ5������ֵ
        while i <= 7:
			print datfile.next(),
			i = i + 1


print '''\n\nThese are all Eigenvalues'''

#�ر�.dat�ļ�
datfile.close()


input('�밴Enter�����ҹرճ���')
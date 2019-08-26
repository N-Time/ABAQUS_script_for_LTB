# -*- coding: cp936 -*-


##������Ϣ
#���빹����Ϣ
instanceName = raw_input('������Instanceʵ������ƣ�')
l = input('�����빹�����ȣ�����һλС������')
nL1, nL2, nL3, nL4, nL5 = map(int,raw_input('''����������
����Ԫ������
�߽�ڵ���a���߽�ڵ���b��
�Ǳ߽�ڵ���ʼ��š��Ǳ߽�ڵ���ĩ���
���Զ��Ż��߿ո��������
''').split())
#������طֲ�������Ϣ
p0, p1, p2, p3 = map(float,raw_input('''����������
�ֲ���������֧��3�ζ���ʽ����ϵ��
���ӵʹ��ʼ����4��ϵ����Ҫ���룩��
''').split())


##����INP�ļ�
inpfile = file('E:\Study\ABAQUS\Temp\Script\undistribution loads.txt','w')


##�����ڵ����б�
#��Ӧ����������Ԫ�������߽�ڵ���a���߽�ڵ���b���Ǳ߽�ڵ���ʼ��š��Ǳ߽�ڵ���ĩ���
def nodesList(elementsSum,a,b,nodesStart,nodesEnd):
    #����ڵ�����
    nodesSum = elementsSum - 2
    #�ڵ���ʼ����ĩ��ţ��Լ�����
    nodesStep = (nodesEnd - nodesStart)/nodesSum
    #��������λ�ýڵ��б�
    nodesList = range(nodesStart,nodesEnd + nodesStep,nodesStep)
    #���ӱ߽�ڵ��ţ�a, b��
    nodesList.insert(0,a)
    nodesList.append(b)
    return nodesList


inpfile.write('�������������ڵ㼯�ϡ�����������\n')
##д��
nodesSet = nodesList(nL1,nL2,nL3,nL4,nL5)
for node in nodesSet:
    inpfile.write('*Nset, nset=DL%d\n %d,\n'%(node,node))
print '�ڵ㼯��д����ϣ�'


##���ش�С
#���طֲ�����
def Dload(z,dx):
    #���طֲ�����
    d = p0 + p1 * z + p2 * z**2 + p3 * z**3
    #�ڵ���ش�С
    f = d * dx
    return f


inpfile.write('\n\n���������������ء�����������\n')
inpfile.write('**\n** LOADS\n**\n')
##д�����
#�ڵ���ػ��ֿ��
dl = l / nL1
for i in range(nL1 + 1):
    load = Dload(i*dl,dl)
    inpfile.write('** Name: DL%d   Type: Concentrated force\n'%nodesSet[i])
    inpfile.write('*Cload\n')
    #Ĭ�Ϻ��ط���Ϊ��
    inpfile.write('%s.DL%d, 2, -%f\n'%(instanceName,nodesSet[i],load))
print '��������д����ϣ�'


##�ر�INP�ļ�
inpfile.close()

input('�밴Enter�����ҹرճ���')

# -*- coding: cp936 -*-


##输入信息
#输入构件信息
instanceName = raw_input('请输入Instance实体的名称：')
l = input('请输入构件长度（保留一位小数）：')
nL1, nL2, nL3, nL4, nL5 = map(int,raw_input('''请依次输入
纵向单元数量、
边界节点编号a、边界节点编号b、
非边界节点起始编号、非边界节点终末编号
（以逗号或者空格隔开）：
''').split())
#输入荷载分布函数信息
p0, p1, p2, p3 = map(float,raw_input('''请依次输入
分布函数（仅支持3次多项式）的系数
（从低次项开始，且4个系数都要输入）：
''').split())


##创建INP文件
inpfile = file('E:\Study\ABAQUS\Temp\Script\undistribution loads.txt','w')


##创建节点结合列表
#对应变量：纵向单元数量、边界节点编号a、边界节点编号b、非边界节点起始编号、非边界节点终末编号
def nodesList(elementsSum,a,b,nodesStart,nodesEnd):
    #纵向节点数量
    nodesSum = elementsSum - 2
    #节点起始与终末编号，以及步长
    nodesStep = (nodesEnd - nodesStart)/nodesSum
    #荷载作用位置节点列表
    nodesList = range(nodesStart,nodesEnd + nodesStep,nodesStep)
    #增加边界节点编号（a, b）
    nodesList.insert(0,a)
    nodesList.append(b)
    return nodesList


inpfile.write('——————节点集合——————\n')
##写入
nodesSet = nodesList(nL1,nL2,nL3,nL4,nL5)
for node in nodesSet:
    inpfile.write('*Nset, nset=DL%d\n %d,\n'%(node,node))
print '节点集合写入完毕！'


##荷载大小
#荷载分布函数
def Dload(z,dx):
    #荷载分布函数
    d = p0 + p1 * z + p2 * z**2 + p3 * z**3
    #节点荷载大小
    f = d * dx
    return f


inpfile.write('\n\n——————荷载——————\n')
inpfile.write('**\n** LOADS\n**\n')
##写入荷载
#节点荷载划分宽度
dl = l / nL1
for i in range(nL1 + 1):
    load = Dload(i*dl,dl)
    inpfile.write('** Name: DL%d   Type: Concentrated force\n'%nodesSet[i])
    inpfile.write('*Cload\n')
    #默认荷载方向为负
    inpfile.write('%s.DL%d, 2, -%f\n'%(instanceName,nodesSet[i],load))
print '荷载作用写入完毕！'


##关闭INP文件
inpfile.close()

input('请按Enter结束且关闭程序')

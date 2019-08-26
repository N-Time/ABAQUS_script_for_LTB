# -*- coding: cp936 -*-
##读取弹性屈曲特征值##


#输入Abaqus中的job名
datname = raw_input('''输入Abaqus中的job名（不包括后缀名）
（默认路径为：E:\Study\ABAQUS\Temp\）：''')


#打开结果.dat文件
datfile = open('E:\Study\ABAQUS\Temp\%s.dat'%datname)


#检索.dat文件
for line in datfile:
    #判断检索至特征值结果行
    if line == ' MODE NO      EIGENVALUE\n':
        i = 1
        print '\n\n MODE NO      EIGENVALUE'
        #输出前5个特征值
        while i <= 7:
			print datfile.next(),
			i = i + 1


print '''\n\nThese are all Eigenvalues'''

#关闭.dat文件
datfile.close()


input('请按Enter结束且关闭程序')
# -*- coding: mbcs -*-
###ACIB v0.4.1 SSB##
###ACIB is AutoConstructe Model of I-Beams for Simple Span Beams
##Update Log
#v0.2 for simple span beams
#适用于单跨梁，包括两端简支、两端固支和悬臂梁
#改用findAt(...)方法寻找对象，而非通用性较弱的getSequenceFromMask(...)方法
#v0.3 for continuous
#连续梁专用版本，仅适用于两跨连续梁
###v0.4 for simple span beams
##适用于单跨梁，包括两端简支、两端固支和悬臂梁（平面内外边界条件相同）。
##增加自动施加均布荷载、跨中集中荷载以及1/3跨两个对称荷载，且荷载作用位置包括上、剪、下。
##在输入窗口中增加必要的文字说明。由于窗口不支持中文，故暂用英文说明。
# Do not delete the following import lines
##导入模块
from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior


##建模函数
def ConstructeIBeams(E,u,b1,b2,s1,s2,tf1,tf2,tw,componentL,BC,ELM,tm,bm,wm,lm):

	##参数说明
	#E = 弹性模量（N/m^2）
	#u = 泊松比
	#s1, s2 = 剪心至上翼缘上缘和下翼缘下缘的距离（mm）
	#b1, b2 = 上翼缘宽度和下翼缘宽度（mm）
	#tf1, tf2, tw = 上翼缘厚度，下翼缘厚度，腹板厚度（mm）
	#componentL = 构件长度（m）
	#A, B = A侧（原点侧）边界条件，B侧（非原点侧）边界条件：
	#边界条件（平面内外）：1、简支：s-s；2、固支：c-c；3、悬臂：cantilever
	
	
	##计算中线节点的坐标(m)
	bf1=b1/2
	bf2=b2/2
	hs1=s1-tf1/2
	hs2=s2-tf2/2

	
	##截面中线轮廓Profile（以剪心为坐标原点）
	s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=5000.0)
	g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
	s.setPrimaryObject(option=STANDALONE)
	#上翼缘中线坐标
	s.Line(point1=(-bf1, hs1), point2=(bf1, hs1))
	s.HorizontalConstraint(entity=g[2], addUndoState=False)
	#腹板中线坐标
	s.Line(point1=(0.0, hs1), point2=(0.0, -hs2))
	s.VerticalConstraint(entity=g[3], addUndoState=False)
	s.PerpendicularConstraint(entity1=g[2], entity2=g[3], addUndoState=False)
	s.CoincidentConstraint(entity1=v[2], entity2=g[2], addUndoState=False)
	s.EqualDistanceConstraint(entity1=v[0], entity2=v[1], midpoint=v[2], 
        addUndoState=False)
	#下翼缘中线坐标
	s.Line(point1=(-bf2, -hs2), point2=(bf2, -hs2))
	s.HorizontalConstraint(entity=g[4], addUndoState=False)
	
	
	##创建部件Part
	p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
	p = mdb.models['Model-1'].parts['Part-1']
	#拉伸截面depth长
	p.BaseShellExtrude(sketch=s, depth=componentL)
	s.unsetPrimaryObject()
	p = mdb.models['Model-1'].parts['Part-1']
	del mdb.models['Model-1'].sketches['__profile__']
	
	
	##创建材料Material（钢材E=2.06e5 Mpa，u=0.3）
	mdb.models['Model-1'].Material(name='steel')
	mdb.models['Model-1'].materials['steel'].Elastic(table=((E, u), 
        ))
	
	
	##创建截面属性Section
	#上翼缘
	mdb.models['Model-1'].HomogeneousShellSection(name='topflange', preIntegrate=OFF, 
        material='steel', thicknessType=UNIFORM, thickness=tf1, 
        thicknessField='', idealization=NO_IDEALIZATION, 
        poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
        useDensity=OFF, integrationRule=SIMPSON, numIntPts=5)
	#腹板
	mdb.models['Model-1'].HomogeneousShellSection(name='web', preIntegrate=OFF, 
        material='steel', thicknessType=UNIFORM, thickness=tw, 
        thicknessField='', idealization=NO_IDEALIZATION, 
        poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
        useDensity=OFF, integrationRule=SIMPSON, numIntPts=5)
	#下翼缘
	mdb.models['Model-1'].HomogeneousShellSection(name='bottomflange', preIntegrate=OFF, 
        material='steel', thicknessType=UNIFORM, thickness=tf2, 
        thicknessField='', idealization=NO_IDEALIZATION, 
        poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
        useDensity=OFF, integrationRule=SIMPSON, numIntPts=5)	
    
	
	##赋予截面
	#上翼缘
	p = mdb.models['Model-1'].parts['Part-1']
	f = p.faces
	faces = f.findAt(((-bf1,hs1,0.0),),((bf1,hs1,0.0),))
	region = p.Set(faces=faces, name='Set-1')
	p = mdb.models['Model-1'].parts['Part-1']
	p.SectionAssignment(region=region, sectionName='topflange', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
	#腹板
	p = mdb.models['Model-1'].parts['Part-1']
	f = p.faces
	faces = f.findAt(((0.0,0.0,0.0),))
	region = p.Set(faces=faces, name='Set-2')
	p = mdb.models['Model-1'].parts['Part-1']
	p.SectionAssignment(region=region, sectionName='web', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
	#下翼缘
	p = mdb.models['Model-1'].parts['Part-1']
	f = p.faces
	faces = f.findAt(((-bf2,-hs2,0.0),),((bf2,-hs2,0.0),))
	region = p.Set(faces=faces, name='Set-3')
	p = mdb.models['Model-1'].parts['Part-1']
	p.SectionAssignment(region=region, sectionName='bottomflange', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
	
	
	##创建1/2, 1/3, 2/3切面
	#1/3
	#建立参考面
	p = mdb.models['Model-1'].parts['Part-1']
	dp13 = p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=componentL/3)
	dp13ID = dp13.id
	#分割构件
	p = mdb.models['Model-1'].parts['Part-1']
	f = p.faces
	pickedFaces = f.findAt(((-bf1,hs1,0.0),),((bf1,hs1,0.0),),
		((0.0,0.0,0.0),),((-bf2,-hs2,0.0),),((bf2,-hs2,0.0),))
	d13 = p.datums
	p.PartitionFaceByDatumPlane(datumPlane=d13[dp13ID], faces=pickedFaces)
	#2/3
	#建立参考面
	p = mdb.models['Model-1'].parts['Part-1']
	dp23 = p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=componentL/3*2)
	dp23ID = dp23.id
	#分割构件
	p = mdb.models['Model-1'].parts['Part-1']
	f = p.faces
	pickedFaces = f.findAt(((-bf1,hs1,componentL),),((bf1,hs1,componentL),),
		((0.0,0.0,componentL),),((-bf2,-hs2,componentL),),((bf2,-hs2,componentL),))
	d23 = p.datums
	p.PartitionFaceByDatumPlane(datumPlane=d23[dp23ID], faces=pickedFaces)
	#1/2
	#建立参考面
	p = mdb.models['Model-1'].parts['Part-1']
	dp12 = p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=componentL/2)
	dp12ID = dp12.id
	#分割构件
	p = mdb.models['Model-1'].parts['Part-1']
	f = p.faces
	pickedFaces = f.findAt(((-bf1,hs1,componentL/2),),((bf1,hs1,componentL/2),),
		((0.0,0.0,componentL/2),),((-bf2,-hs2,componentL/2),),((bf2,-hs2,componentL/2),))
	d12 = p.datums
	p.PartitionFaceByDatumPlane(datumPlane=d12[dp12ID], faces=pickedFaces)
	#shear
	#建立参考面
	p = mdb.models['Model-1'].parts['Part-1']
	dps = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
	dpsID = dps.id
	#分割构件
	p = mdb.models['Model-1'].parts['Part-1']
	f = p.faces
	pickedFaces = f.findAt(((0.0,0.0,0.0),),((0.0,0.0,componentL),),
		((0.0,0.0,componentL/3+0.0001),),((0.0,0.0,componentL/3*2-0.0001),))
	ds = p.datums
	p.PartitionFaceByDatumPlane(datumPlane=ds[dpsID], faces=pickedFaces)
	
	
	##创建边界条件作用域集合（几何型集合）
	#A侧（原点侧）
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	edges = e.findAt(((-bf1+0.0001,hs1,0.0),),((bf1-0.0001,hs1,0.0),),
		((-bf2+0.0001,-hs2,0.0),),((bf2-0.0001,-hs2,0.0),),
		((0.0,0.0001,0.0),),((0.0,-0.0001,0.0),))
	p.Set(edges=edges, name='A')
	#B侧（非原点侧）
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	edges = e.findAt(((-bf1+0.0001,hs1,componentL),),((bf1-0.0001,hs1,componentL),),
		((-bf2+0.0001,-hs2,componentL),),((bf2-0.0001,-hs2,componentL),),
		((0.0,0.0001,componentL),),((0.0,-0.0001,componentL),))
	p.Set(edges=edges, name='B')
	#A侧（约束轴向位移）
	p = mdb.models['Model-1'].parts['Part-1']
	v = p.vertices
	verts = v.findAt(((0.0,0.0,0.0),))
	p.Set(vertices=verts, name='AS')
	
	
	##布置种子Seed以及网格划分Mesh
	htm = tm/2
	hbm = bm/2
	#上翼缘
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	pickedEdges = e.findAt(((-bf1+0.0001,hs1,0.0),),((bf1-0.0001,hs1,0.0),),
		((-bf1+0.0001,hs1,componentL/3),),((bf1-0.0001,hs1,componentL/3),),
		((-bf1+0.0001,hs1,componentL/2),),((bf1-0.0001,hs1,componentL/2),),
		((-bf1+0.0001,hs1,componentL/3*2),),((bf1-0.0001,hs1,componentL/3*2),),
		((-bf1+0.0001,hs1,componentL),),((bf1-0.0001,hs1,componentL),))
	p.seedEdgeByNumber(edges=pickedEdges, number=htm, constraint=FINER)
	#下翼缘
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	pickedEdges = e.findAt(((-bf2+0.0001,-hs2,0.0),),((bf2-0.0001,-hs2,0.0),),
		((-bf2+0.0001,-hs2,componentL/3),),((bf2-0.0001,-hs2,componentL/3),),
		((-bf2+0.0001,-hs2,componentL/2),),((bf2-0.0001,-hs2,componentL/2),),
		((-bf2+0.0001,-hs2,componentL/3*2),),((bf2-0.0001,-hs2,componentL/3*2),),
		((-bf2+0.0001,-hs2,componentL),),((bf2-0.0001,-hs2,componentL),))
	p.seedEdgeByNumber(edges=pickedEdges, number=hbm, constraint=FINER)
	#腹板
	wmt = int((s1-tf1/2)/(s1+s2-tf1/2-tf2/2)*wm+0.5)
	wmb = int((s2-tf2/2)/(s1+s2-tf1/2-tf2/2)*wm+0.5)
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	pickedEdges = e.findAt(((0.0,0.0001,0.0),),((0.0,0.0001,componentL/3),),
		((0.0,0.0001,componentL/2),),((0.0,0.0001,componentL/3*2),),((0.0,0.0001,componentL),))
	p.seedEdgeByNumber(edges=pickedEdges, number=wmt, constraint=FINER)
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	pickedEdges = e.findAt(((0.0,-0.0001,0.0),),((0.0,-0.0001,componentL/3),),
		((0.0,-0.0001,componentL/2),),((0.0,-0.0001,componentL/3*2),),((0.0,-0.0001,componentL),))
	p.seedEdgeByNumber(edges=pickedEdges, number=wmb, constraint=FINER)
	#纵向13
	lm13 = int(lm/3+0.5)
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	pickedEdges = e.findAt(((-bf1,hs1,0.0001),),((bf1,hs1,0.0001),),
		((0.0,hs1,0.0001),),((0.0,-hs2,0.0001),),
		((0.0,0.0,0.0001),),
		((-bf2,-hs2,0.0001),),((bf2,-hs2,0.0001),),
		((-bf1,hs1,componentL-0.0001),),((bf1,hs1,componentL-0.0001),),
		((0.0,hs1,componentL-0.0001),),((0.0,-hs2,componentL-0.0001),),
		((0.0,0.0,componentL-0.0001),),
		((-bf2,-hs2,componentL-0.0001),),((bf2,-hs2,componentL-0.0001),))
	p.seedEdgeByNumber(edges=pickedEdges, number=lm13, constraint=FINER)
	#纵向12和13间
	lm23 = int(lm/6+0.5)
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	pickedEdges = e.findAt(((-bf1,hs1,componentL/3+0.0001),),((bf1,hs1,componentL/3+0.0001),),
		((0.0,hs1,componentL/3+0.0001),),((0.0,-hs2,componentL/3+0.0001),),
		((0.0,0.0,componentL/3+0.0001),),
		((-bf2,-hs2,componentL/3+0.0001),),((bf2,-hs2,componentL/3+0.0001),),
		((-bf1,hs1,componentL/3*2-0.0001),),((bf1,hs1,componentL/3*2-0.0001),),
		((0.0,0.0,componentL/3*2-0.0001),),
		((0.0,hs1,componentL/3*2-0.0001),),((0.0,-hs2,componentL/3*2-0.0001),),
		((-bf2,-hs2,componentL/3*2-0.0001),),((bf2,-hs2,componentL/3*2-0.0001),))
	p.seedEdgeByNumber(edges=pickedEdges, number=lm23, constraint=FINER)
	#根据种子划分网格Mesh
	p = mdb.models['Model-1'].parts['Part-1']
	p.generateMesh()
	
	
	##单元类型Element
	if ELM == 'S8R5':
		#四边形网格单元类型
		elemType1 = mesh.ElemType(elemCode=S8R5, elemLibrary=STANDARD)
		#三角形网格单元类型
		elemType2 = mesh.ElemType(elemCode=STRI65, elemLibrary=STANDARD)
	elif ELM == 'S4R5':
		#四边形网格单元类型
		elemType1 = mesh.ElemType(elemCode=S4R5, elemLibrary=STANDARD)
		#三角形网格单元类型
		elemType2 = mesh.ElemType(elemCode=S3, elemLibrary=STANDARD)
	elif ELM == 'S4R':
		#四边形网格单元类型
		elemType1 = mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD)
		#三角形网格单元类型
		elemType2 = mesh.ElemType(elemCode=S3, elemLibrary=STANDARD)
	else:
		print 'Element type is not defined!'
	p = mdb.models['Model-1'].parts['Part-1']
	faces = p.faces
	pickedRegions =(faces, )
	p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
	
	
	##创建荷载作用点集合
	#2分点荷载集合
	p = mdb.models['Model-1'].parts['Part-1']
	v = p.vertices
	span12top = v.findAt(((0.0,hs1,componentL/2),))
	p.Set(name = 'span12top',vertices = span12top)
	span12shear = v.findAt(((0.0,0.0,componentL/2),))
	p.Set(name = 'span12shear',vertices = span12shear)
	span12bottom = v.findAt(((0.0,-hs2,componentL/2),))
	p.Set(name = 'span12bottom',vertices = span12bottom)
	#3分点对称荷载集合
	p = mdb.models['Model-1'].parts['Part-1']
	v = p.vertices
	span13top = v.findAt(((0.0,hs1,componentL/3),),((0.0,hs1,componentL/3*2),))
	p.Set(name = 'span13top',vertices = span13top)
	span13shear = v.findAt(((0.0,0.0,componentL/3),),((0.0,0.0,componentL/3*2),))
	p.Set(name = 'span13shear',vertices = span13shear)
	span13bottom = v.findAt(((0.0,-hs2,componentL/3),),((0.0,-hs2,componentL/3*2),))
	p.Set(name = 'span13bottom',vertices = span13bottom)
	#均布荷载
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	spanTop = e.findAt(((0.0,hs1,0.0001),),((0.0,hs1,componentL/3+0.0001),),
		((0.0,hs1,componentL/2+0.0001),),((0.0,hs1,componentL/3*2+0.0001),))
	spanShear = e.findAt(((0.0,0.0,0.0001),),((0.0,0.0,componentL/3+0.0001),),
		((0.0,0.0,componentL/2+0.0001),),((0.0,0.0,componentL/3*2+0.0001),))
	spanBottom = e.findAt(((0.0,-hs2,0.0001),),((0.0,-hs2,componentL/3+0.0001),),
		((0.0,-hs2,componentL/2+0.0001),),((0.0,-hs2,componentL/3*2+0.0001),))
	for i in range(len(spanTop)):
		spanTopNodes = spanTop[i].getNodes()
		spanShearNodes = spanShear[i].getNodes()
		spanBottomNodes = spanBottom[i].getNodes()
		p.Set(name = 'spanTop%d'%i,nodes = spanTopNodes)
		p.Set(name = 'spanShear%d'%i,nodes = spanShearNodes)
		p.Set(name = 'spanBottom%d'%i,nodes = spanBottomNodes)
	p.SetByBoolean(name='spanTop', sets=(p.sets['spanTop0'], p.sets['spanTop1'],
		p.sets['spanTop2'],p.sets['spanTop3'],))
	p.SetByBoolean(name='spanShear', sets=(p.sets['spanShear0'], p.sets['spanShear1'],
		p.sets['spanShear2'],p.sets['spanShear3'],))
	p.SetByBoolean(name='spanBottom', sets=(p.sets['spanBottom0'], p.sets['spanBottom1'],
		p.sets['spanBottom2'],p.sets['spanBottom3'],))
	for i in range(len(spanTop)):
		del mdb.models['Model-1'].parts['Part-1'].sets['spanTop%d'%i]
		del mdb.models['Model-1'].parts['Part-1'].sets['spanShear%d'%i]
		del mdb.models['Model-1'].parts['Part-1'].sets['spanBottom%d'%i]
	
	
	##创建实体Instance
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByDefault(CARTESIAN)
	p = mdb.models['Model-1'].parts['Part-1']
	a.Instance(name='Part-1-1', part=p, dependent=ON)
	
	
	##创建分析步Step
	#默认分析前5个模态
	mdb.models['Model-1'].BuckleStep(name='buckling', previous='Initial', 
        numEigen=5, vectors=10, maxIterations=200)
	
	
	##施加边界条件BCs
	#A侧（非原点侧）边界条件
	if BC == 's-s':
		a = mdb.models['Model-1'].rootAssembly
		region = a.instances['Part-1-1'].sets['A']
		mdb.models['Model-1'].DisplacementBC(name='A-Simple', createStepName='Initial', 
			region=region, u1=SET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=SET, 
			amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
			localCsys=None)
		a = mdb.models['Model-1'].rootAssembly
		region = a.instances['Part-1-1'].sets['B']
		mdb.models['Model-1'].DisplacementBC(name='B-Simple', createStepName='Initial', 
			region=region, u1=SET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=SET, 
			amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
			localCsys=None)
		a = mdb.models['Model-1'].rootAssembly
		region = a.instances['Part-1-1'].sets['AS']
		mdb.models['Model-1'].DisplacementBC(name='AS-z', createStepName='Initial', 
			region=region, u1=UNSET, u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
			amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
			localCsys=None)
	elif BC == 'c-c':
		a = mdb.models['Model-1'].rootAssembly
		region = a.instances['Part-1-1'].sets['A']
		mdb.models['Model-1'].EncastreBC(name='A-Clamped', createStepName='Initial', 
			region=region, localCsys=None)
		a = mdb.models['Model-1'].rootAssembly
		region = a.instances['Part-1-1'].sets['B']
		mdb.models['Model-1'].EncastreBC(name='B-Clamped', createStepName='Initial', 
			region=region, localCsys=None)
	elif BC == 'cantilever':
		a = mdb.models['Model-1'].rootAssembly
		region = a.instances['Part-1-1'].sets['A']
		mdb.models['Model-1'].EncastreBC(name='A-Clamped', createStepName='Initial', 
			region=region, localCsys=None)
	
	
	##施加默认单位荷载
	#跨中集中荷载
	a = mdb.models['Model-1'].rootAssembly
	region = a.instances['Part-1-1'].sets['span12top']
	mdb.models['Model-1'].ConcentratedForce(name='span12top', 
		createStepName='buckling', region=region, cf2=-1.0, 
		distributionType=UNIFORM, field='', localCsys=None)
	region = a.instances['Part-1-1'].sets['span12shear']
	mdb.models['Model-1'].ConcentratedForce(name='span12shear', 
		createStepName='buckling', region=region, cf2=-1.0, 
		distributionType=UNIFORM, field='', localCsys=None)
	region = a.instances['Part-1-1'].sets['span12bottom']
	mdb.models['Model-1'].ConcentratedForce(name='span12bottom', 
		createStepName='buckling', region=region, cf2=-1.0, 
		distributionType=UNIFORM, field='', localCsys=None)
	#1/3集中荷载
	region = a.instances['Part-1-1'].sets['span13top']
	mdb.models['Model-1'].ConcentratedForce(name='span13top', 
		createStepName='buckling', region=region, cf2=-1.0, 
		distributionType=UNIFORM, field='', localCsys=None)
	region = a.instances['Part-1-1'].sets['span13shear']
	mdb.models['Model-1'].ConcentratedForce(name='span13shear', 
		createStepName='buckling', region=region, cf2=-1.0, 
		distributionType=UNIFORM, field='', localCsys=None)
	region = a.instances['Part-1-1'].sets['span13bottom']
	mdb.models['Model-1'].ConcentratedForce(name='span13bottom', 
		createStepName='buckling', region=region, cf2=-1.0, 
		distributionType=UNIFORM, field='', localCsys=None)
	#均布荷载
	if ELM == 'S8R5':
		forceAtPerNode = 1.0*(componentL)/(2*lm+1)
	else:
		forceAtPerNode = 1.0*(componentL)/(lm+1)
	region = a.instances['Part-1-1'].sets['spanTop']
	mdb.models['Model-1'].ConcentratedForce(name='spanTop', 
		createStepName='buckling', region=region, cf2=-forceAtPerNode, 
		distributionType=UNIFORM, field='', localCsys=None)
	region = a.instances['Part-1-1'].sets['spanShear']
	mdb.models['Model-1'].ConcentratedForce(name='spanShear', 
		createStepName='buckling', region=region, cf2=-forceAtPerNode, 
		distributionType=UNIFORM, field='', localCsys=None)
	region = a.instances['Part-1-1'].sets['spanBottom']
	mdb.models['Model-1'].ConcentratedForce(name='spanBottom', 
		createStepName='buckling', region=region, cf2=-forceAtPerNode, 
		distributionType=UNIFORM, field='', localCsys=None)
	

	##创建构件显示窗口
	session.viewports['Viewport: 1'].setValues(displayedObject=a)
	session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
        predefinedFields=ON, connectors=ON)
	session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
	
	
	##创建Job
	nwb = int((s1+s2))
	ntf = int((b1))
	nbf = int((b2))
	mdb.Job(name='IB%dX%dX%d'%(nwb,ntf,nbf), model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
        numGPUs=0)
	
	
	##创建INP文件
	mdb.jobs['IB%dX%dX%d'%(nwb,ntf,nbf)].writeInput(consistencyChecking=OFF)


##输入窗口
##单位：力kN；长度mm
modeldata = (('E(MPa)','2.06E5'),('u','0.3'),
	('TopFlange','200'),('BottomFlange','150'),
	('ShearToTop','156.55'),('ShearToBottom','363.45'),
	('TopFlangeThickness','10'),('BottomFlangeThickness','10'),('WebThinckness','8'),
	('Length(m)','12'),
	('BoundaryConduction\n (simple:s-s,clamped:c-c,\n cantilever:cantilever)','s-s'),
	('ElementType\n S4R, S4R5, S8R5','S4R'),
	('TopFlangeMesh','8'),('BottomFlangeMesh','6'),('WebMesh','6'),('LongitudinalMesh','60'))
enterLabel = '''Please enter the data of the model.
The global units are kN and mm without specifying.
The unit of Eigenvalues is N !!!
The default loads includes
1/2 concentrated loads, symmetric 1/3 concentrated loads,
and uniform distributed loads,
subjected at top flange, shear points and bottom flange.'''
E,u,b1,b2,s1,s2,tf1,tf2,tw,componentL,BC,ELM,tm,bm,wm,lm = (
	getInputs(fields=modeldata,label=enterLabel))
dtl = [E,u,b1,b2,s1,s2,tf1,tf2,tw]
dtl = map(float,dtl)
ell = map(int,[tm,bm,wm,lm])
componentL = float(componentL)

			
##调用建模函数
ConstructeIBeams(dtl[0],dtl[1],dtl[2],dtl[3],dtl[4],dtl[5],
	dtl[6],dtl[7],dtl[8],componentL*1000,BC,ELM,ell[0],ell[1],ell[2],ell[3])
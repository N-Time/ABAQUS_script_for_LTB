# -*- coding: mbcs -*-
###ACIB v0.4.1 DCB##
###ACIB is AutoConstructe Model of I-Beams for Doublespans Continuous Beams
##Update Log
#v0.2 for simple span beams
#适用于单跨梁，包括两端简支、两端固支和悬臂梁
#改用findAt(...)方法寻找对象，而非通用性较弱的getSequenceFromMask(...)方法
#v0.3 for continuous
#连续梁专用版本，仅适用于两跨连续梁
###v0.4 for continuous
##适用于双跨连续梁（平面内外皆为连续梁）。
##增加自动施加均布荷载以及跨中集中荷载，且荷载作用位置包括上、剪、下。
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
def ConstructeIBeams(E,u,b1,b2,s1,s2,tf1,tf2,tw,componentL1,componentL2,ELM,tm,bm,wm,lm1,lm2):

	##参数说明
	#E = 弹性模量（N/m^2）
	#u = 泊松比
	#s1, s2 = 剪心至上翼缘上缘和下翼缘下缘的距离（mm）
	#b1, b2 = 上翼缘宽度和下翼缘宽度（mm）
	#tf1, tf2, tw = 上翼缘厚度，下翼缘厚度，腹板厚度（mm）
	#componentL1, componentL2 = 第一跨（靠近原点）构件长度，第二跨构件长度（m）
	#ELM = 单元类型
	#tm, bm, wm = 上翼缘网格数，下翼缘网格数，腹板网格数
	#lm1, lm2 = 连续梁第一跨纵向单元数，第二跨纵向单元数
	
	
	##计算中线节点的坐标
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
	componentL = componentL1 + componentL2
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
	
	
	##分割构件成两跨
	#建立参考面
	p = mdb.models['Model-1'].parts['Part-1']
	dp = p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=componentL1)
	dpID = dp.id
	#分割构件
	p = mdb.models['Model-1'].parts['Part-1']
	f = p.faces
	pickedFaces = f.findAt(((-bf1,hs1,0.0),),((bf1,hs1,0.0),),
		((0.0,0.0,0.0),),((-bf2,-hs2,0.0),),((bf2,-hs2,0.0),))
	d1 = p.datums
	p.PartitionFaceByDatumPlane(datumPlane=d1[dpID], faces=pickedFaces)
	
	
	##切割构件以施加荷载
	#第一跨二分点
	p = mdb.models['Model-1'].parts['Part-1']
	span1dp12 = p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=componentL1/2)
	span1dp12ID = span1dp12.id
	f = p.faces
	pickedFaces = f.findAt(((-bf1,hs1,componentL1/2),),((bf1,hs1,componentL1/2),),
		((0.0,0.0,componentL1/2),),((-bf2,-hs2,componentL1/2),),((bf2,-hs2,componentL1/2),))
	span1d12 = p.datums
	p.PartitionFaceByDatumPlane(datumPlane=span1d12[span1dp12ID], faces=pickedFaces)
	#第二跨二分点
	p = mdb.models['Model-1'].parts['Part-1']
	span2dp12 = p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=componentL1+componentL2/2)
	span2dp12ID = span2dp12.id
	f = p.faces
	pickedFaces = f.findAt(((-bf1,hs1,componentL1+componentL2/2),),
		((bf1,hs1,componentL1+componentL2/2),),((0.0,0.0,componentL1+componentL2/2),),
		((-bf2,-hs2,componentL1+componentL2/2),),((bf2,-hs2,componentL1+componentL2/2),))
	span2d12 = p.datums
	p.PartitionFaceByDatumPlane(datumPlane=span2d12[span2dp12ID], faces=pickedFaces)
	#剪心
	p = mdb.models['Model-1'].parts['Part-1']
	dps = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
	dpsID = dps.id
	f = p.faces
	pickedFaces = f.findAt(((0.0,0.0,0.0),),
		((0.0,0.0,componentL1/2+1.0),),((0.0,0.0,componentL1+1.0),),
		((0.0,0.0,componentL),))
	dps = p.datums
	p.PartitionFaceByDatumPlane(datumPlane=dps[dpsID], faces=pickedFaces)
	
	
	##创建边界条件作用域集合（几何型集合）
	#A侧（原点侧）
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	edges = e.findAt(((-bf1+1.0,hs1,0.0),),((bf1-1.0,hs1,0.0),),
		((-bf2+1.0,-hs2,0.0),),((bf2-1.0,-hs2,0.0),),
		((0.0,1.0,0.0),),((0.0,-1.0,0.0),))
	p.Set(edges=edges, name='A')
	#A侧（约束轴向位移）
	p = mdb.models['Model-1'].parts['Part-1']
	v = p.vertices
	verts = v.findAt(((0.0,0.0,0.0),))
	p.Set(vertices=verts, name='AS')
	#中间支座
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	edges = e.findAt(((-bf1+1.0,hs1,componentL1),),((bf1-1.0,hs1,componentL1),),
		((-bf2+1.0,-hs2,componentL1),),((bf2-1.0,-hs2,componentL1),),
		((0.0,1.0,componentL1),),((0.0,-1.0,componentL1),))
	p.Set(edges=edges, name='B')
	#C侧（非原点侧）
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	edges = e.findAt(((-bf1+1.0,hs1,componentL),),((bf1-1.0,hs1,componentL),),
		((-bf2+1.0,-hs2,componentL),),((bf2-1.0,-hs2,componentL),),
		((0.0,1.0,componentL),),((0.0,-1.0,componentL),))
	p.Set(edges=edges, name='C')
	
	
	##布置种子Seed以及网格划分Mesh
	htm = tm/2
	hbm = bm/2
	#上翼缘
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	pickedEdges = e.findAt(((-bf1+1.0,hs1,0.0),),((bf1-1.0,hs1,0.0),),
		((-bf1+1.0,hs1,componentL1/2),),((bf1-1.0,hs1,componentL1/2),),
		((-bf1+1.0,hs1,componentL1),),((bf1-1.0,hs1,componentL1),),
		((-bf1+1.0,hs1,componentL1+componentL2/2),),((bf1-1.0,hs1,componentL1+componentL2/2),),
		((-bf1+1.0,hs1,componentL),),((bf1-1.0,hs1,componentL),))
	p.seedEdgeByNumber(edges=pickedEdges, number=htm, constraint=FINER)
	#下翼缘
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	pickedEdges = e.findAt(((-bf2+1.0,-hs2,0.0),),((bf2-1.0,-hs2,0.0),),
		((-bf2+1.0,-hs2,componentL1/2),),((bf2-1.0,-hs2,componentL1/2),),
		((-bf2+1.0,-hs2,componentL1),),((bf2-1.0,-hs2,componentL1),),
		((-bf2+1.0,-hs2,componentL1+componentL2/2),),((bf2-1.0,-hs2,componentL1+componentL2/2),),
		((-bf2+1.0,-hs2,componentL),),((bf2-1.0,-hs2,componentL),))
	p.seedEdgeByNumber(edges=pickedEdges, number=hbm, constraint=FINER)
	#腹板
	wmt = int((s1-tf1/2)/(s1+s2-tf1/2-tf2/2)*wm+0.5)
	wmb = int((s2-tf2/2)/(s1+s2-tf1/2-tf2/2)*wm+0.5)
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	pickedEdges = e.findAt(((0.0,1.0,0.0),),((0.0,1.0,componentL1/2),),
		((0.0,1.0,componentL1),),((0.0,1.0,componentL1+componentL2/2),),
		((0.0,1.0,componentL),))
	p.seedEdgeByNumber(edges=pickedEdges, number=wmt, constraint=FINER)
	pickedEdges = e.findAt(((0.0,-1.0,0.0),),((0.0,-1.0,componentL1/2),),
		((0.0,-1.0,componentL1),),((0.0,-1.0,componentL1+componentL2/2),),
		((0.0,-1.0,componentL),))
	p.seedEdgeByNumber(edges=pickedEdges, number=wmb, constraint=FINER)
	#第一跨纵向
	lm112 = int(lm1/2+0.5)
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	pickedEdges = e.findAt(((-bf1,hs1,1.0),),((bf1,hs1,1.0),),
		((0.0,hs1,1.0),),((0.0,-hs2,1.0),),
		((-bf2,-hs2,1.0),),((bf2,-hs2,1.0),),
		((0.0,0.0,1.0),),((0.0,0.0,componentL1/2+1.0),),
		((-bf1,hs1,componentL1/2+1.0),),((bf1,hs1,componentL1/2+1.0),),
		((0.0,hs1,componentL1/2+1.0),),((0.0,-hs2,componentL1/2+1.0),),
		((-bf2,-hs2,componentL1/2+1.0),),((bf2,-hs2,componentL1/2+1.0),))
	p.seedEdgeByNumber(edges=pickedEdges, number=lm112, constraint=FINER)
	#第二跨纵向
	lm212 = int(lm2/2+0.5)
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	pickedEdges = e.findAt(((-bf1,hs1,componentL1+1.0),),((bf1,hs1,componentL1+1.0),),
		((0.0,hs1,componentL1+1.0),),((0.0,-hs2,componentL1+1.0),),
		((-bf2,-hs2,componentL1+1.0),),((bf2,-hs2,componentL1+1.0),),
		((0.0,0.0,componentL1+1.0),),((0.0,0.0,componentL-1.0),),
		((-bf1,hs1,componentL-1.0),),((bf1,hs1,componentL-1.0),),
		((0.0,hs1,componentL-1.0),),((0.0,-hs2,componentL-1.0),),
		((-bf2,-hs2,componentL-1.0),),((bf2,-hs2,componentL-1.0),))
	p.seedEdgeByNumber(edges=pickedEdges, number=lm212, constraint=FINER)
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
	#第一跨2分点荷载集合
	p = mdb.models['Model-1'].parts['Part-1']
	v = p.vertices
	span112top = v.findAt(((0.0,hs1,componentL1/2),))
	p.Set(name = 'span112top',vertices = span112top)
	span112shear = v.findAt(((0.0,0.0,componentL1/2),))
	p.Set(name = 'span112shear',vertices = span112shear)
	span112bottom = v.findAt(((0.0,-hs2,componentL1/2),))
	p.Set(name = 'span112bottom',vertices = span112bottom)
	#第二跨2分点荷载集合
	p = mdb.models['Model-1'].parts['Part-1']
	v = p.vertices
	span212top = v.findAt(((0.0,hs1,componentL1+componentL2/2),))
	p.Set(name = 'span212top',vertices = span212top)
	span212shear = v.findAt(((0.0,0.0,componentL1+componentL2/2),))
	p.Set(name = 'span212shear',vertices = span212shear)
	span212bottom = v.findAt(((0.0,-hs2,componentL1+componentL2/2),))
	p.Set(name = 'span212bottom',vertices = span212bottom)
	#第一跨均布荷载集合
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	span1Top = e.findAt(((0.0,hs1,1.0),),((0.0,hs1,componentL1/2+1.0),))
	span1Shear = e.findAt(((0.0,0.0,1.0),),((0.0,0.0,componentL1/2+1.0),))
	span1Bottom = e.findAt(((0.0,-hs2,1.0),),((0.0,-hs2,componentL1/2+1.0),))
	for i in range(len(span1Top)):
		span1TopNodes = span1Top[i].getNodes()
		span1ShearNodes = span1Shear[i].getNodes()
		span1BottomNodes = span1Bottom[i].getNodes()
		p.Set(name = 'span1Top%d'%i,nodes = span1TopNodes)
		p.Set(name = 'span1Shear%d'%i,nodes = span1ShearNodes)
		p.Set(name = 'span1Bottom%d'%i,nodes = span1BottomNodes)
	p.SetByBoolean(name='span1Top', sets=(p.sets['span1Top0'], p.sets['span1Top1'],))
	p.SetByBoolean(name='span1Shear', sets=(p.sets['span1Shear0'], p.sets['span1Shear1'],))
	p.SetByBoolean(name='span1Bottom', sets=(p.sets['span1Bottom0'], p.sets['span1Bottom1'],))
	for i in range(len(span1Top)):
		del mdb.models['Model-1'].parts['Part-1'].sets['span1Top%d'%i]
		del mdb.models['Model-1'].parts['Part-1'].sets['span1Shear%d'%i]
		del mdb.models['Model-1'].parts['Part-1'].sets['span1Bottom%d'%i]
	#第二跨均布荷载几何
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	span2Top = e.findAt(((0.0,hs1,componentL1+1.0),),((0.0,hs1,componentL-1.0),))
	span2Shear = e.findAt(((0.0,0.0,componentL1+1.0),),((0.0,0.0,componentL-1.0),))
	span2Bottom = e.findAt(((0.0,-hs2,componentL1+1.0),),((0.0,-hs2,componentL-1.0),))
	for i in range(len(span2Top)):
		span2TopNodes = span2Top[i].getNodes()
		span2ShearNodes = span2Shear[i].getNodes()
		span2BottomNodes = span2Bottom[i].getNodes()
		p.Set(name = 'span2Top%d'%i,nodes = span2TopNodes)
		p.Set(name = 'span2Shear%d'%i,nodes = span2ShearNodes)
		p.Set(name = 'span2Bottom%d'%i,nodes = span2BottomNodes)
	p.SetByBoolean(name='span2Top', sets=(p.sets['span2Top0'], p.sets['span2Top1'],))
	p.SetByBoolean(name='span2Shear', sets=(p.sets['span2Shear0'], p.sets['span2Shear1'],))
	p.SetByBoolean(name='span2Bottom', sets=(p.sets['span2Bottom0'], p.sets['span2Bottom1'],))
	for i in range(len(span2Top)):
		del mdb.models['Model-1'].parts['Part-1'].sets['span2Top%d'%i]
		del mdb.models['Model-1'].parts['Part-1'].sets['span2Shear%d'%i]
		del mdb.models['Model-1'].parts['Part-1'].sets['span2Bottom%d'%i]
	
	
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
	#A侧（原点侧）边界条件
	a = mdb.models['Model-1'].rootAssembly
	region = a.instances['Part-1-1'].sets['A']
	mdb.models['Model-1'].DisplacementBC(name='A', createStepName='Initial', 
		region=region, u1=SET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=SET, 
		amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
		localCsys=None)
	a = mdb.models['Model-1'].rootAssembly
	region = a.instances['Part-1-1'].sets['AS']
	mdb.models['Model-1'].DisplacementBC(name='AS-z', createStepName='Initial', 
		region=region, u1=UNSET, u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
		amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
		localCsys=None)
	#B支座（中间支座）边界条件
	a = mdb.models['Model-1'].rootAssembly
	region = a.instances['Part-1-1'].sets['B']
	mdb.models['Model-1'].DisplacementBC(name='B', createStepName='Initial', 
		region=region, u1=SET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=SET, 
		amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
		localCsys=None)
	#C支座（非原点侧）边界条件
	a = mdb.models['Model-1'].rootAssembly
	region = a.instances['Part-1-1'].sets['C']
	mdb.models['Model-1'].DisplacementBC(name='C', createStepName='Initial', 
		region=region, u1=SET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=SET, 
		amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
		localCsys=None)
	
	
	##施加默认单位荷载
	#第一跨跨中集中荷载
	a = mdb.models['Model-1'].rootAssembly
	region = a.instances['Part-1-1'].sets['span112top']
	mdb.models['Model-1'].ConcentratedForce(name='span112top', 
		createStepName='buckling', region=region, cf2=-1.0, 
		distributionType=UNIFORM, field='', localCsys=None)
	region = a.instances['Part-1-1'].sets['span112shear']
	mdb.models['Model-1'].ConcentratedForce(name='span112shear', 
		createStepName='buckling', region=region, cf2=-1.0, 
		distributionType=UNIFORM, field='', localCsys=None)
	region = a.instances['Part-1-1'].sets['span112bottom']
	mdb.models['Model-1'].ConcentratedForce(name='span112bottom', 
		createStepName='buckling', region=region, cf2=-1.0, 
		distributionType=UNIFORM, field='', localCsys=None)
	#第二跨跨中集中荷载
	a = mdb.models['Model-1'].rootAssembly
	region = a.instances['Part-1-1'].sets['span212top']
	mdb.models['Model-1'].ConcentratedForce(name='span212top', 
		createStepName='buckling', region=region, cf2=-1.0, 
		distributionType=UNIFORM, field='', localCsys=None)
	region = a.instances['Part-1-1'].sets['span212shear']
	mdb.models['Model-1'].ConcentratedForce(name='span212shear', 
		createStepName='buckling', region=region, cf2=-1.0, 
		distributionType=UNIFORM, field='', localCsys=None)
	region = a.instances['Part-1-1'].sets['span212bottom']
	mdb.models['Model-1'].ConcentratedForce(name='span212bottom', 
		createStepName='buckling', region=region, cf2=-1.0, 
		distributionType=UNIFORM, field='', localCsys=None)
	#第一跨均布荷载
	if ELM == 'S8R5':
		forceAtPerNode = 1.0*(componentL1/1000.0)/(2*lm1+1.0)
	else:
		forceAtPerNode = 1.0*(componentL1/1000.0)/(lm1+1.0)
	region = a.instances['Part-1-1'].sets['span1Top']
	mdb.models['Model-1'].ConcentratedForce(name='span1Top', 
		createStepName='buckling', region=region, cf2=-forceAtPerNode, 
		distributionType=UNIFORM, field='', localCsys=None)
	region = a.instances['Part-1-1'].sets['span1Shear']
	mdb.models['Model-1'].ConcentratedForce(name='span1Shear', 
		createStepName='buckling', region=region, cf2=-forceAtPerNode, 
		distributionType=UNIFORM, field='', localCsys=None)
	region = a.instances['Part-1-1'].sets['span1Bottom']
	mdb.models['Model-1'].ConcentratedForce(name='span1Bottom', 
		createStepName='buckling', region=region, cf2=-forceAtPerNode, 
		distributionType=UNIFORM, field='', localCsys=None)
	#第二跨均布荷载
	if ELM == 'S8R5':
		forceAtPerNode = 1.0*(componentL1/1000.0)/(2*lm2+1.0)
	else:
		forceAtPerNode = 1.0*(componentL1/1000.0)/(lm2+1.0)
	region = a.instances['Part-1-1'].sets['span2Top']
	mdb.models['Model-1'].ConcentratedForce(name='span2Top', 
		createStepName='buckling', region=region, cf2=-forceAtPerNode, 
		distributionType=UNIFORM, field='', localCsys=None)
	region = a.instances['Part-1-1'].sets['span2Shear']
	mdb.models['Model-1'].ConcentratedForce(name='span2Shear', 
		createStepName='buckling', region=region, cf2=-forceAtPerNode, 
		distributionType=UNIFORM, field='', localCsys=None)
	region = a.instances['Part-1-1'].sets['span2Bottom']
	mdb.models['Model-1'].ConcentratedForce(name='span2Bottom', 
		createStepName='buckling', region=region, cf2=-forceAtPerNode, 
		distributionType=UNIFORM, field='', localCsys=None)


	##创建构件显示窗口
	session.viewports['Viewport: 1'].setValues(displayedObject=a)
	session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
        predefinedFields=ON, connectors=ON)
	session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
	
	
	##创建Job
	nwb = int(s1+s2)
	ntf = int(b1)
	nbf = int(b2)
	mdb.Job(name='IB%dX%dX%d-continuous'%(nwb,ntf,nbf), model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
        numGPUs=0)
	
	
	##创建INP文件
	mdb.jobs['IB%dX%dX%d-continuous'%(nwb,ntf,nbf)].writeInput(consistencyChecking=OFF)


##输入窗口
##单位：力kN；长度mm
modeldata = (('E(kN/mm^2)','2.06E2'),('u','0.3'),
	('TopFlange','200'),('BottomFlange','150'),
	('ShearToTop','156.55'),('ShearToBottom','363.45'),
	('TopFlangeThickness','10'),('BottomFlangeThickness','10'),('WebThinckness','8'),
	('Length1(m)','12'),('Length2(m)','9'),
	('ElementType\n S4R, S4R5, S8R5','S4R5'),
	('TopFlangeMesh','8'),('BottomFlangeMesh','6'),('WebMesh','6'),
	('Longitudinal1Mesh','40'),('Longitudinal2Mesh','30'))
enterLabel = '''Please enter the data of the model.
The global units are kN and mm if without specifying.
The unit of Eigenvalues is kN !!!
The default loads includes
1/2 concentrated loads on both spans,
and uniform distributed loads on both spans,
subjected at top flange, shear points and bottom flange.'''
E,u,b1,b2,s1,s2,tf1,tf2,tw,componentL1,componentL2,ELM,tm,bm,wm,lm1,lm2 = (
	getInputs(fields=modeldata,label=enterLabel))
dtl = [E,u,b1,b2,s1,s2,tf1,tf2,tw]
dtl = map(float,dtl)
ell = map(int,[tm,bm,wm,lm1,lm2])
componentL1 = float(componentL1)*1000
componentL2 = float(componentL2)*1000

			
##调用建模函数
ConstructeIBeams(dtl[0],dtl[1],dtl[2],dtl[3],dtl[4],dtl[5],
	dtl[6],dtl[7],dtl[8],componentL1,componentL2,ELM,ell[0],ell[1],ell[2],ell[3],ell[4])
# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2016 replay file
# Internal Version: 2015_09_25-04.31.09 126547
# Run by Administrator on Fri Aug 30 20:22:09 2019
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
#: Abaqus Error: Execution of "onCaeGraphicsStartup()" in the site directory failed.
#:   File "D:\ProgramFiles\DassaultSystemes\SimulationService\V6R2016x\win_b64\SMA\site\graphicsConfig.env", line 387, in onCaeGraphicsStartup
#:     elif glRenderer[32] == "Mobile Intel(R) 4 Series Express":
#: IndexError: string index out of range
#:
#: This error may have occurred due to a change to the Abaqus Scripting
#: Interface. Please see the Abaqus Scripting Manual for the details of
#: these changes. Also see the "Example environment files" section of
#: the Abaqus Site Guide for up-to-date examples of common tasks in the
#: environment file.
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup

executeOnCaeStartup()
step = mdb.openStep('E:/temp/yanmogangqiu.stp', scaleFromFile=OFF)
mdb.models['Model-1'].PartFromGeometryFile(name='Part-qiu', geometryFile=step,
    combine=False, dimensionality=THREE_D, type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-qiu']
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=20.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=(-2, 2), point2=(2, -2))
p = mdb.models['Model-1'].Part(name='Part-ban', dimensionality=THREE_D,
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-ban']
p.BaseSolidExtrude(sketch=s, depth=2.5)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-ban']
del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].Material(name='Material-ban')
mdb.models['Model-1'].materials['Material-ban'].Density(table=((7.85e-09, ), ))
mdb.models['Model-1'].materials['Material-ban'].Elastic(table=((217000.0, 0.3),
    ))
mdb.models['Model-1'].materials['Material-ban'].Plastic(hardening=JOHNSON_COOK,
    table=((2482.4, 1498.5, 0.19, 0.66, 1698.0, 298.0), ))
mdb.models['Model-1'].materials['Material-ban'].plastic.RateDependent(
    type=JOHNSON_COOK, table=((0.027, 1.0), ))
mdb.models['Model-1'].Material(name='Material-qiu',
    objectToCopy=mdb.models['Model-1'].materials['Material-ban'])
mdb.models['Model-1'].HomogeneousSolidSection(name='Section-ban',
    material='Material-ban', thickness=None)
mdb.models['Model-1'].HomogeneousSolidSection(name='Section-qiu',
    material='Material-qiu', thickness=None)
# session.viewports['Viewport: 1'].view.setValues(nearPlane=4.39276,
    # farPlane=8.16178, width=5.96412, height=2.69063, viewOffsetX=0.305879,
    # viewOffsetY=-0.0922151)
p = mdb.models['Model-1'].parts['Part-ban']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
region = p.Set(cells=cells, name='Set-1')
p = mdb.models['Model-1'].parts['Part-ban']
p.SectionAssignment(region=region, sectionName='Section-ban', offset=0.0,
    offsetType=MIDDLE_SURFACE, offsetField='',
    thicknessAssignment=FROM_SECTION)
p = mdb.models['Model-1'].parts['Part-qiu']
p = mdb.models['Model-1'].parts['Part-qiu']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
region = p.Set(cells=cells, name='Set-1')
p = mdb.models['Model-1'].parts['Part-qiu']
p.SectionAssignment(region=region, sectionName='Section-qiu', offset=0.0,
    offsetType=MIDDLE_SURFACE, offsetField='',
    thicknessAssignment=FROM_SECTION)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-ban']
a.Instance(name='Part-ban-1', part=p, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('Part-ban-1', ), vector=(0.0, 0.0, 20.0))         #part-ban-1装配位置改变
import random
import numpy as np
import xlrd

# class CallingCounter(object):    #函数调用次数统计装饰器
    # def __init__ (self, func):
        # self.func = func
        # self.count = 0

    # def __call__ (self, *args, **kwargs):
        # self.count += 1
        # return self.func(*args, **kwargs)


class excel_read:                  #读取颗粒数据
    def __init__(self, excel_path=r'C:\\Users\\GAO\\Desktop\\shuju.xlsx',encoding='utf-8',index=0):

      self.data=xlrd.open_workbook(excel_path)  ##获取文本对象
      self.table=self.data.sheets()[index]     ###根据index获取某个sheet
      self.rows=self.table.nrows   ##3获取当前sheet页面的总行数,把每一行数据作为list放到 list

    def get_data(self):
        L1=[]                             #存贮速度
        L2=[]                             #存储位置
        for i in range(self.rows):
            col=self.table.row_values(i)  ##获取每一列数据
            #print(col)
            L1.append(col[0:3])
            L2.append(col[3:6])

        return L1,L2

L3,L4=excel_read().get_data()
N=len(L4)                           #单层颗粒数
##  @CallingCounter
def rand_ball():           #随机小球生成函数
    i=0
    while i<N:
        vx=L3[i][0]
        vy=L3[i][2]
        vz=L3[i][1]
        tmp1=np.linalg.norm([vx,vy,vz])     #速度幅值
        SCAL=random.uniform(0,1)
        x1 = L4[i][0]-vx/tmp1*SCAL      #原坐标下x
        y1 = L4[i][2]-vy/tmp1*SCAL      #原坐标下z
        z1 = L4[i][1]-vz/tmp1*SCAL-0.2      #原坐标下y
        p = mdb.models['Model-1'].parts['Part-qiu']
        a = mdb.models['Model-1'].rootAssembly
        a.Instance(name='part-qiu'+str(i+1), part=p, dependent=ON)
        a.translate(instanceList=('part-qiu'+str(i+1),), vector=(x1,y1,z1 ))
        a.ReferencePoint(point=(x1,y1,z1))
        a = mdb.models['Model-1'].rootAssembly
        a.rotate(instanceList=('part-qiu' + str(i+1),), axisPoint=(x1, y1, z1),
        axisDirection=(random.uniform(-1, 1), random.uniform(-1, 1),random.uniform(-1, 1)), angle=random.uniform(0, 360))
        a = mdb.models['Model-1'].rootAssembly
        s1 = a.instances['part-qiu' + str(i+1)].faces
        side1Faces1 = s1.getSequenceFromMask(mask=('[#ffffffff:8 #fff ]', ) )
        a.Surface(side1Faces=side1Faces1, name='Surf-'+str(i+1))
        i=i+1


L_1=rand_ball()        #生成小球

mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial',
    timePeriod=2e-4)

# for i in range(2,n+1):                              #创建多时间步
    # t_1=5e-05
    # mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-'+str(i), previous='Step-'+str(i-1), 
        # timePeriod=t_1, improvedDtMethod=ON)
p = mdb.models['Model-1'].parts['Part-ban']
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.5)
p = mdb.models['Model-1'].parts['Part-ban']
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=1.0)
p = mdb.models['Model-1'].parts['Part-ban']
del p.features['Datum plane-2']
p = mdb.models['Model-1'].parts['Part-ban']
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=1)
p = mdb.models['Model-1'].parts['Part-ban']
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=-1)
p = mdb.models['Model-1'].parts['Part-ban']
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=1)
p = mdb.models['Model-1'].parts['Part-ban']
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=-1)
p = mdb.models['Model-1'].parts['Part-ban']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
d1 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d1[3], cells=pickedCells)
p = mdb.models['Model-1'].parts['Part-ban']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#3 ]', ), )
d2 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d2[8], cells=pickedCells)
p = mdb.models['Model-1'].parts['Part-ban']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#a ]', ), )
d1 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d1[7], cells=pickedCells)
p = mdb.models['Model-1'].parts['Part-ban']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#3f ]', ), )
d2 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d2[6], cells=pickedCells)
p = mdb.models['Model-1'].parts['Part-ban']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#be0 ]', ), )
d1 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d1[5], cells=pickedCells)
mdb.models['Model-1'].ContactProperty('IntProp-1')
mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
    formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF,
    pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((
    0.3, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION,
    fraction=0.005, elasticSlipStiffness=None)
a = mdb.models['Model-1'].rootAssembly
region1=a.surfaces['Surf-1']
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-ban-1'].faces
faces1 = f1.getSequenceFromMask(mask=('[#21200 #2200015 #400 ]', ), )
region2=a.Set(faces=faces1, name='s_Set-3')
mdb.models['Model-1'].SurfaceToSurfaceContactExp(name ='Int1',
    createStepName='Step-1', master = region1, slave = region2,
    mechanicalConstraint=KINEMATIC, sliding=FINITE,
    interactionProperty='IntProp-1', initialClearance=OMIT, datumAxis=None,
    clearanceRegion=None)

v=6
v11=1
f_v=1000                                                 #速度放大系数:m/s转换成mm/s

for i in L3:                                             #预定义速度场加载
    vx11 = i[0]
    vy11 = i[2]
    vz11 = i[1]
    a = mdb.models['Model-1'].rootAssembly
    r1 = a.referencePoints
    refPoints1=(r1[v], )
    region = regionToolset.Region(referencePoints=refPoints1)
    mdb.models['Model-1'].Velocity(name='Predefined Field-'+'_'+str(v11), region=region, 
        field='', distributionType=MAGNITUDE, velocity1=vx11*f_v, velocity2=vy11*f_v, velocity3=vz11*f_v, omega=0.0)
    # if v11<=N*(n-1):
        # mdb.models['Model-1'].boundaryConditions['BC-'+str(i+1)+'_'+str(v11)].deactivate('Step-'+str(i+2))
    v=v+4
    v11=v11+1


a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-ban-1'].faces
faces1 = f1.getSequenceFromMask(mask=('[#0 #80000 ]', ), )
a.Set(faces=faces1, name='Set-2')
j=1
while j<=N:                           #定义接触
    a = mdb.models['Model-1'].rootAssembly
    region1 = a.surfaces['Surf-'+str(j)]
    a = mdb.models['Model-1'].rootAssembly
    region2 = a.sets['s_Set-3']
    mdb.models['Model-1'].SurfaceToSurfaceContactExp(name='Int'+str(j), createStepName='Step-1',
    master=region1, slave=region2,mechanicalConstraint=KINEMATIC, sliding=FINITE,
    interactionProperty='IntProp-1', initialClearance=OMIT,datumAxis=None,clearanceRegion=None)
    j=j+1

k=6
for i in range(N):                          #定义刚体约束
    a = mdb.models['Model-1'].rootAssembly
    region2 = a.instances['part-qiu'+str(i+1)].sets['Set-1']
    a = mdb.models['Model-1'].rootAssembly
    r1 = a.referencePoints
    refPoints1 = (r1[k],)
    region1 = regionToolset.Region(referencePoints=refPoints1)
    mdb.models['Model-1'].RigidBody(name='Constraint-'+str(i+1),
                 refPointRegion=region1, bodyRegion=region2)
    k=k+4


u=1
while u<=N:
    a = mdb.models['Model-1'].rootAssembly
    a.regenerate()
    mdb.models['Model-1'].HistoryOutputRequest(name='H-Output'+str(u+1),createStepName='Step-1',
    variables=ALL, interactions=('Int'+str(u),),sectionPoints=DEFAULT, rebar=EXCLUDE)
    u=u+1
p = mdb.models['Model-1'].parts['Part-qiu']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['Part-qiu']
c = p.cells
pickedRegions = c.getSequenceFromMask(mask=('[#1 ]', ), )
p.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
elemType1 = mesh.ElemType(elemCode=UNKNOWN_HEX, elemLibrary=EXPLICIT)
elemType2 = mesh.ElemType(elemCode=UNKNOWN_WEDGE, elemLibrary=EXPLICIT)
elemType3 = mesh.ElemType(elemCode=C3D10M, elemLibrary=EXPLICIT)
p = mdb.models['Model-1'].parts['Part-qiu']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2,
    elemType3))
p = mdb.models['Model-1'].parts['Part-qiu']
p.seedPart(size=0.04, deviationFactor=0.1, minSizeFactor=0.1)
elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT)
elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT,
    secondOrderAccuracy=OFF, distortionControl=DEFAULT)
p = mdb.models['Model-1'].parts['Part-qiu']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2,
    elemType3))
p = mdb.models['Model-1'].parts['Part-qiu']
p.generateMesh()
p = mdb.models['Model-1'].parts['Part-ban']
p.seedPart(size=0.4, deviationFactor=0.1, minSizeFactor=0.1)                      #全局种子
e = p.edges                                                                       #喷丸撞击区域布种
pickedEdges = e.getSequenceFromMask(mask=('[#4000 #1000000 #102 ]', ), )
p.seedEdgeBySize(edges=pickedEdges, size=0.02, deviationFactor=0.1, 
    constraint=FINER)
e = p.edges                                                                        #侧边布种
pickedEdges1 = e.getSequenceFromMask(mask=('[#208000 #4001000 #42100000 #40 ]', 
    ), )
pickedEdges2 = e.getSequenceFromMask(mask=('[#40000800 #100100 ]', ), )
p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, 
    end2Edges=pickedEdges2, minSize=0.01, maxSize=0.08, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-ban']
p.generateMesh()
elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT,
    kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF,
    hourglassControl=DEFAULT, distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)
p = mdb.models['Model-1'].parts['Part-ban']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#3ffff ]', ), )
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2,
    elemType3))
p = mdb.models['Model-1'].parts['Part-ban']
p.generateMesh()
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-ban-1'].faces
faces1 = f1.getSequenceFromMask(mask=('[#108400 #8010a08 #200 ]', ), )
region = a.Set(faces=faces1, name='Set-25')
mdb.models['Model-1'].EncastreBC(name='BC-0', createStepName='Initial',
    region=region, localCsys=None)
mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=95, 
    memoryUnits=PERCENTAGE, explicitPrecision=SINGLE, 
    nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
    contactPrint=OFF, historyPrint=ON, userSubroutine='', scratch='', 
    resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=6, 
    activateLoadBalancing=False, multiprocessingMode=THREADS, numCpus=6)
leaf = dgo.LeafFromPartInstance(partInstanceName=("PART-BAN-1", ))          #生成显示组
dg = session.DisplayGroup(leaf=leaf, name='DisplayGroup-2')
dg1= session.displayGroups['DisplayGroup-2']
session.viewports['Viewport: 1'].odbDisplay.setValues(visibleDisplayGroups=(
    dg1, ))
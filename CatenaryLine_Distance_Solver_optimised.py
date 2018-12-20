import matplotlib.pylab as plt
from matplotlib.font_manager import  FontProperties
import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
#用于解超越方程
from sympy import Symbol,solveset
import configparser
lineType=['Parabola','Catenary']
#dec计算输出精度，小数点后位数7
dec=4
e=math.e
#电线比载gama  N/mm2*m
gama=0.032719
#导线水平应力delta  N/mm2
delta=31.24
delta_sequence=[5,10,20,30,40,50,60,70]
#L=67
def ch(x):
    return((e**x+e**(-x))/2)
def sh(x):
    return((e**x-e**(-x))/2)

class Point():
    def __init__(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z
    def __str__(self):
        return '({0},{1},{2})'.format(round(self.x,dec),round(self.y,dec),round(self.z,dec))
    def coord(self):
        return [self.x,self.y,self.z]

class CatenaryLine():
    def __init__(self,start,end,gama,delta,formula=lineType[0]): 
        self.formula=formula
#       delta应力
        self.delta=delta
        #gama 比载
        self.gama=gama
        self.h=end.z-start.z
        self.theta_y=math.atan((start.y-end.y)/(end.x-start.x))
        #L悬链线在2D情况下的实际x轴长度
        self.L=(end.x-start.x)/math.cos(self.theta_y)
        #Lv悬链线在xoz平面投影的x轴长度
        self.Lv=(end.x-start.x)
        self.t=self.delta/self.gama
        #beta用于斜抛线方程
        self.beta=math.atan(self.h/self.L)
##        self.Loa2=self.L/2-delta/gama*math.sin(self.beta)
        #Loa真实长度
        if self.formula=="Catenary":
            #利用公式解Loa     
            x=Symbol('x')
            ans=solveset(sh(x)-self.h/2/self.t/sh(self.L/2/self.t),domain=S.Reals) 
            self.Loa=self.L/2-self.t*ans.args[0]         
            # self.Loa=(self.L/2-(self.t)*sh(self.h/2/self.t/sh(self.L/2/self.t)))
        if self.formula=="Parabola":
            self.Loa=self.L/2-self.delta/self.gama*math.sin(self.beta)
        print('Loa: ',self.Loa)
##        print('Loa2:',self.Loa2)
        #Arch_y架构侧点的y坐标，用于calculat_yzs函数
        self.Arch_y=start.y
    def calculate_yzs(self,x):
        ys=list()
        zs=list()
        for i in x:
            ys.append(self.Arch_y-i*math.tan(self.theta_y))
            if self.formula=='Catenary':
            #精确解
                zs.append(self.t*(ch((self.Loa-i/math.cos(self.theta_y))/self.t)-ch(self.Loa/self.t)))
            #斜抛线方程
            if self.formula=="Parabola":
                x1=i/math.cos(self.theta_y)
                zs.append(x1*math.tan(self.beta)-x1*(self.L/math.cos(self.theta_y)-x1)/2/self.t/math.cos(self.beta))
        return (ys,zs)
    
    def calculate_point(self,x):
        #通过给定直角坐标系下的x，计算出直角坐标系下，悬链线上点的位置
        y=self.Arch_y-x*math.tan(self.theta_y)
        if self.formula=='Catenary':
        #精确解
            z=self.t*(ch((self.Loa-x/math.cos(self.theta_y))/self.t)-ch(self.Loa/self.t))
        #斜抛线
        if self.formula=="Parabola":
            x1=x/math.cos(self.theta_y)
            z=x1*math.tan(self.beta)-x1*(self.L/math.cos(self.theta_y)-x1)/2/math.cos(self.beta)/self.t

        #print (x,y,z)
        return (y,z)
    
    def x_lenth(self):
        #返回x轴向长度
        return self.Lv

    def calculate_fm(self):
        if self.formula=="Catenary":
            #悬链线方程过于复杂，最大弧垂采用斜抛线方程计算式
            fm=self.gama*self.L**2/8/self.delta/math.cos(self.beta)
        if self.formula=="Parabola":
            fm=self.gama*self.L**2/8/self.delta/math.cos(self.beta)
        return fm


calculate_sequence=[1,0.2,0.04,0.01]
#计算精度序列
sequence_index=0
def Find_closest_distance(line1,line2,x_start,x_end):
    global sequence_index
    def distance_between_points(point1,point2):
        return math.sqrt((point1.x-point2.x)**2+(point1.y-point2.y)**2+(point1.z-point2.z)**2)
    #print('A',x_start,line1.calculate_point(x_start),x_end,line1.calculate_point(x_end))
    step=calculate_sequence[sequence_index]
    x=np.arange(x_start,x_end,step)
    #初始化最近点和最近距离
    closest_distance=100
    position1=Point(0,0,0)
    position2=Point(0,0,0)
    for i in x:
        for j in x:
            point1=Point(i,*line1.calculate_point(i))
            point2=Point(j,*line2.calculate_point(j))
            dist=distance_between_points(point1,point2)
            if closest_distance>dist:
                closest_distance=dist
                position1=point1
                position2=point2
    if sequence_index<len(calculate_sequence)-1 :
        sequence_index=sequence_index+1
        start=min(position1.x,position2.x)-calculate_sequence[sequence_index]
        if start<0:
            start=0
        end=max(position1.x,position2.x)+calculate_sequence[sequence_index]
        if end>line1.x_lenth():
            end=line1.x_lenth()
        return(Find_closest_distance(line1,line2,start,end))
    else:
        #clear index for next calculation
        sequence_index=0
        #print('dist',distance_between_points(position1,position2))
        return round(closest_distance,4),position1,position2

def read_points():
    def str2float(config):
        assert(type(config)==str)
        p=list(config.strip('[]()').split(','))
        assert(type(p)==list)
        assert(type(p[0])==str)
        p=[float(i) for i in p]
        return p
    config=configparser.ConfigParser()
    config.read('coordinates.txt',encoding='ISO-8859-1')
    Default_config=config['Default']
    print(type(Default_config['p1']))
    p1=str2float(Default_config['p1'])
    p2=str2float(Default_config['p2'])
    p3=str2float(Default_config['p3'])
    p4=str2float(Default_config['p4'])
    p5=str2float(Default_config['p5'])
    p6=str2float(Default_config['p6'])
    return Point(*p1),Point(*p2),Point(*p3),Point(*p4),Point(*p5),Point(*p6)

point1,point2,point3,point4,point5,point6=read_points()
#三相导线对应挂点坐标组合 A相导线
#point1=Point(0,0,10)
#point2=Point(61.1,11.74,18)
#line1=CatenaryLine(point1,point2,gama,delta)

#B相导线
#point3=Point(0,2.2,10)
#point4=Point(61.1,11.84,26.2)
#line2=CatenaryLine(point3,point4,gama,delta)

#C相导线
#point5=Point(0,-2.2,10)
#point6=Point(61.1,11.24,22)
#line3=CatenaryLine(point5,point6,gama,delta)

#根据不同应力（弧垂）来批量计算最短距离
def construct_lines_and_find_distance(p1,p2,p3,p4,p5,p6,gama,delta):
    line1=CatenaryLine(p1,p2,gama,delta)
    line2=CatenaryLine(p3,p4,gama,delta)
    line3=CatenaryLine(p5,p6,gama,delta)
    distanceAtoB,poAtoB1,poAtoB2=Find_closest_distance(line1,line2,0,line1.x_lenth())
    distanceAtoC,poAtoC1,poAtoC2=Find_closest_distance(line1,line3,0,line1.x_lenth())
    distanceBtoC,poBtoC1,poBtoC2=Find_closest_distance(line2,line3,0,line2.x_lenth())
    print('AB相之间距离: ',distanceAtoB,'\n','最近点坐标：',poAtoB1,poAtoB2,'\n')
    print('AC相之间距离: ',distanceAtoC,'\n','最近点坐标：',poAtoC1,poAtoC2,'\n')
    print('BC相之间距离: ',distanceBtoC,'\n','最近点坐标：',poBtoC1,poBtoC2,'\n')
    return [gama,delta,distanceAtoB,poAtoB1,poAtoB2,distanceAtoC,poAtoC1,poAtoC2,distanceBtoC,poBtoC1,poBtoC2]

#批量计算不同应力下的相间距
def batch_cal():
    import csv
    with open('output.csv','w') as f:
        writer=csv.writer(f)
        writer.writerow(['比载','应力','AB间距','A相坐标','B相坐标','AC间距','A相坐标','C相坐标','BC间距','B相坐标','C相坐标'])
        for d in delta_sequence:
            writer.writerow(construct_lines_and_find_distance(point1,point2,point3,point4,point5,point6,gama,d))

def verify(delta):
    line1=CatenaryLine(point1,point2,gama,delta)
    line2=CatenaryLine(point3,point4,gama,delta)
    line3=CatenaryLine(point5,point6,gama,delta)
    #print(line1.calculate_point(0),line1.calculate_point(60))
    import time
    start=time.clock()
    distanceAtoB,poAtoB1,poAtoB2=Find_closest_distance(line1,line2,0,line1.x_lenth())
    distanceAtoC,poAtoC1,poAtoC2=Find_closest_distance(line1,line3,0,line1.x_lenth())
    distanceBtoC,poBtoC1,poBtoC2=Find_closest_distance(line2,line3,0,line2.x_lenth())
    elapsed=(time.clock()-start)
    print('time elapsed:',elapsed)
    print('AB相之间距离: ',distanceAtoB,'\n','最近点坐标：',poAtoB1,poAtoB2)
    print('A相弧垂：',line1.calculate_fm(),'\n')
    print('AC相之间距离: ',distanceAtoC,'\n','最近点坐标：',poAtoC1,poAtoC2)
    print('B相弧垂：',line2.calculate_fm(),'\n')
    print('BC相之间距离: ',distanceBtoC,'\n','最近点坐标：',poBtoC1,poBtoC2)
    print('C相弧垂：',line3.calculate_fm(),'\n')
    print(delta,gama)
    def draw():
        step=0.1
        x=np.arange(0,line1.x_lenth()+step,step)
        import platform
        if platform.system()=='Windows':
            zhfont = FontProperties(fname='C:\Windows\Fonts\simsun.ttc',size=14)
            fig,ax=plt.subplots()
            ax = fig.add_subplot(111, projection='3d')
            plt.title('相间距',fontproperties=zhfont)
            ax.set_xlabel('  单位：m',fontproperties=zhfont)
            ax.set_ylabel('  单位：m',fontproperties=zhfont)
            ax.set_zlabel('  单位：m',fontproperties=zhfont)
        if platform.system()=='Darwin':
            fig,ax=plt.subplots()
            ax = fig.add_subplot(111, projection='3d')
            plt.title('相间距')
            ax.set_xlabel('  单位：m')
            ax.set_ylabel('  单位：m')
            ax.set_zlabel('  单位：m')
        y1,z1=line1.calculate_yzs(x)
        y2,z2=line2.calculate_yzs(x)
        y3,z3=line3.calculate_yzs(x)
        ax.plot(x,y1,z1)
        ax.plot(x,y2,z2)
        ax.plot(x,y3,z3)

        #画出最近点
        a=[poAtoB1.x,poAtoB2.x]
        b=[poAtoB1.y,poAtoB2.y]
        c=[poAtoB1.z,poAtoB2.z]
        ax.plot(a,b,c)

        a=[poAtoC1.x,poAtoC2.x]
        b=[poAtoC1.y,poAtoC2.y]
        c=[poAtoC1.z,poAtoC2.z]
        ax.plot(a,b,c)

        a=[poBtoC1.x,poBtoC2.x]
        b=[poBtoC1.y,poBtoC2.y]
        c=[poBtoC1.z,poBtoC2.z]
        ax.plot(a,b,c)

        #画挂点连线 
        # a=[point1.x,point2.x]
        # b=[point1.y,point2.y]
        # c=[point1.z,point2.z]
        # ax.plot(a,b,c)

        # a=[point4.x,point3.x]
        # b=[point4.y,point3.y]
        # c=[point4.z,point3.z]
        # ax.plot(a,b,c)

        # a=[point5.x,point6.x]
        # b=[point5.y,point6.y]
        # c=[point5.z,point6.z]
        # ax.plot(a,b,c)

        #y轴刻度翻转
        ax.set_ylim(ax.get_ylim()[::-1])
        #print(line1.calculate_point(L),line2.calculate_point(L),line3.calculate_point(L))
        #
        #用于验证理论计算出的挂点高度
        print(z1[-1],z2[-1],z3[-1])
        plt.savefig('line.svg')
        plt.show()
    draw()

#batch_cal()
#用来验证
verify(10)
# for i in delta_sequence:
#     verify(i)

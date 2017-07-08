# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 10:48:20 2017

@author: 郭爽
"""

from tkinter import *
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import numpy as np
from array import array
from scipy import integrate
import tkinter.messagebox as messagebox
import tkinter as tk

global l
global I
global E
global fangdabeishu

l=0
I=0
E=0
fangdabeishu=1000

F=[]
a=[]
M=[]
b=[]
q=[]
x1=[]
x2=[]
unitload_position=[]

#矩形截面计算器
class juxingCalc:
    """docstring for jiheCalc"""
    def __init__(self, daoshuju1):
        self.daoshuju_juxing=daoshuju1
        juxing_window = Tk()
        juxing_window.title('矩形截面计算器')
        juxing_window.geometry('450x230')

        juxing_frame = Frame(juxing_window)
        juxing_frame.pack()

        self.h_label=Label(juxing_frame,text='截面高度h(mm)')
        self.b_label=Label(juxing_frame,text='截面宽度b(mm)')

        self.h_Input=Entry(juxing_frame,width=15,bd=2)
        self.b_Input=Entry(juxing_frame,width=15,bd=2)

        self.juxing_calculate=Button(juxing_frame,text='计算',command=self.juxing_jisuan,bd=3)
        self.juxing_clear=Button(juxing_frame,text='清除',command=self.juxing_qingchu,bd=3)

        self.canvas=Canvas(juxing_frame,bg='white',width=200,height=200)
        self.canvas.create_rectangle(35,30,155,190,fill='blue')
        self.canvas.create_line(35,20,155,20,arrow=tk.BOTH)
        self.canvas.create_line(165,30,165,190,arrow=tk.BOTH)
        self.canvas.create_line(35,25,35,15,width=1)
        self.canvas.create_line(155,25,155,15,width=1)
        self.canvas.create_line(160,30,170,30,width=1)
        self.canvas.create_line(160,190,170,190,width=1)
        self.canvas.create_text((90,13),text='b',anchor=W)
        self.canvas.create_text((170,110),text='h',anchor=W)

        self.canvas.grid(row=0,rowspan=4,column=2,columnspan=2,padx=10,pady=10)

        self.h_label.grid(row=1,column=0,padx=3,sticky=W)
        self.b_label.grid(row=0,column=0,padx=3,sticky=W)

        self.h_Input.grid(row=1,column=1,pady=5,sticky=W)
        self.b_Input.grid(row=0,column=1,pady=5,sticky=W)

        self.juxing_calculate.grid(row=3,column=0,ipadx=20,ipady=5,padx=10,pady=5)
        self.juxing_clear.grid(row=3,column=1,ipadx=20,ipady=5,padx=10,pady=5)

        juxing_window.mainloop()

    def juxing_jisuan(self):
        if (self.b_Input.get()=='' or self.h_Input.get()=='0'):
            messagebox.showinfo('错误提示','请输入完整的截面数据')
        else:
            b=float(self.b_Input.get())/1000
            h=float(self.h_Input.get())/1000
            self.juxing_I=b*h**3/12
            self.daoshuju_juxing.I_Input.insert(0,self.juxing_I)
            
    def juxing_qingchu(self):
        self.h_Input.delete(0,END)
        self.b_Input.delete(0,END)

#圆形截面计算器
class yuanxingCalc:
    def __init__(self,daoshuju2):
        self.daoshuju_yuanxing=daoshuju2
        yuanxing_window = Tk()
        yuanxing_window.title('圆形截面计算器')
        yuanxing_window.geometry('450x230')

        yuanxing_frame = Frame(yuanxing_window)
        yuanxing_frame.pack()

        self.zhijing_label=Label(yuanxing_frame,text='截面直径D(mm)')

        self.zhijing_Input=Entry(yuanxing_frame,width=15,bd=2)

        self.yuanxing_calculate=Button(yuanxing_frame,text='计算',command=self.yuanxing_jisuan,bd=3)
        self.yuanxing_clear=Button(yuanxing_frame,text='清除',command=self.yuanxing_qingchu,bd=3)

        self.canvas=Canvas(yuanxing_frame,bg='white',width=200,height=200)
        self.canvas.create_oval(30,35,170,175,fill='orange')

        self.canvas.create_line(30,20,170,20,arrow=tk.BOTH)
        self.canvas.create_line(30,15,30,25,width=1)
        self.canvas.create_line(170,15,170,25,width=1)
        self.canvas.create_text((95,13),text='D',anchor=W)

        self.canvas.grid(row=0,rowspan=4,column=2,columnspan=2,padx=10,pady=10)

        self.zhijing_label.grid(row=0,rowspan=2,column=0,padx=3,pady=5,sticky=W)
        self.zhijing_Input.grid(row=0,rowspan=2,column=1,pady=5,sticky=W)

        self.yuanxing_calculate.grid(row=2,rowspan=2,ipadx=20,ipady=5,column=0,padx=10,pady=5)
        self.yuanxing_clear.grid(row=2,rowspan=2,column=1,ipadx=20,ipady=5,padx=10,pady=5)

        yuanxing_window.mainloop()

    def yuanxing_jisuan(self):
        if (self.zhijing_Input.get()==''):
            messagebox.showinfo('错误提示','请输入完整的截面数据')
        else:
            zhijing=float(self.zhijing_Input.get())/1000
            self.yuanxing_I=3.1415926/64*(zhijing**4)
            self.daoshuju_yuanxing.I_Input.insert(0,self.yuanxing_I)

    def yuanxing_qingchu(self):
        self.zhijing_Input.delete(0,END)

#薄壁圆梁截面计算器
class yuantongCalc:
    def __init__(self,daoshuju3):
        self.daoshuju_yuantong=daoshuju3
        yuantong_window = Tk()
        yuantong_window.title('薄壁圆梁截面计算')
        yuantong_window.geometry('450x230')

        yuantong_frame = Frame(yuantong_window)
        yuantong_frame.pack()

        self.D_label=Label(yuantong_frame,text='截面外径D(mm)')
        self.d_label=Label(yuantong_frame,text='截面内径d(mm)')

        self.D_Input=Entry(yuantong_frame,width=15,bd=2)
        self.d_Input=Entry(yuantong_frame,width=15,bd=2)

        self.yuantong_calculate=Button(yuantong_frame,text='计算',command=self.yuantong_jisuan,bd=3)
        self.yuantong_clear=Button(yuantong_frame,text='清除',command=self.yuantong_qingchu,bd=3)

        self.canvas=Canvas(yuantong_frame,bg='white',width=200,height=200)
        self.canvas.create_oval(40,40,160,160,fill='yellow')
        self.canvas.create_oval(60,60,140,140,fill='white')
        self.canvas.create_line(40,30,160,30,arrow=tk.BOTH)
        self.canvas.create_line(40,25,40,35,width=1)
        self.canvas.create_line(160,25,160,35,width=1)
        self.canvas.create_line(60,170,140,170,arrow=tk.BOTH)
        self.canvas.create_line(60,165,60,175,width=1)
        self.canvas.create_line(140,165,140,175,width=1)
        self.canvas.create_line(60,170,60,100,dash=(4,4))
        self.canvas.create_line(140,170,140,100,dash=(4,4))
        self.canvas.create_text((95,22),text='D',anchor=W)
        self.canvas.create_text((95,178),text='d',anchor=W)

        self.canvas.grid(row=0,rowspan=4,column=2,columnspan=2,padx=10,pady=10)

        self.d_label.grid(row=1,column=0,padx=3,sticky=W)
        self.D_label.grid(row=0,column=0,padx=3,sticky=W)

        self.d_Input.grid(row=1,column=1,pady=5,sticky=W)
        self.D_Input.grid(row=0,column=1,pady=5,sticky=W)

        self.yuantong_calculate.grid(row=3,column=0,ipadx=20,ipady=5,padx=10,pady=5)
        self.yuantong_clear.grid(row=3,column=1,ipadx=20,ipady=5,padx=10,pady=5)

        yuantong_window.mainloop()

    def yuantong_jisuan(self):
        if (self.D_Input.get()=='' or self.d_Input.get()=='0'):
            messagebox.showinfo('错误提示','请输入完整的截面数据')
        else:
            D=float(self.D_Input.get())/1000
            d=float(self.d_Input.get())/1000
            self.yuantong_I=3.1415926/64*(D**4-d**4)
            self.daoshuju_yuantong.I_Input.insert(0,self.yuantong_I)
    def yuantong_qingchu(self):
        self.D_Input.delete(0,END)
        self.d_Input.delete(0,END)


#绘图弹出窗口
class Application(tk.Tk):

    def __init__(self):

        super().__init__()
        self.wm_title("上方为弯矩图，下方为挠曲线")
        #self.createWidgets()

    def createWidgets(self):


        fig = Figure(figsize=(8,6), dpi=100)
        self.ax1 = fig.add_subplot(211)
        self.ax2 = fig.add_subplot(212)
        self.canvas = FigureCanvasTkAgg(fig, master=self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        toolbar.update()
        footframe = tk.Frame(master=self).pack(side=tk.BOTTOM)
        
        self.draw()

#绘图算法      
    def draw(self):
        global l
        global I
        global E
        global fangdabeishu

#计算两端支座中的约束力        
        def Ra():
            global l
            global I
            global E
            sum_Ra=0
            for i in range(len(F)):
                sum_Ra=sum_Ra+F[i]*(l-a[i])/l
            for j in M:
                sum_Ra=sum_Ra+j/l
            return sum_Ra
        t=Ra()

#定义一阶奇异函数
        def force_Singular_Function(x,force_position):
            if x<=force_position:
                return 0
            else:
                return (x-force_position)

#定义零阶奇异函数
        def moment_Singular_Function(x,moment_position):
            if x<=moment_position:
                return 0
            else:
                return 1

#计算分布载荷单独作用时各点的弯矩
        def distributed_load(x,q,x1,x2):
            global l
            global I
            global E
            Ra=q*(x2-x1)*(2*l-x1-x2)/(2*l)
            Rb=q*(x2-x1)*(x1+x2)/(2*l)
            if x<x1:
                return Ra*x
            elif (x1<=x)  & (x<=x2):
                return Ra*x-q/2*(x-x1)**2
            else:
                return Rb*(l-x)

        t2=moment_Singular_Function
        t1=force_Singular_Function
        t3=distributed_load

#用奇异函数法和叠加法求所有外加载荷引起的总弯矩
        def curve_moment(x):
            global l
            global I
            global E
            curve=t*x
            for i in range(len(F)):
                curve=curve-t1(x,a[i])*F[i]
            for j in range(len(M)):
                curve=curve-t2(x,b[j])*M[j]
            for m in range(len(q)):
                curve=curve+t3(x,q[m],x1[m],x2[m])
            return curve

#单位力法，用莫尔积分计算等号右边的位移向量并存储在一个list中
        tureload=[]

        for x3 in unitload_position:
            def ture_load_Mohr(x):
                Ra3=(l-x3)/l
                Rb3=x3/l
                if x<x3:
                    m3=Ra3*x
                else:
                    m3=Rb3*(l-x)
                m_ture=curve_moment(x)
                return m3*m_ture
            c1,err=integrate.quad(ture_load_Mohr,0,l)
            c=c1/(E*I)
            tureload.append(-c)

#单位力法，用莫尔积分计算正则方程系数矩阵的元素并存储在另一个list中
        unitload=[]
        for u1 in unitload_position:
            for u2 in unitload_position:

                def unit_load_Mohr(x):
                    Ra1=(l-u1)/l
                    Rb1=u1/l
                    if x<u1:
                        m1=Ra1*x
                    else:
                        m1=Rb1*(l-x)
                    Ra2=(l-u2)/l
                    Rb2=u2/l
                    if x<u2:
                        m2=Ra2*x
                    else:
                        m2=Rb2*(l-x)
                    return m1*m2

                o,err=integrate.quad(unit_load_Mohr,0,l)
                h=o/(E*I)
                unitload.append(h)

#将存储系数的list转变为一个n阶矩阵
        coefficient=array("d")
        for i in unitload:
            coefficient.append(i)
        number=int((len(coefficient))**(1/2))
        matrix_of_coefficient=np.frombuffer(coefficient,dtype=np.float).reshape(number,-1)

        A=np.mat(matrix_of_coefficient)
        #print('A\n',A)

        B=np.array(tureload)
        #print('B\n',B)

#解正则方程        
        solution= np.linalg.solve(A, B)
        #print('Solution',solution)

#将多余约束力用叠加法加到总弯矩方程中，得到最终的总弯矩方程
        def sum_bending_moment(x):
            global l
            global I
            global E
            extra=curve_moment(x)
            def yueshu_moment(x,x4):
                Ra1=(l-x4)/l
                Rb1=x4/l
                if x<x4:
                    m5=Ra1*x
                else:
                    m5=Rb1*(l-x)
                return m5

            for i in range(len(solution)):
                extra=extra+solution[i]*yueshu_moment(x,unitload_position[i])
            return extra

#取点拟合弯矩图
        p1=np.linspace(0,l,200)
        wanju=[]
        for j in range(200):
            wanju.append(sum_bending_moment(p1[j]))

#单位力法，用莫尔积分计算梁上各点挠度
        def bending(x):
            global l
            global I
            global E
            def sum_Mohr(z):
                Ra4=(l-x)/l
                Rb4=x/l
                if z<x:
                    m6=Ra4*z
                else:
                    m6=Rb4*(l-z)
                m_sum=sum_bending_moment(z)
                return m6*m_sum
            s2,err=integrate.quad(sum_Mohr,0,l)
            return s2/(E*I)

#挠曲线取点拟合挠曲线
        p=np.linspace(0,l,200)
        v=[]
        for i in range(200):
            v.append(-fangdabeishu*bending(p[i]))

#绘图
        zeroX = np.linspace(0,l,200)
        zeroY = [0]*200

        self.ax1.clear()
        self.ax2.clear()

        self.ax1.plot(p1,wanju,label='$M(x)$',linewidth=1) 
        self.ax1.plot(zeroX,zeroY,color='black')
        self.ax1.set_xlabel('X / m',labelpad=3)
        self.ax1.set_ylabel('N(x) / N*m')

        self.ax2.plot(zeroX,zeroY,color='black')
        self.ax2.plot(p,v,label='$v(x)$',color='red',linewidth=1)
        self.ax2.set_xlabel('X / m')
        self.ax2.set_ylabel('V(x) / m')
        
        self.canvas.show()
        
#主界面
class CalcDemo:
    
    def __init__(self):
        #添加界面上的标签(Label)、按钮(Button)、文本框(Text)、输入框(Entry)、画布(Canvas)等控件，并以网格(grid)方式布局
        window = Tk()
        window.title("连续梁内力分析")
        window.geometry('1050x650')

        frame = Frame(window)
        frame.pack()

        menubar=Menu(window)
        filemenu = Menu(menubar, tearoff=0)
        filemenu.add_command(label='矩形截面',command=self.dakai_juxing)
        filemenu.add_command(label='圆形截面',command=self.dakai_yuanxing)
        filemenu.add_command(label='薄壁圆梁',command=self.dakai_yuantong)
        menubar.add_cascade(label="截面几何计算器",menu=filemenu)
        window.config(menu=menubar)

        self.canvas=Canvas(frame,bg="white",width=600,height=250)
        self.canvas.grid(row=0,rowspan=11,column=4,columnspan=6,padx=10,pady=10)

        self.canvas.create_line(40,125,560,125,fill="black",width=4)
        self.canvas.create_polygon((40,125,50,135,30,135))
        self.canvas.create_polygon((560,125,570,135,550,135))

        llabel=Label(frame,text='梁的长度l(m)')
        Elabel=Label(frame,text='弹性模量E(GPa)')
        Ilabel=Label(frame,text='惯性矩I(m^4)')
        
        fangdalabel=Label(frame,text='挠度放大倍数')

        Flabel=Label(frame,text='集中载荷F(N)')
        alabel=Label(frame,text='作用位置a(m)')
        Mlable=Label(frame,text='集中力偶M(N*m)')
        blabel=Label(frame,text='作用位置b(m)')
        qlabel=Label(frame,text='分布载荷q(N/m)')
        x1label=Label(frame,text='作用位置1(m)')
        x2label=Label(frame,text='作用位置2(m)')

        yueshulabel=Label(frame,text='多余约束位置(m)')
        
        
        self.matrix=Text(frame,width=85,height=26)

        self.l_Input=Entry(frame,width=15,bd=2)
        self.E_Input=Entry(frame,width=15,bd=2)
        self.I_Input=Entry(frame,width=15,bd=2)

        self.F_Input=Entry(frame,width=15,bd=2)
        self.a_Input=Entry(frame,width=15,bd=2)
        self.M_Input=Entry(frame,width=15,bd=2)
        self.b_Input=Entry(frame,width=15,bd=2)
        self.q_Input=Entry(frame,width=15,bd=2)
        self.x1_Input=Entry(frame,width=15,bd=2)
        self.x2_Input=Entry(frame,width=15,bd=2)

        self.yueshu_Input=Entry(frame,width=15,bd=2)
        
        self.fangda_Input=Entry(frame,width=15,bd=2)

        addforce = Button(frame, text = "添加", command = self.force_add,bd=3)
        deleteforce = Button(frame, text="清除", command=self.force_delete,bd=3)
        addmoment = Button(frame, text = "添加", command = self.moment_add,bd=3)
        deletemoment = Button(frame, text="清除", command=self.moment_delete,bd=3)
        adddistributed = Button(frame, text="添加", command=self.distributed_add,bd=3)
        deletedistributed = Button(frame, text="清除", command=self.distributed_delete,bd=3)
        addyueshu=Button(frame,text='添加约束',command=self.yueshu_add,bd=3)
        clearyueshu=Button(frame,text='清除约束',command=self.yueshu_clear,bd=3)
        analysis=Button(frame,text='开始分析',command=self.analyse,bg='orange',bd=3)
        diagram = Button(frame, text="内力图", command=self.new_window,bg='yellow',bd=3)

        llabel.grid(row=0,column=0,sticky=W)
        Elabel.grid(row=1,column=0,sticky=W)
        Ilabel.grid(row=2,column=0,sticky=W)

        Flabel.grid(row=5,column=0,padx=3,sticky=W)
        alabel.grid(row=6,column=0,padx=3,sticky=W)
        Mlable.grid(row=11,column=0,padx=3,sticky=W)
        blabel.grid(row=12,column=0,padx=3,sticky=W)
        qlabel.grid(row=13,column=0,padx=3,sticky=W)
        x1label.grid(row=14,column=0,padx=3,sticky=W)
        x2label.grid(row=15,column=0,padx=3,sticky=W)

        yueshulabel.grid(row=16,column=0,sticky=W)
        
        fangdalabel.grid(row=3,column=0,columnspan=2,sticky=W)

        self.l_Input.grid(row=0,column=1,pady=5,sticky=W)
        self.E_Input.grid(row=1,column=1,pady=5,sticky=W)
        self.I_Input.grid(row=2,column=1,pady=5,sticky=W)

        self.F_Input.grid(row=5,column=1,pady=5,sticky=W)
        self.a_Input.grid(row=6,column=1,pady=5,sticky=W)
        self.M_Input.grid(row=11,column=1,pady=5,sticky=W)
        self.b_Input.grid(row=12,column=1,sticky=W)
        self.q_Input.grid(row=13,column=1,sticky=W)
        self.x1_Input.grid(row=14,column=1,sticky=W)
        self.x2_Input.grid(row=15,column=1,sticky=W)

        self.yueshu_Input.grid(row=16,column=1,pady=5,sticky=W)
        self.fangda_Input.grid(row=3,column=1,pady=5,sticky=W)


        addforce.grid(row=5, column=2,padx=10,pady=3,ipadx=10)
        deleteforce.grid(row=6, column=2,padx=10,pady=3,ipadx=10,ipady=1)
        addmoment.grid(row=11, column=2,padx=10,pady=3,ipadx=10,ipady=1)
        deletemoment.grid(row=12, column=2,padx=10,pady=3,ipadx=10,ipady=1)
        adddistributed.grid(row=13, column=2,padx=10,ipadx=10,ipady=1)
        deletedistributed.grid(row=14, column=2,padx=10,ipadx=10,ipady=1)

        addyueshu.grid(row=17,column=0,padx=20,pady=10,ipadx=10,ipady=1)
        clearyueshu.grid(row=17,column=1,padx=20,pady=10,ipadx=10,ipady=1)
        
        self.matrix.grid(row=11,rowspan=9,column=4,columnspan=4,padx=10,pady=10)
        
        analysis.grid(row=18,column=0,columnspan=1,padx=10,pady=3,ipadx=50,ipady=5)
        diagram.grid(row=18,column=1,columnspan=2,padx=10,pady=3,ipadx=50,ipady=5)      

        window.mainloop()

#添加集中载荷及其作用位置
    def force_add(self):
        jianchachangdu_F=float(self.l_Input.get())
        waijiali=self.F_Input.get()
        lizuoyongweizhi=self.a_Input.get()
        if (waijiali=="" or lizuoyongweizhi==""):
            self.F_Input.delete(0,END)
            self.a_Input.delete(0,END)
            messagebox.showinfo('错误提示','集中载荷与作用位置必须成对添加')
        elif float(jianchachangdu_F)<float(lizuoyongweizhi):
            messagebox.showinfo('错误提示','作用位置超出梁的长度范围')
            self.a_Input.delete(0,END)
        else:
            F.append(float(waijiali))
            a.append(float(lizuoyongweizhi))
            self.huajiantou()
            
#清除集中载荷及其作用位置
    def force_delete(self):
        self.F_Input.delete(0,END)
        F[:]=[]
        self.a_Input.delete(0,END)
        a[:]=[]
        self.clearline()

#输入集中力偶及其作用位置
    def moment_add(self):

        waijialiou=self.M_Input.get()
        jianchachangdu_M=float(self.l_Input.get())
        liouzuoyongweizhi=self.b_Input.get()
        if (waijialiou=="" or liouzuoyongweizhi==""):
            self.M_Input.delete(0,END)
            self.b_Input.delete(0,END)
            messagebox.showinfo('错误提示','集中力偶与作用位置必须成对添加')
        elif float(liouzuoyongweizhi)>jianchachangdu_M:
            messagebox.showinfo('错误提示','作用位置超出梁的长度范围')
            self.b_Input.delete(0,END)
        else:
            M.append(float(waijialiou))
            b.append(float(liouzuoyongweizhi))
            self.huawanju()

#清除集中力偶及其作用位置
    def moment_delete(self):
        self.M_Input.delete(0,END)
        M[:]=[]
        self.b_Input.delete(0,END)
        b[:]=[]
        self.clearwanju()

#添加分布载荷及其作用位置
    def distributed_add(self):
        changdujiancha_q=float(self.l_Input.get())
        waijiafenbuli=self.q_Input.get()
        zuoyongweizhi1=self.x1_Input.get()
        zuoyongweizhi2=self.x2_Input.get()
        if (waijiafenbuli=="" or zuoyongweizhi1=="" or zuoyongweizhi2==""):
            self.q_Input.delete(0,END)
            self.x1_Input.delete(0,END)
            self.x2_Input.delete(0,END)
            messagebox.showinfo('错误提示','分布载荷与作用位置1、2必须成对添加')
        elif (float(zuoyongweizhi1)>float(zuoyongweizhi2)):
            self.x1_Input.delete(0,END)
            self.x2_Input.delete(0,END)
            messagebox.showinfo('错误提示','作用位置1应小于作用位置2')
        elif (float(zuoyongweizhi1)>changdujiancha_q or float(zuoyongweizhi2)>changdujiancha_q):
            self.x1_Input.delete(0,END)
            self.x2_Input.delete(0,END)
            messagebox.showinfo('错误提示','作用范围超出梁的长度范围')
        else:
            q.append(float(waijiafenbuli))
            x1.append(float(zuoyongweizhi1))
            x2.append(float(zuoyongweizhi2))
            self.huajuxing()

#清除分布载荷及其作用位置
    def distributed_delete(self):
        self.q_Input.delete(0,END)
        self.x1_Input.delete(0,END)
        self.x2_Input.delete(0,END)
        q[:]=[]
        x1[:]=[]
        x2[:]=[]
        self.clearjuxing()

#添加多余约束
    def yueshu_add(self):
        yueshuweizhi=self.yueshu_Input.get()
        changdujiancha_yueshu=float(self.l_Input.get())
        if yueshuweizhi=="":
            self.yueshu_clear()
            messagebox.showinfo('错误提示','请输入多余约束位置')
        elif float(yueshuweizhi)>changdujiancha_yueshu:
            self.yueshu_Input.delete(0,END)
            messagebox.showinfo('错误提示','约束位置超出梁的范围')
        elif float(yueshuweizhi) in unitload_position:
            self.yueshu_Input.delete(0,END)
            messagebox.showinfo('错误提示','此处已经添加过支座，同一位置请勿重复添加')
        else:
            unitload_position.append(float(yueshuweizhi))
            self.huazhizuo()

#清除多余约束
    def yueshu_clear(self):
        self.yueshu_Input.delete(0,END)
        unitload_position[:]=[]
        self.clearzhizuo()

#点击“开始分析”按钮：计算并输出正则方程的系数矩阵、右边的位移向量和解向量（多余约束力），算法细节见“绘图算法”
    def analyse(self):
        if len(unitload_position)==0:
            messagebox.showinfo('错误提示','连续梁至少添加一处多余约束')
            return -1

        global l
        global I
        global E
        #print(F)
        #print(a)
        length=self.l_Input.get()
        l=float(length)
        tanxingmoliang=self.E_Input.get()
        E=1000000000*float(tanxingmoliang)
        guanxingju=self.I_Input.get()
        I=float(guanxingju)

        def Ra():
            sum_Ra=0
            for i in range(len(F)):
                sum_Ra=sum_Ra+F[i]*(l-a[i])/l
            for j in M:
                sum_Ra=sum_Ra+j/l
            return sum_Ra
        t=Ra()


        def force_Singular_Function(x,force_position):
            if x<=force_position:
                return 0
            else:
                return (x-force_position)

        def moment_Singular_Function(x,moment_position):
            if x<=moment_position:
                return 0
            else:
                return 1

        def distributed_load(x,q,x1,x2):
            Ra=q*(x2-x1)*(2*l-x1-x2)/(2*l)
            Rb=q*(x2-x1)*(x1+x2)/(2*l)
            if x<x1:
                return Ra*x
            elif (x1<=x)  & (x<=x2):
                return Ra*x-q/2*(x-x1)**2
            else:
                return Rb*(l-x)

        t2=moment_Singular_Function
        t1=force_Singular_Function
        t3=distributed_load

        def curve_moment(x):
            curve=t*x
            for i in range(len(F)):
                curve=curve-t1(x,a[i])*F[i]
            for j in range(len(M)):
                curve=curve-t2(x,b[j])*M[j]
            for m in range(len(q)):
                curve=curve+t3(x,q[m],x1[m],x2[m])
            return curve
        
        tureload=[]

        for x3 in unitload_position:
            def ture_load_Mohr(x):
                Ra3=(l-x3)/l
                Rb3=x3/l
                if x<x3:
                    m3=Ra3*x
                else:
                    m3=Rb3*(l-x)
                m_ture=curve_moment(x)
                return m3*m_ture
            c1,err=integrate.quad(ture_load_Mohr,0,l)
            c=c1/(E*I)
            tureload.append(-c)
    
        unitload=[]
        for u1 in unitload_position:
            for u2 in unitload_position:

                def unit_load_Mohr(x):
                    Ra1=(l-u1)/l
                    Rb1=u1/l
                    if x<u1:
                        m1=Ra1*x
                    else:
                        m1=Rb1*(l-x)
                    Ra2=(l-u2)/l
                    Rb2=u2/l
                    if x<u2:
                        m2=Ra2*x
                    else:
                        m2=Rb2*(l-x)
                    return m1*m2

                o,err=integrate.quad(unit_load_Mohr,0,l)
                h=o/(E*I)
                unitload.append(h)

        coefficient=array("d")
        for i in unitload:
            coefficient.append(i)
        number=int((len(coefficient))**(1/2))
        matrix_of_coefficient=np.frombuffer(coefficient,dtype=np.float).reshape(number,-1)

        A=np.mat(matrix_of_coefficient)
        #print('A\n',A)
        self.matrix.insert(END, "正则方程组系数矩阵：\n")
        self.matrix.insert(END, A)
        self.matrix.insert(END, "\n")

        B=np.array(tureload)
        #print('B\n',B)
        self.matrix.insert(END,"正则方程等号右边的矩阵：\n")
        self.matrix.insert(END,B)
        self.matrix.insert(END,"\n")

        solution= np.linalg.solve(A, B)
        #print('Solution',solution)
        self.matrix.insert(END,"各个支座上的支反力（按照添加顺序排序，向下为正）：\n")
        self.matrix.insert(END,solution)
        self.matrix.insert(END,"\n")

    def dakai_juxing(self):
        juxingCalc(self)

    def dakai_yuanxing(self):
        yuanxingCalc(self)

    def dakai_yuantong(self):
        yuantongCalc(self)

#点击内力图按钮：在弹出的新窗口内显示内力图
    def new_window(self):
        global fangdabeishu
        fangdabeishu=float(self.fangda_Input.get())
        app = Application()
        app.createWidgets()

    def huajiantou(self):
        lidezhengfu=self.F_Input.get()
        chang=self.l_Input.get()
        shijichang=float(chang)
        jiantouweizhi=float(self.a_Input.get())
        if float(lidezhengfu)>0:
            self.canvas.create_line(40+jiantouweizhi/shijichang*520,50,40+jiantouweizhi/shijichang*520,125, fill = "red", arrow = tk.LAST,tag="line")
        else:
            self.canvas.create_line(40+jiantouweizhi/shijichang*520,200,40+jiantouweizhi/shijichang*520,125, fill = "red", arrow = tk.LAST,tag="line")

    def clearline(self):
        self.canvas.delete("line")
        
    def huazhizuo(self):
        chang=self.l_Input.get()
        shijichang=float(chang)
        zhizuoweizhi=float(self.yueshu_Input.get())
        self.canvas.create_polygon((40+zhizuoweizhi/shijichang*520,125,50+zhizuoweizhi/shijichang*520,135,zhizuoweizhi/shijichang*520+30,135),fill='blue',tags = "polygon")
        
    def clearzhizuo(self):
        self.canvas.delete("polygon")
    
    def huawanju(self):
        liouzhengfu=float(self.M_Input.get())
        chang=self.l_Input.get()
        shijichang=float(chang)
        zhexianweizhi=float(self.b_Input.get())
        if liouzhengfu>0:
            self.canvas.create_line(40+zhexianweizhi/shijichang*520,75,40+zhexianweizhi/shijichang*520,175,fill='green',tag="zhexian")
            self.canvas.create_line(40+zhexianweizhi/shijichang*520,175,80+zhexianweizhi/shijichang*520,175,fill='green',arrow=tk.LAST,tag="zhexian")
            self.canvas.create_line(40+zhexianweizhi/shijichang*520,75,0+zhexianweizhi/shijichang*520,75,fil='green',arrow=tk.LAST,tag="zhexian")
        else:
            self.canvas.create_line(40+zhexianweizhi/shijichang*520,75,40+zhexianweizhi/shijichang*520,175,fill='green',tag="zhexian")
            self.canvas.create_line(40+zhexianweizhi/shijichang*520,175,0+zhexianweizhi/shijichang*520,175,fill='green',arrow=tk.LAST,tag="zhexian")
            self.canvas.create_line(40+zhexianweizhi/shijichang*520,75,80+zhexianweizhi/shijichang*520,75,fil='green',arrow=tk.LAST,tag="zhexian")
        
    def clearwanju(self):
        self.canvas.delete("zhexian")
    
    def huajuxing(self):
        chang=self.l_Input.get()
        shijichang=float(chang)
        juxingweizhi1=float(self.x1_Input.get())
        juxingweizhi2=float(self.x2_Input.get())
        self.canvas.create_rectangle(40+juxingweizhi1/shijichang*520,100,40+juxingweizhi2/shijichang*520,123,fill='yellow',tag="rectangle")
        
    def clearjuxing(self):
        self.canvas.delete("rectangle")

CalcDemo()
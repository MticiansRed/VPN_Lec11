import streamlit as st
from PIL import Image
 
menu = st.sidebar.radio('***',
    (
    "Краевая задача", 
    "Аппроксимация расчетной области", 
    "Лагранжевые конечные элементы",
    "Интегральное тождество",
    "Вариационная задача",
    )
)
  
if menu == "Краевая задача":
    r"""
##### Краевая задача

**Уравнение**

$\begin{aligned} 
\mathcal{L} u \equiv - 
 \operatorname{div} \big (k(x) \operatorname{grad} u) + c(x) u = f(x) ,
 \quad x = (x_1, x_2) \in \Omega
\end{aligned}$

с переменными коэффициентами

$\begin{aligned}
  k(x) \geq \kappa > 0,
  \quad c(x) \geq 0
\end{aligned}$

**Смешанные граничные условия**

Граница 

$\begin{aligned}
  \partial \Omega = \Gamma_D \cup \Gamma_N
\end{aligned}$

Условия Дирихле

$\begin{aligned}
  u(x) = \mu(x), 
  \quad x \in \Gamma_D
\end{aligned}$

Условия Неймана

$\begin{aligned}
  k(x) \frac{\partial u}{\partial n} =q(x),
  \quad x \in \Gamma_N
\end{aligned}$

    """    
    
if menu == "Аппроксимация расчетной области":
    r"""
##### Аппроксимация расчетной области

Расчетная область 

$\begin{aligned}
 \Omega, 
 \quad \partial \Omega
\end{aligned}$ 

Сетка  

$\begin{aligned}
   \overline{\omega} = \omega \cup \partial\omega = \{x~ |~ x = \bar{x}_\alpha , \ \bar{x}_\alpha \in \Omega \cup \partial \Omega, \ \alpha  = 0,1,\ldots,m\}
\end{aligned}$

* $\omega$ - множество внутренних узлов ($\bar{x}_\alpha \in \Omega$) 
* $\partial\omega$ - множество граничных узлов ($\bar{x}_\alpha \in \partial\Omega$)

Ячейки

$\begin{aligned}
 \Omega \approx \bigcup_{\beta}^{m} \Omega_\beta
\end{aligned}$
 
    """  
    
    c1, c2 = st.columns([2,1])
    c1.write("$~$")
    image = Image.open("pages/figs/1.png")    
    c1.image(image)   
    
if menu == "Лагранжевые конечные элементы":
    r"""
##### Лагранжевые конечные элементы

**Кусочно линейная аппроксимация**

$y(x) = a_0 + a_1 x_1 + a_2 x_2$ - три узла аппроксимации

Узлы аппроксимации совпадают с узлами сетки 

    """ 
    
    c1, c2 = st.columns([2,1])
    c1.write("$~$")
    image = Image.open("pages/figs/2.png")    
    c1.image(image)   
    
    r"""
$~$

**Аппроксимация полиномом второго порядка**

$y(x) = a_0 + a_1 x_1 + a_2 x_2 + a_3 x_1 x_1 + a_4 x_1 x_2 + a_5 x_2 x_2$ 

Узлы аппроксимации
 
 * три узла сетки
 * три узла на серединах ребер ячейки

    """ 
    
    c1, c2 = st.columns([2,1])
    c1.write("$~$")
    image = Image.open("pages/figs/3.png")    
    c1.image(image)  

if menu == "Интегральное тождество":
    r"""
##### Интегральное тождество

**Теорема о дивергенции**

$\begin{aligned}
 \int_\Omega \operatorname{div} \bm u \, dx = \int_{\partial \Omega}\bm u \cdot \bm n \, d s 
\end{aligned}$

$\bm n ~-~$ внешняя нормаль к $\partial \Omega$

**Пространства**

* $V$ - пространство достаточно гладких функций
* Подпространство $V_D = \{ v \ | \ v \in V, \quad v(x) = \mu(x),
  \quad x \in \Gamma_D \}$ 
* Подпространство $V_0 = \{ v \ | \ v \in V, \quad v(x) = 0,
  \quad x \in \Gamma_D \}$ 

**Интегральное тождество**

$\begin{aligned}
(\mathcal{L} u - f, v) = 0 , 
\quad u \in V_D, \ v \in V_0
\end{aligned}$  

С учетом

$\begin{aligned}
& \operatorname{div} (\varphi \bm w) = \bm w \cdot \operatorname{grad} \varphi + \varphi \operatorname{div} \bm w \\
& \bm w = k(x) \operatorname{grad}u, \quad \varphi = v
\end{aligned}$
  
для $u \in V_D$ и $v \in V_0$  
  
$\begin{aligned}
 (\mathcal{L} u, v) &= - \int_{\Omega} \operatorname{div} (k(x) \operatorname{grad}u) v \, d x + \int_{\Omega}  c(x) u \, v \, dx = \\
 &= - \int_{\Omega} \big (\operatorname{div} (k(x) \operatorname{grad}u \, v ) - k(x) \operatorname{grad} u \operatorname{grad} v \big ) \, d x + \int_{\Omega}  c(x) u \, v \, dx = \\
 & = - {\color {red}  \int_{\partial \Omega} k(x) v \operatorname{grad} u \cdot \bm n \, ds} +  \int_{\Omega} k(x) \operatorname{grad} u \operatorname{grad} v \, d x + \int_{\Omega}  c(x) u \, v \, dx
\end{aligned}$

С учетом $v \in V_0$ и граничных условий Неймана

$\begin{aligned}
 \int_{\partial \Omega} k(x) v \operatorname{grad} u \cdot \bm n \, ds = \int_{\Gamma_N} q(x) v \, ds 
\end{aligned}$

Результат

$\begin{aligned}
\int_{\Omega} k(x) \operatorname{grad} u \operatorname{grad} v \, d x + \int_{\Omega}  c(x) u \, v \, dx = \int_{\Omega} f(x) v \, d x + \int_{\Gamma_N} q(x) v \, ds
\end{aligned}$

    """  
    
if menu == "Вариационная задача":
    r"""
##### Вариационная задача


Билинейная форма 

$\begin{aligned}
 a(u,v) = \int_{\Omega} k(x) \operatorname{grad} u \operatorname{grad} v \, d x + \int_{\Omega}  c(x) u \, v \, dx
\end{aligned}$

Линейная форма 

$\begin{aligned}
 l(v) = \int_{\Omega} f(x) v \, d x + \int_{\Gamma_N} q(x) v \, ds
\end{aligned}$

Вариационная задача: найти  $u \in V_D$ такую, что

$\begin{aligned}
 a(u,v) = l(v) ,
 \quad \forall v \in V_0
\end{aligned}$

    """   
    
if menu == "Тестовая задача":
    r"""
##### Тестовая задача

Уравнение

$\begin{aligned} - 
 \frac{d^2 u}{dx^2} +   c \frac{du}{dx} = 0,
  \quad 0 < x < 1
\end{aligned}$

с коэффициентом $c(x) = c = \mathrm{const}$

Граничные условия 

$\begin{aligned} 
  u(a) = 0,
  \quad u(b) = 1
\end{aligned}$

Особые случаи

+ $|c| \ll 1 ~~- ~$ преобладание диффузии
+ $|c| \gg 1 ~~- ~$ преобладание конвекции

    """ 

if menu == "Ключевые фрагменты кода (FEniCS)":
    r"""
##### Ключевые фрагменты кода (FEniCS)
     
**Граничные условия Дирихле**
    """    
    code = """  
def boundary(x, on_boundary):
    return on_boundary
bc = DirichletBC(V, Expression("x[0]", degree=p+2), boundary)
    """ 
    st.code(code, language="python")   
    r"""    
**Вариационная формулировка задачи**
    """    
    code = """  
f = Expression("0.", degree=p+2)
a = u.dx(0)*v.dx(0)*dx + cc*u.dx(0)*v*dx
L = f*v*dx
    """ 
    st.code(code, language="python")    
    r"""        
**Решение задачи**
    """    
    code = """  
w = Function(V)
solve(a == L, w, bc)    
    """ 
    st.code(code, language="python")       
    
if menu == "Матрица СЛАУ":
    r"""
##### Матрица СЛАУ

**Параметры задачи**

    """

    import matplotlib.pyplot as plt
    import numpy as np 
    from fenics import *
    from scipy.sparse import csr_matrix
    
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Число ячеек $m$")        
    m = c2.slider("", 10, 50, 20)
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Порядок полинома $p$")        
    p = c2.slider("", 1, 3, 1)  
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Коэффициент при конвективном слагаемом $c$")        
    c = c2.slider("", -10, 10, 1, 1)     
     
    mesh = IntervalMesh(m, 0, 1)
    xm = mesh.coordinates()
    ym = np.zeros((m+1), "float") 
    
    V = FunctionSpace(mesh, "CG", p)
    n = V.dim()-1
    
    u = TrialFunction(V)
    v = TestFunction(V)
    
    def boundary(x, on_boundary):
        return on_boundary
    bc = DirichletBC(V, Expression("x[0]", degree=p+2), boundary)
    
    f = Expression("0.", degree=p+2)
    cc = Constant(c)
    a = u.dx(0)*v.dx(0)*dx + cc*u.dx(0)*v*dx
    L = f*v*dx
    
    A, b = assemble_system(a, L, bc)

    mat = as_backend_type(A).mat()
    csr = csr_matrix(mat.getValuesCSR()[::-1], shape=mat.size)
    Ad = csr.toarray()
    
    r"""
**Портрет матрицы**
    """    
          
    fig1 = plt.figure(1)
    plt.spy(csr)        
 
    c1, c2, = st.columns([3,1]) 
    c1.pyplot(fig1) 
    
    r"""
**Элементы матрицы**
    """  
    
    st.table(Ad)
         
if menu == "Влияние конвективного переноса":
    r"""
##### Влияние конвективного переноса

Коэффициент $c$ увеличивается в 10 и 100 раз

    """

    import matplotlib.pyplot as plt
    import numpy as np 
    from fenics import *
    
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Число ячеек $m$")        
    m = c2.slider("", 10, 50, 20)
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Порядок полинома $p$")        
    p = c2.slider("", 1, 3, 1)  
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Коэффициент при конвективном слагаемом $c$")        
    c = c2.slider("", -10, 10, 1, 1)      
          
    fig1 = plt.figure(1)
    
    for kk in range(0, 3): 
    
        mesh = IntervalMesh(m, 0, 1)
        xm = mesh.coordinates()
        ym = np.zeros((m+1), "float") 
        
        V = FunctionSpace(mesh, "CG", p)
        n = V.dim()-1
        
        u = TrialFunction(V)
        v = TestFunction(V)
        
        def boundary(x, on_boundary):
            return on_boundary
        bc = DirichletBC(V, Expression("x[0]", degree=p+2), boundary)
        
        f = Expression("0.", degree=p+2)
        cc = Constant(c)
        a = u.dx(0)*v.dx(0)*dx + cc*u.dx(0)*v*dx
        L = f*v*dx
    
        w = Function(V)
        solve(a == L, w, bc)
        
        N = 500
        xx = np.linspace(0., 1., N) 
        yy = np.linspace(0., 1., N)  
        
        for i in range(0, N): 
            yy[i]  = w(Point(xx[i]))
        s = "$c = $" + str(c)
        plt.plot(xx, yy, label = s) 
        c = c*10 
    
    ss = "$m = $" + str(m) + "$, \ p = $" + str(p) 
    plt.title(ss)
    plt.xlabel('$x$') 
    plt.legend(loc=0)
    plt.grid(True)      
    
    c1, c2, = st.columns([3,1]) 
    c1.pyplot(fig1)  
 

   

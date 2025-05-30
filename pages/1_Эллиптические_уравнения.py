import streamlit as st
 
menu = st.sidebar.radio('***',
    (
    "Уравнение Пуассона", 
    "Операторы векторного анализа", 
    "Формулы векторного анализа",
    "Элементы тензорного анализа",
    "Эллиптическое уравнение второго порядка",
    "Краевые условия",
    )
)

if menu == "Уравнение Пуассона":
    r"""
##### Уравнение Пуассона

**Основные математические модели**

* эллиптические уравнения второго порядка (стационарные)
* параболические уравнения второго порядка (нестационарные)

**Двумерные задачи**

* $x = (x_1, x_2)$
* $\Omega ~-~$ ограниченная двумерная область

**Уравнение Пуассона в декартовых координатах**

$\begin{aligned} -
  \frac{\partial^2 u }{\partial x_1^2} - \frac{\partial^2 u }{\partial x_2^2} = f(x),
  \quad x \in \Omega
\end{aligned}$

**Инвариантная запись**

$\begin{aligned} -
  \operatorname{div} \operatorname{grad}  u = f(x),
  \quad x \in \Omega
\end{aligned}$

или

$\begin{aligned} -
  \triangle  \, u = f(x),
  \quad x \in \Omega
\end{aligned}$

**Оператор Лапласа**

$\begin{aligned} 
  \triangle \equiv \operatorname{div} \operatorname{grad} 
\end{aligned}$

    """  
  
if menu == "Операторы векторного анализа":
    r"""
##### Операторы векторного анализа

**Обозначения**

*  $x = (x_1, \ldots, x_d) ~-~$ точка в $d$-мерной области
*  $u, v ~-~$ скалярные функции
*  $\bm u, \bm v ~-~$ векторные функции

**Градиент**

$\begin{aligned} 
\operatorname{grad} \varphi = \left\{
\frac {\partial \varphi} {\partial x_1}, \dots, \frac {\partial \varphi} {\partial x_d} \right\} 
\end{aligned}$

**Дивергенция**

$\begin{aligned} 
 \operatorname{div} \bm u = \frac {\partial u_1} {\partial x_1} +
 \cdots + \frac {\partial u_d} {\partial x_d}
\end{aligned}$


**Ротор**

3D

$\begin{aligned} 
\operatorname{rot} \bm u =
\left \{ \frac {\partial u_3} {\partial x_2} - \frac {\partial u_2} {\partial x_3}, \ 
\frac {\partial u_1} {\partial x_3} - \frac {\partial u_3} {\partial x_1}, \
\frac {\partial u_2} {\partial x_1} - \frac {\partial u_1} {\partial x_2} \right \}
\end{aligned}$

    """    
    
if menu == "Формулы векторного анализа":
    r"""
##### Формулы векторного анализа

**Обозначения**

*  $\varphi ~-~$ скалярная функция
*  $\bm u, \bm v ~-~$ векторные функции

**Произведение функций**

$\begin{aligned} 
& \operatorname{div} (\varphi \bm u) = \bm u \cdot \operatorname{grad} \varphi + \varphi \operatorname{div} \bm u \\
& \operatorname{rot} (\varphi \bm u) = \operatorname{grad} \varphi \times \bm u + \varphi \operatorname{rot} \bm u \\
& \operatorname{grad} (\bm u \cdot \bm v ) = \bm v \operatorname{grad} \bm u + \bm u \operatorname{grad} \bm v + \bm u \times \operatorname{rot} \bm v + \bm v \times \operatorname{rot} \bm u \\
& \operatorname{div} (\bm u \times \bm v) = \bm v \cdot \operatorname{rot} \bm u - \bm u \cdot \operatorname{rot} \bm v \\
& \operatorname{rot} (\bm u \times \bm v) = \bm u \operatorname{div} \bm v - \bm v \operatorname{div} \bm u + \bm u \operatorname{grad} \bm v - \bm v \operatorname{grad} \bm u
\end{aligned}$

**Повторные операции**

$\begin{aligned} 
& \operatorname{rot} \operatorname{grad} \varphi = 0 \\
& \operatorname{div} \operatorname{rot} \bm u = 0 \\
& \operatorname{rot} \operatorname{rot} \bm u = \operatorname{grad} \operatorname{div} \bm u - \operatorname{div} \operatorname{grad} \bm u
\end{aligned}$

    """  
    
if menu == "Элементы тензорного анализа":
    r"""
##### Элементы тензорного анализа

**Тензорные функции**

* нулевого ранга: $\varphi(x) ~-~$ скаляр
* первого ранга: $\bm u(x) = \{u_1(x), \ldots, u_d(x)\} ~-~$ вектор
* второго ранга: $K(x) = \{K_{\alpha \beta}(x)\}, \quad \alpha, \beta = 1, \ldots, d ~-~$ квадратная матрица

**Оператор Лапласа для векторной функции**

$\begin{aligned} 
\operatorname{div} \operatorname{grad} \bm u =
\{\operatorname{div} \operatorname{grad} u_1, \ldots, \operatorname{div} \operatorname{grad} u_d \}
\end{aligned}$

**Градиент вектора**

$\begin{aligned}
\operatorname{grad} \bm u = \left \{ \frac {\partial u_\alpha} {\partial x_\beta}\right \} =
\begin{pmatrix}
   \displaystyle \frac {\partial u_1} {\partial x_1} & \cdots & \displaystyle \frac {\partial u_1} {\partial x_d}
   \\
   \cdots & \cdots & \cdots
   \\
   \displaystyle \frac {\partial u_d} {\partial x_1} & \cdots & \displaystyle \frac {\partial u_d} {\partial x_d}
\end{pmatrix}
\end{aligned}$

**Дивергенция тензора второго ранга**

$\begin{aligned}
K(x) =
\begin{pmatrix}
   K_{11} (x) & \cdots & K_{1d} (x) \\
   \cdots & \cdots & \cdots \\
   K_{d1} (x) & \cdots & K_{dd} (x)
\end{pmatrix}
\end{aligned}$

$\begin{aligned}
\operatorname{div} K = \left \{\frac {\partial K_{11}} {\partial x_1} + \cdots + \frac {\partial K_{1d}} {\partial x_d} , \ \cdots, \  \frac {\partial K_{d1}} {\partial x_1} + \cdots + \frac {\partial K_{dd}} {\partial x_d} \right \} 
\end{aligned}$

    """   
    
if menu == "Эллиптическое уравнение второго порядка":
    r"""
##### Эллиптическое уравнение второго порядка

**Дивергентная форма**

$\begin{aligned} - 
 \operatorname{div} \big (k(x) \operatorname{grad} u + \bm a u \big ) +
 \bm b \cdot \operatorname{grad} u + c(x) u = f(x) 
\end{aligned}$

с коэффициентами $k(x) > 0, \ c(x) \geq 0$

**Конвективный перенос**

* $\bm a ~-~$ консервативный (дивергентый)
* $\bm b ~-~$ характеристический (недивергентый)

**Анизотропные среды**

$\begin{aligned} 
  k(x) \rightarrow K(x),
  \quad K(x) = 
  \begin{pmatrix}
   K_{11} (x) & \cdots & K_{1d} (x) \\
   \cdots & \cdots & \cdots \\
   K_{d1} (x) & \cdots & K_{dd} (x)
\end{pmatrix}
\end{aligned}$

    """ 

if menu == "Краевые условия":
    r"""
##### Краевые условия

**Внешняя нормаль к границе области**

$\bm n(x) = \{n_1(x), \ldots, n_d(x) \}, \quad x \in \partial \Omega$

**Производная по нормали**

$\begin{aligned} 
 \frac{\partial u}{\partial n} = \operatorname{grad} u \cdot \bm n =
  \frac{\partial u}{\partial x_1} n_1 + \cdots + \frac{\partial u}{\partial x_d} n_d 
\end{aligned}$
     
**Условия третьего рода**

$\begin{aligned} 
 k(x) \frac{\partial u}{\partial n} + \sigma (x) u = \mu(x),
 \quad x \in   \partial \Omega
\end{aligned}$

с коэффициентом $\sigma(x) \geq 0$

**Анизотропные среды**

$\begin{aligned} 
  k(x) \rightarrow K(x),
  \quad K(x) = 
  \begin{pmatrix}
   K_{11} (x) & \cdots & K_{1d} (x) \\
   \cdots & \cdots & \cdots \\
   K_{d1} (x) & \cdots & K_{dd} (x)
\end{pmatrix}
\end{aligned}$

Производная по конормали

$\begin{aligned} 
 & \frac{\partial u}{\partial \nu} = (K(x) \operatorname{grad} u) \cdot \bm n, \\ 
 & K(x) \operatorname{grad} u = \left \{ K_{11} (x) \frac{\partial u}{\partial x_1} + \cdots + K_{1d} (x) \frac{\partial u}{\partial x_d} , \ \ldots, \ K_{d1} (x) \frac{\partial u}{\partial x_1} + \cdots + K_{dd} (x) \frac{\partial u}{\partial x_d} \right \}
\end{aligned}$

    """    


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

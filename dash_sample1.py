# -*- coding: utf-8 -*-
import dash
import sys
import numpy as np
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
#sys.path.append('.')
#sys.path.append('..')
import util
from server import server
import random
import math
import copy


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
external_scripts = ['https://yyrcd-1256568788.cos.na-siliconvalley.myqcloud.com/yyrcd/2020-03-21-iframeResizer.contentWindow.min.js']
app = dash.Dash(name='test1',
                external_stylesheets=external_stylesheets,
                external_scripts=external_scripts)
     #           (server=server,
      #          routes_pathname_prefix='/test1/')

Markdown_text = r"""
## Describing the Free Particle
We start by describing a free particle: a particle that is not under the influence of a potential.   
 As any other particle, the state of a ***Free Particle*** is described with a ket $\left|\psi(x)\right>$. In order to learn about this particle (measure its properties) we must construct *Hermitian Operators*. For example, what is the **momentum operator $\hat P$**?
 
The momentum operator most obey the following property (eigenfunction/eigenvalue equation):
$$\hat P \left| \psi_k(x) \right> =p\left | \psi_k(x)\right>  \tag{1}$$ 
where *p* is an eigenvalue of a *Hermitian operator* and therefore it is a real number.
In the $x$ representation, using the momentum operator as $\hat P =-i\hbar \frac{\partial }{\partial x}$, we can solve equation 1 by proposing a function to represent $\left| \psi_k(x) \right>$ as $\psi_k(x) = c\ e^{ikx}$, where $k$ is a real number.
Let's see if it works:  
$$\hat P \psi_k(x) =p \psi_k(x)$$ 
$$-i\hbar \frac{\partial {c\ e^{ikx}}}{\partial x} =-i\hbar\ c\ ik\ e^{ikx} $$ 
$$\hbar k\ c\ e^{ikx} = \hbar k\ \psi_k(x) \tag{2}$$
with $p=\hbar k$
-------------------
"""
steps=10000
sigma   = 1
epsilon = 2
L=3*sigma
delta=0.5
temp=300
kb=1.38e-23*6.022e23/(1000*4.184)
#start = time.time()
def cutoff(r):
    if r > L:
        r -= L
    elif r < 0:
        r += L
    return r

def energy(r,epsilon):
    energy = 4*epsilon*(((sigma/r)**12)-((sigma/r)**6))
    return energy

def move(diff_e,e,r,temp):
    if diff_e < 0:
        
        e += diff_e
    else:
        rand = random.random()
        if math.exp(-diff_e/(kb*temp)) > rand:
            
            e += diff_e
        else:
            r = origi
            e = pre_e
    return r,e
#initialize system
r = random.uniform(0,L)
e = energy(r,epsilon)
r_data=[r]
r2_data=[r**2]
e_data=[e]
r_ave=[r]
r2_ave=[r**2]
#run mc
acc_step = 0
for step in range(0,steps):
    origi = copy.deepcopy(r)
    pre_e = e
    r += random.uniform(-1,1)*delta
    r = cutoff(r)
    new_e = energy(r,epsilon)
    diff_e = new_e - pre_e
    r,e = move(diff_e,e,r,temp)
    if r != origi:
        acc_step += 1
    r_data.append(r)
    r2_data.append(r**2)
    e_data.append(e)
    r_ave.append(sum(r_data)/len(r_data))
    r2_ave.append(sum(r2_data)/len(r2_data))
#acc_ratio = acc_step/steps
#print(r_data)
Markdown_text = util.convert_latex(Markdown_text)

app.layout = html.Div([
    dcc.Markdown(Markdown_text, dangerously_allow_html=True),
    dcc.Graph(id='graph-with-slider'),
    dcc.Slider(
        id='year-slider',
        min=0,
        max=10000,
        value=3,
        marks={str(x): str(x) for x in np.arange(0, 10001, 1000)},
        step=1000
    )
  #  dcc.Interval(
  #      id='interval-component',
  #      interval=2 * 1000,  # in milliseconds
  #      n_intervals=0
  #  )
])


@app.callback(
    Output('graph-with-slider', 'figure'),
    [Input('year-slider', 'value')])
def update_figure(a):
    x_data = np.arange(0,10001,1)
    print(a)
    y = 1
    x = x_data[a]
    print(y)
    
#    N = 200
#    x = np.linspace(0, 12, N)
#    k = 1
#   # print(sec)
#    w = 1 - 0.3 * sec
#    b = 0
#    y = a * np.sin(k * x + w) + b

    return {
        'data': [dict(
            x=x,
            y=y,
            mode='markers'
  #          opacity=0.7
  #          color='firebrick',
  #          marker={
  #              'size': 100,
  #              'line': {'width': 0.5, 'color': 'red'}
  #         }
        )],
        'layout': dict(
            xaxis={'range': [0, 10000]},
            yaxis={'range': [-2.5,2.5]}
            # margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            # legend={'x': 0, 'y': 1},
            # hovermode='closest',
            # transition={'duration': 500},
        )
    }


if __name__ == '__main__':
    app.run_server(debug=True)

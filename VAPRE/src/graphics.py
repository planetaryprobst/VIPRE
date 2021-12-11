#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#--------------------------------------------------------------------------
#   graphics.py
#--------------------------------------------------------------------------
#
#   Collection of functions that create graphics and visualizations 
#   used by VAPRE.py.
#
#--------------------------------------------------------------------------
#
#   Use:
#   import src.graphics as graphics 
#
#*************************************************************************#
# Language: Python 3 (OSX) using Matlab 2019b
# Author: Alena Probst
# History:
# Version |    Date    |     Name      | Change history
# v1.0    | 06.12.2021 |  A. Probst    | First release
#*************************************************************************#
"""

# GRAPHICS
###################
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from dash_table.Format import Format

# MATH
###################
import math


def generate_table(pointer,data):
    
    dct,vInf,launchDate,ToF,vHyp,arrMass,deltaV,vIDs,tIDs = data
    
    tID = tIDs[pointer]    
    
    params = [
        'Launch date, yr past 2000',
        'Time of flight, yr',
        'Arrival date, yr past 2000',
        f'v<sub>\u221E</sub> @ arrival, EMO2000',
        f'C₃, km²/s²',
        f'\u0394v, km/s',
        'S/C mass @ arrival, kg',
        f'v<sub>\u221E</sub> @ arrival, body-intertial RF',
        'Fly-by sequence',
        f'v<sub>\u221E</sub> ID',
        'trajectory ID',
        ]
    
   # flyby sequence
    lst = dct[tID][14][::2]
    flybys = ''
    
    for item in lst:
        flybys = flybys + str(item) + ' - '    
    
    # v inf
    vHypBodyInert = math.sqrt( (dct[tID][15][0])**2 + (dct[tID][15][1])**2 + (dct[tID][15][2])**2 )
    
    vals = [
        launchDate[pointer]/365,
        ToF[pointer]/365,
        (launchDate[pointer] + ToF[pointer])/365,
        vHyp[pointer],
        dct[tID][11], #C3
        deltaV[pointer],
        arrMass[pointer],
        vHypBodyInert,
        flybys[:-3],
        dct[tID][16],
        tID,
        ]
    
    columns = [
        {"name": "Parameter", "id": "parameter"},
        {"name": "Value", "id": "value",'type': 'numeric',
         'format': Format(
             precision=4, 
             )
         },
        ]
    
    # params = dct['header']
    # vals = dct[tID]
    
    data = [
        {"parameter": param,
         "value": val,
            } for param,val in zip(params,vals)
        ]
    
    return data,columns

def generate_plot_trajOverview(body,data):
    
    dct,vInf,launchDate,ToF,vHyp,arrMass,deltaV,vIDs,tIDs = data

    # strings
    str_vHyp = [f'v<sub>\u221E</sub> = {round(item,2)} km/s' for item in vHyp]
    str_arrMass = [f'm<sub>S/C</sub> = {math.ceil(item)} kg' for item in arrMass]
    str_deltaV = [f'\u0394v = {round(item,2)} km/s' for item in deltaV]
    
    # units
    launchDate = [line/365 for line in launchDate] # in yr after 2000
    ToF = [line/365 for line in ToF] # ToF in yr
    
    fig = make_subplots(rows=1, cols=3, shared_yaxes=True, #shared_xaxes=True,
                          subplot_titles=('v<sub>\u221E</sub> @ arrival, km/s', 'S/C mass @ arrival, kg', '\u0394v, km/s')
                          )
     
    fig.add_trace(
        go.Scatter(x=launchDate, y=ToF, mode = 'markers', 
                    marker = dict(size = 2, color = vHyp, colorscale = 'jet', 
                                  colorbar = dict(x = 0.33, xanchor = 'center'),
                                    showscale = True,),
                    
                    hoverinfo = 'x+y+text',
                    #  hovertext = f'Launch Date: {x}', 
                    hovertext = str_vHyp,
                    # textinfo = 'label',
                    hoverlabel=dict(namelength=0),
                    hovertemplate = 'Launch: 20%{x:.0f} <br>ToF: %{y:.2f)} yr <br> %{hovertext}'
                    ), 
                    row=1, col=1
    )
    
    fig.add_trace(
        go.Scatter(x=launchDate, y=ToF, mode = 'markers',
                    marker = dict(size = 2, color = arrMass, colorscale = 'jet', 
                                  colorbar = dict(x = 0.675, xanchor = 'center'),
                                    showscale = True#, reversescale=True
                                    ),
                    hoverinfo = 'x+y+text',
                    hovertext = str_arrMass,
                    hoverlabel = dict(namelength=0),
                    hovertemplate = 'Launch: 20%{x:.0f} <br>ToF: %{y:.2f)} yr <br> %{hovertext}'
                    ),
        row=1, col=2
    )
    
    fig.add_trace(
        go.Scatter(x=launchDate, y=ToF, mode = 'markers',
                    marker = dict(size = 2, color = deltaV, colorscale = 'jet', 
                                  colorbar = dict(x = 1.00, xanchor = 'left'),
                                    showscale = True),
                    hoverinfo = 'x+y+text',
                    hovertext = str_deltaV,
                    hoverlabel = dict(namelength=0),
                    hovertemplate = 'Launch: 20%{x:.0f} <br>ToF: %{y:.2f)} yr  <br> %{hovertext}'
                    ),
        row=1, col=3
    )
    
    fig.update_layout(
        showlegend = False, 
        title_text=f"Trajectory Data Overview of {body}, launch window 20{math.floor(min(launchDate))} to 20{math.ceil(max(launchDate))}")
    fig.update_layout(plot_bgcolor='white', paper_bgcolor='rgba(0,0,0,0)', font_color='white')
    fig.update_layout(yaxis = dict(range=[math.floor(min(ToF)), math.ceil(max(ToF)) + 0.5]))
    fig.update_xaxes(title_text='Launch Date, yr past 2000', showgrid = True, row=1, col=1,gridcolor='black')
    fig.update_xaxes(title_text='Launch Date, yr past 2000', showgrid = True, row=1, col=2,gridcolor='black')
    fig.update_xaxes(title_text='Launch Date, yr past 2000', showgrid = True, row=1, col=3,gridcolor='black')
    fig.update_yaxes(title_text='Time of Flight, yr', showgrid = True, row=1, col=1,gridcolor='black')
    fig.update_yaxes(showgrid = True, row=1, col=2,gridcolor='black')
    fig.update_yaxes(showgrid = True, row=1, col=3,gridcolor='black')
    
    return fig

def generate_plot_entryConditions(body,arrivalInfo,entryData,hEntry):
    
    
    # unpacking input data
    v_inf,t_arr = arrivalInfo
    x,y,z,FPA,vRel_entry,vRot,safe_rEntry,safe_vEntry = entryData
    
    
    # definition of scatter marker color
    FPA_color = []
    vRel_entry_color = []
    vRot_color = []
    safe_vEntry_color = []
    
    for i in range(len(FPA)):
        
        if FPA[i] is None or math.isnan(FPA[i]):
            FPA_color.append('#3f3f3f')
            vRel_entry_color.append('#3f3f3f')
            vRot_color.append('#3f3f3f')
            safe_vEntry_color.append('#3f3f3f')
            
            # set None to NaN
            FPA[i] = float('NaN')
            vRel_entry[i] = float('NaN')
            vRot[i] = float('NaN')
            safe_rEntry[i] = float('NaN')
            safe_vEntry[i] = float('NaN')
            
        else:
            FPA_color.append(FPA[i])
            vRel_entry_color.append(vRel_entry[i])
            vRot_color.append(vRot[i])
            safe_vEntry_color.append(safe_vEntry[i])
    
    # Initialize figure with 4 3D subplots
    fig = make_subplots(
            rows=2, cols=2,
            subplot_titles = [
                'Flight Path Angle, deg', 'rel. Entry Velocity, km/s',
                'Rotational Vel., m/s', 'Orb. Vel. @ entry, km/s' 
            ],
            specs=[
                [
                    {'type': 'scatter3d'}, {'type': 'scatter3d'}# , {'type': 'scatter3d'}
                ],
                [
                    {'type': 'scatter3d'}, {'type': 'scatter3d'}# , {'type': 'scatter3d'}
                ]
            ],
            vertical_spacing=0.05, horizontal_spacing=0.10,
            column_widths = [600, 600], row_heights = [500, 500],
        )
    
    # adding surfaces to subplots.
    fig.add_trace(
        go.Scatter3d(
            x=x, y=y, z=z, mode = 'markers', #meta = dict(label = {'x': 'X, km', 'y': 'Y, km', 'z':'Z, km'}),
            marker = dict(
                size = 4, color = FPA_color, colorscale = 'jet_r', 
                cmin = -90, cmax = 0, #reversescale=True, 
                colorbar = dict(
                    x = 0.45, xanchor = 'left', 
                    y = 0.75, yanchor = 'middle', len = 0.45,
                )
            ),
            hoverinfo = 'x+y+z+text',
            hovertext = FPA,
            hoverlabel = dict(namelength=0),
            hovertemplate = 'X: %{x:.1f} km <br> Y: %{y:.1f} km <br> Z: %{z:.1f} <br> FPA: %{hovertext:.1f} deg',
        ),
        row=1, col=1,
    )
    
    fig.add_trace( 
        go.Scatter3d(
            x=x, y=y, z=z, mode = 'markers', 
            marker = dict(
                size = 4, color = vRel_entry_color, colorscale = 'jet',  
                cmin = math.floor(min(vRel_entry)*10)/10, cmax = math.ceil(max(vRel_entry)*10)/10,
                colorbar = dict(
                    x = 1.0, xanchor = 'left', y = 0.75, yanchor = 'middle', len = 0.45
                )
            ),
            hoverinfo = 'x+y+z+text',
            hovertext = vRel_entry,
            hoverlabel = dict(namelength=0),
            hovertemplate = 'X: %{x:.1f} km <br> Y: %{y:.1f} km <br> Z: %{z:.1f} <br> v<sub>rel,entry</sub>: %{hovertext:.1f} km/s',
        ),
        row=1, col=2,
    )
    
    fig.add_trace(
        go.Scatter3d(
            x=x, y=y, z=z, mode = 'markers', 
            marker = dict(
                size = 4,color = vRot_color, colorscale = 'jet', 
                cmin = math.floor(min(vRot)), cmax = math.ceil(max(vRot)),
                colorbar = dict(
                    x = 0.45, xanchor = 'left', y = 0.25,
                    yanchor = 'middle', len = 0.45
                )
            ),
            hoverinfo = 'x+y+z+text',
            hovertext = vRot,
            hoverlabel = dict(namelength=0),
            hovertemplate = 'X: %{x:.1f} km <br> Y: %{y:.1f} km <br> Z: %{z:.1f} <br> v<sub>rot</sub>: %{hovertext:.1f} m/s',
        ),
        row=2, col=1,
    )
    
    fig.add_trace(
        go.Scatter3d(
            x=x, y=y, z=z, mode = 'markers', 
            marker = dict(
                size = 4, color = safe_vEntry_color, colorscale = 'jet', 
                cmin = math.floor(min(safe_vEntry)*10)/10, cmax = math.ceil(max(safe_vEntry)*10)/10,
                colorbar = dict(
                    x = 1.0, xanchor = 'left', y = 0.25, 
                    yanchor = 'middle', len = 0.45
                )
            ),
            hoverinfo = 'x+y+z+text',
            hovertext = safe_vEntry,
            hoverlabel = dict(namelength=0),
            hovertemplate = 'X: %{x:.1f} km <br> Y: %{y:.1f} km <br> Z: %{z:.1f} <br> v<sub>entry</sub>: %{hovertext:.1f} km/s',
        ),
        row=2, col=2,
    )
    
    fig.update_layout(
        showlegend = False,
        title_text=f'Arrival @ {body} with v<sub>\u221E</sub> = {v_inf:.2f} km/s on day {t_arr} JD'
    )
    # 
    # fig.update_layout(
    #     plot_bgcolor=colors['background'], paper_bgcolor=colors['background'], 
    #     font_color=colors['text']
    # )
    # fig.update_layout(
    #     height=1000, width = 1200, 
    #     margin = dict(
    #         l = 20, r = 10, t = 100, b = 40
    #     )
    # )
    fig.update_layout(
        xaxis = dict(
            title = 'X, km'
        )
    )    

    return fig

def generate_plot_empty():
    
    # Initialize figure with 4 3D subplots
    fig = make_subplots(
            rows=2, cols=2,
            # specs=[
            #     [
            #         {'type': 'scatter3d'}, {'type': 'scatter3d'}# , {'type': 'scatter3d'}
            #     ],
            #     [
            #         {'type': 'scatter3d'}, {'type': 'scatter3d'}# , {'type': 'scatter3d'}
            #     ]
            # ],
            # vertical_spacing=0.05, horizontal_spacing=0.10,
            # column_widths = [600, 600], row_heights = [500, 500],
        )
    
    fig.update_layout(
        showlegend = False,
        title_text=f'Choose an entry altitude ... '
    )
    # 
    # fig.update_layout(
    #     plot_bgcolor=colors['background'], paper_bgcolor=colors['background'], 
    #     font_color=colors['text']
    # )
    fig.update_layout(
        height=1000, width = 1200, 
        margin = dict(
            l = 20, r = 10, t = 100, b = 40
        )
    )    

    return fig

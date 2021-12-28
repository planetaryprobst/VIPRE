#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#--------------------------------------------------------------------------
#
#   VAPRE - Visualization of Atmospheric PRobe Entry conditions for  
#   different bodies and trajectories
#
#--------------------------------------------------------------------------
#
#   A software tool to visualize the entry conditions for all safe entry 
#   opportunities of a planetary entry probe being released from an 
#   interplanetary trajectory. It is a prototype and proof of concept on    
#   how to display all trajectory options and make the computed entry data  
#   available for filtering and browsing.
#   
#   The data files used as input are the output files from the software  
#   tool IPED.
#
#   VAPRE and IPED are part of the software package VIPRE, a software to 
#   Visualize the Impact of the PRobe Entry location on the science, mission 
#   and spacecraft design. 
#
#--------------------------------------------------------------------------
#   Call:
#   python3 VAPRE.py 
#--------------------------------------------------------------------------
#
#   -----
#   Input
#   -----
#   Output .txt files of IPED, locted in src/data/'body';'/'hEntry' folder.
#
#   ------
#   Output
#   ------
#   browser based app, running on http://127.0.0.1:8050/
#  
#*************************************************************************#
# Language: Python 3 (OSX) using Matlab 2019b
# Author: Alena Probst
# History:
# Version |    Date    |     Name      | Change history
# v1.0    | 06.12.2021 |  A. Probst    | First release
#*************************************************************************#
"""

# DASH
###################

import dashVAPRE
import dash_html_components as html
from dash.dependencies import Input,Output,State
from dash.exceptions import PreventUpdate

import plotly.graph_objects as go
import json

# APP CONTENT
###################

import src.content as var

# APP STRUCTURE
###################

import src.view as view

# APP DATA HANDLING
###################

import src.dataLoading as dL
import src.dataProcessing as dP

# APP GRAPHICS HANDLING
#######################

import src.graphics as graphics


# APP PROPERTIES
###################

app = dash.Dash(__name__)




# APP VISUALIZATION
###################

app.layout = html.Div(
    id = 'vapre-framework',
    className = 'vapre-framework',
    children = [
        
        view.build_header(var.heading,var.subheading),
        
        view.build_body(var.dropdown_body,var.dropdown_hEntry,var.idG1,\
                        var.idG2,var.idG3,var.idT1,var.idT2,var.tag),
 
        
        ]    
    )


# callbacks
@app.callback(
    [Output('store-status','data'),
     Output('store-status-old','data'),
     ],
    [Input('button-update-input', 'n_clicks'),
     ],
    [State('dropdown-menu-body','value'),
     State('dropdown-menu-hEntry','value'),
     State('store-status','data'),
     ]
    )
def update_status(n_clicks,body,hEntry,status):
    
    if n_clicks is None:
        status = {'body': None,'hEntry': None, 'vID': None }
        status_old = {'body': None,'hEntry': None, 'vID': None }
    else:
        status_old = status
        status = {'body': body,'hEntry': hEntry, 'vID': None}
        
    if status == status_old:
        raise PreventUpdate
    else:
        return [status,status_old]    
 
    
@app.callback(
    [Output('store-data','data'),
      Output('store-data-all','data'),
      ],
    [Input('store-status', 'data'),
      ],
    [State('store-status-old','data'),
      State('store-data-all','data'),
      ],
    )
def update_data(status,status_old,data_all):
    if status['body'] is None or status['body'] is status_old['body']:
        raise PreventUpdate
    else:
        body = status['body']
        
        if data_all is None:
            data_all = {body: None}
        elif body not in data_all.keys():
            data_all[body] = None
                    
        if data_all[body] is not None:
            return [data_all[body],data_all]
        else:
            data_all[body] = dL.load_trajectoryData(body)

            return [data_all[body],data_all]


@app.callback(
    [Output('graph-traj-data-overview','figure'),
      Output('graph-traj-data-overview','style'),
      ],
    [Input('store-status', 'data'),
      ],
    [State('store-status-old','data'),
      ],
)
def update_graphTrajOverview(status,status_old):
    
    if status is None or status['body'] is None:
        return [go.Figure(data=[go.Scatter(x=[], y=[])]),{'display': 'none'}]
    elif status['body'] is status_old['body']: 
        raise PreventUpdate
    else:
        data = dL.load_trajectoryData(status['body'])
        
        return [graphics.generate_plot_trajOverview(status['body'],data),{'display': 'block'}]


@app.callback(
    [Output('table-trajectory-details-1','data'),
     Output('table-trajectory-details-1','columns'),
      ],
    [Input('graph-traj-data-overview','clickData')
      ],
    [State('store-status', 'data'),
     State('store-data','data'),
     State('store-clicks','data'),
      ],   
    )
def update_trajDetails1(clickData,status,data,clicks):
    if clickData is None:
        raise PreventUpdate

    if clicks is None:
        clicks = 1
    else:
        clicks = clicks + 1
        
    if (clicks % 2) == 0:
        raise PreventUpdate
    else:
        
        pointer = clickData['points'][0]['pointIndex']
    
        return graphics.generate_table(pointer,data)


@app.callback(
    [Output('table-trajectory-details-2','data'),
     Output('table-trajectory-details-2','columns'),
      ],
    [Input('graph-traj-data-overview','clickData')
      ],
    [State('store-status', 'data'),
     State('store-data','data'),
     State('store-clicks','data'),
      ],   
    )
def update_trajDetails2(clickData,status,data,clicks):
    if clickData is None:
        raise PreventUpdate

    if clicks is None:
        raise PreventUpdate
    else:
        clicks = clicks + 1
        
    if (clicks % 2) == 0:    
        pointer = clickData['points'][0]['pointIndex']
    
        return graphics.generate_table(pointer,data)
    
    else:
        raise PreventUpdate

@app.callback(
    [Output('graph-entry-data','figure'),
      Output('graph-entry-data','style'),
      ],
    [Input('graph-traj-data-overview','clickData'),
     Input('store-status','data'),
      ],
    [State('store-status-old', 'data'),
      State('store-data','data'),
      ],        
    )
def select_trajectory(clickData,status,status_old,data):
    
    if clickData is None or status['hEntry'] is None or status['body'] is None:
        return [go.Figure(data=[go.Scatter(x=[], y=[])]),{'display': 'none'}]
    elif status is status_old:
        raise PreventUpdate
    else:
        
        dct,vInf,launchDate,ToF,vHyp,arrMass,deltaV,vIDs,tIDs = data
        
        pointer = clickData['points'][0]['pointIndex']
        
        vID = vIDs[pointer]
        body = status['body']
        hEntry = status['hEntry']
        
        arrivalInfo = dP.generate_arrivalInfo(vID,vInf)
        entryFile = dL.load_entryFile(body,vID,hEntry)
        entryData = dP.generate_entryData(entryFile)
        
        return [graphics.generate_plot_entryConditions(body,arrivalInfo,entryData,hEntry),{'display': 'inline-block'}]

    
# SERVER 
###################

if __name__ == '__main__':
    app.run_server(debug=True, use_reloader = True)

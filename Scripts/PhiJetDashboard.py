import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from Utils.DashUtils import *
import plotly.graph_objects as go
import numpy as np
from glob import glob
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import plotly.graph_objects as go
import numpy as np
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

out_files = glob("C:/Users/Edwin/Dropbox (GaTech)/DOE/Edwin/May 2019/03-25-2019_PhiJet/OutFiles/*.pickle")
df = pd.concat([pd.read_parquet(file) for file in out_files])
df = df[(df['NO_ppmvd_constraint'] >= 0) & 
        (df['T_constraint'] >= 1900) & 
        (df['CO_ppmvd_constraint'] <= 33)]
# df.to_clipboard()
# available_indicators = df['Indicator Name'].unique()
columns = df.columns

app.layout = html.Div([
    # Crossfiltering and selection div
    html.Div([
        html.Div([ # X-axis selection
            dcc.Dropdown(
                id='crossfilter-xaxis-column',
                options=[{'label': col, 'value': col} for col in columns],
                value='tau_ent_sec (s)'
            ),
            dcc.RadioItems(
                id='crossfilter-xaxis-type',
                options=[{'label': i, 'value': i} for i in ['Linear', 'Log']],
                value='Linear',
                labelStyle={'display': 'inline-block'}
            )
        ],
        id='xaxis-crossfilter-div',
        style={'width': '30%', "margin-left": "auto", "margin-right": "auto", 'display': 'inline-block'}, 
        className="three columns"),

        html.Div([ # Y-axis selection
            dcc.Dropdown(
                id='crossfilter-yaxis-column',
                options=[{'label': i, 'value': i} for i in columns],
                value='NO_ppmvd_constraint'
            ),
            dcc.RadioItems(
                id='crossfilter-yaxis-type',
                options=[{'label': i, 'value': i} for i in ['Linear', 'Log']],
                value='Linear',
                labelStyle={'display': 'inline-block'}
            )
        ], 
        style={'width': '30%', "margin-left": "auto", "margin-right": "auto", 'display': 'inline-block'}, 
        className="three columns"),

        html.Div( # Color dropdown
            dcc.Dropdown(
                id='color-column-dropdown',
                options=[{'label': i, 'value': i} for i in columns] + [{'label': 'None', 'value': 'None'}],
                value='phi_jet_norm'
            ), 
            style={'width': '30%', "margin-left": "auto", "margin-right": "auto", 'display': 'inline-block'},
            className="three columns"       
        )
    ], style={
        'borderBottom': 'thin lightgrey solid',
        'backgroundColor': 'rgb(250, 250, 250)',
        'padding': '10px 5px'
    }),

    # Scatterplot Graph div
    html.Div([
        dcc.Graph(
            id='crossfilter-indicator-scatter',
            hoverData={'points': [{'customdata': 'Japan'}]}
        )
    ], style={'width': '49%', 'display': 'inline-block', 'padding': '0 20'}),
    
    # Timeseries graphs div
    html.Div([
        dcc.Graph(id='x-time-series'),
        dcc.Graph(id='y-time-series'),
    ], style={'display': 'inline-block', 'width': '49%'}),

    # html.Div(dcc.Slider(
    #     id='crossfilter-year--slider',
    #     min=df['Year'].min(),
    #     max=df['Year'].max(),
    #     value=df['Year'].max(),
    #     marks={str(year): str(year) for year in df['Year'].unique()}
    # ), style={'width': '49%', 'padding': '0px 20px 20px 20px'})
])

# Callback function to update the scatterplot's 'figure' attribute based on the chosen x and y filters.
@app.callback(
    dash.dependencies.Output('crossfilter-indicator-scatter', 'figure'),
    [dash.dependencies.Input('crossfilter-xaxis-column', 'value'),
     dash.dependencies.Input('crossfilter-yaxis-column', 'value'),
     dash.dependencies.Input('color-column-dropdown', 'value'),
     dash.dependencies.Input('crossfilter-xaxis-type', 'value'),
     dash.dependencies.Input('crossfilter-yaxis-type', 'value'),])
def update_scatter(xaxis_column_name, yaxis_column_name, color_column_name,
                    xaxis_type, yaxis_type):
    dff = df
    # if (color_column_name != 'None') & (): # can group by color to form multiple traces
        

    return {
        'data': [
            go.Scattergl(
                x=dff[xaxis_column_name],
                y=dff[yaxis_column_name],
                # z=dff['phi_jet_norm'],
                text=[f"T_max = {T_max}" for T_max in dff['T_max']],
                customdata=dff[yaxis_column_name],
                mode='markers',
                marker={
                    'size': 15,
                    'opacity': 0.75,
                    'line': {'width': 0.5, 'color': 'white'},
                    'color': dff[color_column_name].values,
                    'colorscale': 'Viridis', 
                    'colorbar': {'title': color_column_name}
                }
            )
        ],
        'layout': go.Layout(
            xaxis={
                'title': xaxis_column_name,
                'type': 'linear' if xaxis_type == 'Linear' else 'log'
            },
            yaxis={
                'title': yaxis_column_name,
                'type': 'linear' if yaxis_type == 'Linear' else 'log'
            },
            margin={'l': 40, 'b': 30, 't': 10, 'r': 0},
            height=450,
            hovermode='closest'
        )
    }


#helper function used by the update_timeseries callbacks
def create_time_series(dff, axis_type, title):
    return {
        'data': [go.Scatter(
            x=dff['Year'],
            y=dff['Value'],
            mode='lines+markers'
        )],
        'layout': {
            'height': 225,
            'margin': {'l': 20, 'b': 30, 'r': 10, 't': 10},
            'annotations': [{
                'x': 0, 'y': 0.85, 'xanchor': 'left', 'yanchor': 'bottom',
                'xref': 'paper', 'yref': 'paper', 'showarrow': False,
                'align': 'left', 'bgcolor': 'rgba(255, 255, 255, 0.5)',
                'text': title
            }],
            'yaxis': {'type': 'linear' if axis_type == 'Linear' else 'log'},
            'xaxis': {'showgrid': False}
        }
    }

# # 
@app.callback(
    dash.dependencies.Output('x-time-series', 'figure'),
    [dash.dependencies.Input('crossfilter-indicator-scatter', 'hoverData'),
     dash.dependencies.Input('crossfilter-xaxis-column', 'value'),
     dash.dependencies.Input('crossfilter-xaxis-type', 'value')])
def update_y_timeseries(hoverData, xaxis_column_name, axis_type):
    country_name = hoverData['points'][0]['customdata']
    dff = df[df['Country Name'] == country_name]
    dff = dff[dff['Indicator Name'] == xaxis_column_name]
    title = '<b>{}</b><br>{}'.format(country_name, xaxis_column_name)
    return create_time_series(dff, axis_type, title)

# # 
@app.callback(
    dash.dependencies.Output('y-time-series', 'figure'),
    [dash.dependencies.Input('crossfilter-indicator-scatter', 'hoverData'),
     dash.dependencies.Input('crossfilter-yaxis-column', 'value'),
     dash.dependencies.Input('crossfilter-yaxis-type', 'value')])
def update_x_timeseries(hoverData, yaxis_column_name, axis_type):
    dff = df[df['Country Name'] == hoverData['points'][0]['customdata']]
    dff = dff[dff['Indicator Name'] == yaxis_column_name]
    return create_time_series(dff, axis_type, yaxis_column_name)

# 
if __name__ == '__main__':
    app.run_server(debug=True)
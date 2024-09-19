
from dash import Dash, html, dcc, Output, Input
import dash, dash_table
from dash.exceptions import PreventUpdate
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import sqlite3 as sql
from sqlite3 import Error
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import seaborn as sns


app = Dash(__name__)

# Load static data files
gene_pdbs = pd.read_csv("gene_pdbs")
pdb_residual = pd.read_csv("pdb_residual")
gene_names = pd.read_csv("gene_names_list.csv")  # This file contains gene-to-ddg_info filename mapping

# Layout
layout = html.Div(children=[
    html.Br(),
    html.H1(children='Folding Energies'),

    html.Div(children='''
        Use the dropdowns below to select the gene and describe a variant.
    '''),

    # Flex container for histograms
    html.Div([
        html.Div([
            html.Div([
                "Gene: ",
                dcc.Dropdown(
                    options=[{'label': gene, 'value': gene} for gene in gene_pdbs['name_of_gene'].unique() if gene != 'gene_name_value'],
                    id="gene_selected",
                    searchable=True,
                    placeholder="Select a gene...",
                    clearable=True
                ),
            ], style={'padding': 5}),
            dcc.Store(id='ddg_info_store'),  # Store for dynamically loaded ddg_info
            dcc.Store(id='pdb_values'),
            dcc.Graph(id="gene_ddg"),
        ], style={'width': '49%', 'display': 'inline-block', 'padding': 10}),

        html.Div([
            html.Div([
                html.Div([
                    "Residue number: ",
                    dcc.Dropdown(id="residual_selected", searchable=True, placeholder="Select residue number...", clearable=True),
                ], style={'width': '32%', 'display': 'inline-block', 'padding': 5}),
                html.Div([
                    "Mutation from: ",
                    dcc.Dropdown(id="mutfrom_selected", searchable=True, placeholder="Select mutation from...", clearable=True),
                ], style={'width': '32%', 'display': 'inline-block', 'padding': 5}),
                html.Div([
                    "Mutation to: ",
                    dcc.Dropdown(id="mutto_selected", searchable=True, placeholder="Select mutation to...", clearable=True),
                ], style={'width': '32%', 'display': 'inline-block', 'padding': 5}),
            ], style={'display': 'flex', 'justify-content': 'space-between', 'width': '100%'}),
            dcc.Store(id='filtered_ddg_info'),
            dcc.Store(id='median_ddg'),
            dcc.Graph(id="variant_ddg"),
        ], style={'width': '49%', 'display': 'inline-block', 'padding': 10}),
    ], style={'display': 'flex', 'justify-content': 'space-between'}),
    # Full width Markdown text
    html.Div([
        dcc.Markdown(
            id='gene_ddg_markdown',
            style={'width': '100%', 'white-space': 'pre-line', 'padding': '10px', 'box-sizing': 'border-box'}
        ),
    ], style={'padding': 10}),
])
            
             
# Callback to load the appropriate ddg_info file based on gene selection
@app.callback(
    Output('ddg_info_store', 'children'),
    Input('gene_selected', 'value')
)
def load_ddg_info(gene_selected):
    if gene_selected == 'VCP':
        # Special handling for VCP gene
        file_paths = ['ddg_info9a.csv', 'ddg_info9b.csv']
        dfs = [pd.read_csv(file_path) for file_path in file_paths]
        ddg_info = pd.concat(dfs, ignore_index=True)
        return ddg_info.to_dict('records')
    
    else:
        # Find the column in gene_names where the gene_selected is located
        matching_columns = gene_names.columns[gene_names.eq(gene_selected).any()]

        # Use the column name (e.g., 'ddg_info1') to determine the filename
        ddg_info_column = matching_columns[0]
        file_path = f"{ddg_info_column}.csv"  # Construct the filename based on the column name

        # Read the ddG information from the file
        ddg_info = pd.read_csv(file_path)
        
        # Convert the DataFrame to a dictionary for storage in dcc.Store
        return ddg_info.to_dict('records')


@app.callback(
    Output(component_id="residual_selected", component_property="children"),
    Input(component_id="gene_selected", component_property="value")
)

def set_dropdown_options_page1_2a(gene_selected):
    if gene_selected:
        pdb_residual_values = pdb_residual[gene_selected].dropna().astype(int).tolist()
        return [{'label': str(residual), 'value': residual} for residual in pdb_residual_values]
    return []


@app.callback(
    Output(component_id = "mutfrom_selected", component_property = "children"),
    [Input(component_id = "ddg_info_store", component_property = "data"),
     Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value")]
)

def set_dropdown_options_page1_2b(ddg_info_store, gene_selected, residual_selected):
    ddg_info = pd.DataFrame(ddg_info_store) # Convert stored data back to a DataFrame
    filtered_gene_pdbs = gene_pdbs[gene_pdbs['name_of_gene'] == gene_selected]
    pdb_values = filtered_gene_pdbs['pdb'].unique().tolist()
    filtered_ddg_info_var1 = ddg_info[(ddg_info['pdb'].isin(pdb_values)) & (ddg_info['pdb_residual'] == residual_selected)]
    return [{'label': variant, 'value': variant} for variant in filtered_ddg_info_var1['mut_from'].unique() if variant != 'mut_from_value']

@app.callback(
    Output(component_id = "mutto_selected", component_property = "children"),
    [Input(component_id = "ddg_info_store", component_property = "data"),
     Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value")]
)

def set_dropdown_options_page1_2c(ddg_info_store, gene_selected, residual_selected, mutfrom_selected):
    ddg_info = pd.DataFrame(ddg_info_store)
    filtered_gene_pdbs = gene_pdbs[gene_pdbs['name_of_gene'] == gene_selected]
    pdb_values = filtered_gene_pdbs['pdb'].unique().tolist()
    filtered_ddg_info_var2 = ddg_info[(ddg_info['pdb'].isin(pdb_values)) & (ddg_info['pdb_residual'] == residual_selected) & (ddg_info['mut_from'] == mutfrom_selected)]
    return [{'label': variant, 'value': variant} for variant in filtered_ddg_info_var2['mut_to'].unique() if variant != 'mut_to_value']

## Callback for ddg_for_gene
@app.callback(
    Output(component_id = "gene_ddg", component_property = "figure"),
    [Input(component_id = "ddg_info_store", component_property = "data"),
     Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "pdb_values", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value"),
     Input(component_id = "mutto_selected", component_property = "value"),]
)

def get_pdb_values(gene_pdbs, gene_selected):
    filtered_gene_pdbs = gene_pdbs[gene_pdbs['name_of_gene'] == gene_selected]
    pdb_values = filtered_gene_pdbs['pdb'].unique().tolist()
    return pdb_values

# Calculate median of the variant histogram

def calculate_median (ddg_info_store, pdb_values,residual_selected, mutfrom_selected, mutto_selected):
    ddg_info = pd.DataFrame(ddg_info_store)
    if mutfrom_selected and mutto_selected:
        filtered_ddg_info_var3 = ddg_info[(ddg_info['pdb'].isin(pdb_values)) & (ddg_info['pdb_residual'] == residual_selected) & (ddg_info['mut_from'] == mutfrom_selected) & (ddg_info['mut_to'] == mutto_selected)]
        median_ddg = filtered_ddg_info_var3['ddg'].median()
    else:
        median_ddg = None
    return median_ddg

def ddg_for_gene_plot(ddg_info_store, pdb_values, median_ddg):
    ddg_info = pd.DataFrame(ddg_info_store)
    filtered_ddg_info = ddg_info[ddg_info['pdb'].isin(pdb_values)]

    figure = px.histogram(filtered_ddg_info, x='ddg', range_x=[-10, 100], nbins=1000, 
                         title='Histogram of ΔΔG values for selected gene', 
                         labels={'ddg': 'ΔΔG (kcal/mol)', 'count': 'Frequency'})

    if median_ddg is not None:
        figure.add_shape(
            go.layout.Shape(
                type="line",
                x0=median_ddg,
                x1=median_ddg,
                y0=0,
                y1=1,
                xref="x",
                yref="paper",
                line=dict(
                    color="Red",
                    width=2
                )
            )
        )
        figure.add_annotation(
            go.layout.Annotation(
                x=median_ddg,
                y=1,
                xref="x", 
                yref="paper",
                text=f'Variant median: {median_ddg:.2f} kcal/mol',
                showarrow=True,
                arrowhead=2
            )
        )

    return figure

## Callback for ddg_for_variant
@app.callback(
    Output(component_id = "variant_ddg", component_property = "figure"),
    [Input(component_id = "ddg_info_store", component_property = "data"),
     Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "pdb_values", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value"),
     Input(component_id = "mutto_selected", component_property = "value"),]
)

def ddg_for_variant_plot(ddg_info_store, pdb_values, residual_selected=None, mutfrom_selected=None, mutto_selected=None):
    ddg_info = pd.DataFrame(ddg_info_store)
    # Apply the initial filter based on pdb values
    filtered_ddg_info_var3 = ddg_info[(ddg_info['pdb'].isin(pdb_values)) & (ddg_info['pdb_residual'] == residual_selected) & (ddg_info['mut_from'] == mutfrom_selected) & (ddg_info['mut_to'] == mutto_selected)]
   
        # Create the histogram
    figure = px.histogram(filtered_ddg_info_var3, x='ddg', range_x=[-10, 100], nbins=20, title='Histogram of ΔΔG values for selected variant', labels={'ddg': 'ΔΔG (kcal/mol)', 'count': 'Frequency'})
    return figure

##Callback for markdown text
@app.callback(
    Output(component_id = "gene_ddg_markdown", component_property = "children"),
    [Input(component_id = "ddg_info_store", component_property = "data"),
     Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "pdb_values", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value"),
     Input(component_id = "mutto_selected", component_property = "value"),]
)
def calculate_percentile(ddg_info_store, gene_pdbs, gene_selected, ddg_info, residual_selected, mutfrom_selected, mutto_selected):
    ddg_info = pd.DataFrame(ddg_info_store)
    filtered_gene_pdbs = gene_pdbs[gene_pdbs['name_of_gene'] == gene_selected]
    pdb_values = filtered_gene_pdbs['pdb'].unique().tolist()
    filtered_ddg_info = ddg_info[(ddg_info['pdb'].isin(pdb_values))]
    values = filtered_ddg_info['ddg'].values
    filtered_ddg_info_var3 = ddg_info[(ddg_info['pdb'].isin(pdb_values)) & (ddg_info['pdb_residual'] == residual_selected) & (ddg_info['mut_from'] == mutfrom_selected) & (ddg_info['mut_to'] == mutto_selected)]
    median_ddg = filtered_ddg_info_var3['ddg'].median()
    percentile = np.sum(values < median_ddg) / len(values) * 100
    return percentile

def gene_ddg_markdown_text(ddg_info_store, gene_pdbs, gene_selected, ddg_info, residual_selected, mutfrom_selected, mutto_selected, median_ddg):
    ddg_info = pd.DataFrame(ddg_info_store)
    percentile = calculate_percentile(ddg_info_store, gene_pdbs, gene_selected, ddg_info, residual_selected, mutfrom_selected, mutto_selected)
    
    Serrano = "[Serrano](https://www.crg.eu/luis_serrano)"
    Hall = "[Hall, Shorthouse, Alcraft et al. 2023](https://www.nature.com/articles/s42003-023-05136-y)"
    
    if median_ddg is not None:
        if median_ddg > 2.5:
            return (f'A ΔΔG value greater than the {Serrano} value of +2.5 kcal/mol is commonly used as a cut-off for significantly destabilising mutations. '
                    f'Other studies, such as {Hall}, suggest a deleterious value of +0.5 kcal/mol is a threshold for destabilising mutations. '
                    f'The median ΔΔG for the selected variant is {median_ddg:.2f} kcal/mol and in the {percentile:.0f}th percentile. '
                    f'It is greater than the Serrano value of +2.5 kcal/mol and significantly destabilising.')
        elif median_ddg > 0.5:
            return (f'A ΔΔG value greater than the {Serrano} value of +2.5 kcal/mol is commonly used as a cut-off for significantly destabilising mutations. '
                    f'Other studies, such as {Hall}, suggest a deleterious value of +0.5 kcal/mol is a threshold for destabilising mutations. '
                    f'The median ΔΔG for the selected variant is {median_ddg:.2f} kcal/mol and in the {percentile:.0f}th percentile. '
                    f'It is greater than the deleterious value of +0.5 kcal/mol and destabilising.')
        else:
            return (f'A ΔΔG value greater than the {Serrano} value of +2.5 kcal/mol is commonly used as a cut-off for significantly destabilising mutations. '
                    f'Other studies, such as {Hall}, suggest a deleterious value of +0.5 kcal/mol is a threshold for destabilising mutations. '
                    f'The median ΔΔG for the selected variant is {median_ddg:.2f} kcal/mol and in the {percentile:.0f}th percentile. '
                    f'It is not destabilising.')
    return None
        



if __name__ == '__main__':
    app.run_server(debug=True,port=8052)

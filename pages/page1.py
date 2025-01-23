from dash import Dash, html, dcc
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import sqlite3 as sql
import numpy as np

import duckdb

app = Dash(__name__)

# Load static data files
gene_pdbs = pd.read_csv("gene_pdbs")
pdb_residual = pd.read_csv("pdb_residual")
gene_names = pd.read_csv("gene_names_list.csv")  # This file contains gene-to-ddg_info filename mapping
mutfrom_options = pd.read_csv("dropdown_pdb_mut_from.csv", dtype=str)
mutto_options = pd.read_csv("dropdown_pdb_mut_from_to.csv", dtype=str)

# Load ddg info stored in DB
duckdb_con = duckdb.connect('ddg_info/ddg_info.db', read_only=True)

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
            dcc.Loading(
                id="loading-gene-ddg",
                type="default",
                children=dcc.Graph(id="gene_ddg"),
            ),
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
            dcc.Loading(
                id="loading_variant_ddg",
                type="default",
                children=dcc.Graph(id="variant_ddg"),
            ),
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


def set_dropdown_options_page1_2a(gene_selected):
    if gene_selected:
        pdb_residual_values = pdb_residual[gene_selected].dropna().astype(int).tolist()
        return [{'label': str(residual), 'value': residual} for residual in pdb_residual_values]
    return []


def set_dropdown_options_page1_2b(gene_selected,residual_selected):
    if gene_selected and residual_selected:
        column = f"{gene_selected}-{residual_selected}"
        mutfrom_values = mutfrom_options[column].dropna().tolist()
        return [{'label': str(mutfrom), 'value': mutfrom} for mutfrom in mutfrom_values]
    return []
        

def set_dropdown_options_page1_2c(gene_selected, residual_selected, mutfrom_selected):
    if gene_selected and residual_selected and mutfrom_selected:
        column = f"{gene_selected}-{residual_selected}-{mutfrom_selected}"
        mutto_values = mutto_options[column].dropna().tolist()
        return [{'label': str(mutto), 'value': mutto} for mutto in mutto_values]
    return []


def get_pdb_values(gene_pdbs, gene_selected):
    filtered_gene_pdbs = gene_pdbs[gene_pdbs['name_of_gene'] == gene_selected]
    pdb_values = filtered_gene_pdbs['pdb'].unique().tolist()
    return pdb_values

# Calculate median of the variant histogram
def calculate_median(pdb_values, residual_selected, mutfrom_selected, mutto_selected):
    if mutfrom_selected is None or mutto_selected is None:
        return None
    pdb_values_str = ', '.join(f"'{pdb}'" for pdb in pdb_values)
    query = f"""
        SELECT ddg
        FROM ddg_info
        WHERE pdb IN ({pdb_values_str})
        AND pdb_residual = '{residual_selected}'
        AND mut_from = '{mutfrom_selected}'
        AND mut_to = '{mutto_selected}'
    """
    filtered_ddg_info = duckdb_con.execute(query).fetchdf()
    median_ddg = filtered_ddg_info['ddg'].median()
    return median_ddg

def ddg_for_gene_plot(gene_selected, pdb_values, median_ddg):
    if gene_selected is None:
        return go.Figure()

    pdb_values_str = ', '.join(f"'{pdb}'" for pdb in pdb_values)
    query = f"""
        SELECT *
        FROM ddg_info
        WHERE pdb IN ({pdb_values_str})
    """
    filtered_ddg_info = duckdb_con.execute(query).fetchdf()

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

def ddg_for_variant_plot(gene_selected, pdb_values, residual_selected, mutfrom_selected, mutto_selected):
    if gene_selected is None:
        return go.Figure()

    pdb_values_str = ', '.join(f"'{pdb}'" for pdb in pdb_values)
    query = f"""
        SELECT *
        FROM ddg_info
        WHERE pdb IN ({pdb_values_str})
        AND pdb_residual = '{residual_selected}'
        AND mut_from = '{mutfrom_selected}'
        AND mut_to = '{mutto_selected}'
    """
    filtered_ddg_info_var3 = duckdb_con.execute(query).fetchdf()

    # Create the histogram
    figure = px.histogram(filtered_ddg_info_var3, x='ddg', range_x=[-10, 100], nbins=20, title='Histogram of ΔΔG values for selected variant', labels={'ddg': 'ΔΔG (kcal/mol)', 'count': 'Frequency'})
    return figure


##Callback for markdown text
def calculate_percentile(gene_selected, pdb_values, residual_selected, mutfrom_selected, mutto_selected):
    if gene_selected is None:
        return None

    pdb_values_str = ', '.join(f"'{pdb}'" for pdb in pdb_values)
    query_all = f"""
        SELECT ddg
        FROM ddg_info
        WHERE pdb IN ({pdb_values_str})
    """
    # Execute the query and fetch the results
    filtered_ddg_info_all = duckdb_con.execute(query_all).fetchdf()
    values = filtered_ddg_info_all['ddg'].values

    query_variant = f"""
        SELECT ddg
        FROM ddg_info
        WHERE pdb IN ({pdb_values_str})
        AND pdb_residual = '{residual_selected}'
        AND mut_from = '{mutfrom_selected}'
        AND mut_to = '{mutto_selected}'
    """
    filtered_ddg_info_var3 = duckdb_con.execute(query_variant).fetchdf()
    median_ddg = filtered_ddg_info_var3['ddg'].median()

    percentile = np.sum(values < median_ddg) / len(values) * 100
    return percentile

def gene_ddg_markdown_text(median_ddg, percentile):
    
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

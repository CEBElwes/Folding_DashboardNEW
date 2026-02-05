import dash_bootstrap_components as dbc
from dash import html, dcc
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np

import duckdb

# Load static data files
gene_pdbs = pd.read_csv("gene_pdbs")
pdb_residual = pd.read_csv("pdb_residual")
mutfrom_options = pd.read_csv("dropdown_pdb_mut_from.csv", dtype=str)
mutto_options = pd.read_csv("dropdown_pdb_mut_from_to.csv", dtype=str)

# Create in-memory DuckDB connection that queries CSVs directly
duckdb_con = duckdb.connect(':memory:')

# SAMPLE VERSION - Using only 2 CSV files for testing
csv_files = [
    'ddg_info/ddg_info1.csv',
    'ddg_info/ddg_info2.csv',
]

# Create a view that unions all CSV files - DuckDB will read them on-demand
union_query = " UNION ALL ".join([f"SELECT * FROM read_csv_auto('{csv}')" for csv in csv_files])
duckdb_con.execute(f"CREATE VIEW ddg_info AS {union_query}")

# Set memory limits to prevent OOM
duckdb_con.execute("PRAGMA memory_limit='256MB'")
duckdb_con.execute("PRAGMA threads=1")

# Layout
layout = dbc.Container([
    html.Br(),
    # Add banner to show this is sample data
    dbc.Alert(
        "⚠️ DEMO VERSION - This uses sample data (20% of full dataset) for testing. Full version coming soon!",
        color="warning",
        className="text-center mb-3"
    ),
    html.H1('Folding Energies', className='text-center'),
    html.Div(
        'Use the dropdowns below to select the gene and describe a variant.',
        className='text-center mb-4',
    ),

    # Dropdowns
    dbc.Row([
        dbc.Col([
            html.Div("Gene: "),
            dcc.Dropdown(
                options=[{'label': gene, 'value': gene} for gene in gene_pdbs['name_of_gene'].unique()],
                id="gene_selected",
                searchable=True,
                placeholder="Select a gene...",
                clearable=True
            ),
        ], width=3, className='mb-4'),

        dbc.Col([
            html.Div("Residual: "),
            dcc.Dropdown(
                id="residual_selected",
                searchable=True,
                placeholder="Select a residual...",
                clearable=True
            ),
        ], width=3, className='mb-4'),

        dbc.Col([
            html.Div("Mutation From: "),
            dcc.Dropdown(
                id="mutfrom_selected",
                searchable=True,
                placeholder="Select mutation from...",
                clearable=True
            ),
        ], width=3, className='mb-4'),

        dbc.Col([
            html.Div("Mutation To: "),
            dcc.Dropdown(
                id="mutto_selected",
                searchable=True,
                placeholder="Select mutation to...",
                clearable=True
            ),
        ], width=3, className='mb-4'),
    ]),

    # Graphs
    dbc.Row([
        dbc.Col([
            dcc.Loading(
                id="loading-gene-ddg",
                type="default",
                children=dcc.Graph(id="gene_ddg"),
                delay_show=200,
                delay_hide=100,
                show_initially=False,
            ),
        ], width=6, className='mb-4'),

        dbc.Col([
            dcc.Loading(
                id="loading-variant-ddg",
                type="default",
                children=dcc.Graph(id="variant_ddg"),
                delay_show=200,
                delay_hide=100,
                show_initially=False,
            ),
        ], width=6, className='mb-4'),
    ]),

    # Text
    dbc.Row([
        dbc.Col([
            dcc.Markdown(
                id='gene_ddg_markdown',
                style={
                    'width': '100%',
                    'white-space': 'pre-line',
                    'padding': '10px',
                    'box-sizing': 'border-box',
                    },
                ),
        ], width=12, className='mb-4'),
    ]),
], fluid=True)


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
    pdb_values_str = ', '.join(f"'{pdb}'" for pdb in pdb_values)
    query = f"""
        SELECT *
        FROM ddg_info
        WHERE pdb IN ({pdb_values_str})
    """
    filtered_ddg_info = duckdb_con.execute(query).fetchdf()

    figure = px.histogram(
        filtered_ddg_info,
        x='ddg',
        range_x=[-10, 100],
        nbins=1000, 
        title=f'Histogram of ΔΔG values for {gene_selected}', 
        labels={'ddg': 'ΔΔG (kcal/mol)'},
        template="plotly_white",
    )
    figure.update_yaxes(showticklabels=False, title="Frequency")

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

def ddg_for_variant_plot(pdb_values, residual_selected, mutfrom_selected, mutto_selected):
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
    figure = px.histogram(
        filtered_ddg_info_var3, x='ddg',
        range_x=[-10, 100],
        nbins=20,
        title='Histogram of ΔΔG values for selected variant',
        labels={'ddg': 'ΔΔG (kcal/mol)'},
        template="plotly_white",
    )
    figure.update_yaxes(showticklabels=False, title="Frequency")
    return figure


##Callback for markdown text
def calculate_percentile(pdb_values, residual_selected, mutfrom_selected, mutto_selected):
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

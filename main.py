

# Import necessary libraries
from dash import html, dcc
from dash.dependencies import Input, Output, State

# Connect to main app.py file
from app import app

# Connect to your app pages
from pages import page1

# Connect the navbar to the index
from components import navbar

# Define the navbar
nav = navbar.Navbar()

app.config.suppress_callback_exceptions = True

# Define the index page layout
app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    nav,
    html.Div(id='page-content', children=[]),
])

@app.callback(
    Output(component_id = "residual_selected", component_property = "options"),
    Input('gene_selected', 'value'),
    prevent_initial_call=True,
)
def update_dropdown_page1_2a(gene_selected):
    dropdownlist = page1.set_dropdown_options_page1_2a(gene_selected)
    return dropdownlist

@app.callback(
    Output(component_id = "mutfrom_selected", component_property = "options"),
    [Input(component_id = "residual_selected", component_property = "value")],
    [State(component_id = "gene_selected", component_property = "value")],
    prevent_initial_call=True,
)
def update_dropdown_page1_2b(residual_selected, gene_selected):
    dropdownlist = page1.set_dropdown_options_page1_2b(gene_selected, residual_selected)
    return dropdownlist


@app.callback(
    Output(component_id = "mutto_selected", component_property = "options"),
    [Input(component_id = "mutfrom_selected", component_property = "value")],
    [State(component_id = "gene_selected", component_property = "value"),
     State(component_id = "residual_selected", component_property = "value")],
    prevent_initial_call=True,
)
def update_dropdown_page1_2c(mutfrom_selected, gene_selected, residual_selected):
    dropdownlist = page1.set_dropdown_options_page1_2c(gene_selected, residual_selected, mutfrom_selected)
    return dropdownlist


@app.callback(
    [Output(component_id = "gene_ddg", component_property = "figure"),
     Output(component_id = "variant_ddg", component_property = "figure"),
     Output(component_id = "gene_ddg_markdown", component_property = "children")],
    [Input(component_id = "mutto_selected", component_property = "value")],
    [State(component_id = "gene_selected", component_property = "value"),
     State(component_id = "residual_selected", component_property = "value"),
     State(component_id = "mutfrom_selected", component_property = "value")],
    prevent_initial_call=True,
)
def update_graphs_and_markdown(mutto_selected, gene_selected, residual_selected, mutfrom_selected):
    gene_pdbs = page1.gene_pdbs
    pdb_values = page1.get_pdb_values(gene_pdbs, gene_selected)
    median_ddg = page1.calculate_median(pdb_values,residual_selected, mutfrom_selected, mutto_selected)
    percentile = page1.calculate_percentile(gene_selected, pdb_values, residual_selected, mutfrom_selected, mutto_selected)

    gene_figure = page1.ddg_for_gene_plot(gene_selected, pdb_values, median_ddg)
    variant_figure = page1.ddg_for_variant_plot(gene_selected, pdb_values, residual_selected, mutfrom_selected, mutto_selected)
    text = page1.gene_ddg_markdown_text(median_ddg, percentile)

    return [gene_figure, variant_figure, text]


@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/' or pathname == '/page1':
        return page1.layout
    else:  # if redirected to unknown link
        return "404 Page Error! Please choose a link"

# Run the app on localhost:8050
if __name__ == '__main__':
    app.run_server(debug=False, host='localhost', port=8050)



# Import necessary libraries
from dash import html, dcc
from dash.dependencies import Input, Output

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

## CALLBACKS FOR PAGE 1

# Callback to load the appropriate ddg_info file based on gene selection
@app.callback(
    Output('ddg_info_store', 'data'),
    Input('gene_selected', 'value')
)

def update_ddg_info(gene_selected):
    ddg_info = page1.load_ddg_info(gene_selected)
    return ddg_info

## Callbacks for dropdowns

@app.callback(
    Output(component_id = "residual_selected", component_property = "options"),
    [Input(component_id = "gene_selected", component_property = "value")]
)

def update_dropdown_page1_2a(gene_selected):
    dropdownlist = page1.set_dropdown_options_page1_2a(gene_selected)
    return dropdownlist

@app.callback(
    Output(component_id = "mutfrom_selected", component_property = "options"),
    [Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value")]
)

def update_dropdown_page1_2b(gene_selected, residual_selected):
    dropdownlist = page1.set_dropdown_options_page1_2b(gene_selected, residual_selected)
    return dropdownlist
   

@app.callback(
    Output(component_id = "mutto_selected", component_property = "options"),
    [Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value")]
)  

def update_dropdown_page1_2c(gene_selected, residual_selected, mutfrom_selected):
    dropdownlist = page1.set_dropdown_options_page1_2c(gene_selected, residual_selected, mutfrom_selected)
    return dropdownlist

@app.callback(
    Output(component_id = "gene_ddg", component_property = "figure"),
    [Input(component_id = "ddg_info_store", component_property = "data"),
     Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "pdb_values", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value"),
     Input(component_id = "mutto_selected", component_property = "value"),]
)
def update_graph1(ddg_info_store, gene_selected, pdb_values, residual_selected, mutfrom_selected, mutto_selected):
    gene_pdbs = page1.gene_pdbs  
    pdb_values = page1.get_pdb_values(gene_pdbs, gene_selected)
    median_ddg = page1.calculate_median(ddg_info_store, pdb_values,residual_selected, mutfrom_selected, mutto_selected)
    figure = page1.ddg_for_gene_plot(ddg_info_store, pdb_values, median_ddg)
    return figure

## Callback for ddg by variant
@app.callback(
    Output(component_id = "variant_ddg", component_property = "figure"),
    [Input(component_id = "ddg_info_store", component_property = "data"),
     Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "pdb_values", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value"),
     Input(component_id = "mutto_selected", component_property = "value"),]
)
def update_graph2(ddg_info_store, gene_selected, pdb_values, residual_selected, mutfrom_selected, mutto_selected):
    gene_pdbs = page1.gene_pdbs  
    pdb_values = page1.get_pdb_values(gene_pdbs, gene_selected)
    figure = page1.ddg_for_variant_plot(ddg_info_store, pdb_values, residual_selected, mutfrom_selected, mutto_selected)
    return figure

##Callback for markdown text
@app.callback(
    Output(component_id = "gene_ddg_markdown", component_property = "children"),
    [Input(component_id = "ddg_info_store", component_property = "data"),
     Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "pdb_values", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value"),
     Input(component_id = "mutto_selected", component_property = "value"),
     Input(component_id = "median_ddg", component_property = "value"),
     Input(component_id = "percentile", component_property = "value"),]
)

def update_markdown(gene_selected, ddg_info_store, pdb_values, residual_selected, mutfrom_selected, mutto_selected, median_ddg, percentile):
    gene_pdbs = page1.gene_pdbs  
    pdb_values = page1.get_pdb_values(gene_pdbs, gene_selected)
    median_ddg = page1.calculate_median(ddg_info_store, pdb_values, residual_selected, mutfrom_selected, mutto_selected)
    percentile = page1.calculate_percentile(ddg_info_store, pdb_values, residual_selected, mutfrom_selected, mutto_selected)
    text = page1.gene_ddg_markdown_text(median_ddg, percentile)
    return text




@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/' or pathname == '/page1':
        return page1.layout
    else:  # if redirected to unknown link
        return "404 Page Error! Please choose a link"

# Run the app on localhost:8050
if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0', port=80)
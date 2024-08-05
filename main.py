

# Import necessary libraries
from dash import html, dcc
from dash.dependencies import Input, Output

# Connect to main app.py file
from app import app

# Connect to your app pages
from pages import page3

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

## Callbacks for page 3

## Callbacks for dropdowns
@app.callback(
    Output(component_id="gene_selected", component_property="options"),
    [Input(component_id="genes", component_property="value")]
)

@app.callback(
    Output(component_id = "residual_selected", component_property = "options"),
    [Input(component_id = "gene_selected", component_property = "value")]
)

def update_dropdown_page3_2a(gene_selected):
    dropdownlist = page3.set_dropdown_options_page3_2a(gene_selected)
    return dropdownlist

@app.callback(
    Output(component_id = "mutfrom_selected", component_property = "options"),
    [Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value")]
)

def update_dropdown_page3_2b(gene_selected, residual_selected):
    dropdownlist = page3.set_dropdown_options_page3_2b(gene_selected, residual_selected)
    return dropdownlist
   

@app.callback(
    Output(component_id = "mutto_selected", component_property = "options"),
    [Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value")]
)  

def update_dropdown_page3_2c(gene_selected, residual_selected, mutfrom_selected):
    dropdownlist = page3.set_dropdown_options_page3_2c(gene_selected, residual_selected, mutfrom_selected)
    return dropdownlist

@app.callback(
    Output(component_id = "gene_ddg", component_property = "figure"),
    [Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "pdb_values", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value"),
     Input(component_id = "mutto_selected", component_property = "value"),]
)
def update_graph6(gene_selected, pdb_values, residual_selected, mutfrom_selected, mutto_selected):
    gene_pdbs = page3.gene_pdbs  
    ddg_info = page3.ddg_info    
    pdb_values = page3.get_pdb_values(gene_pdbs, gene_selected)
    median_ddg = page3.calculate_median(pdb_values,residual_selected, mutfrom_selected, mutto_selected)
    figure = page3.ddg_for_gene_plot(page3.ddg_info, pdb_values, median_ddg)
    return figure

## Callback for ddg by variant
@app.callback(
    Output(component_id = "variant_ddg", component_property = "figure"),
    [Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "pdb_values", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value"),
     Input(component_id = "mutto_selected", component_property = "value"),]
)
def update_graph7(gene_selected, pdb_values, residual_selected, mutfrom_selected, mutto_selected):
    gene_pdbs = page3.gene_pdbs  
    ddg_info = page3.ddg_info    
    pdb_values = page3.get_pdb_values(gene_pdbs, gene_selected)
    figure = page3.ddg_for_variant_plot(page3.ddg_info, pdb_values, residual_selected, mutfrom_selected, mutto_selected)
    return figure

##Callback for markdown text
@app.callback(
    Output(component_id = "gene_ddg_markdown", component_property = "children"),
    [Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "pdb_values", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value"),
     Input(component_id = "mutto_selected", component_property = "value"),]
)

def update_markdown(gene_selected, pdb_values, residual_selected, mutfrom_selected, mutto_selected):
    gene_pdbs = page3.gene_pdbs
    ddg_info = page3.ddg_info
    pdb_values = page3.get_pdb_values(gene_pdbs, gene_selected)
    median_ddg = page3.calculate_median(pdb_values,residual_selected, mutfrom_selected, mutto_selected)
    percentile = page3.calculate_percentile(gene_pdbs, gene_selected, ddg_info, residual_selected, mutfrom_selected, mutto_selected)
    text = page3.gene_ddg_markdown_text(gene_pdbs, gene_selected, ddg_info, residual_selected, mutfrom_selected, mutto_selected, median_ddg)
    return text





@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/page3':
        return page3.layout
    else:  # if redirected to unknown link
        return "404 Page Error! Please choose a link"

# Run the app on localhost:8050
if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0', port=80)
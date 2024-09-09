# Import necessary libraries
from dash import html
import dash_bootstrap_components as dbc

# Define the navbar structure
def Navbar():

    layout = html.Div([
        dbc.NavbarSimple(
            children=[
                dbc.NavItem(dbc.NavLink("Folding Dashboard", href="/page1")),
            ] ,
            brand="Folding DashBoard",
            brand_href="/page1",
            color="dark",
            dark=True,
        ),
    ])

    return layout

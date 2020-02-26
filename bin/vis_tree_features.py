#!/usr/bin/env python3

import toytree      
import toyplot
import numpy as np
import toyplot.pdf
import toyplot.svg
from io import StringIO

with open('tree_INPUT.newick', 'r') as file:
    newick = file.read();

# Creating tree

tree = toytree.tree(newick, tree_format=0)

# Tree Style

mystyle = {
    "edge_type": 'p',
    "edge_style": {
        "stroke": "black",
        "stroke-width": 1,
    },
    "tip_labels_align": True,
    "tip_labels_colors": "black",
    "node_labels": False,
    "edge_align_style": {
        "stroke": "grey",
        "stroke-width": 0.5,
        "stroke-dasharray": "2"
    }, 
}



canvas, axes = tree.draw(
    width=600,
    height=(tree.ntips * 8.5),
    **mystyle,
    tip_labels_style={"font-size": "8px"}
);


# stretch domain to fit long tip names
axes.x.domain.max = 10

toyplot.svg.render(canvas, "tree.svg")
toyplot.pdf.render(canvas, "tree.pdf")
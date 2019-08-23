# Minimal_IMAP_MCMC

9.8.19
Minimal IMAP MCMC is a rapid Bayesian approach to Causal Structure Learning.  It is written in Python. In addition to calling it from the command line users can also use a Graphical User Interface (GUI) to call the program.
This work was done in collaboration with Tamara Broderick's group and it is based on algorithm developed by Raj Agrawal.  To read more about the algorithm, check out the paper: <https://arxiv.org/pdf/1803.05554.pdf>

In addition to coming with the code for the GUI and the functions for the command line, it also comes pre-loaded with some example data which is run upon startup. 


## GUI Guide and Instructions

To launch the GUI call `bokeh serve --show minIMAPapp` from the folder containing the folder you downloaded.  

### Control Tab and Landing Page:

You will land on the control tab page. It looks like this:

![alt text][con_tab]

Options on the left of the page are for you open data to use in the minimal IMAP MCMC GUI.  Your first option is to open raw data from a csv file format.  Just put in the path to your file, and as long as it is properly formatted (as a csv file) it will open.Your second option is to open the results from a previous MCMC run.  This will save you time if your MCMC takes a long time to run. Just below this is an option to specify a file destination to save the current MCMC run for future use as well. 

In the column to the right of these options are the controls for the MCMC algorithm.  You can change the Conditional Independence lavel used in the Conditional Independence tests.  Number of steps to burn in are the number of graphs to ignore at the beginning of the sampling while you want for the algorithm to move the a more likely space. Number of steps controls the total number of steps, while waiting steps between samples helps gaurantee that samples are less correllated.  Not always necessary.  Gamma controls how much emphasis is put on having sparse graphs through the prior.  It is recommended you leave it at 1, but it too can be changed.  Once you are satisfied with the parameters hit the *Run From CSV* button!

Lastly, we tried to include a couple common questions or problems in the FAQ dropdown menu in the full right column.  If you run into problems, or want mroe information , checkout this dropdown, select and option, and press *Help*!

### Graph View Tab:

After running (or loading) your MCMC algorithm you can move on to viewing the results.  Using this tab you can view the graph of the composite learned graph.  This is an example of what it may look like:

![alt text][graph_plot]

The nodes of the graph represent the variables in your data (for example the concentration of a protein).  If an edge is shown then that means that the algorithm has determined that it is more likely then the threshold you have set for this plot.  If you move your mouse over a node or an edge then a small rectangle will appear near the cursor and tell you the name of the node and the probability of the edge. To set a different threshold, or to change which nodes are plotted, you can use the widgets immediately to the right of the plot, and hit *Plot*:

![alt text][full_wid]

### Edge Matrix Tab:

This plot contains all the same information shown in the graph view tab.  There are certain circumstances where this plot might be more useful though - for example if your learned plot is particularly dense and/or has many nodes. Here is an example of what you can expect this plot to look like:

![alt text][edge_plot] 

The matrix contains all possible edges which could link the node of the row to each of the nodes in the columns.  If the color is just gray then that means that the probability of that edge does not surpass the threshold.  If it is colored then it is greater then the threshold.  The plot cuts the region from the threshold to 100% into 5 piece, and color codes each edge based on what region it falls in so you can get a sense of the most likely edges overall quickly.  In addition, like the full plot tab if you move your mouse over any of the boxes it will tell you the directionality and the probability of the edge. Also like the graph view tab, you can change the threshold and which nodes are in this plot via the panel of widgets which lie to the right of the plot. 

**TODO: Can I make this a pip package?  Add a tab to evaluate convergence of the algorithm.  Add an option to run the algorithm until convergence? Update the GUI to let the suer know how much longer the algorithm will run for?**

[con_tab]: Images/control_tab.PNG "Control Tab for MCMC"
[graph_plot]: Images/graph_view_plot.png "Graph View Plot Example"
[full_wid]: Images/full_widgets.PNG "Controls for the Graph View Plot"
[edge_plot]: Images/edge_mat_plot.png "Edge Matrix View Example"
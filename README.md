# Minimal_IMAP_MCMC

9.8.19
Minimal IMAP MCMC is a rapid Bayesian approach to Causal Structure Learning.  It is written in Python. In addition to calling it from the command line users can also use a Graphical User Interface (GUI) to call the program.
This work was done in collaboration with Tamara Broderick's group and it is based on algorithm developed by Raj Agrawal.  To read more about the algorithm, check out the paper: <https://arxiv.org/pdf/1803.05554.pdf>

In order to use this package, please first download it onto your machine via pip:

'''
pip install minIMAPmcmc
'''

In addition to coming with the code for the GUI and the functions for the command line, it also comes pre-loaded with some example data which is run upon startup. 

To launch the GUI call `bokeh serve --show minIMAPapp`.  You will land on the control tab page. It looks like this:

![alt text][con_tab]

Options on the left of the page are for you open data to use in the minimal IMAP MCMC GUI.  
Your first option is to open raw data from a csv file format.  Just put in the path to your file, and as long as it is properly formatted (as a csv file) it will open.
Your second option is to open the results from a previous MCMC run.  This will save you time if your MCMC takes a long time to run.
Just below this is an option to specify a file destination to save the current MCMC run for future use as well. 

New updates to come!

[con_tab]: Images/control_tab.PNG "Control Tab for MCMC"
# randomWalkOnGraph

Create undirected graph, and set many random walks at different nodes. 

There are 2 ways to create such graph. One way is to randomly add an edge between 2 nodes with some probability. In the second way, we add an edge between a new node and existing nodes such that old nodes with high degrees are more likely to connect to this new node. This second way kind of mimics how the webpages on the internet interact. Then, once the graph is made, the code creates many random walks (starting at different nodes). 

There are 3 plots. Plot 1 (and 2) displays the mean (and variance) distance of all the walks from their starting points with respect to time. Plot 3 displays the degree of the final destination node of all the walks. 

To run, install these libraries. 

# Ubuntu installation instructions

sudo apt-get install python-pip python-dev build-essential

sudo apt-get install python-igraph python-numpy python-scipy python-matplotlib

pip install --upgrade pip

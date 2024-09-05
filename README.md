# Clique community persistence for complex networks
 
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Python3](https://img.shields.io/badge/Python3-blue.svg)
 
## Overview
 
This project is an implementation in Python of the algorithms contained in the paper [“Clique Community Persistence: A Topological Visual Analysis Approach for Complex Networks”](https://ieeexplore.ieee.org/document/8017588), that you can find in [utils.py](./utils.py).
Then there is [a simple example](./toyexample.py) to show how these algorithms work and a full example of their application to ["Les Miserables" co-occurrence network](./lesmiserables_num.gml) in [lesmiserables.py](./lesmiserables.py).
 
## Getting Started
 
### Prerequisites
 
- Python3 (I have tested it with Python3.9.6)

 
### Installation
 
Follow the steps below:
 
1. Clone the repository:
 
    ```bash
    git clone https://github.com/fabertocchi/clique-community-persistence
    ```
 
2. Move into the repository:
 
    ```bash
    cd clique-community-persistence
    ```
 
3. Create a virtual environment, then install the dependencies:

    ```bash
    pip3 install -r requirements.txt
    ```

4. Run the script:

    ```bash
    python3 lesmiserables.py
    ```
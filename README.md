# TE-SAEM v0.1

TE-SAEM is a simulated annealing expectation maximization algorithm designed to quantify transposable elements in bulk RNA-Seq data.
 
Linux Installation:  
    - [Install conda](https://www.anaconda.org)  
    - Create new conda environment: `conda create -n te-saem`  
    - Activate the environment: `conda activate te-saem`
    - Install package with: `conda install te-saem` (currently unpublished)  
    - Fetch genome and TE GTF via `scripts/fetch_resources.py`

Parameterization:  
    - File structure control is done via params.ini (automatically populated with setup.py)  
    - te-saem.py takes three required, and three optional parameters:  
        - "-1" Location of read1 (required)  
        - "-2" Location of read2 (required)  
        - "-o" Location of output count table (required)  
        - "-i" Parameter .ini path (optional)  
        - "-c" Clean intermediate files (optional)  
        - "-t" Number of threads to use (optional)  

Our manuscript with detailed algorithmic implementation and comparative assesment results is available at:  
    - Currently unpublished  

S. Texocaelum  
[Dawei Li Lab](https://dllab.org)  
[TTUHSC](https://www.ttuhsc.edu)

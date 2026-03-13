TE-SAEM is a program designed to quantify transposable elements in bulk RNA-Seq data. It requires a conda environment and utilizes the following tools:
    -- STAR aligner
    -- bedtools
 
Linux Installation:
    -- Install conda
        -- www.anaconda.org
    -- Run setup.py



Parameterization:
    -- File structure control is done via params.ini (automatically populated with setup.py)
    -- te-saem.py takes three required, and three optional parameters:
        -- "-1" Location of read1 (required)
        -- "-2" Location of read2 (required)
        -- "-o" Location of output count table (required)
        -- "-i" Parameter .ini path (optional)
        -- "-c" Clean intermediate files (optional)
        -- "-t" Number of threads to use (optional)



Our manuscript with detailed algorithmic implementation and comparative results from simulation is available at:



To-Do:
    -- Create conda package
    -- More efficient data structures

------------------------------------------------------------------------------------
Transposable Element Quantification via Simulated Annealing Expectation Maximization
    -- S. Texocaelum
------------------------------------------------------------------------------------


SAEM is an algorithm designed to quanitfy transposable elements in bulk RNA sequencing data.


setup.py
    -- Ensures all relevent programs and reference data are installed. It should be run at install, and following any changes to ../.ini
    

align.py
    -- Aligns reads with choice of aligner as indicated in ../.ini
    -- Currently STAR, Bowtie2, and segemehl are supported
    -- Output is saved as "sim_data/alignment/align.sam", with alignment summary saved as "summary.log"

parse_alignment.py
    -- Takes results of alignment and splits FASTAs into: unmapped, uniquely mapped, and multiply mapped reads
    -- These FASTAs are then indexed with and reads are aligned to them
    -- 
te_saem.py
    -- Runs the TE-SAEM algorithm to quanitfy transposable elements.

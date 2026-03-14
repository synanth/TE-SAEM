import sys
import subprocess
import os


if __name__ == '__main__':

    refs_loc = "/".join(sys.argv[1].split("/")[:-2]) + "/"
    species = sys.argv[2]
    fa_loc = refs_loc + "hs1.fa"
    idx_loc = refs_loc + "star_" + species + "_idx/"
    threads = sys.argv[3]
    
    if not os.path.exists(idx_loc):
        if not os.path.exists(fa_loc):
            print("\tGenome fasta not found, downloading from UCSC")
            wget_call = "wget -q --show-progress -P " + refs_loc + " https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz"
            subprocess.run(wget_call, shell=True)
            gunzip_call = "gunzip " + refs_loc + "hs1.fa.gz"
            subprocess.run(gunzip_call, shell=True)
            print("\tGenome fasta downloaded")
        print("\tGenerating genome index")
    
        star_call = "STAR --runThreadN " + threads + " --runMode genomeGenerate --genomeDir " + idx_loc + " --genomeFastaFiles " + fa_loc
        subprocess.run(star_call, shell=True)
    else:
        print("STAR index exists")
        quit()

## clear data ##   
clean_call = "rm -rf _STARtmp"
subprocess.run(clean_call, shell=True)

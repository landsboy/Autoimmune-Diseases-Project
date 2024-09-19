# Autoimmune-Diseases-Project
A computational research examining whether there is a connection between RNA editing and autoimmune diseases.

As part of the research, we went through 1517 known RNA editing sites (edited by the ADAR enzyme) and checked whether there are sites that produce a neoantigen by editing that are shown according to the HLA type 2 common in the population. In contrast, we checked for HLA type 2 common in T1D (patients type 1 diabetes)

# Getting help
If you need help of any kind, feel free to open a new issue.

# Setup
Requires locally installed version of [NetMHCIIpan4.3](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=netMHCIIpan&version=4.3&packageversion=4.3e&platform=Linux)

And requires a local download of the [human genome v.38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
 
1. Clone the repository
```
git clone https://github.com/landsboy/Autoimmune-Diseases-Project.git
```

2. Create conda environment
```
conda env create -f ADP.yml
conda activate ADP  # test successful creation
```
3. Run the following file with the following arguments
```
python tcga_patients.py -p <path/to/yours/netMHCpan4.3> -f <path/to/yours/hg38.fa> -v [diabetes/most_common]
```

# results
There are five options for results:

No1: The editing site does not produce either a weak binder (WB) or a strong binder (SB) because it is located in an intronic region.

No2: The editing site does not produce either WB or SB, despite being located within an exonic region.

No3: The editing site produces either WB or SB, but since the peptide is generated even in the absence of editing, it is not classified as a neo-antigen.

WB: The editing site produces a novel peptide that is likely to be presented by HLA molecules, with a relatively good binding probability.

SB: The editing site produces a novel peptide that is likely to be presented by HLA molecules, with a very high binding probability.

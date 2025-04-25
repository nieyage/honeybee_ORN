# This repository contains the data processing pipeline and code for generating figures used in the article 
"Evolutionary process underlying receptor gene expansion and cellular divergence of olfactory sensory neurons in honeybees"

## Overview

This study employs single-cell multi-omics profiling to investigate the transcriptional regulation of odorant receptor (OR) genes and their contribution to cellular identity divergence in honeybee (*Apis mellifera*) olfactory sensory neurons (OSNs). We reveal how OR gene expansion in social insects aligns with their heightened olfactory demands and uncover molecular mechanisms driving OSN diversification through polycistronic transcription and promoter regulation. These findings provide novel insights into the evolution and adaptation of olfactory systems in social insects.

---
## Data and Code Availability

### Data Access

**Raw FASTQ files and processed data**(including gene expression matrices and chromatin accessibility peak matrices) are deposited in the NCBI GEO database under accession number: **GSE248097**.  

---
### Code Repository
All analysis code for reproducing figures, processing data, and statistical modeling is available on GitHub:  
Key components include:
- Integration of single-cell transcriptomic and chromatin accessibility data
- Classification of OR expression patterns (singular, co-expression with single/multiple active promoters)
- Genomic architecture and regulatory element evolution analysis
- Reproduction scripts for all manuscript figures
## Chemereceptor 
Code used to identify OR, IG, and GR genes via domain scanning in the `Chemoreceptor/` directory. 

## Data Processing Pipeline
The data processing pipeline includes the cellranger and the bam filtering steps
You can find the code for data processing along with relevant documentation in the `Data processing/` directory.

## Figure Generation Code
Code used to generate figures for the paper resides in the `Data visualization/` directory. 

Please note that some code may require specific datasets or dependencies. Refer to the documentation within each folder for detailed information.

---
## Citation

If you use our data or code, please cite our peer-reviewed publication (to be updated).

---
## Contact & Feedback
- **Data/code issues**: Open an issue on the [GitHub Issues page](https://github.com/nieyage/honeybee_ORN/issues).  
- **Project lead**: Yage Nie (GitHub: [@nieyage](https://github.com/nieyage)).  
---
## License
- **Code**: Distributed under the [MIT License](https://opensource.org/licenses/MIT).  
- **Data**: Subject to NCBI GEO's data usage policies. See respective repositories for details.  



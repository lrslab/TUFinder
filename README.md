# TUFinder (Still developing)
An integrated tool for analyzing transcription units (TUs) in prokaryotes using RNA-seq data, with functionalities for identifying transcription start sites (TSS) and transcription termination sites (TTS).
## workflow

```mermaid
graph TD
	A1(RNA-seq)-->|alignment| C1[BAM file]
  A2(Annotation)-->|convert| C2[BED file]
    C1 --> D{extract_TUs.py} 
    C2(BED file) --> D{extract_TUs.py}
    D{extract_TUs.py} -->  E(TU result folder)
    E --> |find_longest_TU_per_gene.py| F(find_longest_TU.csv)
    C2 --> G
    C1 --> |coverage file|G
    F --> G{find_TSS_TTS.py}
    G --> H(TSS_TTS.csv)
    
```


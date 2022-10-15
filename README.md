DaPars(Dynamic analysis of Alternative PolyAdenylation from RNA-seq)
======

Sep, 2022: DaPars was updated to 1.0.0 with the following changes:

   1. Updated to python3
   2. Removed rpy2 and use python for the statistical test
   3. Fixed some minor bugs.
   
**Current version**: 0.9.1

[**Full Documentations**](http://xiazlab.org/dapars_tutorial/html/DaPars.html)
Description
-----
The dynamic usage of the 3’untranslated region (3’UTR) resulting from alternative polyadenylation (APA) is emerging as a pervasive mechanism for regulating mRNA diversity, stability and translation. Though RNA-seq provides the whole-transcriptome information and a lot of tools for analyzing gene/isoform expression are available, very few tool focus on the analysis of 3'UTR from standard RNA-seq. DaPars is the first de novo tool that directly infers the dynamic alternative polyadenylation (APA) usage by comparing standard RNA-seq. Given the annotated gene model, DaPars can infer the de novo proximal APA sites as well as the long and short 3'UTR expression levels. Finally, the dynamic APA usages between two conditions will be identified.



![Flowchart](http://farm6.staticflickr.com/5533/12003068763_87e68075f6.jpg)
![Cancer](http://farm8.staticflickr.com/7459/8858567224_4b0f0214cf.jpg)



Citation
-----
*Please cite the following articles if you use DaPars in your research*:
* Masamha, C.P., Xia, Z., Yang, J., Albrecht, T.R., Li, M., Shyu, A., Li, W., Wagner, E.J. 2014. CFIm25 links Alternative Polyadenylation to Glioblastoma Tumor Suppression. Nature, 510:412-416.

* Xia, Z., Donehower, L.A., Wheeler, D.A., Cooper, T.A., Neilson, J.R., Wagner E.J., Li, W. 2014. Dynamic Analyses of Alternative Polyadenylation from RNA-Seq Reveal 3'-UTR Landscape Across 7 Tumor Types. Nature Communications, 5:5274.

Mailing list
-----------
If you have questions, requests, or bugs to report, please email the [DaPars mailing list](https://groups.google.com/forum/#!forum/DaPars)


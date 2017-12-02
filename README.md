# Allosterome project

## Allosteric Database ASD
- Allosteric Database [ASD](http://mdl.shsmu.edu.cn/ASD/)

### PDB models
- The protein structure that are models have been downloaded under
```
ASD/ASD\_Release\_062015\_SC/
```

### Entry References
- A complete XML list of database entries has been downloaded under
```
ASD/ASD_Release_062015_XF/
```

- The XML entries are suffix numbered, like ASD02030000\_[1,2,3].xml,
where each suffix represents a subumit of the protein.

- The entire set of XML entries is only about 50MB and therefore easily
loadable into R's memory.

### XML entry heandling in R
- Scripts for processing the XML entries are located under
```
ASD/R/
```

- The *as\_list()* function of the *xml2* package converts the XML entry
to an easier-to-handle R list.



Class 11
================

``` r
db <- read.csv("Data Export Summary.csv", row.names = 1)
```

What percentage is protein?

``` r
db$
db$Proteins/db$Total
```

    ## numeric(0)

# Hands-on Section

## Section 3

``` r
library(bio3d)
pdb <- read.pdb("1hsg.pdb")
print(pdb)
```

    ## 
    ##  Call:  read.pdb(file = "1hsg.pdb")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

``` r
plot.bio3d(pdb$atom$b[pdb$calpha], sse=pdb, typ="l", ylab="B-factor")
```

![](class11_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

What type of object is “pdb$atom”?

``` r
str(pdb$atom)
```

    ## 'data.frame':    1686 obs. of  16 variables:
    ##  $ type  : chr  "ATOM" "ATOM" "ATOM" "ATOM" ...
    ##  $ eleno : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ elety : chr  "N" "CA" "C" "O" ...
    ##  $ alt   : chr  NA NA NA NA ...
    ##  $ resid : chr  "PRO" "PRO" "PRO" "PRO" ...
    ##  $ chain : chr  "A" "A" "A" "A" ...
    ##  $ resno : int  1 1 1 1 1 1 1 2 2 2 ...
    ##  $ insert: chr  NA NA NA NA ...
    ##  $ x     : num  29.4 30.3 29.8 28.6 30.5 ...
    ##  $ y     : num  39.7 38.7 38.1 38.3 37.5 ...
    ##  $ z     : num  5.86 5.32 4.02 3.68 6.34 ...
    ##  $ o     : num  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ b     : num  38.1 40.6 42.6 43.4 37.9 ...
    ##  $ segid : chr  NA NA NA NA ...
    ##  $ elesy : chr  "N" "C" "C" "O" ...
    ##  $ charge: chr  NA NA NA NA ...

``` r
# a data frame!
```

Atom selection is done via the function **atom.select()**

``` r
protein.inds <- atom.select(pdb, "protein", value = TRUE)
protein.inds
```

    ## 
    ##  Call:  trim.pdb(pdb = pdb, sele)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1514,  XYZs#: 4542  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 0  (residues: 0)
    ##      Non-protein/nucleic resid values: [ none ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

``` r
lig.inds <- atom.select(pdb, "ligand", value = TRUE)
lig.inds
```

    ## 
    ##  Call:  trim.pdb(pdb = pdb, sele)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 45,  XYZs#: 135  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 45  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

Write a *protein only* PDB file

``` r
write.pdb(protein.inds, file = "1hsg_protein.pdb")
```

Write a *ligand only* PDB file

``` r
write.pdb(lig.inds, file = "1hsg_ligand.pdb")
```

## Section 5

``` r
aa <- get.seq("1ake_A")
```

    ## Warning in get.seq("1ake_A"): Removing existing file: seqs.fasta

``` r
# BLAST search
blast.pdb(aa)
```

    ##  Searching ... please wait (updates every 5 seconds) RID = FDEN6HJ5014 
    ##  ..
    ##  Reporting 97 hits

    ## $hit.tbl
    ##         queryid subjectids identity alignmentlength mismatches gapopens
    ## 1  Query_103391     1AKE_A  100.000             214          0        0
    ## 2  Query_103391     4X8M_A   99.533             214          1        0
    ## 3  Query_103391     4X8H_A   99.533             214          1        0
    ## 4  Query_103391     3HPR_A   99.533             214          1        0
    ## 5  Query_103391     1E4V_A   99.533             214          1        0
    ## 6  Query_103391     5EJE_A   99.065             214          2        0
    ## 7  Query_103391     1E4Y_A   99.065             214          2        0
    ## 8  Query_103391     3X2S_A   98.598             214          3        0
    ## 9  Query_103391     4K46_A   73.239             213         57        0
    ## 10 Query_103391     4NP6_A   72.642             212         58        0
    ## 11 Query_103391     3GMT_A   62.500             216         75        1
    ## 12 Query_103391     4PZL_A   57.346             211         86        2
    ## 13 Query_103391     5G3Y_A   55.505             218         88        2
    ## 14 Query_103391     5G3Z_A   50.459             218         99        2
    ## 15 Query_103391     5G40_A   49.541             218        101        2
    ## 16 Query_103391     5X6J_A   50.000             218         98        3
    ## 17 Query_103391     2C9Y_A   53.723             188         83        1
    ## 18 Query_103391     1S3G_A   49.541             218         99        3
    ## 19 Query_103391     1AK2_A   52.660             188         85        1
    ## 20 Query_103391     3BE4_A   48.611             216        102        3
    ## 21 Query_103391     1AKY_A   46.119             219        108        3
    ## 22 Query_103391     3AKY_A   46.119             219        108        3
    ## 23 Query_103391     3FB4_A   48.165             218        104        2
    ## 24 Query_103391     4QBI_A   47.248             218        106        2
    ## 25 Query_103391     1DVR_A   45.205             219        110        3
    ## 26 Query_103391     3DKV_A   49.772             219         99        3
    ## 27 Query_103391     3DL0_A   48.165             218        104        2
    ## 28 Query_103391     1ZIN_A   45.413             218        110        2
    ## 29 Query_103391     2P3S_A   47.248             218        106        2
    ## 30 Query_103391     2EU8_A   47.248             218        106        2
    ## 31 Query_103391     1P3J_A   47.248             218        106        2
    ## 32 Query_103391     4QBF_A   49.772             219         99        3
    ## 33 Query_103391     2ORI_A   47.248             218        106        2
    ## 34 Query_103391     5X6I_A   46.789             218        107        2
    ## 35 Query_103391     2QAJ_A   47.005             217        106        2
    ## 36 Query_103391     2OO7_A   46.789             218        107        2
    ## 37 Query_103391     2OSB_A   46.789             218        107        2
    ## 38 Query_103391     4MKF_A   46.789             218        107        2
    ## 39 Query_103391     3TLX_A   44.393             214        106        3
    ## 40 Query_103391     4MKH_A   48.624             218        101        3
    ## 41 Query_103391     4QBH_A   45.872             218        109        2
    ## 42 Query_103391     4TYQ_A   48.165             218        102        3
    ## 43 Query_103391     4QBG_B   47.248             218        104        3
    ## 44 Query_103391     4TYP_C   47.248             218        104        3
    ## 45 Query_103391     4JKY_A   44.037             218        103        5
    ## 46 Query_103391     2RGX_A   43.578             218        104        4
    ## 47 Query_103391     4JLO_A   43.578             218        104        5
    ## 48 Query_103391     1ZAK_A   42.326             215        112        3
    ## 49 Query_103391     1ZD8_A   43.915             189         96        3
    ## 50 Query_103391     2AK3_A   44.324             185        101        2
    ## 51 Query_103391     4NTZ_A   38.532             218        119        4
    ## 52 Query_103391     2AR7_A   41.304             184        102        3
    ## 53 Query_103391     3NDP_A   40.761             184        103        3
    ## 54 Query_103391     1P4S_A   39.785             186         77        2
    ## 55 Query_103391     2CDN_A   39.785             186         77        2
    ## 56 Query_103391     3L0P_A   32.735             223        131        7
    ## 57 Query_103391     5X6L_A   35.784             204         98        3
    ## 58 Query_103391     2XB4_A   32.735             223        131        7
    ## 59 Query_103391     5XRU_A   35.294             204         99        3
    ## 60 Query_103391     5YCC_A   35.294             204         99        3
    ## 61 Query_103391     5X6K_A   35.294             204         99        3
    ## 62 Query_103391     5XZ2_A   33.962             212        107        3
    ## 63 Query_103391     5YCF_B   34.906             212        105        3
    ## 64 Query_103391     5YCB_A   34.804             204        100        3
    ## 65 Query_103391     5YCD_A   36.464             181         87        2
    ## 66 Query_103391     3ADK_A   36.066             183         89        3
    ## 67 Query_103391     3UMF_A   33.333             186         92        3
    ## 68 Query_103391     1Z83_A   34.973             183         91        3
    ## 69 Query_103391     3CM0_A   34.434             212        106        5
    ## 70 Query_103391     1UKY_A   27.962             211        117        5
    ## 71 Query_103391     1TEV_A   31.963             219        109        7
    ## 72 Query_103391     2BWJ_A   30.851             188         98        4
    ## 73 Query_103391     1QF9_A   27.907             215        117        5
    ## 74 Query_103391     1RKB_A   40.000              35         21        0
    ## 75 Query_103391     3IIJ_A   40.000              35         21        0
    ## 76 Query_103391     5YBV_B   34.545              55         33        2
    ## 77 Query_103391     5YBV_A   34.545              55         33        2
    ## 78 Query_103391     4HBD_A   34.545              55         33        2
    ## 79 Query_103391     5JZV_A   40.000              35         21        0
    ## 80 Query_103391     5NGH_A   40.000              40         24        0
    ## 81 Query_103391     1Q3T_A   36.842              38         24        0
    ## 82 Query_103391     2GA8_A   53.846              26          8        1
    ## 83 Query_103391     4CVN_A   38.235              34         21        0
    ## 84 Query_103391     4XRU_A   52.632              19          9        0
    ## 85 Query_103391     1XE4_A   31.667              60         35        1
    ## 86 Query_103391     4II9_A   31.667              60         35        1
    ## 87 Query_103391     1XF8_A   43.333              30         17        0
    ## 88 Query_103391     1NE9_A   43.333              30         17        0
    ## 89 Query_103391     3GKR_A   43.333              30         17        0
    ## 90 Query_103391     4XRP_A   52.632              19          9        0
    ## 91 Query_103391     3DC0_A   38.636              44         27        0
    ## 92 Query_103391     5TUE_A   36.364              55         30        1
    ## 93 Query_103391     1J6P_A   40.541              37         21        1
    ## 94 Query_103391     1UA7_A   38.636              44         27        0
    ## 95 Query_103391     1BAG_A   38.636              44         27        0
    ## 96 Query_103391     1P1M_A   40.541              37         21        1
    ## 97 Query_103391     5YWW_A   59.091              22          9        0
    ##    q.start q.end s.start s.end    evalue bitscore positives mlog.evalue
    ## 1        1   214       1   214 9.27e-157    432.0    100.00 359.2790762
    ## 2        1   214       1   214 1.66e-156    432.0    100.00 358.6964569
    ## 3        1   214       1   214 9.18e-156    430.0     99.53 356.9862473
    ## 4        1   214       1   214 1.30e-155    430.0     99.53 356.6383251
    ## 5        1   214       1   214 1.38e-155    430.0     99.53 356.5786059
    ## 6        1   214       1   214 4.13e-155    429.0     99.07 355.4824120
    ## 7        1   214       1   214 2.39e-154    427.0     99.07 353.7268110
    ## 8        1   214       1   214 4.00e-154    426.0     98.60 353.2118100
    ## 9        1   213       1   213 1.07e-115    329.0     84.98 264.7296270
    ## 10       2   213       5   216 5.96e-114    325.0     84.43 260.7096301
    ## 11       2   211      10   225  4.66e-90    265.0     71.30 205.6936429
    ## 12       2   209      26   235  1.11e-86    256.0     74.41 197.9179580
    ## 13       1   214       1   213  1.65e-76    230.0     68.81 174.4956918
    ## 14       1   214       1   213  3.15e-73    221.0     69.27 166.9413093
    ## 15       1   214       1   213  6.46e-71    216.0     68.35 161.6179123
    ## 16       1   213       1   212  1.18e-68    210.0     65.60 156.4102719
    ## 17       1   184      17   204  5.95e-68    209.0     69.68 154.7923951
    ## 18       1   213       1   212  6.12e-68    208.0     65.14 154.7642242
    ## 19       1   184      17   204  1.54e-67    207.0     70.21 153.8414188
    ## 20       2   213       7   217  4.13e-67    206.0     68.06 152.8549238
    ## 21       1   214       5   218  5.83e-67    206.0     65.75 152.5101842
    ## 22       1   214       5   218  7.09e-67    205.0     65.30 152.3145159
    ## 23       1   214       1   213  2.61e-66    204.0     65.14 151.0112659
    ## 24       1   214       1   213  1.26e-64    199.0     65.60 147.1343342
    ## 25       1   214       5   218  2.72e-64    199.0     64.84 146.3648141
    ## 26       1   214       1   213  1.46e-63    197.0     66.67 144.6844244
    ## 27       1   214       1   213  6.45e-63    195.0     66.97 143.1987807
    ## 28       1   214       1   213  7.18e-63    195.0     65.60 143.0915615
    ## 29       1   214       1   213  1.09e-62    194.0     66.97 142.6740981
    ## 30       1   214       1   213  1.14e-62    194.0     66.97 142.6292475
    ## 31       1   214       1   213  1.54e-62    194.0     66.97 142.3284933
    ## 32       1   214       1   213  2.46e-62    194.0     66.21 141.8601144
    ## 33       1   214       1   213  3.85e-62    193.0     66.97 141.4122026
    ## 34       1   214       1   213  7.31e-62    192.0     66.51 140.7710325
    ## 35       1   213       1   212  7.89e-62    192.0     66.82 140.6946796
    ## 36       1   214       1   213  9.18e-62    192.0     66.51 140.5432486
    ## 37       1   214       1   213  1.57e-61    192.0     66.51 140.0066151
    ## 38       1   214       1   213  1.75e-61    191.0     66.51 139.8980749
    ## 39       2   211      31   235  1.10e-60    190.0     64.95 138.0597954
    ## 40       1   213       3   214  1.13e-60    189.0     65.60 138.0328879
    ## 41       1   214       1   213  2.26e-60    189.0     64.68 137.3397408
    ## 42       1   213       1   212  2.77e-60    188.0     65.60 137.1362583
    ## 43       1   213       1   212  4.86e-59    185.0     65.60 134.2714820
    ## 44       1   213       1   212  9.13e-59    184.0     65.14 133.6409548
    ## 45       1   214       1   203  4.19e-56    177.0     66.51 127.5120645
    ## 46       1   214       1   203  4.80e-56    177.0     65.60 127.3761493
    ## 47       1   214       1   203  1.43e-55    176.0     66.51 126.2845057
    ## 48       1   214       6   209  2.88e-54    173.0     63.72 123.2818047
    ## 49       1   185       8   190  2.24e-50    164.0     64.55 114.3227788
    ## 50       1   185       7   189  3.08e-50    163.0     65.41 114.0043251
    ## 51       1   213       6   213  1.10e-46    154.0     62.39 105.8236041
    ## 52       1   182      28   207  3.88e-45    150.0     64.13 102.2604940
    ## 53       1   182       6   185  3.25e-44    148.0     63.59 100.1350891
    ## 54       1   182       1   155  5.39e-39    133.0     56.99  88.1162732
    ## 55       1   182      21   175  9.70e-39    133.0     56.99  87.5286927
    ## 56       1   209       1   218  2.07e-31    115.0     54.26  70.6525893
    ## 57       3   205      13   184  2.53e-31    114.0     52.45  70.4519186
    ## 58       1   209       1   218  2.74e-31    114.0     54.26  70.3721800
    ## 59       3   205      11   182  3.04e-31    113.0     52.45  70.2682804
    ## 60       3   205      11   182  3.57e-31    113.0     52.45  70.1075723
    ## 61       3   205      13   184  4.23e-31    113.0     52.45  69.9379359
    ## 62       3   213      13   192  4.99e-31    113.0     52.83  69.7727020
    ## 63       3   213      11   190  5.09e-31    113.0     51.89  69.7528601
    ## 64       3   205      11   182  5.48e-31    113.0     52.45  69.6790328
    ## 65       3   182      11   164  1.45e-30    112.0     54.70  68.7059892
    ## 66       3   184      12   167  1.74e-30    111.0     54.10  68.5236677
    ## 67       3   185      32   188  9.21e-30    110.0     56.45  66.8572629
    ## 68       3   184      12   167  1.08e-29    109.0     54.64  66.6980067
    ## 69       3   214       7   185  5.21e-29    107.0     50.94  65.1243878
    ## 70       3   210      18   196  5.63e-25     97.8     53.55  55.8365179
    ## 71       3   213       6   192  3.99e-24     95.5     49.77  53.8782510
    ## 72       3   187      15   173  4.12e-22     90.1     49.47  49.2410189
    ## 73       3   213       9   189  1.60e-20     85.9     47.44  45.5816982
    ## 74       2    36       6    40  4.30e-01     31.6     57.14   0.8439701
    ## 75       2    36      13    47  5.30e-01     31.6     57.14   0.6348783
    ## 76     153   204      59   113  6.80e-01     31.6     49.09   0.3856625
    ## 77     153   204      56   110  7.20e-01     31.2     49.09   0.3285041
    ## 78     153   204      74   128  7.90e-01     31.2     49.09   0.2357223
    ## 79       2    36      39    73  8.00e-01     31.2     57.14   0.2231436
    ## 80     134   173     102   141  8.50e-01     30.8     55.00   0.1625189
    ## 81       1    38      17    54  2.80e+00     29.6     57.89  -1.0296194
    ## 82       3    24      27    52  4.80e+00     28.9     69.23  -1.5686159
    ## 83       1    34      12    45  6.10e+00     28.5     52.94  -1.8082888
    ## 84       3    21      14    32  6.50e+00     28.5     84.21  -1.8718022
    ## 85     110   169      80   133  6.90e+00     28.5     45.00  -1.9315214
    ## 86     110   169      81   134  7.00e+00     28.5     45.00  -1.9459101
    ## 87     110   139      80   109  7.10e+00     28.5     63.33  -1.9600948
    ## 88     110   139      80   109  7.10e+00     28.5     63.33  -1.9600948
    ## 89     110   139      81   110  7.10e+00     28.5     63.33  -1.9600948
    ## 90       3    21      14    32  7.40e+00     28.5     84.21  -2.0014800
    ## 91      27    70       1    44  7.40e+00     28.5     45.45  -2.0014800
    ## 92     116   170      82   131  7.40e+00     28.5     43.64  -2.0014800
    ## 93      98   133     345   381  7.70e+00     28.5     54.05  -2.0412203
    ## 94      27    70       1    44  7.80e+00     28.5     45.45  -2.0541237
    ## 95      27    70       4    47  8.20e+00     28.5     45.45  -2.1041342
    ## 96      98   133     333   369  9.50e+00     28.1     54.05  -2.2512918
    ## 97       3    24     251   272  9.80e+00     28.1     68.18  -2.2823824
    ##    pdb.id    acc
    ## 1  1AKE_A 1AKE_A
    ## 2  4X8M_A 4X8M_A
    ## 3  4X8H_A 4X8H_A
    ## 4  3HPR_A 3HPR_A
    ## 5  1E4V_A 1E4V_A
    ## 6  5EJE_A 5EJE_A
    ## 7  1E4Y_A 1E4Y_A
    ## 8  3X2S_A 3X2S_A
    ## 9  4K46_A 4K46_A
    ## 10 4NP6_A 4NP6_A
    ## 11 3GMT_A 3GMT_A
    ## 12 4PZL_A 4PZL_A
    ## 13 5G3Y_A 5G3Y_A
    ## 14 5G3Z_A 5G3Z_A
    ## 15 5G40_A 5G40_A
    ## 16 5X6J_A 5X6J_A
    ## 17 2C9Y_A 2C9Y_A
    ## 18 1S3G_A 1S3G_A
    ## 19 1AK2_A 1AK2_A
    ## 20 3BE4_A 3BE4_A
    ## 21 1AKY_A 1AKY_A
    ## 22 3AKY_A 3AKY_A
    ## 23 3FB4_A 3FB4_A
    ## 24 4QBI_A 4QBI_A
    ## 25 1DVR_A 1DVR_A
    ## 26 3DKV_A 3DKV_A
    ## 27 3DL0_A 3DL0_A
    ## 28 1ZIN_A 1ZIN_A
    ## 29 2P3S_A 2P3S_A
    ## 30 2EU8_A 2EU8_A
    ## 31 1P3J_A 1P3J_A
    ## 32 4QBF_A 4QBF_A
    ## 33 2ORI_A 2ORI_A
    ## 34 5X6I_A 5X6I_A
    ## 35 2QAJ_A 2QAJ_A
    ## 36 2OO7_A 2OO7_A
    ## 37 2OSB_A 2OSB_A
    ## 38 4MKF_A 4MKF_A
    ## 39 3TLX_A 3TLX_A
    ## 40 4MKH_A 4MKH_A
    ## 41 4QBH_A 4QBH_A
    ## 42 4TYQ_A 4TYQ_A
    ## 43 4QBG_B 4QBG_B
    ## 44 4TYP_C 4TYP_C
    ## 45 4JKY_A 4JKY_A
    ## 46 2RGX_A 2RGX_A
    ## 47 4JLO_A 4JLO_A
    ## 48 1ZAK_A 1ZAK_A
    ## 49 1ZD8_A 1ZD8_A
    ## 50 2AK3_A 2AK3_A
    ## 51 4NTZ_A 4NTZ_A
    ## 52 2AR7_A 2AR7_A
    ## 53 3NDP_A 3NDP_A
    ## 54 1P4S_A 1P4S_A
    ## 55 2CDN_A 2CDN_A
    ## 56 3L0P_A 3L0P_A
    ## 57 5X6L_A 5X6L_A
    ## 58 2XB4_A 2XB4_A
    ## 59 5XRU_A 5XRU_A
    ## 60 5YCC_A 5YCC_A
    ## 61 5X6K_A 5X6K_A
    ## 62 5XZ2_A 5XZ2_A
    ## 63 5YCF_B 5YCF_B
    ## 64 5YCB_A 5YCB_A
    ## 65 5YCD_A 5YCD_A
    ## 66 3ADK_A 3ADK_A
    ## 67 3UMF_A 3UMF_A
    ## 68 1Z83_A 1Z83_A
    ## 69 3CM0_A 3CM0_A
    ## 70 1UKY_A 1UKY_A
    ## 71 1TEV_A 1TEV_A
    ## 72 2BWJ_A 2BWJ_A
    ## 73 1QF9_A 1QF9_A
    ## 74 1RKB_A 1RKB_A
    ## 75 3IIJ_A 3IIJ_A
    ## 76 5YBV_B 5YBV_B
    ## 77 5YBV_A 5YBV_A
    ## 78 4HBD_A 4HBD_A
    ## 79 5JZV_A 5JZV_A
    ## 80 5NGH_A 5NGH_A
    ## 81 1Q3T_A 1Q3T_A
    ## 82 2GA8_A 2GA8_A
    ## 83 4CVN_A 4CVN_A
    ## 84 4XRU_A 4XRU_A
    ## 85 1XE4_A 1XE4_A
    ## 86 4II9_A 4II9_A
    ## 87 1XF8_A 1XF8_A
    ## 88 1NE9_A 1NE9_A
    ## 89 3GKR_A 3GKR_A
    ## 90 4XRP_A 4XRP_A
    ## 91 3DC0_A 3DC0_A
    ## 92 5TUE_A 5TUE_A
    ## 93 1J6P_A 1J6P_A
    ## 94 1UA7_A 1UA7_A
    ## 95 1BAG_A 1BAG_A
    ## 96 1P1M_A 1P1M_A
    ## 97 5YWW_A 5YWW_A
    ## 
    ## $raw
    ##         queryid subjectids identity alignmentlength mismatches gapopens
    ## 1  Query_103391     1AKE_A  100.000             214          0        0
    ## 2  Query_103391     4X8M_A   99.533             214          1        0
    ## 3  Query_103391     4X8H_A   99.533             214          1        0
    ## 4  Query_103391     3HPR_A   99.533             214          1        0
    ## 5  Query_103391     1E4V_A   99.533             214          1        0
    ## 6  Query_103391     5EJE_A   99.065             214          2        0
    ## 7  Query_103391     1E4Y_A   99.065             214          2        0
    ## 8  Query_103391     3X2S_A   98.598             214          3        0
    ## 9  Query_103391     4K46_A   73.239             213         57        0
    ## 10 Query_103391     4NP6_A   72.642             212         58        0
    ## 11 Query_103391     3GMT_A   62.500             216         75        1
    ## 12 Query_103391     4PZL_A   57.346             211         86        2
    ## 13 Query_103391     5G3Y_A   55.505             218         88        2
    ## 14 Query_103391     5G3Z_A   50.459             218         99        2
    ## 15 Query_103391     5G40_A   49.541             218        101        2
    ## 16 Query_103391     5X6J_A   50.000             218         98        3
    ## 17 Query_103391     2C9Y_A   53.723             188         83        1
    ## 18 Query_103391     1S3G_A   49.541             218         99        3
    ## 19 Query_103391     1AK2_A   52.660             188         85        1
    ## 20 Query_103391     3BE4_A   48.611             216        102        3
    ## 21 Query_103391     1AKY_A   46.119             219        108        3
    ## 22 Query_103391     3AKY_A   46.119             219        108        3
    ## 23 Query_103391     3FB4_A   48.165             218        104        2
    ## 24 Query_103391     4QBI_A   47.248             218        106        2
    ## 25 Query_103391     1DVR_A   45.205             219        110        3
    ## 26 Query_103391     3DKV_A   49.772             219         99        3
    ## 27 Query_103391     3DL0_A   48.165             218        104        2
    ## 28 Query_103391     1ZIN_A   45.413             218        110        2
    ## 29 Query_103391     2P3S_A   47.248             218        106        2
    ## 30 Query_103391     2EU8_A   47.248             218        106        2
    ## 31 Query_103391     1P3J_A   47.248             218        106        2
    ## 32 Query_103391     4QBF_A   49.772             219         99        3
    ## 33 Query_103391     2ORI_A   47.248             218        106        2
    ## 34 Query_103391     5X6I_A   46.789             218        107        2
    ## 35 Query_103391     2QAJ_A   47.005             217        106        2
    ## 36 Query_103391     2OO7_A   46.789             218        107        2
    ## 37 Query_103391     2OSB_A   46.789             218        107        2
    ## 38 Query_103391     4MKF_A   46.789             218        107        2
    ## 39 Query_103391     3TLX_A   44.393             214        106        3
    ## 40 Query_103391     4MKH_A   48.624             218        101        3
    ## 41 Query_103391     4QBH_A   45.872             218        109        2
    ## 42 Query_103391     4TYQ_A   48.165             218        102        3
    ## 43 Query_103391     4QBG_B   47.248             218        104        3
    ## 44 Query_103391     4TYP_C   47.248             218        104        3
    ## 45 Query_103391     4JKY_A   44.037             218        103        5
    ## 46 Query_103391     2RGX_A   43.578             218        104        4
    ## 47 Query_103391     4JLO_A   43.578             218        104        5
    ## 48 Query_103391     1ZAK_A   42.326             215        112        3
    ## 49 Query_103391     1ZD8_A   43.915             189         96        3
    ## 50 Query_103391     2AK3_A   44.324             185        101        2
    ## 51 Query_103391     4NTZ_A   38.532             218        119        4
    ## 52 Query_103391     2AR7_A   41.304             184        102        3
    ## 53 Query_103391     3NDP_A   40.761             184        103        3
    ## 54 Query_103391     1P4S_A   39.785             186         77        2
    ## 55 Query_103391     2CDN_A   39.785             186         77        2
    ## 56 Query_103391     3L0P_A   32.735             223        131        7
    ## 57 Query_103391     5X6L_A   35.784             204         98        3
    ## 58 Query_103391     2XB4_A   32.735             223        131        7
    ## 59 Query_103391     5XRU_A   35.294             204         99        3
    ## 60 Query_103391     5YCC_A   35.294             204         99        3
    ## 61 Query_103391     5X6K_A   35.294             204         99        3
    ## 62 Query_103391     5XZ2_A   33.962             212        107        3
    ## 63 Query_103391     5YCF_B   34.906             212        105        3
    ## 64 Query_103391     5YCB_A   34.804             204        100        3
    ## 65 Query_103391     5YCD_A   36.464             181         87        2
    ## 66 Query_103391     3ADK_A   36.066             183         89        3
    ## 67 Query_103391     3UMF_A   33.333             186         92        3
    ## 68 Query_103391     1Z83_A   34.973             183         91        3
    ## 69 Query_103391     3CM0_A   34.434             212        106        5
    ## 70 Query_103391     1UKY_A   27.962             211        117        5
    ## 71 Query_103391     1TEV_A   31.963             219        109        7
    ## 72 Query_103391     2BWJ_A   30.851             188         98        4
    ## 73 Query_103391     1QF9_A   27.907             215        117        5
    ## 74 Query_103391     1RKB_A   40.000              35         21        0
    ## 75 Query_103391     3IIJ_A   40.000              35         21        0
    ## 76 Query_103391     5YBV_B   34.545              55         33        2
    ## 77 Query_103391     5YBV_A   34.545              55         33        2
    ## 78 Query_103391     4HBD_A   34.545              55         33        2
    ## 79 Query_103391     5JZV_A   40.000              35         21        0
    ## 80 Query_103391     5NGH_A   40.000              40         24        0
    ## 81 Query_103391     1Q3T_A   36.842              38         24        0
    ## 82 Query_103391     2GA8_A   53.846              26          8        1
    ## 83 Query_103391     4CVN_A   38.235              34         21        0
    ## 84 Query_103391     4XRU_A   52.632              19          9        0
    ## 85 Query_103391     1XE4_A   31.667              60         35        1
    ## 86 Query_103391     4II9_A   31.667              60         35        1
    ## 87 Query_103391     1XF8_A   43.333              30         17        0
    ## 88 Query_103391     1NE9_A   43.333              30         17        0
    ## 89 Query_103391     3GKR_A   43.333              30         17        0
    ## 90 Query_103391     4XRP_A   52.632              19          9        0
    ## 91 Query_103391     3DC0_A   38.636              44         27        0
    ## 92 Query_103391     5TUE_A   36.364              55         30        1
    ## 93 Query_103391     1J6P_A   40.541              37         21        1
    ## 94 Query_103391     1UA7_A   38.636              44         27        0
    ## 95 Query_103391     1BAG_A   38.636              44         27        0
    ## 96 Query_103391     1P1M_A   40.541              37         21        1
    ## 97 Query_103391     5YWW_A   59.091              22          9        0
    ##    q.start q.end s.start s.end    evalue bitscore positives
    ## 1        1   214       1   214 9.27e-157    432.0    100.00
    ## 2        1   214       1   214 1.66e-156    432.0    100.00
    ## 3        1   214       1   214 9.18e-156    430.0     99.53
    ## 4        1   214       1   214 1.30e-155    430.0     99.53
    ## 5        1   214       1   214 1.38e-155    430.0     99.53
    ## 6        1   214       1   214 4.13e-155    429.0     99.07
    ## 7        1   214       1   214 2.39e-154    427.0     99.07
    ## 8        1   214       1   214 4.00e-154    426.0     98.60
    ## 9        1   213       1   213 1.07e-115    329.0     84.98
    ## 10       2   213       5   216 5.96e-114    325.0     84.43
    ## 11       2   211      10   225  4.66e-90    265.0     71.30
    ## 12       2   209      26   235  1.11e-86    256.0     74.41
    ## 13       1   214       1   213  1.65e-76    230.0     68.81
    ## 14       1   214       1   213  3.15e-73    221.0     69.27
    ## 15       1   214       1   213  6.46e-71    216.0     68.35
    ## 16       1   213       1   212  1.18e-68    210.0     65.60
    ## 17       1   184      17   204  5.95e-68    209.0     69.68
    ## 18       1   213       1   212  6.12e-68    208.0     65.14
    ## 19       1   184      17   204  1.54e-67    207.0     70.21
    ## 20       2   213       7   217  4.13e-67    206.0     68.06
    ## 21       1   214       5   218  5.83e-67    206.0     65.75
    ## 22       1   214       5   218  7.09e-67    205.0     65.30
    ## 23       1   214       1   213  2.61e-66    204.0     65.14
    ## 24       1   214       1   213  1.26e-64    199.0     65.60
    ## 25       1   214       5   218  2.72e-64    199.0     64.84
    ## 26       1   214       1   213  1.46e-63    197.0     66.67
    ## 27       1   214       1   213  6.45e-63    195.0     66.97
    ## 28       1   214       1   213  7.18e-63    195.0     65.60
    ## 29       1   214       1   213  1.09e-62    194.0     66.97
    ## 30       1   214       1   213  1.14e-62    194.0     66.97
    ## 31       1   214       1   213  1.54e-62    194.0     66.97
    ## 32       1   214       1   213  2.46e-62    194.0     66.21
    ## 33       1   214       1   213  3.85e-62    193.0     66.97
    ## 34       1   214       1   213  7.31e-62    192.0     66.51
    ## 35       1   213       1   212  7.89e-62    192.0     66.82
    ## 36       1   214       1   213  9.18e-62    192.0     66.51
    ## 37       1   214       1   213  1.57e-61    192.0     66.51
    ## 38       1   214       1   213  1.75e-61    191.0     66.51
    ## 39       2   211      31   235  1.10e-60    190.0     64.95
    ## 40       1   213       3   214  1.13e-60    189.0     65.60
    ## 41       1   214       1   213  2.26e-60    189.0     64.68
    ## 42       1   213       1   212  2.77e-60    188.0     65.60
    ## 43       1   213       1   212  4.86e-59    185.0     65.60
    ## 44       1   213       1   212  9.13e-59    184.0     65.14
    ## 45       1   214       1   203  4.19e-56    177.0     66.51
    ## 46       1   214       1   203  4.80e-56    177.0     65.60
    ## 47       1   214       1   203  1.43e-55    176.0     66.51
    ## 48       1   214       6   209  2.88e-54    173.0     63.72
    ## 49       1   185       8   190  2.24e-50    164.0     64.55
    ## 50       1   185       7   189  3.08e-50    163.0     65.41
    ## 51       1   213       6   213  1.10e-46    154.0     62.39
    ## 52       1   182      28   207  3.88e-45    150.0     64.13
    ## 53       1   182       6   185  3.25e-44    148.0     63.59
    ## 54       1   182       1   155  5.39e-39    133.0     56.99
    ## 55       1   182      21   175  9.70e-39    133.0     56.99
    ## 56       1   209       1   218  2.07e-31    115.0     54.26
    ## 57       3   205      13   184  2.53e-31    114.0     52.45
    ## 58       1   209       1   218  2.74e-31    114.0     54.26
    ## 59       3   205      11   182  3.04e-31    113.0     52.45
    ## 60       3   205      11   182  3.57e-31    113.0     52.45
    ## 61       3   205      13   184  4.23e-31    113.0     52.45
    ## 62       3   213      13   192  4.99e-31    113.0     52.83
    ## 63       3   213      11   190  5.09e-31    113.0     51.89
    ## 64       3   205      11   182  5.48e-31    113.0     52.45
    ## 65       3   182      11   164  1.45e-30    112.0     54.70
    ## 66       3   184      12   167  1.74e-30    111.0     54.10
    ## 67       3   185      32   188  9.21e-30    110.0     56.45
    ## 68       3   184      12   167  1.08e-29    109.0     54.64
    ## 69       3   214       7   185  5.21e-29    107.0     50.94
    ## 70       3   210      18   196  5.63e-25     97.8     53.55
    ## 71       3   213       6   192  3.99e-24     95.5     49.77
    ## 72       3   187      15   173  4.12e-22     90.1     49.47
    ## 73       3   213       9   189  1.60e-20     85.9     47.44
    ## 74       2    36       6    40  4.30e-01     31.6     57.14
    ## 75       2    36      13    47  5.30e-01     31.6     57.14
    ## 76     153   204      59   113  6.80e-01     31.6     49.09
    ## 77     153   204      56   110  7.20e-01     31.2     49.09
    ## 78     153   204      74   128  7.90e-01     31.2     49.09
    ## 79       2    36      39    73  8.00e-01     31.2     57.14
    ## 80     134   173     102   141  8.50e-01     30.8     55.00
    ## 81       1    38      17    54  2.80e+00     29.6     57.89
    ## 82       3    24      27    52  4.80e+00     28.9     69.23
    ## 83       1    34      12    45  6.10e+00     28.5     52.94
    ## 84       3    21      14    32  6.50e+00     28.5     84.21
    ## 85     110   169      80   133  6.90e+00     28.5     45.00
    ## 86     110   169      81   134  7.00e+00     28.5     45.00
    ## 87     110   139      80   109  7.10e+00     28.5     63.33
    ## 88     110   139      80   109  7.10e+00     28.5     63.33
    ## 89     110   139      81   110  7.10e+00     28.5     63.33
    ## 90       3    21      14    32  7.40e+00     28.5     84.21
    ## 91      27    70       1    44  7.40e+00     28.5     45.45
    ## 92     116   170      82   131  7.40e+00     28.5     43.64
    ## 93      98   133     345   381  7.70e+00     28.5     54.05
    ## 94      27    70       1    44  7.80e+00     28.5     45.45
    ## 95      27    70       4    47  8.20e+00     28.5     45.45
    ## 96      98   133     333   369  9.50e+00     28.1     54.05
    ## 97       3    24     251   272  9.80e+00     28.1     68.18
    ## 
    ## $url
    ##                                                                                                                                                        FDEN6HJ5014 
    ## "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=Alignment&ALIGNMENT_VIEW=Tabular&RESULTS_FILE=on&FORMAT_TYPE=CSV&ALIGNMENTS=20000&RID=FDEN6HJ5014" 
    ## 
    ## attr(,"class")
    ## [1] "blast"
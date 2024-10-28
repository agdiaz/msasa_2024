from io import StringIO

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sequences = [
    'BBS12014',
    'BB11013',
    'BB11004',
    'BBS20002',
    'BBS20009',
    'BBS20007',
    'BBS30016',
    'BBS30017',
    'BBS30001',
    'BBS30015',
    'BBS30003',
    'BB40010',
    'BB40014',
    'BBS50005',
    'BOX132',
    'BOX212',
    'BOX122',
    'BBA0142',
    'BBA0117',
    'BBA0030',
    'BBA0011',
    'BBA0192',
    'BBA0065',
]
dimension_length = [
    52,
    101,
    410,
    49,
    219,
    433,
    66,
    287,
    190,
    650,
    514,
    214,
    609,
    806,
    159,
    625,
    1041,
    78,
    77,
    165,
    251,
    799,
    786
]
dimension_quantity = [
    9,
    5,
    4,
    20,
    29,
    23,
    37,
    15,
    116,
    21,
    142,
    9,
    9,
    11,
    13,
    8,
    174,
    19,
    460,
    19,
    102,
    4,
    100
]
group_length = [
    'CL',
    'ML',
    'LL',
    'CL',
    'ML',
    'LL',
    'CL',
    'ML',
    'ML',
    'LL',
    'LL',
    'ML',
    'LL',
    'LL',
    'ML',
    'LL',
    'LL',
    'CL',
    'CL',
    'ML',
    'ML',
    'LL',
    'LL'
]
group_quantity = [
    'BC',
    'BC',
    'BC',
    'BC',
    'BC',
    'BC',
    'BC',
    'BC',
    'AC',
    'BC',
    'AC',
    'BC',
    'BC',
    'BC',
    'BC',
    'BC',
    'AC',
    'BC',
    'AC',
    'BC',
    'AC',
    'BC',
    'AC'
]
study_case = [
    'C1',
    'C2',
    'C3',
    'C4',
    'C5',
    'C6',
    'C7',
    'C8',
    'C9',
    'C10',
    'C11',
    'C12',
    'C13',
    'C14',
    'C15',
    'C16',
    'C17',
    'C18',
    'C19',
    'C20',
    'C21',
    'C22',
    'C23',
]

data = """
Caso de Estudio;Software;Overlap Score
BB11004;MAFFT;0.608653
BB11004;MUSCLE;0.607178
BB11004;TCOFFE;0.604228
BB11004;KAlign;0.546214
BB11004;Clustal;0.506391
BB11004;MSASA Similitud Gonnet92 Estricto;0.040315
BB11004;MSASA Similitud Blosum62 Estricto;0.039823
BB11004;MSASA Similitud Gonnet92 Libre;0.029990
BB11004;MSASA Coincidencias Estricto;0.024090
BB11004;MSASA Global Libre;0.023599
BB11004;MSASA Similitud PAM250 Libre;0.016716
BB11004;MSASA Similitud Blosum62 Libre;0.014258
BB11004;MSASA Global Estricto;0.011799
BB11004;MSASA Local Estricto;0.010816
BB11004;MSASA Coincidencias Libre;0.010324
BB11004;MSASA Similitud PAM250 Estricto;0.010324
BB11004;MSASA Identidad Estricto;0.008850
BB11004;MSASA Identidad Libre;0.006391
BB11004;MSASA Local Libre;0.017699
BB11013;KAlign;0.301075
BB11013;MUSCLE;0.163082
BB11013;MSASA Similitud Blosum62 Libre;0.137993
BB11013;TCOFFE;0.130824
BB11013;MAFFT;0.089606
BB11013;MSASA Coincidencias Estricto;0.089606
BB11013;MSASA Similitud PAM250 Estricto;0.089606
BB11013;MSASA Similitud PAM250 Libre;0.078853
BB11013;MSASA Local Estricto;0.069892
BB11013;MSASA Similitud Gonnet92 Libre;0.064516
BB11013;Clustal;0.048387
BB11013;MSASA Similitud Gonnet92 Estricto;0.046595
BB11013;MSASA Global Estricto;0.030466
BB11013;MSASA Coincidencias Libre;0.026882
BB11013;MSASA Global Libre;0.017921
BB11013;MSASA Identidad Estricto;0.016129
BB11013;MSASA Identidad Libre;0.010753
BB11013;MSASA Similitud Blosum62 Estricto;0.000000
BB11013;MSASA Local Libre;0.035842
BB40010;Clustal;0.815350
BB40010;MAFFT;0.801672
BB40010;TCOFFE;0.798632
BB40010;MUSCLE;0.798252
BB40010;KAlign;0.780775
BB40010;MSASA Similitud PAM250 Estricto;0.445289
BB40010;MSASA Similitud PAM250 Libre;0.417933
BB40010;MSASA Similitud Gonnet92 Estricto;0.402736
BB40010;MSASA Similitud Gonnet92 Libre;0.335106
BB40010;MSASA Similitud Blosum62 Libre;0.321809
BB40010;MSASA Similitud Blosum62 Estricto;0.306231
BB40010;MSASA Local Estricto;0.147796
BB40010;MSASA Global Libre;0.140578
BB40010;MSASA Global Estricto;0.136778
BB40010;MSASA Coincidencias Libre;0.053191
BB40010;MSASA Coincidencias Estricto;0.039514
BB40010;MSASA Identidad Libre;0.024696
BB40010;MSASA Identidad Estricto;0.006459
BB40010;MSASA Local Libre;0.155015
BB40014;TCOFFE;0.638698
BB40014;Clustal;0.638218
BB40014;MUSCLE;0.632859
BB40014;MAFFT;0.630139
BB40014;KAlign;0.605263
BB40014;MSASA Similitud Gonnet92 Libre;0.114622
BB40014;MSASA Similitud Gonnet92 Estricto;0.113742
BB40014;MSASA Similitud PAM250 Estricto;0.094865
BB40014;MSASA Coincidencias Estricto;0.046553
BB40014;MSASA Coincidencias Libre;0.037994
BB40014;MSASA Local Estricto;0.035514
BB40014;MSASA Similitud Blosum62 Libre;0.026636
BB40014;MSASA Similitud Blosum62 Estricto;0.024636
BB40014;MSASA Similitud PAM250 Libre;0.015677
BB40014;MSASA Global Estricto;0.013518
BB40014;MSASA Identidad Libre;0.011358
BB40014;MSASA Global Libre;0.007199
BB40014;MSASA Identidad Estricto;0.005839
BB40014;MSASA Local Libre;0.015518
BBA0011;MAFFT;0.983106
BBA0011;Clustal;0.980758
BBA0011;TCOFFE;0.976574
BBA0011;MUSCLE;0.974694
BBA0011;KAlign;0.965540
BBA0011;MSASA Similitud Gonnet92 Libre;0.732607
BBA0011;MSASA Similitud Blosum62 Libre;0.718393
BBA0011;MSASA Similitud Blosum62 Estricto;0.708311
BBA0011;MSASA Global Estricto;0.696062
BBA0011;MSASA Global Libre;0.691707
BBA0011;MSASA Coincidencias Libre;0.362936
BBA0011;MSASA Similitud PAM250 Estricto;0.258648
BBA0011;MSASA Similitud Gonnet92 Estricto;0.256632
BBA0011;MSASA Local Estricto;0.256252
BBA0011;MSASA Similitud PAM250 Libre;0.239002
BBA0011;MSASA Coincidencias Estricto;0.137981
BBA0011;MSASA Identidad Libre;0.051429
BBA0011;MSASA Identidad Estricto;0.047051
BBA0011;MSASA Local Libre;0.726552
BBA0030;MAFFT;0.944087
BBA0030;Clustal;0.900249
BBA0030;TCOFFE;0.893949
BBA0030;MUSCLE;0.890668
BBA0030;KAlign;0.869449
BBA0030;MSASA Similitud Blosum62 Estricto;0.476397
BBA0030;MSASA Similitud Gonnet92 Estricto;0.440478
BBA0030;MSASA Similitud Gonnet92 Libre;0.410640
BBA0030;MSASA Similitud PAM250 Libre;0.408321
BBA0030;MSASA Similitud Blosum62 Libre;0.401365
BBA0030;MSASA Similitud PAM250 Estricto;0.353283
BBA0030;MSASA Local Estricto;0.295271
BBA0030;MSASA Global Estricto;0.272564
BBA0030;MSASA Global Libre;0.260795
BBA0030;MSASA Coincidencias Estricto;0.225270
BBA0030;MSASA Coincidencias Libre;0.190532
BBA0030;MSASA Identidad Libre;0.057925
BBA0030;MSASA Identidad Estricto;0.055475
BBA0030;MSASA Local Libre;0.286302
BBA0065;MUSCLE;0.822006
BBA0065;MAFFT;0.814270
BBA0065;TCOFFE;0.809057
BBA0065;Clustal;0.794041
BBA0065;KAlign;0.774565
BBA0065;MSASA Similitud Blosum62 Libre;0.045674
BBA0065;MSASA Similitud Blosum62 Estricto;0.045483
BBA0065;MSASA Global Libre;0.045305
BBA0065;MSASA Similitud PAM250 Estricto;0.045005
BBA0065;MSASA Similitud Gonnet92 Libre;0.043476
BBA0065;MSASA Similitud Gonnet92 Estricto;0.042798
BBA0065;MSASA Similitud PAM250 Libre;0.041655
BBA0065;MSASA Local Estricto;0.033671
BBA0065;MSASA Global Estricto;0.033199
BBA0065;MSASA Coincidencias Estricto;0.021629
BBA0065;MSASA Coincidencias Libre;0.020526
BBA0065;MSASA Identidad Libre;0.015323
BBA0065;MSASA Identidad Estricto;0.014853
BBA0065;MSASA Local Libre;0.044180
BBA0117;MUSCLE;0.825744
BBA0117;Clustal;0.807104
BBA0117;TCOFFE;0.800397
BBA0117;MAFFT;0.783769
BBA0117;KAlign;0.730688
BBA0117;MSASA Similitud Gonnet92 Estricto;0.237514
BBA0117;MSASA Global Libre;0.233606
BBA0117;MSASA Similitud Gonnet92 Libre;0.231653
BBA0117;MSASA Similitud PAM250 Estricto;0.223820
BBA0117;MSASA Coincidencias Estricto;0.145380
BBA0117;MSASA Similitud PAM250 Libre;0.112337
BBA0117;MSASA Similitud Blosum62 Libre;0.109768
BBA0117;MSASA Similitud Blosum62 Estricto;0.108270
BBA0117;MSASA Local Estricto;0.099063
BBA0117;MSASA Global Estricto;0.097678
BBA0117;MSASA Identidad Libre;0.072166
BBA0117;MSASA Coincidencias Libre;0.060362
BBA0117;MSASA Identidad Estricto;0.051067
BBA0117;MSASA Local Libre;0.094460
BBA0142;Clustal;0.915546
BBA0142;KAlign;0.914347
BBA0142;MAFFT;0.913490
BBA0142;TCOFFE;0.908437
BBA0142;MUSCLE;0.889936
BBA0142;MSASA Similitud PAM250 Estricto;0.596745
BBA0142;MSASA Similitud Blosum62 Libre;0.587152
BBA0142;MSASA Similitud Gonnet92 Estricto;0.583298
BBA0142;MSASA Similitud PAM250 Libre;0.569764
BBA0142;MSASA Local Estricto;0.528737
BBA0142;MSASA Global Libre;0.501071
BBA0142;MSASA Similitud Blosum62 Estricto;0.486510
BBA0142;MSASA Global Estricto;0.417730
BBA0142;MSASA Similitud Gonnet92 Libre;0.397859
BBA0142;MSASA Coincidencias Estricto;0.186552
BBA0142;MSASA Coincidencias Libre;0.151949
BBA0142;MSASA Identidad Libre;0.112805
BBA0142;MSASA Identidad Estricto;0.106467
BBA0142;MSASA Local Libre;0.521884
BBA0192;Clustal;0.824113
BBA0192;KAlign;0.823591
BBA0192;MUSCLE;0.802975
BBA0192;TCOFFE;0.800887
BBA0192;MAFFT;0.797495
BBA0192;MSASA Similitud PAM250 Libre;0.190762
BBA0192;MSASA Coincidencias Estricto;0.110647
BBA0192;MSASA Global Estricto;0.051670
BBA0192;MSASA Local Estricto;0.024791
BBA0192;MSASA Coincidencias Libre;0.014614
BBA0192;MSASA Global Libre;0.012265
BBA0192;MSASA Identidad Estricto;0.011221
BBA0192;MSASA Identidad Libre;0.007046
BBA0192;MSASA Similitud Gonnet92 Libre;0.005741
BBA0192;MSASA Similitud Gonnet92 Estricto;0.005741
BBA0192;MSASA Similitud Blosum62 Libre;0.003914
BBA0192;MSASA Similitud Blosum62 Estricto;0.003653
BBA0192;MSASA Similitud PAM250 Estricto;0.000000
BBA0192;MSASA Local Libre;0.017745
BBS12014;MUSCLE;1.000000
BBS12014;MAFFT;1.000000
BBS12014;TCOFFE;1.000000
BBS12014;KAlign;1.000000
BBS12014;Clustal;0.995465
BBS12014;MSASA Similitud Blosum62 Estricto;0.931973
BBS12014;MSASA Similitud Gonnet92 Libre;0.931973
BBS12014;MSASA Similitud Gonnet92 Estricto;0.931973
BBS12014;MSASA Similitud PAM250 Estricto;0.845805
BBS12014;MSASA Similitud Blosum62 Libre;0.841270
BBS12014;MSASA Similitud PAM250 Libre;0.836735
BBS12014;MSASA Local Estricto;0.769841
BBS12014;MSASA Global Estricto;0.722222
BBS12014;MSASA Global Libre;0.587868
BBS12014;MSASA Coincidencias Estricto;0.308957
BBS12014;MSASA Coincidencias Libre;0.210884
BBS12014;MSASA Identidad Estricto;0.174036
BBS12014;MSASA Identidad Libre;0.082200
BBS12014;MSASA Local Libre;0.760771
BBS20002;Clustal;0.908806
BBS20002;MAFFT;0.902919
BBS20002;MUSCLE;0.898033
BBS20002;TCOFFE;0.882876
BBS20002;KAlign;0.869848
BBS20002;MSASA Similitud Gonnet92 Libre;0.663786
BBS20002;MSASA Similitud Blosum62 Estricto;0.651635
BBS20002;MSASA Similitud Blosum62 Libre;0.643868
BBS20002;MSASA Local Estricto;0.613804
BBS20002;MSASA Similitud Gonnet92 Estricto;0.566704
BBS20002;MSASA Similitud PAM250 Libre;0.535763
BBS20002;MSASA Global Estricto;0.524364
BBS20002;MSASA Similitud PAM250 Estricto;0.519604
BBS20002;MSASA Global Libre;0.518351
BBS20002;MSASA Coincidencias Estricto;0.291369
BBS20002;MSASA Coincidencias Libre;0.212451
BBS20002;MSASA Identidad Libre;0.092071
BBS20002;MSASA Identidad Estricto;0.086809
BBS20002;MSASA Local Libre;0.495177
BBS20007;MUSCLE;0.805245
BBS20007;MAFFT;0.796680
BBS20007;TCOFFE;0.791916
BBS20007;KAlign;0.767280
BBS20007;Clustal;0.765809
BBS20007;MSASA Similitud PAM250 Estricto;0.272484
BBS20007;MSASA Local Estricto;0.230602
BBS20007;MSASA Similitud Gonnet92 Estricto;0.142255
BBS20007;MSASA Similitud PAM250 Libre;0.132664
BBS20007;MSASA Similitud Gonnet92 Libre;0.126555
BBS20007;MSASA Similitud Blosum62 Libre;0.123040
BBS20007;MSASA Similitud Blosum62 Estricto;0.115915
BBS20007;MSASA Global Estricto;0.101327
BBS20007;MSASA Global Libre;0.099897
BBS20007;MSASA Coincidencias Estricto;0.087574
BBS20007;MSASA Coincidencias Libre;0.066136
BBS20007;MSASA Identidad Libre;0.049165
BBS20007;MSASA Identidad Estricto;0.043120
BBS20007;MSASA Local Libre;0.178928
BBS20009;MUSCLE;0.966470
BBS20009;TCOFFE;0.965124
BBS20009;Clustal;0.962209
BBS20009;KAlign;0.959212
BBS20009;MAFFT;0.951801
BBS20009;MSASA Similitud Gonnet92 Libre;0.706358
BBS20009;MSASA Global Estricto;0.704485
BBS20009;MSASA Similitud Gonnet92 Estricto;0.689781
BBS20009;MSASA Similitud PAM250 Libre;0.667092
BBS20009;MSASA Similitud Blosum62 Estricto;0.663381
BBS20009;MSASA Global Libre;0.353852
BBS20009;MSASA Coincidencias Estricto;0.350656
BBS20009;MSASA Similitud PAM250 Estricto;0.325767
BBS20009;MSASA Similitud Blosum62 Libre;0.320873
BBS20009;MSASA Local Estricto;0.275835
BBS20009;MSASA Coincidencias Libre;0.146809
BBS20009;MSASA Identidad Estricto;0.083063
BBS20009;MSASA Identidad Libre;0.081377
BBS20009;MSASA Local Libre;0.312912
BBS30001;MAFFT;0.874342
BBS30001;Clustal;0.872335
BBS30001;MUSCLE;0.869704
BBS30001;TCOFFE;0.867090
BBS30001;KAlign;0.865158
BBS30001;MSASA Global Estricto;0.432586
BBS30001;MSASA Similitud PAM250 Estricto;0.428673
BBS30001;MSASA Global Libre;0.334708
BBS30001;MSASA Similitud Blosum62 Estricto;0.332612
BBS30001;MSASA Similitud PAM250 Libre;0.331899
BBS30001;MSASA Similitud Blosum62 Libre;0.326572
BBS30001;MSASA Similitud Gonnet92 Libre;0.321051
BBS30001;MSASA Similitud Gonnet92 Estricto;0.313756
BBS30001;MSASA Local Estricto;0.286886
BBS30001;MSASA Coincidencias Estricto;0.188506
BBS30001;MSASA Coincidencias Libre;0.183356
BBS30001;MSASA Identidad Libre;0.065256
BBS30001;MSASA Identidad Estricto;0.051884
BBS30001;MSASA Local Libre;0.431160
BBS30003;MUSCLE;0.653171
BBS30003;TCOFFE;0.643967
BBS30003;MAFFT;0.612653
BBS30003;KAlign;0.577945
BBS30003;Clustal;0.559501
BBS30003;MSASA Similitud Gonnet92 Libre;0.141259
BBS30003;MSASA Similitud Blosum62 Estricto;0.139970
BBS30003;MSASA Similitud Blosum62 Libre;0.139265
BBS30003;MSASA Similitud PAM250 Libre;0.139118
BBS30003;MSASA Local Estricto;0.130352
BBS30003;MSASA Similitud Gonnet92 Estricto;0.118935
BBS30003;MSASA Similitud PAM250 Estricto;0.118400
BBS30003;MSASA Global Estricto;0.118257
BBS30003;MSASA Global Libre;0.117655
BBS30003;MSASA Coincidencias Estricto;0.076545
BBS30003;MSASA Coincidencias Libre;0.054879
BBS30003;MSASA Identidad Estricto;0.036027
BBS30003;MSASA Identidad Libre;0.033738
BBS30003;MSASA Local Libre;0.135955
BBS30015;TCOFFE;0.847631
BBS30015;MUSCLE;0.845834
BBS30015;Clustal;0.740336
BBS30015;KAlign;0.697469
BBS30015;MAFFT;0.662252
BBS30015;MSASA Similitud PAM250 Libre;0.384106
BBS30015;MSASA Similitud Blosum62 Libre;0.382823
BBS30015;MSASA Similitud Blosum62 Estricto;0.378870
BBS30015;MSASA Similitud Gonnet92 Estricto;0.367678
BBS30015;MSASA Local Estricto;0.342779
BBS30015;MSASA Similitud PAM250 Estricto;0.319832
BBS30015;MSASA Similitud Gonnet92 Libre;0.317521
BBS30015;MSASA Global Estricto;0.314903
BBS30015;MSASA Global Libre;0.285487
BBS30015;MSASA Coincidencias Estricto;0.162277
BBS30015;MSASA Coincidencias Libre;0.061553
BBS30015;MSASA Identidad Libre;0.049849
BBS30015;MSASA Identidad Estricto;0.039786
BBS30015;MSASA Local Libre;0.212331
BBS30016;MUSCLE;0.897875
BBS30016;MAFFT;0.884002
BBS30016;TCOFFE;0.724826
BBS30016;Clustal;0.717993
BBS30016;KAlign;0.635093
BBS30016;MSASA Similitud Gonnet92 Libre;0.575860
BBS30016;MSASA Similitud PAM250 Libre;0.573652
BBS30016;MSASA Similitud Blosum62 Libre;0.566196
BBS30016;MSASA Similitud PAM250 Estricto;0.531695
BBS30016;MSASA Similitud Gonnet92 Estricto;0.524706
BBS30016;MSASA Global Libre;0.509223
BBS30016;MSASA Global Estricto;0.505196
BBS30016;MSASA Similitud Blosum62 Estricto;0.485919
BBS30016;MSASA Local Estricto;0.472566
BBS30016;MSASA Coincidencias Libre;0.438870
BBS30016;MSASA Coincidencias Estricto;0.368726
BBS30016;MSASA Identidad Estricto;0.116232
BBS30016;MSASA Identidad Libre;0.086849
BBS30016;MSASA Local Libre;0.543542
BBS30017;MUSCLE;0.736176
BBS30017;TCOFFE;0.731534
BBS30017;MAFFT;0.678720
BBS30017;Clustal;0.616378
BBS30017;KAlign;0.573418
BBS30017;MSASA Similitud Gonnet92 Libre;0.205147
BBS30017;MSASA Similitud PAM250 Estricto;0.173060
BBS30017;MSASA Similitud PAM250 Libre;0.160762
BBS30017;MSASA Similitud Blosum62 Libre;0.158645
BBS30017;MSASA Similitud Gonnet92 Estricto;0.155509
BBS30017;MSASA Similitud Blosum62 Estricto;0.149849
BBS30017;MSASA Global Estricto;0.099316
BBS30017;MSASA Global Libre;0.093941
BBS30017;MSASA Local Estricto;0.075006
BBS30017;MSASA Coincidencias Estricto;0.054605
BBS30017;MSASA Identidad Estricto;0.044181
BBS30017;MSASA Coincidencias Libre;0.036566
BBS30017;MSASA Identidad Libre;0.036485
BBS30017;MSASA Local Libre;0.091661
BBS50005;MUSCLE;0.943200
BBS50005;TCOFFE;0.943027
BBS50005;MAFFT;0.929583
BBS50005;Clustal;0.910064
BBS50005;KAlign;0.905775
BBS50005;MSASA Global Estricto;0.193977
BBS50005;MSASA Similitud PAM250 Libre;0.192768
BBS50005;MSASA Similitud PAM250 Estricto;0.182145
BBS50005;MSASA Similitud Gonnet92 Libre;0.163116
BBS50005;MSASA Similitud Gonnet92 Estricto;0.160957
BBS50005;MSASA Local Estricto;0.105337
BBS50005;MSASA Similitud Blosum62 Estricto;0.070820
BBS50005;MSASA Similitud Blosum62 Libre;0.067336
BBS50005;MSASA Global Libre;0.055389
BBS50005;MSASA Coincidencias Libre;0.041743
BBS50005;MSASA Identidad Libre;0.037799
BBS50005;MSASA Identidad Estricto;0.027781
BBS50005;MSASA Coincidencias Estricto;0.018799
BBS50005;MSASA Local Libre;0.102718
BOX122;MUSCLE;0.762914
BOX122;TCOFFE;0.754686
BOX122;Clustal;0.733084
BOX122;MAFFT;0.728942
BOX122;KAlign;0.727589
BOX122;MSASA Similitud PAM250 Libre;0.125504
BOX122;MSASA Similitud PAM250 Estricto;0.116972
BOX122;MSASA Local Estricto;0.115589
BOX122;MSASA Global Libre;0.111590
BOX122;MSASA Coincidencias Estricto;0.075304
BOX122;MSASA Similitud Gonnet92 Libre;0.072423
BOX122;MSASA Similitud Blosum62 Libre;0.072154
BOX122;MSASA Global Estricto;0.069889
BOX122;MSASA Similitud Gonnet92 Estricto;0.069103
BOX122;MSASA Similitud Blosum62 Estricto;0.068469
BOX122;MSASA Identidad Libre;0.045770
BOX122;MSASA Identidad Estricto;0.044573
BOX122;MSASA Coincidencias Libre;0.041548
BOX122;MSASA Local Libre;0.118737
BOX132;MUSCLE;0.852800
BOX132;MAFFT;0.840263
BOX132;TCOFFE;0.837443
BOX132;Clustal;0.806623
BOX132;KAlign;0.795863
BOX132;MSASA Similitud Blosum62 Estricto;0.304430
BOX132;MSASA Similitud Gonnet92 Estricto;0.263686
BOX132;MSASA Similitud PAM250 Libre;0.229628
BOX132;MSASA Similitud PAM250 Estricto;0.204242
BOX132;MSASA Similitud Gonnet92 Libre;0.195884
BOX132;MSASA Global Libre;0.142081
BOX132;MSASA Similitud Blosum62 Libre;0.135395
BOX132;MSASA Global Estricto;0.102486
BOX132;MSASA Local Estricto;0.101024
BOX132;MSASA Coincidencias Estricto;0.081070
BOX132;MSASA Coincidencias Libre;0.064563
BOX132;MSASA Identidad Libre;0.036147
BOX132;MSASA Identidad Estricto;0.026327
BOX132;MSASA Local Libre;0.152215
BOX212;MUSCLE;0.336735
BOX212;Clustal;0.330771
BOX212;TCOFFE;0.321892
BOX212;MAFFT;0.318844
BOX212;KAlign;0.285317
BOX212;MSASA Similitud Gonnet92 Libre;0.022131
BOX212;MSASA Similitud Gonnet92 Estricto;0.010734
BOX212;MSASA Coincidencias Estricto;0.010602
BOX212;MSASA Coincidencias Libre;0.009011
BOX212;MSASA Global Estricto;0.007156
BOX212;MSASA Similitud PAM250 Libre;0.006891
BOX212;MSASA Similitud Blosum62 Estricto;0.005698
BOX212;MSASA Similitud Blosum62 Libre;0.004241
BOX212;MSASA Global Libre;0.003578
BOX212;MSASA Identidad Libre;0.003313
BOX212;MSASA Local Estricto;0.003180
BOX212;MSASA Identidad Estricto;0.001988
BOX212;MSASA Similitud PAM250 Estricto;0.000000
BOX212;MSASA Local Libre;0.006228
"""

# Use StringIO to simulate reading from a file
data_io = StringIO(data)

# Create the DataFrame
df_scores = pd.read_csv(data_io, sep=';')  # Specify the separator as semicolon

# Create DataFrame from additional data
df_additional = pd.DataFrame({
    'Input Alignment': sequences,
    'Dimension Length': dimension_length,
    'Dimension Quantity': dimension_quantity,
    'Group Length': group_length,
    'Group Quantity': group_quantity,
    'Study Case': study_case
})

df_merged = pd.merge(df_scores, df_additional, left_on="Caso de Estudio", right_on='Input Alignment', how='left')

# Assume df_merged is the merged DataFrame containing all relevant data
softwares = [
    'MSASA Identidad Libre',
    'MSASA Identidad Estricto',
    'MSASA Coincidencias Libre',
    'MSASA Coincidencias Estricto',
    'MSASA Similitud Blosum62 Libre',
    'MSASA Similitud Blosum62 Estricto',
    'MSASA Similitud PAM250 Libre',
    'MSASA Similitud PAM250 Estricto',
    'MSASA Similitud Gonnet92 Libre',
    'MSASA Similitud Gonnet92 Estricto',
    'MSASA Global Libre',
    'MSASA Global Estricto',
    'MSASA Local Libre',
    'MSASA Local Estricto',
]

# Create a figure and axis objects for subplots
fig, axs = plt.subplots(nrows=7, ncols=2, figsize=(14, 28), constrained_layout=True)
axs = axs.flatten()  # Flatten the array to index it more easily

# Define a color map and norm for the new category
cmap = mcolors.ListedColormap(['blue', 'red'])  # Blue for BC, Red for AC
norm = mcolors.BoundaryNorm([0, 100, np.max(df_merged['Dimension Quantity'])], cmap.N)

# Loop through each software and plot data
for i, software in enumerate(softwares):
    # Filter data for the specific software
    data = df_merged[df_merged['Software'] == software]

    # # Create scatter plot
    # scatter = axs[i].scatter(data['Dimension Length'], data['Overlap Score'], c=data['Dimension Quantity'], cmap='viridis', label=f'{software} by Dimension Quantity')

    # # Add a color bar
    # cbar = fig.colorbar(scatter, ax=axs[i])
    # cbar.set_label('Cantidad de secuencias')
    # Categorize 'Dimension Quantity' into 'BC' and 'AC'
    data['Category'] = np.where(data['Dimension Quantity'] < 100, 'BC', 'AC')

    # Create scatter plot
    scatter = axs[i].scatter(data['Dimension Length'], data['Overlap Score'], c=data['Dimension Quantity'], cmap=cmap, norm=norm, label=f'{software} by Category')

    # Add a color bar
    cbar = fig.colorbar(scatter, ax=axs[i], ticks=[50, 300])
    cbar.set_label('Cantidad de secuencias')
    cbar.set_ticklabels(['BC', 'AC'])

    # Setting the title and labels
    axs[i].set_title(software)
    axs[i].set_xlabel('Longitud de las secuencias')
    axs[i].set_ylabel('MUMSA Overlap Score')

    # Highlight the groups in zones
    axs[i].axvspan(0, 100, color='red', alpha=0.25, label='CL: 0-100')
    axs[i].axvspan(100, 400, color='blue', alpha=0.25, label='ML: 100-400')
    axs[i].axvspan(400, 1100, color='green', alpha=0.25, label='LL: 400+')

    # Set the y-axis limits to be the same for all subplots
    axs[i].set_ylim(0.0, 1.0)

# Adjust layout to not overlap plots
# plt.tight_layout()

# Show plot
# plt.show()
plt.savefig("figura_actualizada.png", dpi=300)

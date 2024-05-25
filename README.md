
<a target="_blank" href="https://colab.research.google.com/github/https://colab.research.google.com/github/kangmg/rmsd_ipython/blob/master/notebooks/usage_tutorials.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

## rmsd_ipython

`rmsd_ipython` is a modified version of [rmsd](https://github.com/charnley/rmsd). `rmsd_ipython.py` enables to use rmsd in Ipython notebook enviroments but with restricted fuctionalities.

Please refer to the [rmsd](https://github.com/charnley/rmsd) for more details about rmsd.

## Usage
```python
from rmsd_ipython import voting_RMSD

xyz1 = """6

Cl    0.000315113091   -4.851256508727    0.000000000000
Si    0.000008122666   -0.081479435428    0.000000000000
H   -1.409996918298   -0.481349820968    0.000000000000
H    0.705178765477   -0.481189816756   -1.221101893621
H    0.705178765477   -0.481189816756    1.221101893621
Br   -0.000147111879    2.196922164251    0.000000000000"""

xyz2 = """6

Cl    0.000308701706   -4.739114565961    0.000000000000
Si    0.000011888584   -0.119054586111    0.000000000000
H   -1.411266251779   -0.514082066748    0.000000000000
H    0.705205731485   -0.514578303328   -1.222045485239
H    0.705205731485   -0.514578303328    1.222045485239
Br   -0.000130084759    2.161823224660    0.000000000000"""

RMSD_dic = voting_RMSD(xyz1, xyz2, round_digit=7)

RMSD_with_H = RMSD_dic["with_H"]
RMSD_without_H = RMSD_dic["without_H"]

print(f"RMSD_with_H     : {RMSD_with_H}")
print(f"RMSD_without_H  : {RMSD_without_H}")
```


Here there is the new python code relative to the paper "Estimating the vertical ionization potential of single-stranded DNA molecules" (M.Rooman and F. Pucci, submitted).

## Install as a library

```bash
pip install "git+ssh://git@github.com/Benoitdw/vIPer.git"
```

OR

```bash
pip install "git+https://git@github.com/Benoitdw/vIPer.git"
```

## RUN

```python
from viper import Viper

# All these call will return 7.5826 as is it single strand
Viper().compute_score("GAC")
Viper().compute_score("GAC", double_strand=False)
Viper().compute_single_strand_score("GAC")


# All these call will return 7.5826 as is it double strand
Viper().compute_score("GAC")
Viper().compute_score("GAC", double_strand=True)
Viper().compute_double_strand_score("GAC")
```

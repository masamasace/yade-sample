############## 基本的なモジュールのインポート ##############
from __future__ import print_function
from yade import pack, plot, O, Matrix3, utils
from yade.wrapper import FrictMat, ForceResetter, \
    InsertionSortCollider, InteractionLoop, NewtonIntegrator, PyRunner
import numpy as np
import datetime
from pathlib import Path
import os
import sys
import json

# search .gz file from temp directory

root_dir = Path("temp")
gz_files = list(root_dir.glob("*.gz"))

if len(gz_files) == 0:
    raise FileNotFoundError("No .gz files found in the temp directory.")
elif len(gz_files) > 1:
    raise ValueError("Multiple .gz files found in the temp directory.")
else:
    gz_file = gz_files[0]

# load simulation data from .gz file
O.load(str(gz_file))

# check some parameters
print(O.cell.velGrad)
print(utils.getStress())
print(O.cell.trsf)


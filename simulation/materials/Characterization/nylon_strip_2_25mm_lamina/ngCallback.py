# Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.
#
# Modified on 11/15/19 by Connor McCann to log more details about progress

import json
import time
import warnings
import datetime
from typing import Any, Union, List, Dict
from pathlib import Path
import numpy as np
import nevergrad.optimization.base as base
from nevergrad.parametrization import parameter as p
import ngFunctions
import os 
import matplotlib.pyplot as plt

class ngLogger:
    """Logs parameter and run information throughout into a file during
    optimization.
    Parameters
    ----------
    filepath: str or pathlib.Path
        the path to dump data to
    Usage
    -----
    logger = ParametersLogger(filepath)
    optimizer.register_callback("tell",  logger)
    optimizer.minimize()
    list_of_dict_of_data = logger.load()
    Note
    ----
    - arrays are converted to lists
    - this class will eventually contain display methods
    """

    def __init__(self,parameterlist, filepath: Union[str, Path], delete_existing_file: bool = False,delete_folder=False) -> None:
        self._session = datetime.datetime.now().strftime("%y-%m-%d %H:%M")
        self._filepath = Path(filepath)
        self.delete_folder = delete_folder
        if self._filepath.exists() and delete_existing_file:
            self._filepath.unlink()  # missing_ok argument added in python 3.8
        self.parameterlist = parameterlist
        self.loss_list = []
        # If the result is not in the top 5, then delete result
        self.topN = 5
        self.WorkerNumber_list = []
        for i in range(5): self.WorkerNumber_list.append(0)

    def __call__(self, optimizer: base.Optimizer, candidate: p.Parameter, value: float) -> None:
        X = []
        for p in self.parameterlist:
            X.append(candidate.kwargs[p])
        Xp = np.hstack([X])
        WorkerNumber = hash(tuple(Xp))
        self.loss_list.append(value)
        self.loss_list = sorted(self.loss_list)
        print('WorkerNumber: '+str(WorkerNumber)+' completed, loss: '+str(value))
        if self.delete_folder and value != 5e+20:
            if self.loss_list.index(value) >=self.topN:      # Check if result is in the top 5, if not, delete
                owd = os.getcwd()
                WorkerFolder='Files_'+str(WorkerNumber)
                ngFunctions.CleanWorkerDirectory(WorkerFolder,owd,WorkerNumber)
            else:
                indx = self.loss_list.index(value)
                if self.WorkerNumber_list[indx] == 0:
                    self.WorkerNumber_list[indx] = WorkerNumber
                else:
                    owd = os.getcwd()
                    WorkerFolder='Files_'+str(self.WorkerNumber_list[indx])
                    ngFunctions.CleanWorkerDirectory(WorkerFolder,owd,self.WorkerNumber_list[indx])
                    self.WorkerNumber_list[indx]= WorkerNumber
        data = {"#num-ask": optimizer.num_ask,
                "#num-tell": optimizer.num_tell,
                "#loss": value,
                "#job_id": str(WorkerNumber),
                "#arg": X}

        #params = dict(candidate.kwargs/10)
        #params.update({f"#arg{k}": arg for k, arg in enumerate(candidate.args)})
        #data.update({k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in params.items()})
        try:  # avoid bugging as much as possible
            with self._filepath.open("a") as f:
                f.write(json.dumps(data) + "\n")
        except Exception:  # pylint: disable=broad-except
            warnings.warn("Failing to json data")

    def load(self) -> List[Dict[str, Any]]:
        """Loads data from the log file
        """
        data: List[Dict[str, Any]] = []
        if self._filepath.exists():
            with self._filepath.open("r") as f:
                for line in f.readlines():
                    data.append(json.loads(line))
        return data

    def load_flattened(self, max_list_elements: int = 24) -> List[Dict[str, Any]]:
        """Loads data from the log file, and splits lists (arrays) into multiple arguments
        Parameters
        ----------
        max_list_elements: int
            Maximum number of elements displayed from the array, each element is given a
            unique id of type list_name#i1_i2_...
        """
        data = self.load()
        flat_data: List[Dict[str, Any]] = []
        for element in data:
            list_keys = {key for key, val in element.items() if isinstance(val, list)}
            flat_data.append({key: val for key, val in element.items() if key not in list_keys})
            for key in list_keys:
                for k, (indices, value) in enumerate(np.ndenumerate(element[key])):
                    if k >= max_list_elements:
                        break
                    flat_data[-1][key + "#" + "_".join(str(i) for i in indices)] = value
        return flat_data
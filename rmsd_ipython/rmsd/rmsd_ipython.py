LICENSE = """
BSD 2-Clause License
--------------------

Copyright (c) 2024, Kang mingi <kangmg@korea.ac.kr>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
"""

PROJECT_PAGE  = "https://github.com/kangmg/rmsd_ipython"


from .rmsd import centroid, check_reflections , get_coordinates_xyz_from_string # utilities
from .rmsd import quaternion_rmsd, kabsch_rmsd, rmsd # rmsd
from .rmsd import reorder_hungarian, reorder_inertia_hungarian, reorder_brute, reorder_distance # reorder
import copy
import pandas as pd
import numpy as np
import warnings
import seaborn as sns
import itertools
import matplotlib.pyplot as plt


# available methods option
AVAILABLE_ROTATION_METHODS = ['kabsch', 'quaternion', None]
AVAILABLE_REORDER_METHODS = ['hungarian', 'inertia-hungarian', 'brute', 'distance', None]


rmsd_options = lambda : print("""
____________.____________________________________________________________________________\n""",
" \033[0;34mParameters\033[0m ","|","                            \033[0;34mSupported methods\033[0m\n"
"""____________.____________________________________________________________________________
  rotation  |                      None, 'kabsch', 'quaternion'
____________.____________________________________________________________________________
  reorder   |          None, 'hungarian', 'inertia-hungarian', 'brute', 'distance'
____________.____________________________________________________________________________

""", sep="")


def get_rmsd(P_str: str, Q_str: str, rotation_method: str = None, 
             reorder_method: str = None, ignore_hydrogen: bool = False)->float:
    '''
    Description
    -----------

    Parameters
    ----------
      - P_str(str) : xyz format string
      - Q_str(str) : xyz format string
      - rotation_method:(str) : Rotation method. Default is None.
      - reorder_method(str) : Reorder method. Default is None.
      - ignore_hydrogen(bool) : whether exclude hydrogen atoms

    Returns
    -------
      - rmsd(float)

    Usage
    -----
    >>> H2_pre = """2
    ...
    ...H 0 0 1
    ...H 0 0 0"""

    >>> H2_opt = """2
    ...
    ...H 0 0 0.74
    ...H 0 0 0.00"""

    >>> get_rmsd(H2_pre, H2_opt, "kabsch", "hungarian")
    
    0.13
    '''

    P_all_atoms, P_all = get_coordinates_xyz_from_string(P_str, return_atoms_as_int=True)
    Q_all_atoms, Q_all = get_coordinates_xyz_from_string(Q_str, return_atoms_as_int=True)

    # Set local view
    p_view = None
    q_view = None

    if ignore_hydrogen:
        p_view = np.where(P_all_atoms != 1)
        q_view = np.where(Q_all_atoms != 1)

    # Set local view
    if p_view is None:
        P_coord = copy.deepcopy(P_all)
        Q_coord = copy.deepcopy(Q_all)
        P_atoms = copy.deepcopy(P_all_atoms)
        Q_atoms = copy.deepcopy(Q_all_atoms)

    else:
        P_coord = copy.deepcopy(P_all[p_view])
        Q_coord = copy.deepcopy(Q_all[q_view])
        P_atoms = copy.deepcopy(P_all_atoms[p_view])
        Q_atoms = copy.deepcopy(Q_all_atoms[q_view])

    P_size = P_atoms.shape[0]
    Q_size = Q_atoms.shape[0]

    # size check
    if not P_size == Q_size:
        raise ValueError(f"Structures not same size, first xyz : second xyz = {P_size} : {Q_size}")

    # Recenter to centroid
    P_cent = centroid(P_coord)
    Q_cent = centroid(Q_coord)
    P_coord -= P_cent
    Q_coord -= Q_cent

    # set rotation method
    if rotation_method == 'kabsch':
        rmsd_method = kabsch_rmsd
    elif rotation_method == 'quaternion':
        rmsd_method = quaternion_rmsd
    else:
        rmsd_method = rmsd

    # set reorder method
    if reorder_method == 'hungarian':
        reorder_method = reorder_hungarian
    elif reorder_method == 'inertia-hungarian':
        reorder_method = reorder_inertia_hungarian
    elif reorder_method == 'brute':
        reorder_method = reorder_brute
    elif reorder_method == 'distance':
        reorder_method = reorder_distance

    if reorder_method is not None:
        Q_review = reorder_method(P_atoms, Q_atoms, P_coord, Q_coord)
    else:
        Q_review = None

    # If there is a reorder, then apply before print
    if Q_review is not None:

        Q_atoms = Q_atoms[Q_review]
        Q_coord = Q_coord[Q_review]

        assert all(
            P_atoms == Q_atoms
        ), "error: Structure not aligned. Please submit bug report at http://github.com/charnley/rmsd"

    result_rmsd = rmsd_method(P_coord, Q_coord)

    return result_rmsd


def RMSD_TABLE(P_str:str, Q_str:str, exclude_None_option:bool=True, return_None:bool=True, round_digit:str|int="full"):
  '''
  Description
  -----------
  Generate RMSD table comparing two molecular geometry.

  Parameters
  ----------
    - P_str (str) : XYZ string representing the first molecular structure.
    - Q_str (str) : XYZ string representing the second molecular structure.
    - exclude_None_option (bool) : Whether to exclude None options in rotation and reorder methods.
    - return_None (bool) : Whether to return the RMSD tables or display them.
    - round_digit (str|int) : N for round(RMSD, N). 

  Returns
  -------
  - if `return_None=True`:
    - None
  - if `return_None=False`:
    - (include_Hydrogen, exclude_Hydrogen) : tuple[pd.io.formats.style.Styler, pd.io.formats.style.Styler]

  Usage
  -----

  >>> F2_pre = """2
  ...
  ...F 0 0 1
  ...F 0 0 0"""

  >>> F2_opt = """2
  ...
  ...F 0 0 0.74
  ...F 0 0 0.00"""
  
  >>> RMSD_TABLE(F2_opt, F2_pre, exclude_None_option=False, return_None=True)
  '''
  warnings.filterwarnings("ignore")

  # inner function for styling
  def RMSD_highlight(df: pd.DataFrame, RMSD: float, caption: str) -> pd.DataFrame:
    styled_df = df.style.apply(
      lambda row: ['opacity: 0.5' if row['RMSD'] != RMSD else None for _ in row], axis=1)
    styled_df = styled_df.set_caption(f"<span style='font-size:150%;'>{caption}</span>")
    return styled_df

  # include_None_option
  if exclude_None_option:
    ROTATION_METHODS = AVAILABLE_ROTATION_METHODS.copy()
    ROTATION_METHODS.remove(None)
    REORDER_METHODS = AVAILABLE_REORDER_METHODS.copy()
    REORDER_METHODS.remove(None)
  else:
    ROTATION_METHODS = AVAILABLE_ROTATION_METHODS
    REORDER_METHODS = AVAILABLE_REORDER_METHODS
  
  # define df framework
  include_Hydrogen = pd.DataFrame(columns=["Rotation", "Reorder", "ignore_H", "RMSD"])
  exclude_Hydrogen = pd.DataFrame(columns=["Rotation", "Reorder", "ignore_H", "RMSD"])
  
  # calculate rmsd
  for rotation_option in ROTATION_METHODS:
    for reorder_option in REORDER_METHODS:
      for ignore_H in [True, False]:
        RMSD = get_rmsd(P_str, Q_str, rotation_method=rotation_option, reorder_method = reorder_option, ignore_hydrogen=ignore_H)
        if round_digit == "full":
          pass
        elif (type(round_digit) == int) and (round_digit > 0):
          RMSD = round(RMSD, round_digit)
        else:
           raise ValueError(f"round_digit expected positive integer, but got {round_digit}")
        new_row = pd.DataFrame([{"Rotation": rotation_option, "Reorder": reorder_option, "ignore_H": ignore_H, "RMSD": RMSD}], columns=["Rotation", "Reorder", "ignore_H", "RMSD"])
        if ignore_H == True:
          exclude_Hydrogen = pd.concat([exclude_Hydrogen, new_row], ignore_index=True)
        elif ignore_H == False:
          include_Hydrogen = pd.concat([include_Hydrogen, new_row], ignore_index=True)

  # voting
  exld_H_RMSD = exclude_Hydrogen["RMSD"].value_counts().idxmax()
  icld_H_RMSD = include_Hydrogen["RMSD"].value_counts().idxmax()
  
  # highlightling
  include_Hydrogen = RMSD_highlight(include_Hydrogen, icld_H_RMSD, "RMSD Table (including Hydrogens)")
  exclude_Hydrogen = RMSD_highlight(exclude_Hydrogen, exld_H_RMSD, "RMSD Table (excluding Hydrogens)")

  if return_None:
    display(include_Hydrogen)
    display(exclude_Hydrogen)
  else:
    return include_Hydrogen, exclude_Hydrogen


def voting_RMSD(P_str:str, Q_str:str, exclude_None_option:bool=True, round_digit:str|int="full"):
  """
  Description
  -----------
  Calculating RMSD values with all available options and selecting the most frequently occurring one.

  Parameters
  ----------
    - P_str (str) : XYZ string representing the first molecular structure.
    - Q_str (str) : XYZ string representing the second molecular structure.
    - exclude_None_option (bool) : Whether to exclude None options in rotation and reorder methods.
    - round_digit (str|int) : N for round(RMSD, N). 

  Returns
  -------
    - voted RMSD(float) 
  """
  icld_df, ecld_df = RMSD_TABLE(P_str, Q_str, exclude_None_option, return_None=False, round_digit=round_digit)
  icld_H_voting_RMSD = icld_df.data["RMSD"].value_counts().idxmax()
  ecld_H_voting_RMSD = ecld_df.data["RMSD"].value_counts().idxmax()
  return {"with_H" : icld_H_voting_RMSD, "without_H" : ecld_H_voting_RMSD}


def rmsd_heatmap(xyz_set:dict, **kwds):
  """
  Description
  -----------
  Generates a heatmap of Root Mean Square Deviation (RMSD) values for a given set of molecular coordinates.
    
  Parameters:
  ----------
  - xyz_set (dict) : A dictionary where keys are label (e.g., method name) and values are xyz format string.
  - **kwds (dict) : plot keyword arguments
    - title (str) : title of the heatmap
    - cmap (str) : Colormap used for the heatmap
    - vmin (float) : colorbar min
    - vmax (float) : colorbar max
    - xylabel (str) : xlabel, ylabel
  """
  # parsing plot kwards
  title = kwds.get("title", "Molecular RMSD Heatmap")
  cmap = kwds.get("cmap", "Blues")
  vmin = kwds.get("vmin", None)
  vmax = kwds.get("vmax", None)
  xylabel = kwds.get("xylabel", None)
  
  # seperate xyz data and labels
  labels = xyz_set.keys()
  coordinates = xyz_set.values()
  
  # calculate RMSD and makes RMSD matrix
  xyz_pair = itertools.product(coordinates, repeat=2)
  RMSD_result = np.array(list(get_rmsd(xyz1, xyz2) for xyz1, xyz2 in xyz_pair))
  NumOfXYZ = len(labels)
  heatmap_matrix = RMSD_result.reshape(NumOfXYZ, NumOfXYZ)
  
  # plot heatmap
  mask = np.tril(np.ones_like(heatmap_matrix, dtype=bool))
  ax = sns.heatmap(heatmap_matrix, annot=True, linewidth=.5, cmap=cmap, mask=~mask, xticklabels=labels, yticklabels=labels, vmin=vmin, vmax=vmax)
  ax.set_xlabel(xylabel)
  ax.set_ylabel(xylabel)
  ax.set_title(title)
  plt.show()
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit.Chem import AllChem, MACCSkeys, Draw
from rdkit import Chem, DataStructs
from rdkit.Avalon import pyAvalonTools as fpAvalon


class tanimoto_similarity:
    """
    Calculate Tanimoto coefficient using different type of molecular fingerprints
    For non-coding users, please visit this website to use: https://huggingface.co/spaces/MedAILab/TanimotoSimilarities

    Parameters
    ------
    smiles_list : list
        list contain smiles string.
    fps: str
        type of molecular fingerprints using to calculate Tanimoto coefficient: ecfp2, ecfp4, ecfp6, maccs, avalon, fcfp2, fcfp4, fcfp6, rdk5, rdk6, rdk7
    id_list: str
        name of identity list
    Returns
    --------
    fig1: figure
        similarity matrix; triangle heatmap and square heatmap.
    fig2: figure
        pair of best similarity
    Example
    --------
    ```python
    import pandas as pd
    from similarities import tanimoto_similarity

    # 1. Load data
    smi = pd.read_csv('test_smiles.smi', header = None, sep="\t")
    smiles_list=smi.iloc[:,0].tolist() #convert dataframe to smiles list

    # 2. Call class tanimoto_similarity
    simi = tanimoto_similarity(smiles_list = smiles_list, fps ='ecfp2')

    # 3. Visualize triangle similarity matrix
    triangle_heatmap = simi.visualize_triangle()
    display(triangle_heatmap)

    # 4. Visualize two molecules having best similarity
    pair_visualize = simi.pair_best_similarity()
    display(pair_visualize)
    ```

      
    """
    def __init__(self,smiles_list:list, fps:str, id_list:list = None):
        self.smiles_list = smiles_list # identify smiles_list in this class
        mol_list = []

        # Check if smiles valid, if not raise error and break
        for smiles in smiles_list:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    raise ValueError('SMILES is not valid!')
                else:
                    mol_list.append(mol)
            except ValueError as e:
                print(e)
                break
        self.mol_list = mol_list # identify mol_list in this class
        self.fps = fps.lower() # make sure that the character is lowercase

        # Check lenghth of id_list and smiles_list must be the same
        if id_list is not None:
          if len(id_list) != len(self.smiles_list):
            print('Error: ID list and SMILES list have different lengths.')
          else:
            self.id_list = id_list
        else:
            self.id_list = list(range(len(smiles_list)))
        sns.set(font_scale=1)
        sns.set(style ='darkgrid')
        self.table=pd.DataFrame() # create the blank table to gather similarity
        # calculate molecular fingerprints, it would be depend on self.fps attribute above
        self.calculate()

    def calculate(self):
        """calculate fingerprints and similarity matrix
        
        Parameters:
        --------
        self.fps: str
            type of molecular fingerprints
        Returns:
        --------
        self.triangle: numpy
            similarity matrix, using triangle to reduce computational cost
        """

        # 1. Calculate molecular fingerprint
        if self.fps == "ecfp4":
          fps= [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits = 2048) for mol in self.mol_list]
        if self.fps == "ecfp2":
          fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 1, nBits = 1024) for mol in self.mol_list]
        if self.fps == "ecfp6":
          fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits = 4096) for mol in self.mol_list]
        if self.fps == "maccs":
          fps = [MACCSkeys.GenMACCSKeys(mol) for mol in self.mol_list]
        if self.fps =="avalon":
          fps = [fpAvalon.GetAvalonFP(mol, 1024)  for mol in self.mol_list]
        if self.fps == "fcfp2":
          fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 1, nBits=1024, useFeatures=True) for mol in self.mol_list]
        if self.fps == "fcfp4":
          fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048, useFeatures=True) for mol in self.mol_list]
        if self.fps == "fcfp6":
          fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=4096, useFeatures=True) for mol in self.mol_list]
        if self.fps == "rdk5":
          fps = [Chem.RDKFingerprint(mol, maxPath=5, fpSize=2048, nBitsPerHash=2) for mol in self.mol_list]
        if self.fps =="rdk6":
          fps = [Chem.RDKFingerprint(mol, maxPath=6, fpSize=2048, nBitsPerHash=2) for mol in self.mol_list]
        if self.fps == "rdk7":
          fps = [Chem.RDKFingerprint(mol, maxPath=7, fpSize=2048, nBitsPerHash=2) for mol in self.mol_list]

        # 2. Calculate similarity matrix
        for index, i in enumerate(fps):
            for jndex, j in enumerate(fps):
                similarity=DataStructs.FingerprintSimilarity(i,j, metric=DataStructs.TanimotoSimilarity)
                self.table.loc[self.id_list[index],self.id_list[jndex]]=similarity
        corr = self.table.corr()
        self.mask = np.tril(np.ones_like(corr, dtype=bool))
        self.triangle = np.tril(self.table.to_numpy(),-1)
        self.max_sim = self.triangle.max()

    def visualize_triangle(self,plot_setting:dict = None,save_to:str = None):
        # generating the plot
        if plot_setting is not None:
            fig,ax = plt.subplots(**plot_setting)
        else:
            fig,ax = plt.subplots(figsize=(25,25))
        ax.set_title('Heatmap of Tanimoto Similarities', fontsize=24) 

        ax = sns.heatmap(self.table, annot = True, annot_kws={"fontsize":10}, center=0,
                    square=True,  linewidths=.7, cbar_kws={"shrink": .5}, mask = self.mask)

        if save_to is not None:
            plt.savefig(save_to,dpi=250)
        plt.close()
        return fig


    def visualize_square(self,plot_setting:dict = None,save_to:str = None):
        # generating the plot
        if plot_setting is not None:
            fig,ax = plt.subplots(**plot_setting)
        else:
            fig,ax = plt.subplots(figsize=(25,25))
        ax.set_title('Heatmap of Tanimoto Similarities', fontsize=24) 

        ax = sns.heatmap(self.table, annot = True, annot_kws={"fontsize":10}, center=0,
                    square=True,  linewidths=.7, cbar_kws={"shrink": .5}, vmin = 0, vmax = 1)

        if save_to is not None:
            plt.savefig(save_to,dpi=250)
        plt.close()
        return fig
      
    def max_sim_idx(self):
      
      id_smiles_pairs = []
      r,c = np.where(self.triangle == self.max_sim)
      for i,j in zip(r,c):
        id_smiles_pairs.append({self.id_list[i]:self.smiles_list[i], self.id_list[j]:self.smiles_list[j]})
      return id_smiles_pairs

    def pair_best_similarity(self):
      """Identify pair of molecules having best similarity
     
        Parameters:
        --------
        self.pair: list
            list contains smiles and id of two molecules having best similarity
        Returns:
        --------
        fig: figure
            pair of molecules having best similarity and their ID
      """
      pairs = self.max_sim_idx()
      for pair in pairs:
        id_ls = list(pair.keys())
        sm_ls = list(pair.values())
        mol_ls = [Chem.MolFromSmiles(sm_ls[i]) for i in range(len(sm_ls))]
      print('Best similarity values=', np.round(self.max_sim,3))

      return Draw.MolsToGridImage(
          [mol_ls[i] for i in range(len(mol_ls))],
          legends=[str(i) for i in id_ls],
          molsPerRow=2, subImgSize=(300, 300))
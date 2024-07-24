from typing import Dict, List, Optional

import numpy as np
from .scoring.peptriever import PeptrieverScorer
from modlamp.descriptors import PeptideDescriptor
from softnanotools.logger import Logger

logger = Logger(__name__)

class LossFn:
    """Loss function class to compute energy of a peptide

    Over time, more loss functions can be added to this class as a new method

    These methods must be called with the following naming convention:
    {model_name}_energy

    Parameters
    ----------
    loss_fn : Dict
        Dictionary containing loss function parameters

    consider look at points. 
    Gaussian around 0 
    """
    def __init__(
            self,  
            target_proteins: List[str], 
            non_targets: List[str],
            name_to_seq: Dict[str, str],
            concentrations: Dict[str, float],
            loss_fn: Optional[Dict] = None,
        ):

        if loss_fn is None:
            self.loss_fn = {
                #window either 5 or 7:
                "janin": {"window": 5},
                "peptriever": {"mid_point": 3.7, "cutoff": 4.0},
                #the weights need to be a formula that is dependent on concentration of the protein relative to other proteins 
                "target_weights": {"janin": 0.5, "peptriever": 0.5},
                "non_target_weights": {"janin": 0.5, "peptriever": 0.5},
                "alpha": {"adjustment factor": None},
            }
        else:
            self.loss_fn = loss_fn

        self.concentrations = concentrations
        self.name_to_seq = name_to_seq
        self.models = self.initialise_models()

        self.target_proteins = target_proteins
        self.non_targets = non_targets

        self.proteins = target_proteins + non_targets

        #debugging protein sequence insertion:
        print("Target proteins provided:", target_proteins)
        print("Non-target proteins provided:", non_targets)

    def initialise_models(self):
        """Add more if statements here for each unique model that needs to be initialised""" 
        models = {}
        if "peptriever" in self.loss_fn:
            models["peptriever"] = PeptrieverScorer()
        if "janin" in self.loss_fn:
            assert "window" in self.loss_fn["janin"], "Need to provide window parameter"
           
        return models

    #def peptriever_energy(self, peptide):
    def peptriever_energy(self, proteins: List[str] = None) -> float:
        """
        Compute Peptriever scores

        Parameters
        ----------
        peptide : str
            Peptide sequence to compute energy for (to be scored)
        """
        #1. calls the DNN model to score the peptide against all proteins [NO ISSUE]
        peptriever_scorer: PeptrieverScorer = self.models["peptriever"]
        #scores = list(peptriever_scorer.score_all_vs_all(self.proteins, [peptide], score_fn=["cosine_similarity"]).flatten())
        scores = list(peptriever_scorer.score_all_vs_all(self.proteins, [proteins], score_fn=["cosine_similarity"]).flatten())
        assert len(scores) == (len(self.target_proteins) + len(self.non_targets)), "Number of scores should be equal to number of proteins"

        #2. convert scores to energies by calling the class method below
        energies = self.peptriever_score_to_energy(scores)

        # #Outputting Peptriever Energies 
        # target_es: List[float] = [energies[i] for i in range(len(self.target_proteins))]
        # nontarget_es: List[float] = [energies[i] for i in range(len(self.non_targets), len(energies))]   
        # target_e = max(target_es)
        # nontarget_e = min(nontarget_es)
        # #GAP = OUTPUT
        # output = target_e - nontarget_e
        
        return energies

    def peptriever_score_to_energy(self, scores: List[float]) -> List[float]:
        """
        Equation to convert peptriever score to energy

        Notes
        -----
        energy = (score - mid_point_score)^2
        """
        for score in scores:
            assert score <= 1.0, AssertionError(f'{score} Max expected score is 1.0')

        alpha_adj = 5
        energies = [( 1.0 - score) ** 2 * alpha_adj for score in scores]

        #print(f"this is peptriever energy: {energies}")
        return energies 

    #def janin_score(self, peptide) -> float:
    def janin_score(self, proteins: List[str] = None) -> float:
        """Compute Janin hydrophobicity scores
        ----------
        Changed it to a loop to iterate through the proteins and store it in a list* (remove loop to revert back changes)
        """
        scores = []
        for _ in proteins:
            descriptor = PeptideDescriptor(seqs=[proteins], scalename='janin')
            descriptor.calculate_moment(modality="mean")
            value = descriptor.descriptor
            scores.append(value[0][0])

        #print(f"There are {len(scores)} scores in the list: {scores}")
        return scores
    
    #def janin_energy(self, peptide) -> float:
    def janin_energy(self, proteins: List[str] = None) -> float:
        """Compute Janin hydrophobicity scores and convert to energy"""
        scores = self.janin_score(proteins)
        energies = []
        #Convert this score to energy:
        mid_point = self.loss_fn["janin"]["window"]
        #print(f"mid_point is: {mid_point}")
        
        alpha_adj = 0.2
        #alpha_adj = 1.5
        target_janin = -0.30

        for score in scores:
            hydrophobicity = (score - target_janin) ** 2 / alpha_adj
            energies.append(hydrophobicity) 

        energies = [ (score - target_janin) ** 2 / alpha_adj for score in scores]

        #hydrophobicity = (score - mid_point) ** 2 / alpha_adj 
        #print (f"hydrophobicity equation output: {hydrophobicity}")

        return energies
    
    def apply_weights (self, concentrations: Dict[str, float], scores: Dict[str, float]) -> float:
        """Modular equation to compute weights for both targets and non_targets
        ---------
        Compute the weighted score for proteins.
        Parameters:
        ----------
        concentrations : Dict[str, float]
            Concentration of proteins in mol/L.
        scores : Dict[str, float]
            DNN scores of sequences.

        Returns:
        -------
        float
        Weighted score.
        """
        # label = list(concentrations.keys())
        # protein_to_seq = self.name_to_seq

        sum_Ci_Si = 0
        sum_Ci = 0

        for sequence, score in scores.items():
            # Ensure that the sequence exists in the protein-to-sequence mapping
            protein = None
            for prot, seq in self.name_to_seq.items():
                if seq == sequence:
                    protein = prot
                    break
            
            if protein and protein in concentrations:
                concentration = concentrations[protein]
                sum_Ci_Si += concentration * score
                sum_Ci += concentration

        assert sum_Ci != 0, "Sum of concentrations should not be zero"
        W_i = sum_Ci_Si / sum_Ci
        #print(f"this is the apply weights score: {W_i}")
        return W_i

    #def compute_energy(self, peptide: str) -> float:
    def compute_energy(self, proteins: List[str]) -> float:
        """Compute energy of the input peptide
        ------
        """      
        #Create empty dictionaries to store values for energy and weights for the E_final equation: 
        energy_components = {"target": {}, "non_target": {}}
        weighted_components = {"target": {}, "non_target": {}}

        #Obtains the energies of the peptide; inclusive of targets and non_targets 
        peptriever_energies = self.peptriever_energy(proteins)
        janin_energies = self.janin_energy(proteins) 

        #ORGANIZE DATA (Si data is obtained from here)
        target_peptriever_es = peptriever_energies[:len(self.target_proteins)]
        non_targets_peptriever_es = peptriever_energies[len(self.target_proteins):]
        target_janin = janin_energies[:len(self.target_proteins)]
        non_target_janin = janin_energies[len(self.target_proteins):]  

        #print(f"buckets of NT data - peptriever: {non_targets_peptriever_es} and janin: {non_target_janin}")
        
        #Peptriever #GAP ANALYSIS       
        target_es: List[float] = target_peptriever_es
        nontarget_es: List[float] = non_targets_peptriever_es
        target_e = max(target_es)
        nontarget_e = min(nontarget_es)

        energy_components["target"]["peptriever"] = target_e
        energy_components["target"]["janin"] = target_janin
        energy_components["non_target"]["peptriever"] = nontarget_e
        energy_components["non_target"]["janin"] = non_target_janin

        #target_weights = self.loss_fn.get('target_weights')

        #Unpack the scores for targets and non-targets:
        S_ntjanin: Dict[str,float] = {protein: score for protein, score in zip(self.non_targets, non_target_janin)}
        S_ntpeptriever: Dict[str, float] = {protein: score for protein, score in zip(self.non_targets, non_targets_peptriever_es)}
        S_targ_janin: Dict[str, float] = {protein: score for protein, score in zip(self.target_proteins, target_janin)}
        S_targ_peptriever: Dict[str, float] = {protein: score for protein, score in zip(self.target_proteins, target_peptriever_es)}
        
        #NEW APPROACH
        # 1. weighted_score_X: Dict[str, float] = self.apply_weights(weights=self.concentration, scores=S_nt<X>)
        # 2. weighted_components["target"][X] = sum(weighted_score_X.values())

        #1A. non_targets
        w_score_janin: float = self.apply_weights(self.concentrations, S_ntjanin)
        w_score_peptriever: float  = self.apply_weights(self.concentrations, S_ntpeptriever)
        
        #1B. targets
        w_targ_score_janin: float = self.apply_weights(self.concentrations, S_targ_janin)   
        w_targ_score_peptriever: float = self.apply_weights(self.concentrations, S_targ_peptriever) 
 
        #print(f"the weighted scores of X should be a dictionary: janin = {w_score_janin} and peptriever = {w_score_peptriever}")

        #2. weighted_components
        weighted_components["target"]["peptriever"] = w_targ_score_peptriever
        weighted_components["target"]["janin"] = w_targ_score_janin
        weighted_components["non_target"]["peptriever"] = w_score_peptriever
        weighted_components["non_target"]["janin"] = w_score_janin

        print(f"weighted components before minimization: {weighted_components}")

        final_energy = (
            sum(weighted_components["target"].values()) -
            sum(weighted_components["non_target"].values())
        )
        print(f"final energy: {final_energy}")

        return final_energy, energy_components


     # final_energy = sum(
        #     weighted_components["target.values()"] * weighted_components["peptriever"],
        #     weighted_components["target.values()"] * weighted_components["janin"]
        #     ) - sum(
        #     weighted_components["non_target.values()"] * weighted_components["peptriever"],
        #     weighted_components["non_target.values()"] * weighted_components["janin"]
        #     )

        # #Alpha adjustment if necessary:
        # alpha = self.loss_fn.get('alpha', 1.0)
        # final_energy *= alpha 
        #[energies[i] for i in range(len(self.target_proteins))]

     # #Store the energies in a list for further processing.
        # total_peptriever_energies = [peptriever_energies] if isinstance(peptriever_energies, float) else list(peptriever_energies)
        # total_janin_energies = [janin_energies] if isinstance(janin_energies, float) else list(janin_energies)

        # target_janin = total_janin_energies[:len(self.target_proteins)]
        # non_target_janin = total_janin_energies[len(self.target_proteins):]

        # #Use a loop to iterate for the len of target sequences and the loop to itrate thel en of dictionary basket of non targets. 
        # if isinstance(peptriever_energies, float) and isinstance(janin_energies, float):
        #     total_peptriever_data = [peptriever_energies]
        #     total_janin_data = [janin_energies]
        # else:
        #     print(f"beware of data type operations - debugging required")
        #     total_peptriever_data = list(peptriever_energies) if not isinstance(peptriever_energies, list) else peptriever_energies
        #     total_janin_data = list(janin_energies) if not isinstance(janin_energies, list) else janin_energies

        # #Approach 2: split total data into its respective components. 
        # target_peptriever = total_peptriever_data[:len(self.target_proteins)]
        # target_janin = total_janin_data[:len(self.target_proteins)]

        # non_target_start_index = len(self.target_proteins)
        # non_target_peptriever = total_peptriever_data[non_target_start_index:]
        # non_target_janin = total_janin_data[non_target_start_index:]

        # #concentrations is a dictionary. self.concentration.values() may be useful:
        # sum_Ci_Si_janin = sum(self.concentrations * S_ntjanin for _ in self.concentrations)
        # sum_Ci_janin = sum(self.concentrations)
        # W_nt_janin = ((sum_Ci_Si_janin) / sum_Ci_janin)

        # sum_Ci_Si_peptriever = sum(self.concentrations * S_ntpeptriever for _ in self.concentrations)
        # sum_Ci_peptriever = sum(self.concentrations)
        # W_nt_peptriever = ((sum_Ci_Si_peptriever) / sum_Ci_peptriever)

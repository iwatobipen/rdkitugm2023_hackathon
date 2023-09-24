import logging
import knime.extension as knext
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import pandas as pd
LOGGER = logging.getLogger(__name__)


@knext.node(name="Tautmer Enumerator", node_type=knext.NodeType.LEARNER, icon_path="demo.png", category="/")
@knext.input_table(name="SMLES Column", description="please provide vaild smiles!")
@knext.output_table(name="Output Data", description="List of all Tautomers")
class EnumerateTautomerNodeNode((knext.PythonNode)):
    """Enumerate all possible tautomers
    Input SMILES should be Vaild!!!!
    """
    def enumerate_tautomer(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        enum = rdMolStandardize.TautomerEnumerator()
        tautomers = enum.Enumerate(mol)
        uniq_smi = set([Chem.MolToSmiles(m) for m in tautomers])
        return [Chem.MolFromSmiles(smi) for smi in uniq_smi]

    mol_id_column = knext.ColumnParameter(label="Select the column of molecule ID", 
                                            description="Choose the column from the input table containing the MOLID",
                                            )
    molecule_column = knext.ColumnParameter(label="Select the column containing SMILES", 
                                            description="Choose the column from the input table containing the SMILES") 
    def configure(self, configure_context, input_schema_1):
    # def configure(self, configure_context, input_schema_1, input_schema_2):  ### Tutorial step 11: Uncomment to configure the new port (and comment out the previous configure header)
        return input_schema_1

        ### Tutorial step 12: Uncomment the following to adjust to the changes we do in this step in the execute method (and comment out the previous return statement)
        # return input_schema_1.append(knext.Column(knext.double(), "column2"))
        ### Tutorial step 13: Uncomment to set a warning for the configuration, which will be shown in the workflow
        # configure_context.set_warning("This is a warning during configuration")

 
    def execute(self, exec_context, input_1):
    # def execute(self, exec_context, input_1, input_2):  ### Tutorial step 11: Uncomment to accept the new port (and comment out the previous execute header)

        res_df = {"original_mol_id":[], "tautomer_idx":[], "original_smiles":[], "tautomer":[]}
        input_1_pandas = input_1.to_pandas()
        mol_ids = input_1_pandas[self.mol_id_column].to_list()
        for i, smi in enumerate(input_1_pandas[self.molecule_column].to_list()):
            res = self.enumerate_tautomer(smi)
            for j, tautomer in enumerate(res):
                res_df["original_mol_id"].append(mol_ids[i])
                res_df["original_smiles"].append(smi)
                Chem.SanitizeMol(tautomer)
                try:
                    Chem.Kekulize(tautomer)
                except:
                    pass
                res_df["tautomer_idx"].append(j)
                res_df["tautomer"].append(tautomer)
        output_df = pd.DataFrame(res_df)
        return knext.Table.from_pandas(output_df)
import numpy as np
import pandas as pd


class StanfordRNA():
    def __init__(self, labels: pd.DataFrame, sequences: pd.DataFrame):
        self.labels = labels
        self.sequences = sequences
        
        self.labels = self.get_xyz_dataframe(self.labels)
        
        self.labels['ohe'] = self.labels['resname'].apply(self.one_hot_encode)
        self.sequences['ohe'] = self.sequences['sequence'].apply(self.one_hot_encode)

    @staticmethod
    def get_xyz_dataframe(df):
        """
        Extracts and groups x, y, z coordinates by unique sequence ID from a DataFrame.

        Arguments:
            df(pd.DataFrame): A DataFrame containing an 'ID' column and at least three numerical columns 
                        representing x, y, and z coordinates.
            

        Returns:
            pd.DataFrame: A DataFrame with two columns:
                    - 'Unique_ID': The extracted unique sequence ID.
                    - 'Coordinates': A list of [x, y, z] coordinate triplets for each unique ID.
        """

        # Extract unique ID (first 4 characters of 'ID')
        string_literal = r'^([^_]+_[^_]+)' # returns all values before the second underscore

        df["Unique_ID"] = df['ID'].str.extract(string_literal, expand=False)

        # Create a dictionary to store results
        data = []
        
        # Iterate over each unique prefix and extract flattened x, y, z coordinates

        unique_ids = df["Unique_ID"].unique()

        for unique_id in unique_ids:
            filtered_df = df[df["Unique_ID"] == unique_id]
            xyz_list = filtered_df.iloc[:, 3:6].values.tolist() # gets the x,y,z coordinates of each sequence
            data.append([unique_id, xyz_list])

        # Convert to a DataFrame
        xyz_df = pd.DataFrame(data, columns=["Unique_ID", "Coordinates"])

        return xyz_df
    
    @staticmethod
    def one_hot_encode(seq):
        """
        Converts an RNA sequence into a one-hot encoded representation.

        Parameters:
            seq (str): A string consisting of RNA nucleotide characters ('A', 'C', 'G', 'U').

        Returns:
            np.ndarray: A 2D NumPy array where each nucleotide is represented as a one-hot encoded vector.
        """

        mapping = dict(zip("ACGU", range(4))) # pairs nucleotide with index   
        seq2 = [mapping[i] for i in seq] # assigns a numerical value for each nucleotide in the sequence
        return np.eye(4, dtype=int)[seq2] # returns an array of 3 zeros and a 1 based on the assigned value
    
    
    def remove_invalid_rows(self):
        """
        Removes rows that have any empty coordinates or invalid characters
        This function must be called after get_xyz_dataframe()

        Arguments:
            df (pd.Dataframe): input data frame of sequences and x,y,z coordinates

        Returns:
            cleaned dataframe with it's rows removed
        
        """

        self.labels = self.self.labels

        condition_1 = self.sequences['sequence'].str.contains('-')
        condition_2 = self.sequences['sequence'].str.contains('X')

        self.sequences = self.sequences[~(condition_1) & ~(condition_2)]

        self.labels = self.labels[self.labels['Unique_ID'].isin(self.sequences['target_id'])]

        self.labels = self.labels.groupby('Unique_ID').filter(lambda x: (x != '').all().all())

    

x = pd.read_csv('stanford-rna-3d-folding/train_labels.csv')

y = pd.read_csv('stanford-rna-3d-folding/train_sequences.csv')


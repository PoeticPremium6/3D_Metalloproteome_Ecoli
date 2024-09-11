###############
#To start this workflow, we should obtain a list of PDB IDs from our favorite model organism
#This list can then be used to query all of our proteins
#I named mine like 'rcsb_pdb_ids_Ecoli.csv'
import pandas as pd
import os
from urllib.request import urlopen
from urllib.error import HTTPError
from biopandas.pdb import PandasPdb

def download_pdb(pdb_id, download_dir="."):
    """
    Download the PDB file with the given PDB ID.

    Parameters:
    - pdb_id: The ID of the PDB to download.
    - download_dir: The directory to which the PDB file should be saved. Defaults to the current directory.
    """
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    try:
        response = urlopen(url)
        data = response.read()

        with open(os.path.join(download_dir, f"{pdb_id}.pdb"), "wb") as pdb_file:
            pdb_file.write(data)
        print(f"{pdb_id}.pdb downloaded successfully.")

    except HTTPError as e:
        # If PDB download fails, try downloading PDBx/mmCIF format
        url_cif = f"https://files.rcsb.org/download/{pdb_id}.cif"
        try:
            response = urlopen(url_cif)
            data = response.read()

            with open(os.path.join(download_dir, f"{pdb_id}.cif"), "wb") as cif_file:
                cif_file.write(data)
            print(f"{pdb_id}.cif downloaded successfully.")

        except HTTPError as e_cif:
            print(f"Error downloading both {pdb_id}.pdb and {pdb_id}.cif: {e_cif}")
# Path to your CSV file
csv_file = "C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\rcsb_pdb_ids_Ecoli.csv"
df = pd.read_csv(csv_file)

# Extract PDB Names
pdb_names = df['PDB_Names'].tolist()

download_directory = "C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\PDB_Total"
for pdb_name in pdb_names:
    download_pdb(pdb_name, download_directory)

print("Download process completed!")
####################################################
#Let's convert the CIF Files into PDB for seamless processing
import os
import gemmi

pdb_directory = "/Users/josspa/GPS-M/PDB_Total"

def rename_long_chain_ids(structure):
    """Rename chains that have identifiers too long for the PDB format."""
    for model in structure:
        for chain in model:
            if len(chain.name) > 1:
                chain.name = chain.name[0]  # Take only the first character
def cif_to_pdb(cif_path, pdb_path):
    """Converts a CIF file to PDB format using gemmi."""
    # Load the CIF file
    structure = gemmi.read_structure(cif_path)

    # Rename long chain identifiers
    rename_long_chain_ids(structure)

    # Write it out as a PDB file
    structure.write_pdb(pdb_path)
def main():
    files = [f for f in os.listdir(pdb_directory) if f.endswith('.cif')]

    for file in files:
        cif_path = os.path.join(pdb_directory, file)
        pdb_path = os.path.join(pdb_directory, file.replace('.cif', '.pdb'))
        cif_to_pdb(cif_path, pdb_path)
        print(f"Converted {file} to PDB format.")

if __name__ == "__main__":
    main()

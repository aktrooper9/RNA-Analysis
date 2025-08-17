import os
import pandas as pd

input_dir = "GSE135092_RAW"
output_file = "GSE135092_merged_raw_counts.tsv"
annotation_file = "/Users/galaxy/Documents/Version 2/Human.GRCh38.p13.annot.tsv"

# Load annotation file: EnsemblGeneID <-> GeneID
annot = pd.read_csv(annotation_file, sep='\t', dtype=str)
annot = annot[['EnsemblGeneID', 'GeneID']].dropna()
annot = annot[(annot['EnsemblGeneID'] != '') & (annot['GeneID'] != '')]  # Remove empty EnsemblGeneID or GeneID

# Get all .tsv files, sorted by sample number
files = sorted(
    [f for f in os.listdir(input_dir) if f.endswith('.tsv')],
    key=lambda x: int(x.split('_sample')[-1].split('.tsv')[0])
)

merged_df = None
expected_id_refs = None
error_files = []

for filename in files:
    file_path = os.path.join(input_dir, filename)
    sample_name = filename.split('_sample')[0]
    try:
        df = pd.read_csv(file_path, sep='\t', skiprows=3, usecols=['ID_REF', 'count'])
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        error_files.append(filename)
        continue

    # Check that ID_REFs match across files
    if expected_id_refs is None:
        expected_id_refs = df['ID_REF'].tolist()
    else:
        if not df['ID_REF'].tolist() == expected_id_refs:
            print(f"ID_REF mismatch in file: {filename}")
            error_files.append(filename)
            continue

    df = df.rename(columns={'count': sample_name})
    if merged_df is None:
        merged_df = df
    else:
        merged_df = pd.merge(merged_df, df, on='ID_REF', how='outer')

# Map ID_REF to GeneID and filter for matches
merged_df = pd.merge(annot, merged_df, left_on='EnsemblGeneID', right_on='ID_REF', how='inner')
merged_df = merged_df.drop(columns=['EnsemblGeneID', 'ID_REF'])
merged_df = merged_df.rename(columns={'GeneID': 'GENEID'})
merged_df = merged_df[['GENEID'] + [col for col in merged_df.columns if col != 'GENEID']]

# Fill missing values with 'null'
merged_df = merged_df.fillna('null')

# Save the merged dataframe
merged_df.to_csv(output_file, sep='\t', index=False)

print(f"Combined {len(files) - len(error_files)} files into {output_file}")

if error_files:
    print("The following files had errors or mismatched ID_REFs:")
    for ef in error_files:
        print(f" - {ef}")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import Counter

# List of peptides (including the initial 15 and the new ones added)
peptides = [
    "TKKAFFSQRPKKKWLCGVMTCECLQHCPWL",
    "HKKFYQTMSCKKKCGRPFWLSQRNDPGLQF",
    "HKKFRPYSQWKKKMYFVCVCGMGPALGGWL",
    "HKKPTHVQNHKKKCNCGSAMLFCVDPTTRN",
    "HKKFTMYECFKKKRNWLCWRPFDPAHCGCV",
    "HKKAPFMGHFKKKTMFAISCCAFWVAGWCT",
    "HKKGPRGLQNKKKWLFCGAWGWSCIGFNHA",
    "HKKDCPQMLDKKKPWYFFQGFCGWMMWLWL",
    "TKKPMGAWPTKKKCVRPMCDCTIFCPSHFN",
    "SRKCDPNRQFKKKCTVVKHVTCYPRNWAYA",
    "TKKCFTMFQNKKKPRWVSCQRNSYYPDPCG",
    "TKKHTMTMWVKKKCFCGFVFWYRNMGPNWA",
    "TKKPWVFNAWKKKWLRNCPRPFVWCGDPTC",
    "HKKFFPWPGHKKKAYFTMCGRNMYSHWLWF",
    "TKKPARMGRFKKKWVFRNCSNGWHVTTQCT",
    "TKKCPRWVSHKKKCFWCTCCARFAWLCTSS",
    "TKKPNPTMMGKKKMMCSQASVTCWMGPNCP",
    "TKKSDGCPTMKKKCCGWLCVTPECMQWVEF",
    "HKKCFNWPCHKKKTMFSFTCTACMGWWVWF",
    "HKKYDCTRQCKKKWFTCCTCVYCLMWDPSH"
]

# Recalculate the position counts for the updated peptide list
position_counts = {i: Counter() for i in range(30)}

# Count the amino acids at each position
for peptide in peptides:
    for i, amino_acid in enumerate(peptide):
        position_counts[i][amino_acid] += 1

# Convert the counts to a dataframe for easier plotting
df = pd.DataFrame(position_counts).fillna(0)

# Define specific hex colors for each group of amino acids
colors_dict = {
    'A': '#1f77b4', 'V': '#1f77b4', 'L': '#1f77b4', 'I': '#1f77b4', 'P': '#1f77b4', 'F': '#1f77b4', 'W': '#1f77b4', 'M': '#1f77b4', 'G': '#1f77b4',  # Nonpolar (Blue)
    'S': '#2ca02c', 'T': '#2ca02c', 'C': '#2ca02c', 'Y': '#2ca02c', 'N': '#2ca02c', 'Q': '#2ca02c',  # Polar (Green)
    'K': '#d62728', 'R': '#ff0000', 'H': '#ff7f0e',  # Positively charged (Red and Orange)
    'D': '#9467bd', 'E': '#8c564b'  # Negatively charged (Purple)
}

# Plotting the histogram with the updated peptide list and fixed hex colors
plt.figure(figsize=(15, 7))

for amino_acid in df.index:
    plt.bar(df.columns, df.loc[amino_acid], label=amino_acid, color=colors_dict.get(amino_acid, '#7f7f7f'), alpha=0.7)

# Update x-axis to 1-30
plt.xlabel('Position in Peptide')
plt.ylabel('Amino Acid Prevalence')
plt.title('Amino Acid Frequency at Each Position in Peptide Sequence')

# Grouped legend with hex colors
plt.legend(title='Amino Acid', bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2)
plt.xticks(np.arange(0, 30, step=1))  # Ensure the x-axis matches positions correctly
plt.savefig("Amino_Acid_Frequency.png")

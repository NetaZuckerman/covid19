import os
from collections import Counter

import pandas as pd

output_df = pd.DataFrame(columns=["sample", "position", "reference", "depth", "A", "G", "C", "T"])

for pileup_file in os.listdir("pileup"):
    pileup = pd.read_csv("pileup/"+pileup_file, delimiter="\t", header=None, names=["name", "position", "reference",
                                                                                    "depth", "pile", "quel"])

    for index, row in pileup.iterrows():
        base = {"A": 0, "T": 0, "C": 0, "G": 0, "deletion": 0}

        pos_count = Counter(map(str.upper, row["pile"]))
        base[row["reference"]] = pos_count[","] + pos_count["."]
        if row["depth"] == base[row["reference"]]:  # no mutations at this position
            continue

        output_df = output_df.append({"sample": pileup_file.split(".")[0],
                                      "position": row["position"],
                                      "reference": row["reference"],
                                      "depth": row["depth"],
                                      "A": base["A"] + pos_count["A"],
                                      "G": base["G"] + pos_count["G"],
                                      "C": base["C"] + pos_count["C"],
                                      "T": base["T"] + pos_count["T"],
                                      "deletion": pos_count["*"]
                                      }, ignore_index=True)

output_df.to_csv("results/quasispecies.csv", index=False)

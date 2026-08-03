#!/usr/bin/env python

import matplotlib.pyplot as plt


file = "tables/R4.tsv"
new_image_name = "DistributionR4.png"



#Record table and chart graph
x = []
y = []


with open(file, "r") as fh:
    for line in fh:
        pos, avg =line.strip().split()
        x.append(int(pos))
        y.append(float(avg))


plt.bar(x, y, color="cornflowerblue")
plt.title(f"Quality score distribution per nucleotide (R1)")
plt.xlabel("Nucleotide Position")
plt.ylabel("Average Quality Score")
plt.savefig(new_image_name)
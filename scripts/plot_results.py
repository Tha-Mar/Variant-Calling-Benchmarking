import matplotlib.pyplot as plt
import numpy as np
import os

# check that the output directory exists
os.makedirs("results/plots", exist_ok=True)

tools = ["bcftools", "freebayes"]

# Metrics as percentages 
precision = [29.0, 0.26]     # %
recall    = [98.0, 58.0]     # %
f1        = [45.0, 0.5]      # %

runtime_minutes = [ #got this from the logs 
    6 + 26/60,      # bcftools: 6 min 26 sec
    107 + 21/60     # freebayes: 1 hr 47 min 21 sec
]

max_ram_mb = [ #got this from the logs 
    111,   # bcftools
    197    # freebayes
]

def add_value_labels(ax):
    for p in ax.patches:
        height = p.get_height()
        ax.text(
            p.get_x() + p.get_width() / 2,
            height,
            f"{height:.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=0
        )

# precision
plt.figure(figsize=(6,4))
ax = plt.gca()
ax.bar(tools, precision)
ax.set_ylabel("Precision (%)")
ax.set_ylim(0, 100)
ax.set_title("Precision per Caller")
add_value_labels(ax)
plt.tight_layout()
plt.savefig("results/plots/precision.png")
plt.close()

# recall
plt.figure(figsize=(6,4))
ax = plt.gca()
ax.bar(tools, recall)
ax.set_ylabel("Recall (%)")
ax.set_ylim(0, 100)
ax.set_title("Recall per Caller")
add_value_labels(ax)
plt.tight_layout()
plt.savefig("results/plots/recall.png")
plt.close()

# f1
plt.figure(figsize=(6,4))
ax = plt.gca()
ax.bar(tools, f1)
ax.set_ylabel("F1 Score (%)")
ax.set_ylim(0, 100)
ax.set_title("F1 Score per Caller")
add_value_labels(ax)
plt.tight_layout()
plt.savefig("results/plots/f1.png")
plt.close()

# runtime
plt.figure(figsize=(6,4))
ax = plt.gca()
ax.bar(tools, runtime_minutes)
ax.set_ylabel("Runtime (minutes)")
ax.set_title("Runtime per Caller")
add_value_labels(ax)
plt.tight_layout()
plt.savefig("results/plots/runtime.png")
plt.close()

# memory usage
plt.figure(figsize=(6,4))
ax = plt.gca()
ax.bar(tools, max_ram_mb)
ax.set_ylabel("Max RAM (MB)")
ax.set_title("Memory Usage per Caller")
add_value_labels(ax)
plt.tight_layout()
plt.savefig("results/plots/memory.png")
plt.close()

# CPU usage extracted from time logs
tools = ["bcftools", "freebayes"]
cpu_percent = [128, 99]  

plt.figure(figsize=(6,4))
plt.bar(tools, cpu_percent, color=["skyblue", "salmon"])
plt.ylabel("CPU Usage (%)")
plt.title("CPU Utilization of Variant Callers")

for i, v in enumerate(cpu_percent):
    plt.text(i, v + 2, str(v) + "%", ha="center")

plt.savefig("results/plots/cpu_usage.png")
plt.close()

print("Plots saved in results/plots/")


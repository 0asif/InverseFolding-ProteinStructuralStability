import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch.nn.functional as F

def get_df_entropies(infile, AA20):
    # Assume logits is a PyTorch tensor of shape [L, 20]
    # where L = number of residues, 20 = amino acid classes
    # (from ProteinMPNN.forward() or model.sample())

    df = pd.read_csv(infile)
    df = df.drop('Unnamed: 0', axis=1)

    logits = torch.from_numpy(df.values)

    # --- Temperature scaling and probability normalization ---
    temperature = 0.1
    probs = F.softmax(logits / temperature, dim=-1)  # shape [L, 20]

    # --- Per-residue confidence metrics ---
    entropy = (-probs * torch.log(probs + 1e-12)).sum(dim=-1)           # in nats
    perplexity = torch.exp(entropy)                                     # exp(H)
    top1_prob, top1_idx = torch.max(probs, dim=-1)                      # confidence per residue
    logp = torch.log(top1_prob)                                         # log-prob of predicted AA

    # --- Build dataframe for analysis ---
    seq = "".join(AA20[i] for i in top1_idx)
    df = pd.DataFrame({
        "position": np.arange(1, len(seq) + 1),
        "top1_AA": [AA20[i] for i in top1_idx],
        "top1_prob": top1_prob.detach().numpy(),
        "entropy": entropy.detach().numpy(),
        "perplexity": perplexity.detach().numpy(),
        "log_prob": logp.detach().numpy(),
    })

    return [df, probs]

def main():
    infile = '../data/output/mutation_frequencies_trasposed.csv'
    ofolder = '../figures/'
    AA20 = list("SEMDKPTAQRVNGHLIYFWC")
    data = get_df_entropies(infile, AA20)
    df = data[0]
    probs = data[1]

    # Entropy & Probability Plots

    plt.figure(figsize=(8,4))
    plt.plot(df["position"], df["entropy"], color="tab:blue")
    plt.title("Per-residue Entropy (ProteinMPNN)", fontsize=20)
    plt.xlabel("Residue position", fontsize=20)
    plt.ylabel("Entropy (nats)", fontsize=20)
    plt.tight_layout()
    plt.savefig(ofolder + 'entropy.tiff',dpi=150)

    plt.figure(figsize=(8,4))
    plt.imshow(probs.detach().numpy().T, aspect="auto", interpolation="nearest")
    plt.yticks(np.arange(20), AA20)
    plt.xlabel("Residue position", fontsize=20)
    plt.ylabel("Amino acid", fontsize=20)
    plt.title("Residue-wise Probability Heatmap", fontsize=20)
    plt.colorbar(label="Probability")
    plt.tight_layout()
    plt.savefig(ofolder + 'probs.tiff',dpi=150)

if __name__ == "__main__":
    main()


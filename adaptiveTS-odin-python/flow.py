import argparse
import pandas as pd
import numpy as np
import torch
from nflows.distributions import StandardNormal
from nflows.flows.base import Flow
from nflows.transforms import CompositeTransform, MaskedAffineAutoregressiveTransform

def logit(x, eps=1e-6):
    r'''
    Maps data from [0,1] to the real line. Used so flow sampling does not have to be aware of bounds

    Parameters
    ----------
    x: ndarray_like
        Data
    eps: float
        clipping bounds for numerical stability
    
    Returns
    -------
    log(x/1-x)
    '''
    x = np.clip(x, eps, 1 - eps)
    return np.log(x / (1 - x))

def inv_logit(z):
    r'''
    Inverse logit transformation. Map z from real line to [0,1]

    Parameters
    ----------
    z: ndarray_like

    Returns:

    1 / (1 + exp(-z))
    '''
    return 1 / (1 + np.exp(-z))

def train_and_sample(df, nsamples, seeds):
    r'''
    Train normalizing flow and sample from it

    Parameters
    ----------
    df: pd.DataFrame
        dataframe of observations. must contain calibration parameters, weights, and seeds
    nsamples: int
        number of samples from normalizing flow
    
    Returns
    -------
    out_df: pd.DataFrame
        dataframe with sampled outputs (and seeds)

    '''
    seed = np.random.randint(low=1,high=100_000)
    torch.manual_seed(seed)
    # Load samples
    weights = df["weight"].values
    r_data = pd.Categorical(df["r"].astype(int),categories=seeds)
    x_columns = [c for c in df.columns if c not in ["r", "weight"]]

    # Logit-transform x (maps [0,1] -> R)
    x_data = df[x_columns].to_numpy()
    x_data = logit(x_data)

    # One-hot encode r
    r_onehot = pd.get_dummies(r_data).to_numpy()

    # Combine features
    x_cond = np.hstack([x_data, r_onehot])
    x_tensor = torch.tensor(x_cond, dtype=torch.float64)
    w_tensor = torch.tensor(weights / np.sum(weights), dtype=torch.float64)

    # Define flow
    n_features = x_tensor.shape[1]
    transform = CompositeTransform([
        MaskedAffineAutoregressiveTransform(
            features=n_features,
            hidden_features=64,
            num_blocks=2
        ) for _ in range(3)
    ])
    base_dist = StandardNormal([n_features])
    flow = Flow(transform, base_dist)
    optimizer = torch.optim.Adam(flow.parameters(), lr=1e-3)

    # Train
    for epoch in range(100):
        optimizer.zero_grad()
        log_prob = flow.log_prob(inputs=x_tensor)
        loss = -(w_tensor * log_prob).sum()
        loss.backward()
        optimizer.step()

    # Sample from the flow
    with torch.no_grad():
        samples = flow.sample(nsamples).numpy()

    # Split x and r back
    n_x = len(x_columns)
    x_samples = inv_logit(samples[:, :n_x])  # map back to [0,1]

    r_cols = [col for col in pd.get_dummies(r_data).columns]
    r_indices = samples[:, n_x:].argmax(axis=1)
    r_decoded = [r_cols[i] for i in r_indices]

    # Save output
    out_df = pd.DataFrame(x_samples, columns=x_columns)
    out_df["r"] = r_decoded
    #out_df.to_csv(output_csv, index=False)
    #print(f"Generated samples saved to {output_csv}")
    return out_df
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--nsamples", type=int, required=True)
    args = parser.parse_args()

    train_and_sample(args.input, args.output, args.nsamples)
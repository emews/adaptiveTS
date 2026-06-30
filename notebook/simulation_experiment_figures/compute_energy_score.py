#!/usr/bin/env python3
"""
Compute the multivariate Energy Score per method from exp1_search_space.csv
based only on predictive samples (p1, p2), given a ground truth (p1_true, p2_true).

Energy score (sample estimator):
ES_hat = (1/n) * mean_i ||x_i - y|| - (1/n^2) * sum_{i<j} ||x_i - x_j||

Where:
- x_i are the (p1, p2) samples for a method (across replicates/rows)
- y is the ground truth vector (p1_true, p2_true)

By default, groups by 'method'. You can override with --group-by (comma-separated),
e.g. --group-by method,time or method,iter to compute energy scores within finer groupings.
"""

import argparse
import sys
from typing import List

import numpy as np
import pandas as pd


def pairwise_distance_sum_euclidean(X: np.ndarray) -> float:
    """
    Return sum of Euclidean distances over all unique pairs i<j in X.
    Uses scipy.pdist if available; otherwise falls back to a memory-efficient numpy loop.
    """
    try:
        from scipy.spatial.distance import pdist  # type: ignore
        return float(pdist(X, metric="euclidean").sum())
    except Exception:
        n = X.shape[0]
        total = 0.0
        # Memory-efficient: iterate rows, accumulate distances to subsequent rows only (i<j)
        for i in range(n - 1):
            diffs = X[i] - X[(i + 1):]
            dists = np.sqrt(np.sum(diffs * diffs, axis=1))
            total += float(dists.sum())
        return total


def energy_score_for_group(X: np.ndarray, y: np.ndarray) -> dict:
    """
    Compute energy score components for a group of samples X against ground truth y.
    Returns a dict with energy_score, n, term1_mean, term2_mean.
    """
    n = X.shape[0]
    if n == 0:
        return {"energy_score": np.nan, "n": 0, "term1_mean": np.nan, "term2_mean": np.nan}

    # First term: mean distance to ground truth
    term1 = float(np.linalg.norm(X - y, axis=1).mean())

    # Second term: mean pairwise distance component
    if n > 1:
        pdist_sum = pairwise_distance_sum_euclidean(X)  # sum over i<j
        term2 = float(pdist_sum / (n * n))  # equals (1/(2n^2)) * sum_{i,j} d_ij
    else:
        term2 = 0.0

    es = term1 - term2
    return {"energy_score": es, "n": n, "term1_mean": term1, "term2_mean": term2}


def parse_group_by(arg: str) -> List[str]:
    if not arg:
        return ["method"]
    cols = [c.strip() for c in arg.split(",") if c.strip()]
    return cols or ["method"]


def main():
    parser = argparse.ArgumentParser(description="Compute Energy Score per group from (p1, p2) samples.")
    parser.add_argument(
        "--input",
        default="notebook/simulation_experiment_figures/exp1_search_space.csv",
        help="Path to exp1_search_space.csv",
    )
    parser.add_argument(
        "--output",
        default="notebook/simulation_experiment_figures/energy_scores_by_method.csv",
        help="Path to write the energy scores CSV",
    )
    parser.add_argument(
        "--p1_true",
        type=float,
        required=True,
        help="Ground truth p1 value",
    )
    parser.add_argument(
        "--p2_true",
        type=float,
        required=True,
        help="Ground truth p2 value",
    )
    parser.add_argument(
        "--group-by",
        default="method",
        help="Comma-separated columns to group by (default: method). Examples: 'method,time', 'method,iter'",
    )
    parser.add_argument(
        "--filter-methods",
        default="",
        help="Optional comma-separated list of methods to include (others will be excluded).",
    )

    args = parser.parse_args()
    group_cols = parse_group_by(args.group_by)
    methods_filter = [m.strip() for m in args.filter_methods.split(",") if m.strip()]

    # Load CSV
    try:
        df = pd.read_csv(args.input)
    except Exception as e:
        print(f"ERROR: Failed to read input CSV '{args.input}': {e}", file=sys.stderr)
        sys.exit(1)

    # Validate required columns
    required_cols = {"p1", "p2"}
    missing = required_cols - set(df.columns)
    if missing:
        print(f"ERROR: Input CSV missing required columns: {sorted(missing)}", file=sys.stderr)
        sys.exit(1)

    # Validate group-by columns
    missing_group = set(group_cols) - set(df.columns)
    if missing_group:
        print(f"ERROR: Group-by columns not found in CSV: {sorted(missing_group)}", file=sys.stderr)
        print(f"Available columns: {sorted(df.columns)}", file=sys.stderr)
        sys.exit(1)

    # Optional filter by method(s)
    if methods_filter:
        if "method" not in df.columns:
            print("ERROR: --filter-methods provided but 'method' column not in CSV.", file=sys.stderr)
            sys.exit(1)
        df = df[df["method"].isin(methods_filter)]
        if df.empty:
            print("ERROR: No rows left after filtering by methods.", file=sys.stderr)
            sys.exit(1)

    # Ground truth vector
    y = np.array([args.p1_true, args.p2_true], dtype=np.float64)

    # Group and compute energy scores
    results = []
    for keys, group in df.groupby(group_cols, sort=False):
        # keys can be scalar (single group col) or tuple (multiple)
        X = group[["p1", "p2"]].to_numpy(dtype=np.float64)
        stats = energy_score_for_group(X, y)

        # Build result row with group identifiers
        if isinstance(keys, tuple):
            row = {col: val for col, val in zip(group_cols, keys)}
        else:
            row = {group_cols[0]: keys}

        row.update(stats)
        results.append(row)

    out_df = pd.DataFrame(results)

    # Sort by energy score ascending (lower is better)
    if "energy_score" in out_df.columns:
        out_df = out_df.sort_values(by="energy_score", ascending=True, kind="mergesort")

    # Write output
    try:
        out_df.to_csv(args.output, index=False)
    except Exception as e:
        print(f"ERROR: Failed to write output CSV '{args.output}': {e}", file=sys.stderr)
        sys.exit(1)

    # Print brief summary
    print(f"Wrote energy scores for {len(out_df)} group(s) to: {args.output}")
    if "method" in out_df.columns:
        top = out_df.head(5)[["method", "energy_score", "n"]]
        print("Top (lowest energy score) methods:")
        print(top.to_string(index=False))


if __name__ == "__main__":
    main()
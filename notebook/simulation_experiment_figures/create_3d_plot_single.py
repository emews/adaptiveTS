"""
Single 3D plot utility for the simulation experiment figures.

Usage (from a notebook):
    from notebook.simulation_experiment_figures.create_3d_plot_single import plot_single_3d
    fig = plot_single_3d(dfx, method="fHet", outfile="fig1.png")

Notes:
- dfx is expected to have columns: 'method', 'rep' (or a custom rep_col), 'p1', 'p2'.
- Optional columns:
    - 'iter' (rows with iter > 1 are kept unless method == 'Initial Design')
    - 'num_replicates' (used to color points by replicate threshold)
- Exporting to PNG/SVG requires the 'kaleido' package: `pip install kaleido`
"""

import numpy as np
import plotly.graph_objects as go


def plot_single_3d(
    dfx,
    method: str = "fHet",
    rep_col: str = "rep",
    p1_col: str = "p1",
    p2_col: str = "p2",
    point_color: str = "tomato",
    replicate_threshold: int = 20,
    ground_truth: dict | None = None,
    plane_x: float = 10,
    grid_n: int = 10,
    width: int = 600,
    height: int = 720,
    dot_size: int = 5,
    eye_dict: dict | None = None,
    margin_dict: dict | None = None,
    outfile: str | None = "fig1.png",
    show: bool = True,
):
    """
    Create a single 3D plot (Scatter3d + Surface) with Plotly.

    Parameters
    ----------
    dfx : pandas.DataFrame
        DataFrame containing at least columns [method, rep_col, p1_col, p2_col].
    method : str
        The method label to filter dfx (e.g., 'fHet').
    rep_col, p1_col, p2_col : str
        Column names for x (replicates), y (p1), and z (p2).
    replicate_threshold : int
        Threshold used to color points if 'num_replicates' exists in dfx.
    ground_truth : dict | None
        Dict with keys {'p1', 'p2', 'rep'}. Defaults to {'p1': 0.45, 'p2': 0.35, 'rep': 50}.
    plane_x : float
        Constant x-position for the reference surface plane.
    grid_n : int
        Resolution of the plane grid in y/z.
    width, height : int
        Figure size in pixels.
    outfile : str | None
        File path to export image (png/svg). Set to None to skip export.
    show : bool
        Whether to call fig.show().

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The constructed figure.
    """
    if ground_truth is None:
        ground_truth = {"p1": 0.45, "p2": 0.35, "rep": 50}

    # Filter to the requested method
    dfm = dfx[(dfx["method"] == method)].copy()
    if method != "Initial Design" and "iter" in dfm.columns:
        dfm = dfm[dfm["iter"] > 1]

    if eye_dict is None:
        eye_dict=dict(  x=1,    # out to the right
                        y=-2.2, # pulled forward
                        z=0.9)  # lifted up
    else:
        eye_dict=eye_dict

    if margin_dict is None:
        margin_dict=dict(l=0, r=0, t=30, b=0)
    else:
        margin_dict = margin_dict
    
    colors = point_color
    # # Build colors based on num_replicates if present
    # if "num_replicates" in dfm.columns:
    #     colors = [point_color if rr < replicate_threshold else "tomato" for rr in dfm["num_replicates"]]
    # else:
    #     colors = "tomato"

    # Determine x-axis range
    if method.endswith("Het") and rep_col in dfx.columns and dfx[rep_col].size > 0:
        try:
            range_max = float(dfx[rep_col].max())
        except Exception:
            range_max = 55.0
    else:
        range_max = 55.0

    # Ground-truth line data
    gt_x = np.linspace(0, range_max, max(int(range_max) + 1, 2))
    gt_y = np.full_like(gt_x, float(ground_truth["p1"]), dtype=float)
    gt_z = np.full_like(gt_x, float(ground_truth["p2"]), dtype=float)

    # Reference plane at x = plane_x over y,z in [0,1]
    p1_lin = np.linspace(0, 1, grid_n)
    p2_lin = np.linspace(0, 1, grid_n)
    P1, P2 = np.meshgrid(p1_lin, p2_lin)
    X = np.full_like(P1, float(plane_x), dtype=float)

    fig = go.Figure()

    # Scatter points
    fig.add_trace(
        go.Scatter3d(
            x=dfm[rep_col],
            y=dfm[p1_col],
            z=dfm[p2_col],
            mode="markers",
            marker=dict(size=dot_size, color=colors, opacity=1),
            showlegend=False,
        )
    )

    # Ground-truth point
    fig.add_trace(
        go.Scatter3d(
            x=[ground_truth["rep"]],
            y=[ground_truth["p1"]],
            z=[ground_truth["p2"]],
            mode="markers",
            marker=dict(size=4, color="black"),
            showlegend=False,
        )
    )

    # Ground-truth line
    fig.add_trace(
        go.Scatter3d(
            x=gt_x,
            y=gt_y,
            z=gt_z,
            mode="lines",
            marker=dict(size=2, color="black"),
            showlegend=False,
        )
    )

    # Plane at x = plane_x
    fig.add_trace(
        go.Surface(
            x=X,
            y=P1,
            z=P2,
            surfacecolor=np.zeros_like(P1),
            colorscale=[[0, "cornflowerblue"], [1, "cornflowerblue"]],
            opacity=0.5,
            showscale=False,
            showlegend=False,
        )
    )

    camera = dict(
        eye=eye_dict,
        up=dict(x=0, y=0, z=1),
    )

    fig.update_layout(
        width=width,
        height=height,
        margin=margin_dict,
        paper_bgcolor="white",
        font = dict(size=24),
        scene=dict(
            bgcolor="white",
            aspectmode="cube",
            # x-axis (replicates), descending
            xaxis=dict(
                range=[range_max, -1],
                title="replicate",
                showgrid=True,
                zeroline=False,
                showline=True,
                mirror=True,
                linewidth=2,
                linecolor="black",
                showticklabels=True,
                backgroundcolor="white",
                tickfont = dict(size=14),
            ),
            # y-axis (p1)
            yaxis=dict(
                title="β",
                range=[0, 1],
                showgrid=False,
                zeroline=False,
                showline=True,
                mirror=True,
                linewidth=2,
                linecolor="black",
                showticklabels=False,
                backgroundcolor="light gray",
            ),
            # z-axis (p2)
            zaxis=dict(
                title="γ",
                range=[0, 1],
                showgrid=False,
                zeroline=False,
                showline=True,
                mirror=True,
                linewidth=2,
                linecolor="black",
                showticklabels=False,
                backgroundcolor="light gray",
            ),
            camera=camera,
        ),
        annotations=[
            dict(
                text="",
                x=0.5,
                y=0.85,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=24),
            )
        ],
    )

    if show:
        fig.show()

    if outfile:
        try:
            fig.write_image(outfile, scale = 2)
        except Exception as e:
            # Writing images requires 'kaleido' to be installed
            print(f"Warning: write_image failed: {e}. Install 'kaleido' to enable image export.")

    return fig


if __name__ == "__main__":
    # This module defines plot_single_3d(dfx, ...).
    # Import and call it from your notebook:
    #   from notebook.simulation_experiment_figures.create_3d_plot_single import plot_single_3d
    #   fig = plot_single_3d(dfx, method="fHet", outfile="fig1.png")
    print("Module ready: import plot_single_3d(dfx, ...) from your notebook.")

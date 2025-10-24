# Adaptive Grid Thompson Sampling for Trajectory-Oriented Optimization

This repository contains the source code to reproduce the analysis presented in an upcoming paper on Adaptive Grid Thompson Sampling for Trajectory-Oriented Optimization (TOO).

## Project Overview

This project implements and evaluates an Adaptive Grid Thompson Sampling algorithm for Trajectory-Oriented Optimization. The primary surrogate model used is a Common Random Number Gaussian Process (GP).

The repository includes:
*   A full bake-off experiment comparing different GP surrogates within the Thompson Sampling framework, using a compartmental model as the simulation model.
*   An application of the Adaptive Thompson Sampling algorithm to CityCOVID, a detailed agent-based model for simulating the spread of COVID-19 in the city of Chicago.

The code is designed to be used with the EMEWS (Extreme-scale Model Exploration with Swift) framework to manage large-scale computational experiments.

## Paper

arxiv: [https://arxiv.org/abs/2510.18099](https://arxiv.org/abs/2510.18099)

## Directory Structure

*   `emews/`: Contains files related to the EMEWS framework, including Swift scripts, R scripts for running the model, and data files.
*   `experiments/`: Contains configuration files for experiments.
*   `notebook/`: Contains Jupyter and R Markdown notebooks for analysis and visualization of results.
*   `plots/`: Contains plots generated from the analysis of experimental results.
*   `script/`: Contains various R and Python scripts for analysis and helper functions.


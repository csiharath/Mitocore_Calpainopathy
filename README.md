# Mitocore_calpainopathy <!-- omit in toc -->

## Table of content <!-- omit in toc -->

- [Introduction](#introduction)
- [This project was implemented in python using COBRApy methods such as:](#this-project-was-implemented-in-python-using-cobrapy-methods-such-as)
- [Installation](#installation)
- [Usage](#usage)
  - [Notebooks](#notebooks)
  - [Library](#library)
- [Model and data description](#model-and-data-description)
  - [MitoCore](#mitocore)
- [Results](#results)


## Introduction

Energy metabolism is a complex process involving the conversion and use of energy in the body, and can undergo modifications depending on conditions (age, activity, health...). Neuromuscular diseases are a group of disorders affecting the motor unit composed of a motor neuron, its axon and all the muscle fibers it innervates, and, as a side consequence, affecting the metabolism of muscle cells, responsible of contractions. Modelling metabolic networks aids in understanding changes happening due to physical activity, and to take account of the known regulations caused by different factors, such as diseases.This project aims to implement a method to build a constraint-based computer model, illustrating altered energy metabolism in neuromuscular diseases, taking the example of calpainopathy.


This project was implemented in python using COBRApy methods such as:
- 

## Installation 

First, you need to have conda/mamba/anaconda

Then, you can follow the code below:

```bash
ENV_NAME="mitoenv"
# Renme it if you prefer

conda env create -n $ENV_NAME -f requirements.yml 
# This will create a new env with all depedencies specified in 'requirements.yml'
# Use 'mamba' instead of 'conda' or append '--solver libmamba', this will speed up the process

conda run -n $ENV_NAME python -m ipykernel install --user --name $ENV_NAME
# This will install the env for Jupyter

conda activate $ENV_NAME
```

## Usage
### Notebooks

- Atleration_model.ipynb

Calculates new bounds from known regulations and muscle activity of healthy cells, and creates a new model representing the studied neuromuscular disease.

```bash
jupyter notebook code/Alteration_model.ipynb
```

- Ateration_model_study.ipynb
  
Studies different set of paramaters in order to find possible compensation.

```bash
jupyter notebook code/Alteration_model_study.ipynb
```

- comparison_control.ipynb

Compares control and alterations models (with the same uptakes parameters, but could be change according results found with `Alteration_model.ipynb`)
```bash
jupyter notebook code/comparison_control.ipynb
```

### Library

- metabolic_analysis_library.py
  
Library of functions used to build and analyse healthy and deseased models, specifically during physical exercise.
This library contains functions to run FBA or different method of analysis (pFBA, FVA) on models during physical exercise, build a dysregulated model, build data visualisation, and plot results.
It can be used in a python scritp or notebook by using the following code:

```python
import metabolic_analysis_library as mal
```

## Model and data description

### MitoCore

This study was carried out using the [Mitocore](https://doi.org/10.1186/s12918-017-0500-7) model, which was originally build to represent the metabolism of a human heart cell mitochondrion and therefore does not correspond exactly to the model needed. MitoCore includes 481 reactions covering all parts of central metabolism.

## Results

Here, comparison of healthy and calpainopathy models reveals a compromised aerobic metabolism in diseased cells, but compensatory pathways like glycolysis and increased amino acids inputs offer alternatives to counter the decreased ATP production.

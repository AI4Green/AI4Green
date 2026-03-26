## Understanding Reaction Predictions

### What are Reaction SMILES?

SMILES (Simplified Molecular Input Line Entry System) is a chemical notation language specifically designed for fast and flexible input, storage, and retrieval of molecular structures on computer systems (Weininger, 1988).

SMILES encodes a molecule’s structure as a simple ASCII string using a few basic rules that represent atoms, bonds, branches, ring closures, and disconnected structures. Reaction SMILES are a way to represent chemical reactions using text.

Reactants, reagents, and products are written using SMILES notation. A '.' separates reactants while '>>' separates reactants from the products. ChemDraw can be used to generate SMILES.

<pre>CCO.C=O>>CCOC=O</pre>

### Prediction Models

Machine learning (ML) is a subset of artificial intelligence where algorithms learn patterns from data and improve their predictions or decisions without being explicitly programmed (Hey et al., 2020). ML techniques, especially deep learning, have made breakthroughs in fields like image recognition and natural language processing by leveraging vast amounts of data and computational power. In scientific contexts, machine learning is being increasingly used to automate data analysis and enable new discoveries from complex, large-scale experimental datasets.

Our system uses machine learning models trained on chemical reaction data derived from the United States Patent and Trademark Office (USPTO) database to predict the most likely products of reactions. The model takes reactant SMILES as input and returns 5 Product SMILES and their corresponding confidence scores as output, which are then translated to a full reaction scheme and displayed as 2D structures. The Top-N5 accuracy metric was used, and it was observed that the first prediction is often accurate. However, this is not always the case, which highlights the importance of validating predictions against both the literature and Reaxy’s. This approach can enhance understanding of reaction mechanisms and support exploratory data analysis during the planning phase.

### Understanding Scores

ML models can be evaluated using Top 10 or Top 5 accuracy. This means that the correct prediction will be found within the first 10 or 5 predicted outputs. Each predicted product is given a score, indicating the confidence of the prediction. A higher score means the model is more confident that the predicted product is correct.

For Product Prediction Model, the highest score is -0.01 and the lowest is -10.

### How to Read Results

- __High score:__ Very high confidence prediction.
- __Medium score:__ Likely correct, but double-check.
- __Low score:__ Be cautious; alternative products may exist.

### References

1. Hey, T., Butler, K., Jackson, S., & Thiyagalingam, J. (2020). Machine learning and big scientific data. Philosophical Transactions of the Royal Society A, 378(2166), 20190054. https://doi.org/10.1098/rsta.2019.0054

1. Weininger, D. (1988). SMILES, a chemical language and information system. 1. Introduction to methodology and encoding rules. Journal of Chemical Information and Computer Sciences, 28(1), 31–36. https://doi.org/10.1021/ci00057a005

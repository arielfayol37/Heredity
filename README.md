# Gene Inference using Bayesian Network

This repository contains Python code to make inferences about the gene distribution and trait expression in a population based on a Bayesian Network model. The model takes into account the probability of different gene distributions and trait expressions based on genetic inheritance and mutations.

## Getting Started

1. Clone the repository to your local machine.
2. Make sure you have Python installed (Python 3.x).
3. Open the `heredity.py` file in your preferred code editor.

## Understanding the Data

The data directory contains sample data sets in CSV format. Each CSV file represents a family, and the columns are defined as follows:

- `name`: The name of the individual.
- `mother`: The mother's name of the individual (if available).
- `father`: The father's name of the individual (if available).
- `trait`: Indicates whether the individual exhibits the trait (hearing impairment) or not (1 for True, 0 for False, and empty if unknown).

## Model Definitions

The `PROBS` dictionary in the `heredity.py` file contains probability distributions for various events:

- `PROBS["gene"]`: Represents the unconditional probability distribution over the gene (i.e., the probability if we know nothing about that person's parents).
- `PROBS["trait"]`: Represents the conditional probability that a person exhibits a trait based on the gene they possess.
- `PROBS["mutation"]`: Represents the probability that a gene mutates from being the gene in question to not being that gene, and vice versa.

## Functions

### `joint_probability(people, one_gene, two_genes, have_trait)`

Compute the joint probability of gene distributions and trait expressions for a given family.

Parameters:
- `people`: A dictionary of people as described in the "Understanding" section. The keys represent names, and the values are dictionaries containing mother, father, and trait information.
- `one_gene`: A set of people with one copy of the gene in the current joint distribution.
- `two_genes`: A set of people with two copies of the gene in the current joint distribution.
- `have_trait`: A set of people with the trait in the current joint distribution.

Returns:
- The joint probability of all of those events taking place.

### `update(probabilities, one_gene, two_genes, have_trait, p)`

Add the newly computed joint probability to the existing probability distributions in probabilities.

Parameters:
- `probabilities`: A dictionary of people as described in the "Understanding" section. Each person is mapped to a "gene" distribution and a "trait" distribution.
- `one_gene`: A set of people with one copy of the gene in the current joint distribution.
- `two_genes`: A set of people with two copies of the gene in the current joint distribution.
- `have_trait`: A set of people with the trait in the current joint distribution.
- `p`: The probability of the joint distribution.

Returns:
- None. The function updates the probabilities dictionary.

### `normalize(probabilities)`

Normalize a dictionary of probabilities such that each probability distribution is normalized (i.e., sums to 1, with relative proportions the same).

Parameters:
- `probabilities`: A dictionary of people as described in the "Understanding" section. Each person is mapped to a "gene" distribution and a "trait" distribution.

Returns:
- None. The function updates the probabilities dictionary.

## Example Usage

1. Open one of the sample data sets in the `data` directory (e.g., `data/family0.csv`) to see the family information.
2. Run the program with appropriate inputs to make inferences about the gene distribution and trait expression.


Let's make gene inferences with the Bayesian Network model!

import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    zero_gene = set([name for name in people if name not in one_gene and name not in two_genes])
    one_gene_p = {} # key: value --> name: probability of having one copy of the gene
    two_genes_p = {} # key: value --> name: probability of having two copies of the gene
    zero_gene_p = {}
    
    # Computing the probability of person having one copy of the gene.
    for name in one_gene:
        # A person has one copy of the gene either from the mother AND NOT the father, 
        # or from the father AND NOT the mother. In case we have no info about parents, use
        # unconditional probability.
        father, mother = people[name]["father"], people[name]["mother"]
        if father == None or mother == None: # TODO: and should give the same result 
                                             # as or per the problem description
            one_gene_p[name] = PROBS["gene"][1]
            continue
        
        poppa_p, momma_p = person_prob(father, one_gene, two_genes), person_prob(mother, one_gene, two_genes)
        one_gene_p[name] = poppa_p * (1 - momma_p)  +  momma_p * (1 - poppa_p)
        

    # Computing the probability of person having two copies of the gene.
    for name in two_genes:
        # A person has two copies of the gene, one from BOTH parents.
        # Note, the parents may have no copies but the copy is mutated.
        father, mother = people[name]["father"], people[name]["mother"]
        if father == None or mother == None: # TODO: and should give the same result 
                                             # as or per the problem description
            two_genes_p[name] = PROBS["gene"][2]
            continue
        poppa_p, momma_p = person_prob(father, one_gene, two_genes), person_prob(mother, one_gene, two_genes)
        two_genes_p[name] = poppa_p * momma_p

    # Computing the probability of person having zero copies of the gene.
    for name in zero_gene:
        # A person has zero copy of the gene if BOTH parents gave zero copies.
        # This can happen if parent has zero copy of gene, or via mutation.
        father, mother = people[name]["father"], people[name]["mother"]
        if father == None or mother == None: # TODO: and should give the same result 
                                             # as or per the problem description
            zero_gene_p[name] = PROBS["gene"][0]
            continue
        poppa_p, momma_p = person_prob(father, one_gene, two_genes), person_prob(mother, one_gene, two_genes)
        zero_gene_p[name] = (1-poppa_p) * (1 - momma_p)

    # Computing the join probability of all the events. 
    joint_p = 1
    for index, prob_dict in enumerate([zero_gene_p, one_gene_p, two_genes_p]): # !Important: Must be in this order!
        joint_p *= trait_prob(prob_dict, index, have_trait)

    return joint_p

def person_prob(person_name, one_gene_names, two_genes_names):
    """
    Helper function for joint probability.
    It takes person_name (father or mother), and the context (one_gene people and two gene)
    and returns the probability that one copy is gotten from person_name if gives_copy == True else
    returns that no copy is gotten from person_name
    """
   
    if person_name in two_genes_names:
        # gives_copy == True: Probability that the name in the iterator in the function above this one
        # gets a copy from person_name, given that the person_name has two copies 
        person_p = 0.99 
                        
    elif person_name in one_gene_names:
        person_p = 0.5 * 0.99 
    else: 
        person_p = PROBS["mutation"]
    return person_p

def trait_prob(prob_dict, num_copies, have_trait):
    """
    Helper function for joint probability.
    Return value: joint_probability over prob_dict only. 
        Should be multiply with other joint probability for all probability dicts.
    Takes a dictionary of the form name:probability, num_copies (associated with prob_dict), and have_trait(iterable) names.
    E.g one_gene_p is a dictionary with names and the probability they have 1 copy. So the function can be called this way:
    trait_prob(one_gene_p, 1, have_trait)
    """
    product_ = 1
    for name in prob_dict:
        if name in have_trait:
            product_ *= prob_dict[name] * PROBS['trait'][num_copies][True]
        else:
            product_ *= prob_dict[name] * PROBS['trait'][num_copies][False]  
    return product_  

def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    for person in probabilities:
        if person in have_trait:
            probabilities[person]["trait"][True] += p
        else:
            probabilities[person]["trait"][False] += p

        if person in one_gene:
            probabilities[person]["gene"][1] += p
              
        elif person in two_genes:
            probabilities[person]["gene"][2] += p
        else: # No copy of the gene
            probabilities[person]["gene"][0] += p


def get_dist_sum(dist):
    """
    Takes a dictionary(dist) and returns the sum of its values
    """
    sum_ = 0
    for item in dist:
        sum_ += dist[item]
    return sum_

def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for person_name in probabilities:
        for dist_name in probabilities[person_name]: # dist_name == distribution_name (e.g "gene")
            sum_ = get_dist_sum(probabilities[person_name][dist_name])
            for item in probabilities[person_name][dist_name]:
                probabilities[person_name][dist_name][item] /= sum_



if __name__ == "__main__":
    main()



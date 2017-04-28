__author__ = 'Aaron'
def gini(vector_of_n):
    h, inequality = 0., 0.
    for n in vector_of_n:
        h += n
        inequality += h - n / 2.
    equality = h * len(vector_of_n) / 2.
    return (equality - inequality) / equality

from symnmf import *
from analysis import *
import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def main():
    k=0
    file_name=""
    elements = []
    num_of_elements=0
    d=0

    if len(sys.argv) == 4:
        if sys.argv[1].isdigit():
            k = int(sys.argv[1])
        else:
            print("Invalid number of clusters!")
            exit(0)
        file_name = str(sys.argv[3])
        goal = str(sys.argv[2])
    else:
        print("An Error Has Occurred")
        exit(0)

    file = pd.read_csv(file_name, header=None)

    elements = file.values
    elements = elements[:,0:2]

    num_of_elements=len(elements)
    d=len(elements[0])

    if k<1 or k>num_of_elements:
        print("An Error Has Occurred")
        exit(0)

    if goal=="sym":
        printMatrix(mf.sym(elements.tolist(), num_of_elements, d), num_of_elements, num_of_elements)
    elif goal=="ddg":
        printMatrix(mf.ddg(elements.tolist(), num_of_elements, d), num_of_elements, num_of_elements)
    elif goal=="norm":
        printMatrix(mf.norm(elements.tolist(), num_of_elements, d), num_of_elements, num_of_elements)
    elif goal=="symnmf":
        W = mf.norm(elements.tolist(), num_of_elements, d)
        H = initialize_H(W, num_of_elements, k)
        new_H = mf.symnmf(H.tolist(), W, k, num_of_elements)
        printMatrix(mf.symnmf(H.tolist(), W, k, num_of_elements), num_of_elements, k)
    else:
        print("An Error Has Occurred")
        exit(0)

    symnmf_clusters = symnmfClusterAssign(new_H, num_of_elements)

    f = plt.figure(1, figsize=[10,10])
    x = [d[0] for d in elements]
    y = [d[1] for d in elements]
    plt.scatter(x,y, c=symnmf_clusters)

    plt.xticks(np.arange(1, 11, step=1))
    plt.savefig("elbow.png")
    plt.show()

if __name__ == "__main__":
    main()



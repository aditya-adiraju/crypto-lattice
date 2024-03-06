from gram_schmidt import gram_schmidt

"""
An implementation of the LLL latice-reduction algorithm
INPUT: 
    -- Original Lattice Basis 
OUTPUT:
    -- New LLL-reduced Lattice basis
"""

def LLL(basis, DELTA=0.75):

    LLL_basis = [False for i in range(len(basis))]

    while True:
        lovasz_condition = True
        for i in range(1, len(basis)):
            # from j = 0 to i - 1
            for j in range(i):
                # reduce basis[i] by basis[i-1]
                gs_basis = gram_schmidt(basis[:j+1])
                u_ij = proj_length(basis[i], gs_basis[j])
                basis[i] = basis[i] - (round(u_ij)*basis[j])

            lovasz_condition = basis[i].dot_product(basis[i]) >= (DELTA - (proj_length(basis[i], basis[i-1]))^2)*basis[i-1].dot_product(basis[i-1])
            if not lovasz_condition:
                tmp = basis[i]
                basis[i] = basis[i-1]
                basis[i-1] = tmp
                break
        
        if lovasz_condition:
            return basis

"""
returns the length of the projection of v1 onto v2
"""
def proj_length(v1, v2):
    return v2.dot_product(v1)/v2.dot_product(v2)
    




            




    


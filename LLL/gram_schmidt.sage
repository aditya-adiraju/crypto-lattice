"""
An implementation of the Gram Schmidt Orthogonalization process.

INPUTS:
    -- basis : list[vector] all must be same dimension
OUTPUT:
    -- orthogonal_basis: list[vector]
"""

def gram_schmidt(basis):
    # STEP 0: Verify wheter they are all the same dimension
    assert len(basis) > 0
    dim = len(basis[0])
    assert all(len(v) == dim for v in basis)

    orthogonal_basis = []
    # Step 1
    orthogonal_basis.append(basis[0])

    # Step 2
    for i in range(1, len(basis)):
        # Find \Sigma_{j = 1}^{i -1}((b_i \dot b_j*)/(b_j* \dot b_j*)) * b_j*
        orthogonal_projection = vector([0 for k in range(dim)])
        for j in range(0, i):
            orthogonal_projection += (basis[i].dot_product(basis[j]) / basis[j].dot_product(basis[j]) ) * basis[j]

        # Find b_i^*
        orthogonal_basis.append(basis[i] - orthogonal_projection)

    return orthogonal_basis
        

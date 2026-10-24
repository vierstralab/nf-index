from typing import List


def get_matrices_mapping_by_types(matrices: List[str], matrices_types):
    matrices_mapping = {}
    for matrix in matrices:
        for matrix_type in matrices_types:
            if matrix.endswith(f'.{matrix_type}.npy') or matrix.endswith(f'{matrix_type}.matrix.npy'):
                matrices_mapping[matrix_type] = matrix
                break
        else:
            raise ValueError(f"Matrix {matrix} not recognized. Expected one of {matrices_types}")
    return matrices_mapping

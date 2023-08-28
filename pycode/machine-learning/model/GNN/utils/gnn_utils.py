import tensorflow as tf
from functools import partial 
from spektral.utils.convolution import _incidence_matrix_single
import os 
import numpy as np 



def incidence_matrix(adjacency):
    """
    Creates the corresponding incidence matrices for graphs with particular
    adjacency matrices.
    :param adjacency: The binary adjacency matrices. Should have shape
        ([batch], n_nodes, n_nodes).
    :return: The computed incidence matrices. It will have a shape of
        ([batch], n_nodes, n_edges).
    """
    adjacency = tf.convert_to_tensor(adjacency, dtype=tf.float32)
    added_batch = False
    if tf.size(tf.shape(adjacency)) == 2:
        # Add the extra batch dimension if needed.
        adjacency = tf.expand_dims(adjacency, axis=0)
        added_batch = True

    # Compute the maximum number of edges. We will pad everything in the
    # batch to this dimension
    ######## I changed this line because my graph has directed edges and a unsymmetrical adjacency matrix
    #adjacency_upper = _triangular_adjacency(adjacency)   
    num_edges = tf.math.count_nonzero(adjacency, axis=(1, 2))
    max_num_edges = tf.reduce_max(num_edges)

    # Compute all the transformation matrices.
    make_single_matrix = partial(_incidence_matrix_single, num_edges=max_num_edges)
    transformation_matrices = tf.map_fn(
        make_single_matrix,
        adjacency,
        fn_output_signature=tf.TensorSpec(shape=[None, None], dtype=tf.float32),
    )

    if added_batch:
        # Remove the extra batch dimension before returning.
        transformation_matrices = transformation_matrices[0]
    return transformation_matrices

def getBaselineMatrix():
    """! loads the baselinematrix
    """

    baseline_contact_matrix0 = os.path.join(
        "./data/contacts/baseline_home.txt")
    baseline_contact_matrix1 = os.path.join(
        "./data/contacts/baseline_school_pf_eig.txt")
    baseline_contact_matrix2 = os.path.join(
        "./data/contacts/baseline_work.txt")
    baseline_contact_matrix3 = os.path.join(
        "./data/contacts/baseline_other.txt")

    baseline = np.loadtxt(baseline_contact_matrix0) \
        + np.loadtxt(baseline_contact_matrix1) + \
        np.loadtxt(baseline_contact_matrix2) + \
        np.loadtxt(baseline_contact_matrix3)

    return baseline


#def local_normalization():

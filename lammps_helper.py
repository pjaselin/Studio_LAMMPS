import numpy as np

def c_to_np(carray_input):
    '''
    Description: Converts a LAMMPS C array to a NumPy array (only for reading RDF data)
    Input:
        carray_input: lmp variable with C array data
    
    Output:
        result: NumpPy array
    '''
    # create an empty NumPy array
    result = np.array(np.zeros((50,2)))
    
    # loop through C array and enter into NumPy array
    for i in range(0,50):
        result[i,0] = carray_input[i][0]
        result[i,1] = carray_input[i][1]
    
    # return NumPy array
    return(result)

__author__ = "github.com/pjaselin"
__date__ = "7/30/2018"
__version__ = "1.2"
__license__ = "GPLv3"

'''
lammps_helper.py

Description: Python file of helper functions to handle C-level LAMMPS-Python
wrapper.

This file is part of Studio LAMMPS.

Studio LAMMPS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Studio LAMMPS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Studio LAMMPS.  If not, see <https://www.gnu.org/licenses/>.
'''

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

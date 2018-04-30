import numpy as np

def c_to_np(carray_input):
  result = np.array(np.zeros((50,2)))
  for i in range(0,50):
    result[i,0] = carray_input[i][0]
    result[i,1] = carray_input[i][1]
  
  return(result)

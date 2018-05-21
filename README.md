# ID3-implementation
Implementation of ID3, in MATLAB.


## Build the ID3:
No builds necessary, just open MATLAB interpreter in the root folder of this project.
But the data to train the tree must have an specific format... Basically:

- Import the data to train the tree into a numerical matrix.
- Make sure the "class" to classify must be on the first column of the matrix of data, so the rest of the columns should be features.

## How to use it

To generate the tree just create a new  object `id3`
```
  data = csvread('fisher.csv')
  example = id3
  example.createTree(data)

```

### Additional notes:

The Tree is generated, the find method to test the model in to be implemented... soon.

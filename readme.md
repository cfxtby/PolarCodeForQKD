
# The construction of polar code by upgrading and degrading channels

**From:** 
This algorithm is from the paper "How To Construct Polar Code", is construction of polar code by upgrading and degrading channels.

**Advantage:** I have implemented the algorithm with C++. It can be used with the frame length more than 1Mb. Meanwhile, the implementation supports to calculate the final value with multiple processes so as to save the calculation time.

**Instructions:** 
The code have several instructions: 
-N: the block length $$n$$ is $$n=2^N$$ 
-mil: $$\mu$$ 
-p: the error bits rate is $$0.01p$$ 

**Output files:** 
The output file contians two columns: the indexes of each channel and the estimated probability of error under maximum-likelihood decision.

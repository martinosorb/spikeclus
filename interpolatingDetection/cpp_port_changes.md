# C++ (cython) port changes:
* 3 files (.pyx, .cpp, .h)    
* GUI removed
* GUI output redirected to console
* Input variables are now hardcoded (provisional)
* StreamWriter/FileStream are now a single fstreams using std::iostream (<<)

* Notation: Math.max/min functions -> max/min
* Notation: declaration [][] -> **
* Notation: decimal -> standard float (equivalent?)
* Notation: access [,] -> [][]
* Notation: foreach (type x in y) -> for (;;;) (c++11 "for (x : X)" cannot be used here)
* Notation(c++11): Array.Sort(x) -> std::sort(std::begin(x), std::end(x))
* Notation(optional): ! -> not 

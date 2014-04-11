=======
HALOGEN
=======

HALOGEN is a C program for creating realistic synthetic dark matter halo catalogs 
in a flash. 

It uses a statistical approach for generating the masses of the halos from a 
parameterised halo mass distribution, and placing them down in a spatial 
distribution with cheap 2LPT field as scaffolding.


Installation
------------
HALOGEN has no external code dependencies, and therefore, compiling is as simple
as typing ``make`` in the top-level directory of the desired version. 

The program will be compiled into the same top-level directory as ``halogen``.

To optionally install, save ``halogen`` somewhere on the system path, eg. 
in ``/usr/bin/``.


Usage
-----
To run HALOGEN, type ``halogen <input.file>``.
  
An example input file is found in the input/ subdirectory of each version.

Notably, several files must be provided to the input:

* 2LPT file: A gadget file, which is output of the 2LPT run. This could also
  be the output of a full NBODY run with gadget.
  
* Output file: A file in which the output halo masses and positions will be written.

* alpha-file: A file containing a table of M and alpha, which encodes the 
  relative biasing of halos of different masses. For constant bias, this should
  be one row long.
  
* hmf file: A file specifying the tabulated mass function with columns m,n(>m).
  This file can be generated in any way, but is most easily done with the 
  native support at hmf.icrar.org
  
  
Acknowledgments
---------------
This code will soon have a paper to cite. For now, just let us know you are
using the code.

Authors
-------
Santiago Avila Perez: santiagoavilaperez@gmail.com
Steven Murray: steven.murray@uwa.edu.au 

 
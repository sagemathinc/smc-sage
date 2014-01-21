r"""
Free module bases.

The class :class:`FreeModuleBasis` implements bases over a free module `M`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version

EXAMPLES:

"""
#******************************************************************************
#       Copyright (C) 2014 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2014 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.sage_object import SageObject

class FreeModuleBasis(SageObject):
    r""" 
    Basis of a free module over a commutative ring.
    
    INPUT:
    
    - ``fmodule`` -- free module `M` (must be an instance of 
      :class:`GenFreeModule`)
    - ``symbol`` -- a letter (of a few letters) to denote a generic element of
      the basis
    - ``latex_symbol`` -- (default: None) symbol to denote a generic element of
      the basis; if None, the value of ``symbol`` is used. 

    EXAMPLES:

    """
    def __init__(self, fmodule, symbol, latex_symbol=None):
        from free_module_tensor import FreeModuleVector
        self.fmodule = fmodule
        self.name = "(" + \
          ",".join([symbol + "_" + str(i) for i in self.fmodule.irange()]) +")"
        if latex_symbol is None:
            latex_symbol = symbol
        self.latex_name = r" ,\left(" + \
          ",".join([latex_symbol + "_" + str(i) 
                    for i in self.fmodule.irange()]) + r"\right)"
        # The basis is added to the module list of bases; moreover the first 
        # defined basis is considered as the default one
        for other in self.fmodule.known_bases:
            if repr(self) == repr(other):
                raise ValueError("The " + str(self) + " already exist on the " +
                                 str(self.fmodule))
        self.fmodule.known_bases.append(self)
        if self.fmodule.def_basis is None:
            self.fmodule.def_basis = self
        # The individual vectors:
        vl = list()
        for i in self.fmodule.irange():
            v_name = symbol + "_" + str(i)
            v_symb = latex_symbol + "_" + str(i)
            v = FreeModuleVector(self.fmodule, v_name, v_symb)
            for j in self.fmodule.irange():
                v.set_comp(self)[j] = 0
            v.set_comp(self)[i] = 1
            vl.append(v)
        self.vec = tuple(vl)

    def _repr_(self):
        r"""
        String representation of the object.
        """
        return "basis " + self.name + " on the " + str(self.fmodule)

    def _latex_(self):
        r"""
        LaTeX representation of the object.
        """
        return self.latex_name

    def __hash__(self):
        r"""
        Hash function (since instances of :class:`FreeModuleBasis` are used as
        dictionary keys).
        """
        return id(self)

    def __eq__(self, other):
        r"""
        Comparison operator
        """
        return other is self
        
    def __getitem__(self, index):
        r"""
        Returns the basis element corresponding to a given index.
        
        INPUT:
        
        - ``index`` -- the index of the basis element

        """
        n = self.fmodule.rank
        si = self.fmodule.sindex
        i = index - si
        if i < 0 or i > n-1:
            raise ValueError("Index out of range: " +
                              str(i+si) + " not in [" + str(si) + "," +
                              str(n-1+si) + "]")
        return self.vec[i]

    def __len__(self):
        r"""
        Return the basis length, i.e. the rank of the free module.
        """
        return self.fmodule.rank


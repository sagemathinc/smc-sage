r"""
Homomorphisms between free modules of finite rank

The class :class:`FiniteRankFreeModuleMorphism` implements homomorphisms
between two free modules of finite rank over the same commutative ring. 

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version

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

from sage.categories.morphism import Morphism
from sage.categories.homset import Hom
from finite_rank_free_module import FiniteRankFreeModule

class FiniteRankFreeModuleMorphism(Morphism):
    r"""
    Homomorphism between free modules of finite rank
    
    This is an *element* class, whose *parent* class is
    :class:`~sage.tensor.modules.free_module_homset.FreeModuleHomset`.

    INPUT:
    
    - ``parent`` -- hom-set Hom(M,N) to which the homomorphism belongs
    - ``name`` -- (string; default: None) name given to the homomorphism
    - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
      to the homomorphism; if None, ``name`` will be used. 
    
    """

    def __init__(self, parent, matrix_rep, basis1=None, basis2=None, 
                 name=None, latex_name=None):
        r"""
        TESTS::
        """
        from sage.matrix.constructor import matrix
        Morphism.__init__(self, parent)
        fmodule1 = parent.domain()
        fmodule2 = parent.codomain()
        if basis1 is None:
            basis1 = fmodule1.default_basis()
        elif basis1 not in fmodule1.bases():
            raise TypeError(str(basis1) + " is not a basis on the " + \
                            str(fmodule1) + ".")
        if basis2 is None:
            basis2 = fmodule2.default_basis()
        elif basis2 not in fmodule2.bases():
            raise TypeError(str(basis2) + " is not a basis on the " + \
                            str(fmodule2) + ".")
        ring = fmodule1.base_ring()
        n1 = fmodule1.rank()
        n2 = fmodule2.rank()
        self._matrices = {(basis1, basis2): matrix(ring, n2, n1, matrix_rep)}
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
    

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
    
    - ``fmodule1`` -- domain of the homomorphism; must be
      an instance of 
      class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``fmodule2`` -- codomain of the homomorphism; must be
      an instance of 
      class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``name`` -- (string; default: None) name given to the homomorphism
    - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
      to the homomorphism; if None, ``name`` will be used. 
    
    """

    def __init__(self, fmodule1, fmodule2, name=None, latex_name=None):
        r"""
        TESTS::
        """
        if not isinstance(fmodule1, FiniteRankFreeModule):
            raise TypeError("The argument fmodule1 must be a free module of" + 
                            " finite rank.")
        if not isinstance(fmodule2, FiniteRankFreeModule):
            raise TypeError("The argument fmodule2 must be a free module of" + 
                            " finite rank.")
        parent = Hom(fmodule1, fmodule2) 
        Morphism.__init__(self, parent)
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
    

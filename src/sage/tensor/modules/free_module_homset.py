r"""
Set of homomorphisms between free modules of finite rank

The class :class:`FreeModuleHomset` implements sets (actually free modules) of
homomorphisms between two free modules of finite rank over the same 
commutative ring. 

Given two free modules `M` and `N` of respective ranks `m` and `n` over a 
commutative ring `R`, the set `\mathrm{Hom}(M,N)` of homomorphisms 
`M\rightarrow N` is a free module of rank `mn` over `R`. 

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

from sage.categories.homset import Homset
from free_module_morphism import FiniteRankFreeModuleMorphism

class FreeModuleHomset(Homset):
    r"""
    Set of homomorphisms between free modules of finite rank
    
    Given two free modules `M` and `N` of respective ranks `m` and `n` over a 
    commutative ring `R`, the class :class:`FreeModuleHomset` implements the 
    set `\mathrm{Hom}(M,N)` of homomorphisms `M\rightarrow N`. It is a 
    *parent* class, whose *elements* are instances of 
    class:`~sage.tensor.modules.free_module_morphism.FiniteRankFreeModuleMorphism`

    INPUT:
    
    - ``fmodule1`` -- free module `M` (domain of the homomorphisms); must be
      an instance of 
      class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``fmodule2`` -- free module `N` (codomain of the homomorphisms); must be
      an instance of 
      class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``name`` -- (string; default: None) name given to the hom-set; if None, 
      Hom(M,N) will be used
    - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
      hom-set; if None, `\mathrm{Hom}(M,N)` will be used

    """

    Element = FiniteRankFreeModuleMorphism

    def __init__(self, fmodule1, fmodule2, name=None, latex_name=None):
        r"""
        TESTS::
        """
        Homset.__init__(self, fmodule1, fmodule2)
        if name is None:
            self._name = "Hom(" + fmodule1._name + "," + fmodule2._name + ")"
        else:
            self._name = name
        if latex_name is None:
            self._latex_name = \
                    r"\mathrm{Hom}\left(" + fmodule1._latex_name + "," + \
                    fmodule2._latex_name + "\right)"
        else:
            self._latex_name = latex_name
        
    def _latex_(self):
        r"""
        LaTeX representation of the object.
        
        EXAMPLES::
        
        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name
        

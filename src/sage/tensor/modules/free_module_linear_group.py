"""
General linear group of a free module

AUTHORS:

- Eric Gourgoulhon (2015): initial version

"""
#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.groups import Groups
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.tensor.modules.free_module_automorphism import FreeModuleAutomorphism


class FreeModuleLinearGroup(UniqueRepresentation, Parent):
    r"""
    General linear group of a free module.

    INPUT:

    - ``fmodule`` -- free module `M` of finite rank (must be an instance of
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`)

    """
    
    Element = FreeModuleAutomorphism
    
    def __init__(self, fmodule):
        r"""
        See :class:`FreeModuleLinearGroup` for documentation and examples.

        TESTS::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: from sage.tensor.modules.free_module_linear_group import FreeModuleLinearGroup
            sage: GL = FreeModuleLinearGroup(M) ; GL
            General linear group of the Rank-3 free module M over the Integer Ring
            sage: GL.category()
            Category of groups
            sage: TestSuite(GL).run()

        """
        if not isinstance(fmodule, FiniteRankFreeModule):
            raise TypeError("{} is not a free module of finite rank".format(
                            fmodule))
        Parent.__init__(self, category=Groups())
        self._fmodule = fmodule
        self._one = None # to be set by self.one()
 
    #### Parent methods ####

    def _element_constructor_(self, comp=[], basis=None, name=None,
                              latex_name=None):
        r"""
        Construct a free module automorphism.

        """
        if comp == 1:
            return self.one()
        # standard construction
        resu = self.element_class(self._fmodule, name=name,
                                  latex_name=latex_name)
        if comp:
            resu.set_comp(basis)[:] = comp
        return resu


    def _an_element_(self):
        r"""
        Construct some specific free module automorphism.

        """
        resu = self.element_class(self._fmodule)
        if self._fmodule._def_basis is not None:
            comp = resu.set_comp()
            for i in self._fmodule.irange():
                if i%2 == 0:
                    comp[[i,i]] = self._fmodule._ring.one()
                else:
                    comp[[i,i]] = -(self._fmodule._ring.one())
        return resu
        
    #### End of parent methods ####

    #### Monoid methods ####

    def one(self):
        r"""
        Return the group unit element.
        """
        if self._one is None:
            self._one = self.element_class(self._fmodule, is_identity=True)
            if self._fmodule.bases():
                self._one.components(self._fmodule.bases()[0]) # initializes
                                                # the components in some basis
        return self._one

    #### End of monoid methods ####

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        """
        return "General linear group of the {}".format(self._fmodule)

    def _latex_(self):
        r"""
        Return a string representation of ``self``.

        """
        from sage.misc.latex import latex
        return r"\mathrm{GL}\left(" + latex(self._fmodule) + r"\right)"


    def base_module(self):
        r"""
        Return the free module on which ``self`` is constructed.

        OUTPUT:

        - instance of :class:`FiniteRankFreeModule` representing the free
          module of which ``self`` is the general linear group

        """
        return self._fmodule

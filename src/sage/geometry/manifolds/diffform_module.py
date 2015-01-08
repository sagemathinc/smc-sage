r"""
Differential form modules

The set of differential forms along an open subset `U` of some manifold `S`
with values in a open subset `V` of a manifold `M` (possibly `S=M` and `U=V`)
is a module over the algebra `C^\infty(U)` of differentiable scalar fields
on `U`. It is a free module iff `V` is parallelizable.
Accordingly, two classes are devoted to differential form modules:

- :class:`DiffFormModule` for differential forms with values in a generic (in
  practice, not parallelizable) open set `V`
- :class:`DiffFormFreeModule` for differential forms with values in a
  parallelizable open set `V`

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
from sage.categories.modules import Modules
from sage.tensor.modules.ext_pow_free_module import ExtPowerFreeModule
from diffform import DiffForm, DiffFormParal

class DiffFormModule(UniqueRepresentation, Parent):
    r"""
    Module of differential forms of a given degree `p` along an open subset `U`
    of some manifold `S` with values in a open subset `V` of a manifold `M`.

    This is a module over `C^\infty(U)`, the ring (algebra) of differentiable
    scalar fields on `U`.

    The standard case of differential forms *on* a manifold corresponds to
    `U=V` (and hence `S=M`). Another common case is `\Phi` being an
    immersion.

    If `V` is parallelizable, the class :class:`DiffFormFreeModule` should
    be used instead.

    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.tensor.geometry.manifolds.DiffForm`.

    INPUT:

    - ``vector_field_module`` -- module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` associated with the mapping `\Phi:\; U \rightarrow V`.
    - ``degree`` -- positive integer; the degree `p` of the differential forms

    """

    Element = DiffForm

    def __init__(self, vector_field_module, degree):
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        name = "/\^{}(".format(degree) + domain._name
        latex_name = r"\Lambda^{" + str(degree) + r"}\left(" + \
                     domain._latex_name
        if dest_map is domain._identity_map:
            name += ")"
            latex_name += r"\right)"
        else:
            name += "," + dest_map._name + ")"
            latex_name += "," + dest_map._latex_name + r"\right)"
        self._vmodule = vector_field_module
        self._degree = degree
        self._name = name
        self._latex_name = latex_name
        # the member self._ring is created for efficiency (to avoid calls to
        # self.base_ring()):
        self._ring = domain.scalar_field_algebra()
        Parent.__init__(self, base=self._ring, category=Modules(self._ring))
        self._domain = domain
        self._dest_map = dest_map
        self._ambient_domain = vector_field_module._ambient_domain
        # NB: self._zero_element is not constructed here, since no element
        # can be constructed here, to avoid some infinite recursion.

    #### Parent methods

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None):
        r"""
        Construct a differential form.

        """
        if comp == 0:
            if not hasattr(self, '_zero_element'):
                self._zero_element = self._element_constructor_(name='zero',
                                                                latex_name='0')
                for frame in self._domain._frames:
                    if self._dest_map.restrict(frame._domain) == \
                                                               frame._dest_map:
                        self._zero_element.add_comp(frame)
                        # (since new components are initialized to zero)
            return self._zero_element
        if isinstance(comp, DiffForm):
            if self._degree == comp._tensor_type[1] and \
               self._domain.is_subset(comp._domain) and \
               self._ambient_domain.is_subset(comp._ambient_domain):
                return comp.restrict(self._domain)
            else:
                raise TypeError("Cannot coerce the {} ".format(comp) +
                                "to a differential form in {}".format(self))
        resu = self.element_class(self._vmodule, self._degree, name=name,
                                  latex_name=latex_name)
        if comp != []:
            resu.set_comp(frame)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unamed) differential form.

        """
        resu = self.element_class(self._vmodule, self._degree)
        #!# a zero element is returned
        return resu

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.
        
        """
        if isinstance(other, (DiffFormModule, DiffFormFreeModule)):
            return self._degree == other._degree and \
                   self._domain.is_subset(other._domain) and \
                   self._ambient_domain.is_subset(other._ambient_domain)
        else:
            return False

    #### End of Parent methods

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        """
        description = "Module "
        if self._name is not None:
            description += self._name + " "
        description += "of {}-forms ".format(self._degree)
        if self._dest_map is self._domain._identity_map:
            description += "on the {}".format(self._domain)
        else:
            description += "along the {} mapped into the {}".format(
                                            self._domain, self._ambient_domain)
        return description

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.
        
        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def base_module(self):
        r"""
        Return the vector field module on which ``self`` is constructed.

        OUTPUT:

        - instance of class
          :class:`~sage.geometry.manifolds.vectorfield_module.VectorFieldModule`
          representing the module on which the tensor module is defined.

        """
        return self._vmodule


#******************************************************************************

class DiffFormFreeModule(ExtPowerFreeModule):
    r"""
    Module of differential forms of a given degree `p` along an open subset `U`
    of some manifold `S` with values in a parallelizable open subset `V` of
    a manifold `M`.

    Since `V` is parallelizable, the module is a free module over `C^\infty(U)`,
    the ring (algebra) of differentiable scalar fields on `U`.

    The standard case of differential forms *on* a manifold corresponds to
    `U=V` (and hence `S=M`). Another common case is `\Phi` being an
    immersion.

    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.tensor.geometry.manifolds.DiffFormParal`.

    INPUT:

    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` associated with the mapping `\Phi:\; U \rightarrow V`.
    - ``degree`` -- positive integer; the degree `p` of the differential forms

    """

    Element = DiffFormParal

    def __init__(self, vector_field_module, degree):
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        name = "/\^{}(".format(degree) + domain._name
        latex_name = r"\Lambda^{" + str(degree) + r"}\left(" + \
                     domain._latex_name
        if dest_map is domain._identity_map:
            name += ")"
            latex_name += r"\right)"
        else:
            name += "," + dest_map._name + ")"
            latex_name += "," + dest_map._latex_name + r"\right)"
        ExtPowerFreeModule.__init__(self, vector_field_module, degree,
                                    name=name, latex_name=latex_name)
        self._domain = domain
        self._dest_map = dest_map
        self._ambient_domain = vector_field_module._ambient_domain

    #### Parent methods

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None, sym=None, antisym=None):
        r"""
        Construct a differential form.

        """
        if comp == 0:
            return self._zero_element
        if isinstance(comp, DiffForm):
            if self._degree == comp._tensor_type[1] and \
               self._domain.is_subset(comp._domain) and \
               self._ambient_domain.is_subset(comp._ambient_domain):
                return comp.restrict(self._domain)
            else:
                raise TypeError("Cannot coerce the {} ".format(comp) +
                                "to a differential form in {}".format(self))
        resu = self.element_class(self._fmodule, self._degree, name=name,
                                  latex_name=latex_name)
        if comp != []:
            resu.set_comp(frame)[:] = comp
        return resu

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.
        
        """
        if isinstance(other, (DiffFormModule, DiffFormFreeModule)):
            return self._degree == other._degree and \
                   self._domain.is_subset(other._domain) and \
                   self._ambient_domain.is_subset(other._ambient_domain)
        else:
            return False

    #### End of Parent methods

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        """
        description = "Free module "
        if self._name is not None:
            description += self._name + " "
        description += "of {}-forms ".format(self._degree)
        if self._dest_map is self._domain._identity_map:
            description += "on the {}".format(self._domain)
        else:
            description += "along the {} mapped into the {}".format(
                                            self._domain, self._ambient_domain)
        return description











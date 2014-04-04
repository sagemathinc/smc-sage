r"""
Vector field module. 


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

from sage.tensor.modules.finite_free_module import FiniteFreeModule
from scalarfield import ScalarField
from vectorfield import VectorFieldParal

class VectorFieldFreeModule(FiniteFreeModule):
    r"""
    Module of vector fields along an open subset `U` of some manifold `S`
    with values in a parallelizable open subset `V` of a manifold `M`. 
    
    Given a differential mapping

    .. MATH::

        \Phi:\ U\subset S \longrightarrow V\subset \mathcal{M}
    
    the module `\mathcal{X}(U,\Phi)` is the set of all vector fields of 
    the type

    .. MATH::

        v:\ U  \longrightarrow TM
        
    such that 

    .. MATH::

        \forall p \in U,\ v(p) \in T_{\Phi(p)}M
        
    
    Since `V` is parallelizable, the `\mathcal{X}(U,\Phi)` is a free module 
    over `C^\infty(U)`, the ring of differentiable scalar fields on `U`.
    Its rank is the dimension of `M`. 
    
    The standard case of vector fields *on* a manifold corresponds to `S=M`, 
    `U=V` and `\Phi = \mathrm{Id}`. 

    Another common case is `\Phi` being an immersion.
    
    INPUT:
    
    - ``domain`` -- open subset `U` on which the vector fields are defined
    - ``dest_map`` -- (default: None) destination map `\Phi:\ U \rightarrow V` 
      (type: :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`); 
      if none is provided, the identity is assumed (case of vector fields *on* 
      `U`)
    
    """
    
    Element = VectorFieldParal

    def __init__(self, domain, dest_map=None):
        self.domain = domain
        name = "X(" + domain.name
        latex_name = r"\mathcal{X}\left(" + domain.latex_name
        if dest_map is None:
            self.dest_map = None
            self.ambient_domain = domain
            name += ")" 
            latex_name += r"\right)"
        else:
            self.dest_map = dest_map
            self.ambient_domain = dest_map.codomain
            name += "," + self.dest_map.name + ")" 
            latex_name += "," + self.dest_map.latex_name + r"\right)" 
        manif = self.ambient_domain.manifold
        FiniteFreeModule.__init__(self, domain.scalar_field_ring(), 
                                  manif.dim, name=name, latex_name=latex_name, 
                                  start_index=manif.sindex,
                                  output_formatter=ScalarField.function_chart)

    #### Methods to be redefined by derived classes of FiniteFreeModule ####

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "free module "
        if self.name is not None:
            description += self.name + " "
        description += "of vector fields "
        if self.dest_map is None:
            description += "on the " + str(self.domain)
        else:
            description += "along the " + str(self.domain) + \
                           " mapped into the " + str(self.ambient_domain)
        return description

    def tensor_module(self, k, l):
        r"""
        Return the free module of all tensors of type (k,l) defined on 
        ``self``. 
        
        INPUT: 
        
        - ``k`` -- (non-negative integer) the contravariant rank, the tensor type 
          being (k,l)
        - ``l`` -- (non-negative integer) the covariant rank, the tensor type 
          being (k,l)
        
        OUTPUT:

        - instance of 
          :class:`TensorFieldFreeModule` 
          representing the free module of type-`(k,l)` tensors on the 
          free module ``self``. 
        
        EXAMPLES:
        
        
        """
        from tensorfield_module import TensorFieldFreeModule
        if (k,l) not in self._tensor_modules:
            self._tensor_modules[(k,l)] = TensorFieldFreeModule(self, (k,l))
        return self._tensor_modules[(k,l)]

    def basis(self, symbol=None, latex_symbol=None, from_frame=None):
        r""" 
        Define a basis (vector frame) of the free module.
        
        If the basis specified by the given symbol already exists, it is
        simply returned.
        If no argument is provided the module's default basis is returned. 
        
        INPUT:
        
        - ``symbol`` -- (string; default: None) a letter (of a few letters) to 
          denote a generic element of the basis; if None and ``from_frame=None`` 
          the module's default basis is returned.
        - ``latex_symbol`` -- (string; default: None) symbol to denote a 
          generic element of the basis; if None, the value of ``symbol`` is 
          used. 
        - ``from_frame`` -- (default: None) vector frame `\tilde e` on the 
          codomain `V` of the destination map `\Phi` of ``self``; the returned
          basis `e` is then such that 
          `\forall p \in U, e(p) = \tilde e(\Phi(p))`

        OUTPUT:
        
        - instance of 
          :class:`~sage.geometry.manifolds.vector.VectorFrame` 
          representing a basis on ``self``.
        
        EXAMPLES:
            
        """
        from vectorframe import VectorFrame
        if symbol is None:
            return self.default_basis()
        else:
            for other in self.known_bases:
                if symbol == other.symbol:
                    return other
            return VectorFrame(self.domain, symbol=symbol, 
                               latex_symbol=latex_symbol, dest_map=self.dest_map)

    def tensor(self, tensor_type, name=None, latex_name=None, sym=None, 
               antisym=None):
        r"""
        Construct a tensor on the free module. 

        The tensor is actually a tensor field on the domain of ``self``. 
        
        INPUT:
        
        - ``tensor_type`` -- pair (k,l) with k being the contravariant rank 
          and l the covariant rank
        - ``name`` -- (string; default: None) name given to the tensor
        - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
          tensor; if none is provided, the LaTeX symbol is set to ``name``
        - ``sym`` -- (default: None) a symmetry or a list of symmetries among 
          the tensor arguments: each symmetry is described by a tuple 
          containing the positions of the involved arguments, with the 
          convention position=0 for the first argument. For instance:

          * sym=(0,1) for a symmetry between the 1st and 2nd arguments 
          * sym=[(0,2),(1,3,4)] for a symmetry between the 1st and 3rd
            arguments and a symmetry between the 2nd, 4th and 5th arguments.

        - ``antisym`` -- (default: None) antisymmetry or list of antisymmetries 
          among the arguments, with the same convention as for ``sym``. 
          
        OUTPUT:
        
        - instance of 
          :class:`~sage.geometry.geometry.tensorfield.TensorFieldParal` 
          representing the tensor defined on ``self`` with the provided 
          characteristics.
          
        EXAMPLES:
                    
        See :class:`~sage.geometry.geometry.tensorfield.TensorFieldParal` 
        for more examples and documentation.
                
        """
        from tensorfield import TensorFieldParal
        from rank2field import EndomorphismFieldParal, SymBilinFormFieldParal
        from diffform import DiffFormParal, OneFormParal
        if tensor_type==(1,0):
            return VectorFieldParal(self, name=name, latex_name=latex_name)
        elif tensor_type==(0,1):
            return OneFormParal(self, name=name, latex_name=latex_name)
        elif tensor_type==(1,1):
            return EndomorphismFieldParal(self, name=name, 
                                                         latex_name=latex_name)
        elif tensor_type==(0,2) and sym==(0,1):
            return SymBilinFormFieldParal(self, name=name, 
                                                         latex_name=latex_name)
        elif tensor_type[0]==0 and tensor_type[1]>1 and antisym is not None:
            if len(antisym)==tensor_type[1]:
                return DiffFormParal(self, tensor_type[1], name=name, 
                                                         latex_name=latex_name)
            else:
                return TensorFieldParal(self, tensor_type, name=name, 
                                        latex_name=latex_name, sym=sym, 
                                        antisym=antisym)
        else:
            return TensorFieldParal(self, tensor_type, name=name, 
                                    latex_name=latex_name, sym=sym, 
                                    antisym=antisym) 

    def tensor_from_comp(self, tensor_type, comp, name=None, latex_name=None):
        r"""
        Construct a tensor on the free module from a set of components.
        
        The tensor is actually a tensor field on the domain of ``self``. 
        The tensor symmetries are deduced from those of the components.
        
        INPUT:
        
        - ``tensor_type`` -- pair (k,l) with k being the contravariant rank and l 
          the covariant rank
        - ``comp`` -- instance of :class:`~sage.tensor.modules.comp.Components` 
          representing the tensor components in a given basis
        - ``name`` -- (string; default: None) name given to the tensor
        - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the tensor; 
          if none is provided, the LaTeX symbol is set to ``name``
          
        OUTPUT:
        
        - instance of 
          :class:`~sage.geometry.manifolds.tensorfield.TensorFieldParal` 
          representing the tensor defined on ``self`` with the provided 
          characteristics.
          
        EXAMPLES:
                        
        """
        from tensorfield import TensorFieldParal
        from rank2field import EndomorphismFieldParal, SymBilinFormFieldParal
        from diffform import DiffFormParal, OneFormParal
        from sage.tensor.modules.comp import CompWithSym, CompFullySym, \
                                                               CompFullyAntiSym
        #
        # 0/ Compatibility checks:
        if comp.ring is not self.ring:
             raise TypeError("The components are not defined on the same" + 
                            " ring as the module.")           
        if comp.frame not in self.known_bases:
            raise TypeError("The components are not defined on a basis of" + 
                            " the module.")
        if comp.nid != tensor_type[0] + tensor_type[1]:
            raise TypeError("Number of component indices not compatible with "+
                            " the tensor type.")
        #
        # 1/ Construction of the tensor:
        if tensor_type == (1,0):
            resu = VectorFieldParal(self, name=name, latex_name=latex_name)
        elif tensor_type == (0,1):
            resu = OneFormParal(self, name=name, latex_name=latex_name)
        elif tensor_type == (1,1):
            resu = EndomorphismFieldParal(self, name=name, 
                                          latex_name=latex_name)
        elif tensor_type == (0,2) and isinstance(comp, CompFullySym):
            resu = SymBilinFormFieldParal(self, name=name, 
                                          latex_name=latex_name)
        elif tensor_type[0] == 0 and tensor_type[1] > 1 and \
                                        isinstance(comp, CompFullyAntiSym):
            resu = DiffFormParal(self, tensor_type[1], name=name, 
                                                         latex_name=latex_name)
        else:
            resu = TensorFieldParal(self, tensor_type, name=name, 
                                    latex_name=latex_name) 
            # Tensor symmetries deduced from those of comp:
            if isinstance(comp, CompWithSym):
                resu.sym = comp.sym
                resu.antisym = comp.antisym
        #
        # 2/ Tensor components set to comp:
        resu.components[comp.frame] = comp
        #
        return resu

    def alternating_form(self, degree, name=None, latex_name=None):
        r"""
        Construct an alternating form on the free module. 
        
        INPUT:
    
        - ``degree`` -- the degree of the alternating form (i.e. its tensor rank)
        - ``name`` -- (string; default: None) name given to the alternating 
          form
        - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
          alternating form; if none is provided, the LaTeX symbol is set to 
          ``name``
          
        OUTPUT:
        
        - instance of 
          :class:`~sage.geometry.manifolds.diffform.DiffFormParal` 
          (``degree`` > 1) or 
          :class:`~sage.geometry.manifolds.diffform.OneFormParal` 
          (``degree`` = 1)

        EXAMPLES:
        
        
        See 
        :class:`~sage.geometry.manifolds.diffform.DiffFormParal` 
        for further documentation. 

        """
        from diffform import DiffFormParal, OneFormParal
        if degree == 1:
            return OneFormParal(self, name=name, latex_name=latex_name)
        else:
            return DiffFormParal(self, degree, name=name, 
                                 latex_name=latex_name)


    def linear_form(self, name=None, latex_name=None):
        r"""
        Construct a linear form on the free module. 
        
        A linear form on the vector free module ``self`` is actually a field
        of linear forms (i.e. a 1-form) along the open subset `U` on which 
        ``self`` is defined.

        INPUT:
    
        - ``name`` -- (string; default: None) name given to the linear 
          form
        - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
          linear form; if none is provided, the LaTeX symbol is set to 
          ``name``
          
        OUTPUT:
        
        - instance of 
          :class:`~sage.geometry.manifolds.diffform.OneFormParal` 

        EXAMPLES:
        

        """
        from diffform import OneFormParal
        return OneFormParal(self, name=name, latex_name=latex_name)

    def endomorphism(self, name=None, latex_name=None):
        r"""
        Construct an endomorphism on the free module ``self``.
        
        An endomorphism on the vector free module ``self`` is actually a field
        of tangent-space endomorphisms along the open subset `U` on which 
        ``self`` is defined.
        
        INPUT:
    
        - ``name`` -- (string; default: None) name given to the endomorphism
        - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
          endomorphism; if none is provided, the LaTeX symbol is set to 
          ``name``
          
        OUTPUT:
        
        - instance of 
          :class:`~sage.geometry.manifolds.rank2field.EndomorphismFieldParal`
          
        EXAMPLES:

   
        See 
        :class:`sage.geometry.manifolds.rank2field.EndomorphismFieldParal` 
        for further documentation. 

        """
        from rank2field import EndomorphismFieldParal
        return EndomorphismFieldParal(self, name=name, latex_name=latex_name)

    def automorphism(self, name=None, latex_name=None):
        r"""
        Construct an automorphism on the free module ``self``.
        
        An automorphism on the vector free module ``self`` is actually a field
        of tangent-space automorphisms along the open subset `U` on which 
        ``self`` is defined.
        
        INPUT:
    
        - ``name`` -- (string; default: None) name given to the automorphism
        - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
          automorphism; if none is provided, the LaTeX symbol is set to 
          ``name``
          
        OUTPUT:
        
        - instance of 
          :class:`~sage.geometry.manifolds.rank2field.AutomorphismFieldParal`
          
        EXAMPLES:

   
        See 
        :class:`sage.geometry.manifolds.rank2field.AutomorphismFieldParal` 
        for further documentation. 

        """
        from rank2field import AutomorphismFieldParal
        return AutomorphismFieldParal(self, name=name, latex_name=latex_name)

    def identity_map(self, name='Id', latex_name=None):
        r"""
        Construct the identity map on the free module ``self``. 
        
        The identity map on the vector free module ``self`` is actually a field
        of tangent-space identity maps along the open subset `U` on which 
        ``self`` is defined.
        
        INPUT:
    
        - ``name`` -- (string; default: 'Id') name given to the identity map
        - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
          identity map; if none is provided, the LaTeX symbol is set to 
          ``name``
          
        OUTPUT:
        
        - instance of 
          :class:`~sage.geometry.manifolds.rank2field.IdentityMapParal`
          
        EXAMPLES:


        See
        :class:`~sage.geometry.manifolds.rank2field.IdentityMapParal` 
        for further documentation. 
 
        """
        from rank2field import IdentityMapParal
        return IdentityMapParal(self, name=name, latex_name=latex_name)

    def sym_bilinear_form(self, name=None, latex_name=None):
        r"""
        Construct a symmetric bilinear form on the free module ``self``.
        
        A symmetric bilinear form on the vector free module ``self`` is 
        actually a field of tangent-space symmetric bilinear forms along 
        the open subset `U` on which ``self`` is defined.
        
        INPUT:
    
        - ``name`` -- (string; default: None) name given to the automorphism
        - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
          automorphism; if none is provided, the LaTeX symbol is set to 
          ``name``
          
        OUTPUT:
        
        - instance of 
          :class:`~sage.geometry.manifolds.rank2field.SymBilinFormFieldParal`
          
        EXAMPLES:

   
        See 
        :class:`sage.geometry.manifolds.rank2field.SymBilinFormFieldParal` 
        for further documentation. 

        """
        from rank2field import SymBilinFormFieldParal
        return SymBilinFormFieldParal(self, name=name, latex_name=latex_name)

    #### End of methods to be redefined by derived classes of FiniteFreeModule ####

